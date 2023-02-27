import argparse
import json
from multiprocessing.sharedctypes import Value
import os,io,time
import pickle
import numpy as np
import sys

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from scripts.data.preprocessing.cleanStructs import extractBackbone
from terminator.data.data import BinderScoringIterableDataset,ComplexScoringDataset
from terminator.models.TERMinator import TERMinator
from terminator.utils.model.loop_utils import run_epoch,_to_dev
from terminator.utils.model.loss_fn import construct_loss_fn
from terminator.utils.common import int_to_3lt_AA,AA_to_int,int_to_AA

# from scripts.score_binders.hbond_utils import bbHBond

# pylint: disable=unspecified-encoding

def get_binder_names(pathToListFile):
    if pathToListFile == '':
        return None
    with open(pathToListFile,'r') as file:
        names = set([line.strip() for line in file])
    return names

# Functions for computing the score of a binder 
class interfaceScorer:
    def __init__(self, target_data, binder_complex_data, model, dev='cuda:0'):
        self.model = model
        self.dev = dev

        # Set model to eval mode
        self.model.eval()
        torch.set_grad_enabled(False)
    
        self.score_file = None
        self.structure_score_file = None
        self.seq_file = None
        self.pep_seq_prob_file = None

        self.set_new_target(target_data)
        self.set_new_binder_and_complex(binder_complex_data)

    def set_new_target(self,target_data):
        # target is always a single structure
        self.target_data = target_data
        self.target_netout = None
        self.unbound_target_aa_probs = None
        self.unbound_target_seq_probs = None
        # self.new_seq = None

    def set_new_binder_and_complex(self,binder_complex_data):
        '''Note: assumes that the target structure is unchanged'''
        # binder and complex data are lists 
        self.binder_complex_data = binder_complex_data
        self.binder_netout = None
        self.complex_netout = None
        # self.new_seq = None
        # self.unbound_binder_aa_probs = None
        # self.unbound_binder_seq_probs = None
        # self.bound_aa_probs = None
        # self.bound_seq_probs = None
        # self.unbound_seq_probs = None

        self.binder_chain_ids = [chain_id for idx,chain_id in enumerate(self.binder_complex_data['binder_chain_id']) if idx % 2 == 0]

    def get_pdb_name(self,idx):
        return self.complex_netout[idx]['id']

    def generate_etabs(self):
        # run COORDinator to generate etabs
        if self.target_netout == None:
            self.target_netout = self.run_model(self.target_data, self.dev)[0]
        binder_complex_netout = self.run_model(self.binder_complex_data, self.dev) #alternating netout: 1) binder, 2) complex
        print(f"Ran model to generate network output for {len(binder_complex_netout)} structures")
        self.binder_netout = [x for idx,x in enumerate(binder_complex_netout) if idx % 2 == 0]
        self.complex_netout = [x for idx,x in enumerate(binder_complex_netout) if idx % 2 == 1]
        print(f"binder_netout = {len(self.binder_netout)}, complex_netout = {len(self.complex_netout)}")

        # get the residue info
        (self.protein_res,
        self.peptide_res,
        self.prot_res_info,
        self.pep_res_info) = self.get_res_idx([x['res_info'] for x in self.complex_netout],self.binder_chain_ids)
        self.prot_res_range = [(x.min(),x.max()+1) for x in self.protein_res]
        self.pep_res_range =  [(x.min(),x.max()+1) for x in self.peptide_res]
        # print('prot_res',self.protein_res)
        # print('pep_res',self.peptide_res)
        
    def run_model(self, data, dev):
        """ Run :code:`model` on data from :code:`dataloader`

        Note: this function uses the masks to gather 

        Args
        ----
        model : terminator.model.TERMinator.TERMinator
            An instance of TERMinator
        data : dict 
            Input data from IterableDataset
        dev : str, default="cuda:0"
            What device to compute on

        Returns
        -------
        dump : list of dicts
        """
        # Move data to the GPU
        _to_dev(data, dev)
        max_seq_len = max(data['seq_lens'].tolist())
        try:
            etab, E_idx, sscore = self.model(data, max_seq_len)
            print('etab: ',etab.shape)
            print('E_idx:',E_idx.shape)
            print('sscore: ',sscore.shape)
        except Exception as e:
            print(data['ids'])
            raise e
        dump_list = []
        for idx in range(etab.shape[0]):
            l, k = etab[idx].shape[0:2]

            # etab and E_idx were padded to fit in the batch, need to remove the extra positions
            l_crop = data['seq_lens'][idx].item()
            # print('l',l)
            # print('l_crop',l_crop)
            # print(idx,'etab_shape',etab[idx].shape)
            # print(idx,'E_idx_shape',E_idx[idx].shape)
            # print('etab after slice',etab[idx][:l_crop,:min(l_crop,k)].shape)
            # print('E_idx after slice',E_idx[idx][:l_crop,:min(l_crop,k)].shape)

            # print('E_idx before slice',E_idx)

            dump = {
                'id': data['ids'][idx],
                # 'etab': etab[idx].view(l, k, 20, 20).cpu(),
                # 'E_idx': E_idx[idx].cpu(),
                'etab': etab[idx][:l_crop,:min(l_crop,k)].view(l_crop, min(l_crop,k), 20, 20).cpu(),
                'E_idx': E_idx[idx][:l_crop,:min(l_crop,k)].cpu(),
                'sscore': sscore[idx][:l_crop].cpu(),
                'seq':data['seqs'][idx][:l_crop].cpu(),
                'res_info':data['res_info'][idx]
            }
            # print('id',dump['id'])
            # print('etab',dump['etab'].shape)
            # print('E_idx',dump['E_idx'].shape)
            # print('seq',dump['seq'].shape)
            # print('res_info',len(dump['res_info']))
            dump_list.append(dump)
        return dump_list

    @staticmethod
    def get_res_idx(res_info_list, peptide_chain=None):
        protein_res_list = list()
        peptide_res_list = list()
        prot_res_info_list = list()
        pep_res_info_list = list()
        if peptide_chain == None:
            peptide_chain = [''] * len(res_info_list)
        for idx,res_info in enumerate(res_info_list):
            protein_res = list()
            peptide_res = list()
            prot_res_info = list()
            pep_res_info = list()
            for i,(chain,resnum) in enumerate(res_info):
                if chain == peptide_chain[idx]:
                    peptide_res.append(i)
                    pep_res_info.append(res_info[i])
                else:
                    protein_res.append(i)
                    prot_res_info.append(res_info[i])
            if peptide_chain[idx] != '' and len(peptide_res) == 0:
                raise ValueError('No residues found matching peptide chain: '+peptide_chain[idx])
            protein_res_list.append(torch.tensor(protein_res))
            peptide_res_list.append(torch.tensor(peptide_res))
            prot_res_info_list.append(prot_res_info)
            pep_res_info_list.append(pep_res_info)
        return protein_res_list, peptide_res_list, prot_res_info_list, pep_res_info_list

    def get_binder_scores(self,seq_mode,score_mode):
        self.bound_aa_probs, self.new_seq = [], []
        self.unbound_binder_aa_probs = []
        self.unbound_binder_seq_probs = []
        self.bound_seq_probs, self.unbound_seq_probs = [], []
        self.binder_scores = []
        self.unbound_structure_scores = []
        self.unbound_structure_score_targettotal = None
        self.bound_structure_scores = []

        if self.target_netout == None or self.complex_netout == None:
            self.generate_etabs()

        for idx in range(len(self.complex_netout)):
            print(f"scoring {self.get_pdb_name(idx)} in {seq_mode} and {score_mode} mode")
            
            # get the unbound target amino acid probabilities
            if self.unbound_target_aa_probs == None:
                self.unbound_target_aa_probs,_,_ = interfaceScorer.get_aa_prob(self.target_netout['etab'],self.target_netout['E_idx'],self.target_netout['seq'],torch.arange(0,self.target_netout['etab'].shape[0]),torch.tensor([]),'native_aa',score_mode)

            # get the bound amino acid probabilities
            bound_aa_probs,new_seq,_ = interfaceScorer.get_aa_prob(self.complex_netout[idx]['etab'],self.complex_netout[idx]['E_idx'],self.complex_netout[idx]['seq'],self.protein_res[idx],self.peptide_res[idx],seq_mode,score_mode)
            self.bound_aa_probs.append(bound_aa_probs)
            self.new_seq.append(new_seq)

            # print('target',interfaceScorer.tensorToSeq(self.target_netout['seq']))
            # print('nat_complex',interfaceScorer.tensorToSeq(self.complex_netout[idx]['seq']))
            # print('new_complex',interfaceScorer.tensorToSeq(new_seq))
            # print('binder_seq',interfaceScorer.tensorToSeq(self.new_seq[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]]))

            # get the unbound binder amino acid probabilities (we do this last because we need the consensus sequence from the previous step)
            unbound_binder_aa_probs,_,_ = interfaceScorer.get_aa_prob(self.binder_netout[idx]['etab'],self.binder_netout[idx]['E_idx'],self.new_seq[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]],torch.arange(0,self.binder_netout[idx]['etab'].shape[0]),torch.tensor([]),'native_aa',score_mode)
            self.unbound_binder_aa_probs.append(unbound_binder_aa_probs)

            # get the sequence probability for the target residues
            if self.unbound_target_seq_probs == None:
                self.unbound_target_seq_probs = interfaceScorer.get_seq_prob(self.unbound_target_aa_probs,self.target_netout['seq'])
            unbound_binder_seq_probs = interfaceScorer.get_seq_prob(self.unbound_binder_aa_probs[idx],self.new_seq[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]])
            self.unbound_binder_seq_probs.append(unbound_binder_seq_probs)

            bound_seq_probs = interfaceScorer.get_seq_prob(self.bound_aa_probs[idx],self.new_seq[idx])
            unbound_seq_probs = torch.cat((self.unbound_target_seq_probs,self.unbound_binder_seq_probs[idx]),dim=0)
            self.bound_seq_probs.append(bound_seq_probs)
            self.unbound_seq_probs.append(unbound_seq_probs)

            binder_scores = -torch.log(bound_seq_probs/unbound_seq_probs)
            self.binder_scores.append(binder_scores)

            # get the structure score from the bound peptide + protein
            if self.unbound_structure_score_targettotal == None:
                self.unbound_structure_score_targettotal = self.target_netout['sscore'].sum()
            self.unbound_structure_scores.append(torch.cat((self.target_netout['sscore'],self.binder_netout[idx]['sscore']),dim=0))
            self.bound_structure_scores.append(self.complex_netout[idx]['sscore'])

        return self.binder_scores, self.bound_structure_scores
    
    @staticmethod
    def get_aa_prob(etab, E_idx, seq, fixed_res, variable_res, seq_mode, score_mode):
        """ Get the probability of each protein amino acid

        NOTE: assumes that peptide residues are in a contiguous sequence

        Args
        ----
        etab : torch.tensor
            The self and pair energies for the given structure
            L x k x 20 x 20
        E_idx : torch.tensor 
            The k-NN indices which specify the interacting partner residue for each pair energy
            L x k
        seq : torch.tensor
            The int-encoded sequence of the structure
            L
        fixed_res : torch.tensor
            The indices of the residues in the structure with fixed identities (e.g., target protein residues)
            L - p
        variable_res : torch.tensor
            The indices of the residues in the structure with variable identities (e.g., the peptide residues)
            p
        mode : str
            Controls how the energies from the peptide residues are used to calculate the probabilities of target amino acids
            omit_peptide: pair energies from the peptide residues are not used in the calculation
            omit_target: pair energies from target are not used
            native_aa: the peptide residues are assumed to have their native identity (according to whatever sequence is provided as `seq`)
            -removed mean_energy: each pair energy is averaged over all possible residue identities at a variable residue
            consensus_aa: the lowest energy residue identity is selected for each variable residue (ignoring pair energies with other peptide residues) and this sequence assumed in the following calculations
            -removed marginal_prob_e: The pair energies are converted to a joint probability distribution and marginalized over the peptide residue identity. In practice this has a similar output to `mean_energy`

        Returns
        -------
        # nat_aa_probs: torch.tensor
        #     Native amino acid probabilities
        aa_probs: torch.tensor L x 20
            All amino acid probabilities
        new_seq: torch.tensor L
            The sequence that was used to calculate the probability (can differ from the input depending on mode)
        # native_seq: string
        #     The native peptide sequence
        # consensus_seq: string
        #     The peptide sequence after consensus design
        """

        L, k, aa1, aa2 = etab.shape
        p = variable_res.shape[0]

        new_seq = seq.detach().clone()

        # get the self energies for each position
        selfE = etab[:,0:1].squeeze(1) # L x 20 x 20
        selfE = torch.diagonal(selfE,offset=0,dim1=-2,dim2=-1) # L x 20

        # # select the self energies from protein residues
        # sel_res_idx = fixed_res.view(*fixed_res.shape,1).expand(-1,20) # L - p x 20
        # selfE_filt = torch.gather(selfE,dim=0,index=sel_res_idx) # L - p x 20

        # get all pair energies for each position
        pairE = etab[:,1:] # L x k - 1 x 20 x 20
        E_idx_pair = E_idx[:,1:k] # L x k - 1

        # replace the peptide sequence, if necessary
        if seq_mode == 'native_aa':
            # nothing to do here
            pass
        # elif mode == 'native_aa':
        #     # assume the peptide residues have the native sequence
        #     pairE_sum = pairE_u.sum(axis=-2) # L x 20
        #     sel_res_idx = fixed_res.view(*fixed_res.shape,1).expand(-1,20) # L - p x 20
        #     pairE_sum_filt = torch.gather(pairE_sum,dim=0,index=sel_res_idx) # L - p x 20
        elif seq_mode == 'consensus_aa':
            pepE_minaa = interfaceScorer.get_consensus_aa(etab,E_idx,new_seq,fixed_res,variable_res)
            
            # replace the sequence of the peptide
            p_min_idx,p_max_idx = (variable_res.min().item(),variable_res.max().item()) if variable_res.numel() != 0 else (0,-1)
            new_seq[p_min_idx:p_max_idx+1] = pepE_minaa # L
            
        # # "soft" scoring functions are a different kind of beast as they do not assume a specific sequence, they should probably be refactored into a distinct function
        # elif mode == 'mean_energy':
        #     # set energies from peptide to the average value
            
        #     # set non-selected residues to 0 with mask
        #     mask_nonsel = sum(E_idx_pair==res for res in fixed_res).bool() # L x 29
        #     pairE_u_masked = pairE_u*(mask_nonsel.view(*mask_nonsel.shape,1).expand(-1,-1,20))
        #     pairE_u_masked_sum = pairE_u_masked.sum(axis=-2) # L x 20
            
        #     # set selected residues to 0 with mask
        #     pairE_mean = pairE.mean(axis=-1) # L x k-1 x 20
        #     mask_sel = ~mask_nonsel
        #     pairE_mean_invmasked = pairE_mean*(mask_sel.view(*mask_sel.shape,1).expand(-1,-1,20)) # L x k-1 X 20
        #     pairE_mean_invmasked_sum = pairE_mean_invmasked.sum(axis=-2) # L x 20
            
        #     # sum the pair energies given rj_u for protein and averaged over rj for peptide
        #     sel_res_idx = fixed_res.view(*fixed_res.shape,1).expand(-1,20) # L - p x 20
        #     pairE_sum_filt = torch.gather(pairE_u_masked_sum,dim=0,index=sel_res_idx) # L - p x 20
        #     pairE_sum_filt +=  torch.gather(pairE_mean_invmasked_sum,dim=0,index=sel_res_idx) # L - p x 20
        # elif mode == 'marginal_prob_e':
        #     # set pair energies from non-selected residues to the energy computed from the marginal probability
            
        #     # set non-selected residues to 0 with mask
        #     mask_nonsel = sum(E_idx_pair==res for res in fixed_res).bool() # L x 29
        #     pairE_u_masked = pairE_u*(mask_nonsel.view(*mask_nonsel.shape,1).expand(-1,-1,20))
        #     pairE_u_masked_sum = pairE_u_masked.sum(axis=-2) # L x 20
            
        #     # compute the marginal probabilities, then convert back to energy
        #     pairE_marg = -torch.log((-pairE).view(L,k-1,aa1*aa2).softmax(axis=-1).view(L,k-1,aa1,aa2).sum(axis=-1)) # L x k-1 x 20
            
        #     # shift the energies to match the previous min_energy
        #     min_pairE,_ = pairE.view(L,k-1,aa1*aa2).min(axis=-1) # L x k-1 x 1
        #     min_pairE_marg,_ = pairE_marg.min(axis=-1) # L x k-1 x 1
        #     diff = min_pairE_marg - min_pairE
        #     pairE_marg = pairE_marg - diff.view(L,k-1,1)
            
        #     # set selected residues to 0 with mask
        #     mask_sel = ~mask_nonsel
        #     pairE_marg_invmasked = pairE_marg*(mask_sel.view(*mask_sel.shape,1).expand(-1,-1,20)) # L x k-1 X 20
        #     pairE_marg_invmasked_sum = pairE_marg_invmasked.sum(axis=-2) # L x 20
            
        #     # sum the pair energies given rj_u for protein and averaged over rj for peptide
        #     sel_res_idx = fixed_res.view(*fixed_res.shape,1).expand(-1,20) # L - p x 20
        #     pairE_sum_filt = torch.gather(pairE_u_masked_sum,dim=0,index=sel_res_idx) # L - p x 20
        #     pairE_sum_filt += torch.gather(pairE_marg_invmasked_sum,dim=0,index=sel_res_idx) # L - p x 20
        else:
            raise ValueError('seq mode: ',seq_mode,' not implemented')

        # get E_p(R_i=m|R_j=u), the energy of R_i=m given R_j=u
        seq_expand = new_seq.view(*new_seq.shape,1,1,1).expand(-1,k-1,20,-1) # L x 29 x 20 x 1
        # if (ignore_target_pair_e):
        #     print('Ignoring pair energies in amino acid probability calculation')
        #     E_idx_pair_jn = torch.zeros(*E_idx_pair.shape,20,1,dtype=E_idx_pair.dtype)
        # else:
        E_idx_pair_jn = E_idx_pair.view(*E_idx_pair.shape,1,1).expand(-1,-1,20,1) # L x 29 x 20 x 1
        rj_seq = torch.gather(seq_expand,dim=-4,index=E_idx_pair_jn) # L x 29 x 20 x 1
        pairE_u = torch.gather(pairE,dim=-1,index=rj_seq).squeeze(-1) # L x 29 x 20

        if score_mode == 'omit_peptide':
            # ignore pair energies involving a peptide residue
            mask_pep = sum(E_idx_pair==res for res in fixed_res).bool() # L x 29
            pairE_u_masked = pairE_u*(mask_pep.view(*mask_pep.shape,1).expand(-1,-1,20))
            pairE_u_sum = pairE_u_masked.sum(axis=-2) # L x 20
        elif score_mode == 'omit_target':
            # ignore pair energies involving a target residue
            mask_target = ~sum(E_idx_pair==res for res in fixed_res).bool() # L x 29
            pairE_u_masked = pairE_u*(mask_target.view(*mask_target.shape,1).expand(-1,-1,20))
            pairE_u_sum = pairE_u_masked.sum(axis=-2) # L x 20
        elif score_mode == 'interface_only':
            # only consider pair energies between residues across the interface
            E_idx_pair_target = E_idx_pair[0:L-p]
            E_idx_pair_peptide = E_idx_pair[L-p:]
            
            # mask out pair energies between a target residue and another target residue
            mask_target = ~sum(E_idx_pair_target==res for res in fixed_res).bool() # L - p x 29

            # mask out pair energies between a peptide residue and another peptide residue
            mask_peptide = ~sum(E_idx_pair_peptide==res for res in variable_res).bool() if variable_res.numel() != 0 else torch.tensor([]) # p x 29

            # concatenate masks and apply
            samechain_mask = torch.concat((mask_target,mask_peptide),dim=0)
            pairE_u_masked = pairE_u*(samechain_mask.view(*samechain_mask.shape,1).expand(-1,-1,20))
            pairE_u_sum = pairE_u_masked.sum(axis=-2) # L x 20
        elif score_mode == "default":
            pairE_u_sum = pairE_u.sum(axis=-2) # L x 20
        else:
            raise ValueError('score mode: ',score_mode,' not implemented')
            
        # sum the energies at each position
        positionE = selfE + pairE_u_sum  # L x 20

        # sum the energies over all positions
        totalE = torch.sum(torch.gather(positionE,axis=-1,index=seq.view(-1,1)))
        print('Energy of native sequence is: ',totalE.item())
        if seq_mode == 'consensus_aa':
            totalE = torch.sum(torch.gather(positionE,axis=-1,index=new_seq.view(-1,1)))
            print('Energy of new sequence is: ',totalE.item())

        # convert to probabilities
        aa_probs = torch.softmax(-positionE,axis=-1)

        return aa_probs,new_seq,pairE_u

    @staticmethod
    def get_consensus_aa(etab, E_idx, seq, protein_res, peptide_res):
        # for each peptide residue, find the lowest energy amino acid ignoring pair energies involving other peptide residues
        # this should serve as a fast-to-compute approximation of the lowest energy peptide sequence
        L, k, aa1, aa2 = etab.shape
        # p_min_idx,p_max_idx = (peptide_res.min().item(),peptide_res.max().item()) if peptide_res.numel() != 0 else (0,-1)
        
        selfE = etab[:,0:1].squeeze(1) # L x 20 x 20
        selfE = torch.diagonal(selfE,offset=0,dim1=-2,dim2=-1) # L x 20
        pairE = etab[:,1:] # L x 29 x 20 x 20
        E_idx_pair = E_idx[:,1:] # L x 29

        # get E_p(R_i=m|R_j=u), the energy of R_i=m given R_j=u
        seq_expand = seq.view(*seq.shape,1,1,1).expand(-1,k-1,aa1,-1) # L x 29 x 20 x 1
        E_idx_pair_jn = E_idx_pair.view(*E_idx_pair.shape,1,1).expand(-1,-1,aa1,1) # L x 29 x 20 x 1
        rj_seq = torch.gather(seq_expand,dim=-4,index=E_idx_pair_jn) # L x 29 x 20 x 1
        pairE_u = torch.gather(pairE,dim=-1,index=rj_seq).squeeze(-1) # L x 29 x 20

        # set E_p(r_i=m,r_u=u) to 0 where r_u in the peptide
        mask_nonsel = sum(E_idx_pair==res for res in protein_res).bool() # L x 29
        pairE_u_masked = pairE_u*(mask_nonsel.view(*mask_nonsel.shape,1).expand(-1,-1,aa1))
        pairE_u_masked_sum = pairE_u_masked.sum(axis=-2) # L x 20

        # get the minimum energy amino acid for each residue in the peptide
        # self energies
        peptide_res_idx = peptide_res.view(*peptide_res.shape,1).expand(-1,aa1) # p x 20
        pep_selfE = torch.gather(selfE,dim=0,index=peptide_res_idx) # p x 20

        # pair energies
        pep_pairE_u_masked_sum = torch.gather(pairE_u_masked_sum,dim=0,index=peptide_res_idx) # p x 20

        # sum and take the min energy aa at each position
        pepE = pep_selfE + pep_pairE_u_masked_sum # p x 20
        pepE_minaa = pepE.argmin(axis=-1) # p
        
        # nat_seq = interfaceScorer.tensorToSeq(seq[p+1])
        # consensus_seq = interfaceScorer.tensorToSeq(pepE_minaa)
        # seq_id = interfaceScorer.seqID(nat_seq,consensus_seq)
        # print('native pep seq:',nat_seq)
        # print('minE pep seq: ',consensus_seq)
        # print('seq ID: ',seq_id)
        return pepE_minaa

    @staticmethod
    def get_seq_prob(aa_probs, seq, sel_res_idx = None):
        # gather native aa from the protein and get native probability 
        seq_expand = seq.view(*seq.shape,1) # L x 1
        if sel_res_idx is not None:
            sel_res_idx_expand = sel_res_idx.view(-1,1).expand(-1,20)
            aa_probs_filt = torch.gather(aa_probs,dim=0,index=sel_res_idx_expand) # len(sel_res_idx) x 20
            seq_filt = torch.gather(seq_expand,dim=0,index=sel_res_idx.view(-1,1)) # len(sel_res_idx) x 1
            seq_prob = torch.gather(aa_probs_filt,dim=-1,index=seq_filt) # len(sel_res_idx) x 1
        else:
            seq_prob = torch.gather(aa_probs,dim=-1,index=seq_expand) # L x 1
        return seq_prob.squeeze()

    ### File output methods
    def write_residue_scores(self, score_mode: str):
        if (self.score_file == None):
            self.score_file = open('residueBinderScores.csv','w')
            self.score_file.write('name,score_mode,res_idx,chain_id,resnum,seqstruct_score,unbound_nat_aa_prob,bound_nat_aa_prob,unbound_struct_score,bound_struct_score\n')
        for idx in range(len(self.complex_netout)):
            assert len(self.binder_scores[idx]) == len(self.unbound_structure_scores[idx]) == len(self.prot_res_info[idx]) + len(self.pep_res_info[idx]) == len(self.bound_seq_probs[idx]) == len(self.unbound_seq_probs[idx]), print(len(self.binder_scores[idx]),len(self.bound_structure_scores[idx]),len(self.prot_res_info[idx]),len(self.pep_res_info[idx]))
            for i,((chain_id,res_num),score,unbound_prob,bound_prob,unbound_struct_score,bound_struct_score) in enumerate(zip(self.prot_res_info[idx]+self.pep_res_info[idx],self.binder_scores[idx],torch.cat((self.unbound_target_seq_probs,self.unbound_binder_seq_probs[idx]),dim=0),self.bound_seq_probs[idx],self.unbound_structure_scores[idx],self.bound_structure_scores[idx])):
                line = f"{self.complex_netout[idx]['id']},{score_mode},{i},{chain_id},{res_num},{score.item()},{unbound_prob.item()},{bound_prob.item()},{unbound_struct_score.item()},{bound_struct_score.item()}"
                self.score_file.write(line+'\n')

    def write_structure_scores(self, score_mode: str):
        if (self.structure_score_file == None):
            self.structure_score_file = open('structureBinderScores.csv','w')
            self.structure_score_file.write('name,score_mode,N_target_res,N_peptide_res,seqstruct_score,struct_score\n')
        for idx in range(len(self.complex_netout)):
            N_target_res = len(self.protein_res[idx])
            N_peptide_res = len(self.peptide_res[idx])
            seqstruct_score = torch.sum(self.binder_scores[idx]).item()
            struct_score = (torch.sum(self.bound_structure_scores[idx]) - self.unbound_structure_score_targettotal).item() # subtract unbound target for comparison
            line = f"{self.complex_netout[idx]['id']},{score_mode},{N_target_res},{N_peptide_res},{seqstruct_score},{struct_score}"
            self.structure_score_file.write(line+'\n')

    def write_pep_seqs(self):
        if self.seq_file == None:
            self.seq_file = open('pepSeqs.csv','w')
            self.seq_file.write('name,nat_seq,consensus_seq,seq_id\n')
        for idx in range(len(self.complex_netout)):
            native_seq = self.tensorToSeq(self.complex_netout[idx]['seq'][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]])
            consensus_seq = self.tensorToSeq(self.new_seq[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]])
            assert len(native_seq) == len(consensus_seq)
            self.seq_file.write(f"{self.complex_netout[idx]['id']},{native_seq},{consensus_seq},{self.seqID(native_seq,consensus_seq)}"+'\n')

    def write_pep_seq_probs(self):
        if self.pep_seq_prob_file == None:
            self.pep_seq_prob_file = open('pepSeqProb.csv','w')
            self.pep_seq_prob_file.write('name,pep_res_idx,res_num,chain_id,aa_prob,aa\n')
        for idx in range(len(self.complex_netout)):
            pep_seq = self.tensorToSeq(self.new_seq[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]])
            assert len(self.pep_res_info[idx]) == self.bound_seq_probs[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]].shape[0] == len(pep_seq)
            for i,((chain_id,res_num),aa_prob,seq_aa) in enumerate(zip(self.pep_res_info[idx],self.bound_seq_probs[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]],pep_seq)):
                line = f"{self.complex_netout[idx]['id']},{i},{chain_id},{res_num},{aa_prob.item()},{seq_aa}"
                self.pep_seq_prob_file.write(line+'\n')

    @staticmethod
    def get_energy(energy_per_res, seq):
        # energy_per_res: L x 20
        # seq: L
        return torch.sum(torch.gather(energy_per_res,dim=-1,index=seq)) # L x 1 -> 1

    @staticmethod
    def seqID(seq1,seq2):
        return 100*np.array([x==y for x,y in zip(seq1,seq2)]).mean()

    @staticmethod
    def tensorToSeq(tensor):
        return ''.join([int_to_AA[x.item()] for x in tensor])

    def close_files(self):
        # not a great solution, I know
        if self.score_file is not None:
            self.score_file.close()
        if self.structure_score_file is not None:
            self.structure_score_file.close()
        if self.seq_file is not None:
            self.seq_file.close()
        if self.pep_seq_prob_file is not None:
            self.pep_seq_prob_file.close()

    def pickle_network_output_target(self, path_prefix:str):
        target_netout_pickle = {
            'etab':self.target_netout['etab'].numpy(),
            'E_idx':self.target_netout['E_idx'].numpy(),
            'sscore':self.target_netout['sscore'].numpy(),
            'seq':self.target_netout['seq'].numpy(),
            'protein_res':self.protein_res[0].numpy(),
        }
        with open(path_prefix+self.target_netout['id']+'_target.pickle','wb') as file:
            pickle.dump(target_netout_pickle, file, pickle.HIGHEST_PROTOCOL)

    def pickle_network_output_binder(self, path_prefix:str):
        for idx in range(len(self.complex_netout)):
            netout_pickle = {
                'etab':self.binder_netout[idx]['etab'].numpy(),
                'E_idx':self.binder_netout[idx]['E_idx'].numpy(),
                'sscore':self.binder_netout[idx]['sscore'].numpy(),
                'seq':self.binder_netout[idx]['seq'].numpy(),
                'peptide_res':self.peptide_res[idx].numpy(),
                'final_seq':self.new_seq[idx][self.pep_res_range[idx][0]:self.pep_res_range[idx][1]].numpy()
            }
            with open(path_prefix+self.binder_netout[idx]['id']+'_binder.pickle','wb') as file:
                pickle.dump(netout_pickle, file, pickle.HIGHEST_PROTOCOL)

    def pickle_network_output_complex(self, path_prefix:str):
        for idx in range(len(self.complex_netout)):
            complex_netout_pickle = {
                'etab':self.complex_netout[idx]['etab'].numpy(),
                'E_idx':self.complex_netout[idx]['E_idx'].numpy(),
                'sscore':self.complex_netout[idx]['sscore'].numpy(),
                'seq':self.complex_netout[idx]['seq'].numpy(),
                'protein_res':self.protein_res[idx].numpy(),
                'peptide_res':self.peptide_res[idx].numpy(),
                'final_seq':self.new_seq[idx].numpy()
            }
            with open(path_prefix+self.complex_netout[idx]['id']+'.pickle','wb') as file:
                pickle.dump(complex_netout_pickle, file, pickle.HIGHEST_PROTOCOL)    
