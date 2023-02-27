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

# pylint: disable=unspecified-encoding

# Functions for computing the score of a binder 
'''
TODO: allow the user to provide another structure that models the "unfolded" structure of the protein
'''
class singlechainScorer:
    def __init__(self, structure_data, model, dev='cuda:0'):
        self.model = model
        self.dev = dev

        # Set model to eval mode
        self.model.eval()
        torch.set_grad_enabled(False)
    
        self.score_file = None
        self.seq_file = None
        self.pep_seq_prob_file = None

        self.set_new_structure(structure_data)

    def set_new_structure(self,structure_data):
        # target is always a single structure
        self.structure_data = structure_data
        self.structure_netout = None
        self.structure_aa_probs = None

    def get_pdb_name(self,idx):
        return self.structure_netout[idx]['id']

    def generate_etabs(self):
        # run COORDinator to generate etabs
        if self.structure_netout == None:
            self.structure_netout = self.run_model(self.structure_data, self.dev)
        print(f"Ran model to generate network output for {len(self.structure_netout)} structures")

        # get the residue info
        (self.protein_res,
        self.prot_res_info) = self.get_res_idx([x['res_info'] for x in self.structure_netout])
        self.prot_res_range = [(x.min(),x.max()+1) for x in self.protein_res]
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
    def get_res_idx(res_info_list):
        protein_res_list = list()
        prot_res_info_list = list()
        for idx,res_info in enumerate(res_info_list):
            protein_res = list()
            prot_res_info = list()
            for i,(chain,resnum) in enumerate(res_info):
                protein_res.append(i)
                prot_res_info.append(res_info[i])
            protein_res_list.append(torch.tensor(protein_res))
            prot_res_info_list.append(prot_res_info)
        return protein_res_list, prot_res_info_list

    def get_scores(self):
        self.structure_aa_probs = []
        self.structure_seq_probs = []
        self.seqstruct_scores = []
        self.structure_scores = []

        if self.structure_netout == None:
            self.generate_etabs()

        for idx in range(len(self.structure_netout)):            
            # get the unbound target amino acid probabilities
            self.structure_aa_probs.append(singlechainScorer.get_aa_prob(self.structure_netout[idx]['etab'],self.structure_netout[idx]['E_idx'],self.structure_netout[idx]['seq']))

            # get the sequence probability for the target residues
            self.structure_seq_probs.append(singlechainScorer.get_seq_prob(self.structure_aa_probs[idx],self.structure_netout[idx]['seq']))

            seqstruct_scores = -torch.log(self.structure_seq_probs[idx]/0.05) # assume uniform background
            self.seqstruct_scores.append(seqstruct_scores)

            # get the structure score from the bound peptide
            self.structure_scores.append(self.structure_netout[idx]['sscore'])

        return self.seqstruct_scores, self.structure_scores
    
    @staticmethod
    def get_aa_prob(etab, E_idx, seq):
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

        # get E_p(R_i=m|R_j=u), the energy of R_i=m given R_j=u
        seq_expand = new_seq.view(*new_seq.shape,1,1,1).expand(-1,k-1,20,-1) # L x 29 x 20 x 1
        E_idx_pair_jn = E_idx_pair.view(*E_idx_pair.shape,1,1).expand(-1,-1,20,1) # L x 29 x 20 x 1
        rj_seq = torch.gather(seq_expand,dim=-4,index=E_idx_pair_jn) # L x 29 x 20 x 1
        pairE_u = torch.gather(pairE,dim=-1,index=rj_seq).squeeze(-1) # L x 29 x 20
        pairE_u_sum = pairE_u.sum(axis=-2) # L x 20
            
        # sum the energies at each position
        positionE = selfE + pairE_u_sum  # L x 20

        # sum the energies over all positions
        totalE = torch.sum(torch.gather(positionE,axis=-1,index=seq.view(-1,1)))
        print('Energy of native sequence is: ',totalE.item())

        # convert to probabilities
        aa_probs = torch.softmax(-positionE,axis=-1)

        return aa_probs

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
    def write_scores(self):
        if (self.score_file == None):
            self.score_file = open('singlechainScores.csv','w')
            self.score_file.write('name,res_idx,chain_id,resnum,seqstruct_score,folded_nat_aa_prob,structure_score\n')
        for idx in range(len(self.structure_netout)):
            assert len(self.seqstruct_scores[idx]) == len(self.structure_scores[idx]) == len(self.prot_res_info[idx]) == len(self.structure_seq_probs[idx]), print(len(self.seqstruct_scores[idx]),len(self.structure_scores[idx]),len(self.prot_res_info[idx]),len(self.pep_res_info[idx]))
            for i,((chain_id,res_num),seqstruct_score,structure_score,aa_fold_prob) in enumerate(zip(self.prot_res_info[idx],self.seqstruct_scores[idx],self.structure_scores[idx],self.structure_seq_probs[idx])):
                line = f"{self.structure_netout[idx]['id']},{i},{chain_id},{res_num},{seqstruct_score.item()},{aa_fold_prob.item()},{structure_score}"
                self.score_file.write(line+'\n')

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
        if self.seq_file is not None:
            self.seq_file.close()
        if self.pep_seq_prob_file is not None:
            self.pep_seq_prob_file.close()

    def pickle_network_output_structure(self, path_prefix:str):
        structure_netout_pickle = {
            'etab':self.target_netout['etab'].numpy(),
            'E_idx':self.target_netout['E_idx'].numpy(),
            'sscore':self.target_netout['sscore'].numpy(),
            'seq':self.target_netout['seq'].numpy(),
            'protein_res':self.protein_res[0].numpy(),
        }
        with open(path_prefix+self.target_netout['id']+'_target.pickle','wb') as file:
            pickle.dump(structure_netout_pickle, file, pickle.HIGHEST_PROTOCOL)
