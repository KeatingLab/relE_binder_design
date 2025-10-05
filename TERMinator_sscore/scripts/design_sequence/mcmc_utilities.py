import time
import regex as re
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
import sys
import os
import numexpr as ne

from terminator.utils.common import AA_to_int,int_to_AA

'''
Note regarding sequence complexity: this metric is maximized when the amino acid usage for all types is even.
This means that a sequence that has low local sequence complexity ("AAAAEEEECCCC") could still have a high 
complexity. As a result, I found that when optimizing a subset of positions, it's important to calculate complexity
over those positions only. Otherwise, complexity can be increased by simply making mutations that are rare inside of 
the fixed region (e.g. HHHHH)
'''

def stirlingApprox(N):
    # return N*np.log(N) - N
    # more accurate version
    return 0.5*np.log(2*np.pi*N)+N*np.log(N/np.e)

def sequenceComplexity(sequence: np.array):
    # Use stirlings approximation
    assert sequence.ndim == 1
    L = sequence.size
    values, counts = np.unique(sequence,return_counts=True)
    log_N_possible_sequences = stirlingApprox(L)
    log_N_nonunique_rearrangements = 0
    for count in counts:
        log_N_nonunique_rearrangements += stirlingApprox(count)
    return (log_N_possible_sequences - log_N_nonunique_rearrangements)/L

def get_merge_dups_mask(E_idx):
    N = E_idx.shape[0]
    tens_place = np.expand_dims(np.arange(N), (-1))
    # tens_place = tens_place.unsqueeze(0).unsqueeze(-1)
    min_val = np.minimum(E_idx, tens_place)
    max_val = np.maximum(E_idx, tens_place)
    edge_indices = min_val*N + max_val
    edge_indices = edge_indices.flatten()
    uniq, unique_inv = np.unique(edge_indices, return_inverse=True)
    num_edges = len(uniq)
    return unique_inv, num_edges

class mcmcEneryTable:
    """Efficient computation of energy given a sequence
    """

    def __init__(self,etab,E_idx):
        """ Reformat the potts model
        
        Parameters
        ----------
        etab : np.ndarray
            The self and pair energies
            size: L x k x 20 x 20
        E_idx : np.ndarray
            The k-NN graph providing the indices of the residues 
            size: L x k
        """
        # Store dimensions and verify
        self.L, self.k, aa1, aa2, = etab.shape
        assert aa1 == aa2 == 20, "Energy table should be formatted to handle 20 amino acids"

        mapping, n_edges = get_merge_dups_mask(E_idx)

        etab_r = etab.reshape(self.L * self.k, aa1, aa2)
        etab_dedup = np.zeros((n_edges, aa1, aa2))
        seen_mappings = set()
        for i, i_map in enumerate(mapping):
            if i_map in seen_mappings:
                continue
            etab_dedup[i_map] = etab_r[i]
            seen_mappings.add(i_map)
        self.etab_dedup = etab_dedup
        self.mapping = mapping

        self.E_idx = E_idx

        # create indexing arrays
        self.res_i = np.expand_dims(np.arange(0,self.L),axis=1)
        self.res_j = np.expand_dims(np.arange(0,self.k-1),axis=0)
        self.res_ij_mapped = np.arange(np.max(mapping)+1)


        _, indices = np.unique(mapping, return_index=True) 
        pair_id_to_i = dict(zip(mapping[indices], np.floor(indices / self.k)))  
        self.pair_id_to_i = [int(id) for id in pair_id_to_i.values()]
        pair_id_to_j = dict(zip(mapping[indices], indices))  
        self.pair_id_to_j = E_idx.flatten()[list(pair_id_to_j.values())]
        self.pair_E_idx = [[]]*self.L
        self.pair_E_idx_t = [[]]*self.L
        self.pair_ij_mapped = [[]]*self.L
        self.pair_ij_mapped_t = [[]]*self.L
        self.self_ij_mapped = [[]]*self.L
        for i_res in range(self.L):
            for iid, (id_to_i, id_to_j) in enumerate(zip(self.pair_id_to_i, self.pair_id_to_j)):
                if i_res == id_to_i and i_res == id_to_j:
                    self.self_ij_mapped[i_res] = self.self_ij_mapped[i_res] + [iid]
                elif i_res == id_to_i:
                    self.pair_E_idx[i_res] = self.pair_E_idx[i_res] + [id_to_j]
                    self.pair_ij_mapped[i_res] = self.pair_ij_mapped[i_res] + [iid]
                elif i_res == id_to_j:
                    self.pair_E_idx_t[i_res] = self.pair_E_idx_t[i_res] + [id_to_i]
                    self.pair_ij_mapped_t[i_res] = self.pair_ij_mapped_t[i_res] + [iid]

        # divide table into self and pair E
        self.self_E = np.diagonal(etab[:,0],axis1=-2,axis2=-1)
        self.pair_E = etab[:,1:,]
        # self.pair_E_idx = self.E_idx[:,1:]

    def getEnergy(self,sequence):
        """ Get energy given sequence
        
        Parameters
        ----------
        sequence : np.ndarray
            size: L,
        """
        assert self.self_E.shape[0] == sequence.shape[0]

        # # self energy
        # energy = self.self_E[np.arange(0,self.L),sequence].sum()

        # pair energy
        res_seq_i = sequence[self.pair_id_to_i]
        res_seq_j = sequence[self.pair_id_to_j]

        energy = self.etab_dedup[self.res_ij_mapped,res_seq_i,res_seq_j].sum()
        # energy += self.pair_E[self.res_i,self.res_j,res_seq_i,res_seq_j].sum()
        return energy
    
    def getResidueEnergy(self,sequence, position, mutation):
        """ Score energy of specific mutation

        Parameters
        ----------
        sequence : np.ndarray
            size: L,
        position : int
            index of position of mutation
        mutation : int
            new residue at position
        """
        self_id = self.self_ij_mapped[position]
        E = self.etab_dedup[self_id, mutation, mutation].sum() - self.etab_dedup[self_id, sequence[position], sequence[position]].sum()

        res_seq_j = sequence[self.pair_E_idx[position]]
        pair_ids = self.pair_ij_mapped[position]
        E += self.etab_dedup[pair_ids, mutation, res_seq_j].sum() - self.etab_dedup[pair_ids, sequence[position], res_seq_j].sum()

        res_seq_i = sequence[self.pair_E_idx_t[position]]
        pair_ids_t = self.pair_ij_mapped_t[position]
        E += self.etab_dedup[pair_ids_t, res_seq_i, mutation].sum() - self.etab_dedup[pair_ids_t, res_seq_i, sequence[position]].sum()
        
        # res_seq_j = sequence[self.pair_E_idx[position]]
        # print('res_seq_j: ', res_seq_j)
        # print(self.pair_E_idx)
        # print(self.pair_E_idx[position])
        # print(res_seq_j)
        # print(position)
        # raise ValueError
        # E = self.self_E[position][mutation] + self.pair_E[position,self.res_j[0],mutation,res_seq_j].sum() - self.self_E[position][sequence[position]] - self.pair_E[position,self.res_j[0],sequence[position],res_seq_j].sum()
        return E
    
    def check_update(self, cur_seq, update_seq, pos, mut):
        res_seq_j = cur_seq[self.pair_E_idx[pos]]
        print("\t", res_seq_j)

        res_seq_i = cur_seq[self.pair_id_to_i]
        res_seq_j = cur_seq[self.pair_id_to_j]

        # print("\t", res_seq_i, "  ", res_seq_j)

        res_seq_iu = update_seq[self.pair_id_to_i]
        res_seq_ju = update_seq[self.pair_id_to_j]

        # print("\t", res_seq_iu, "  ", res_seq_ju)

        E_cur = self.etab_dedup[self.res_ij_mapped,res_seq_i,res_seq_j]
        E_u = self.etab_dedup[self.res_ij_mapped,res_seq_iu,res_seq_ju]
        dc = 0
        for i, (ec, eu) in enumerate(zip(E_cur, E_u)):
            if ec != eu:
                print("\t", i, self.res_ij_mapped[i], res_seq_i[i], res_seq_j[i], res_seq_iu[i], res_seq_ju[i], self.pair_id_to_i[i], self.pair_id_to_j[i], ec, eu)
                dc += 1
        print("\t", dc)

        print("\t", cur_seq)
        print("\t", update_seq)

        res_seq_j = cur_seq[self.pair_E_idx[pos]]
        pair_ids = self.pair_ij_mapped[pos]        # E_p = self.pair_E[pos,self.res_j[0],mut,res_seq_j]
        print('t', len(res_seq_j), res_seq_j, pair_ids, self.pair_E_idx[pos], mut)
        print('\t', self.etab_dedup[pair_ids, mut, res_seq_j])
        # print('\t', E_p)

        E = self.self_E[pos][cur_seq[pos]] 
        # E_p = self.pair_E[pos,self.res_j[0],cur_seq[pos],res_seq_j]
        print('\t', self.etab_dedup[pair_ids, cur_seq[pos], res_seq_j])
        # print('\t', E_p)


@dataclass
class residueGroup:
    """Defines a set of residues
    
    Note the following logic:
    - If no chain is set, then no residues are selected
    - If a chain is set, but no residues are selected, then ALL residues of that chain are selected
    - If a chain and residues are selected, then residues matching those numbers within the chain are selected
    """
    chain_id: str = str()
    res_num: set = field(default_factory=set)

    def is_res_selected(self,chain_id,res_num):
        if chain_id == '' or chain_id != self.chain_id:
            return False
        else:
            if len(self.res_num) == 0:
                return True
            elif res_num in self.res_num:
                return True
            else:
                return False

@dataclass
class mcmcParams:
    '''
    n_it : int
        The number of iterations per cycle
    n_cycles : int
        The number of cycles to perform (each cycle is reinitialized)
    Ti : int
        The initial temperature
    Tf : int
        The final temperature
    seed : int
        default = None
        If set, will ensure reproducibility
    complexity_weight : float
        If > 0, the complexity of the sequence will be computed and added to the energy
    '''
    n_it: int = 100000
    n_cycles: int = 10
    Ti: float = 1
    Tf: float = 0.01
    seed: float = -1
    complexity_weight: float = -1.0
    convergence: int = 5

    def printParams(self):
        print(f"n_it = {self.n_it}")
        print(f"n_cycles = {self.n_cycles}")
        print(f"Ti = {self.Ti}")
        print(f"Tf = {self.Tf}")
        print(f"seed = {self.seed}")
        print(f"complexity_weight = {self.complexity_weight}")
        print(f"convergence = {self.convergence}")

class simulatedAnnealingMCMC:
    """A class for optimizing an energy table"""

    def __init__(self):
        self.etab: mcmcEneryTable = None
        self.residue_info: list = None #[(chain_id,res_num)] size: L
        self.constraint_sequence: np.array = None # size: L
        self.variable_pos: np.array = None # size: V < L
        self.aa_size = 20

        # store info during optimization
        self.cycle_store = []
        self.init_seq_store = []
        self.best_seq_store = []
        self.init_seq_e_store = []
        self.best_seq_e_store = []
        self.best_seq_complexity_store = [] # note: this is only calculated over the variable region
        self.top_100_seq_store = []
        self.top_100_E_store = []
        self.top_100_marker = -1

    def loadEnergyTable(self, etab: np.array, E_idx: np.array, residue_info=None):
        if residue_info is not None:
            assert len(residue_info) == etab.shape[0], f"Number of residues does not match: {len(residue_info)} vs {etab.shape[0]}"
        else:
            residue_info = [('A',x+1) for x in range(etab.shape[0])]
        self.residue_info = residue_info
        self.etab = mcmcEneryTable(etab,E_idx)
        self.variable_pos = np.array([i for i in range(len(self.residue_info))])

    def setConstraints(self, constraint_sequence, selected_var_res: list):
        """ Provide a sequence as a constraint and selected residues that will be allowed to mutate
        
        Parameters
        ----------
        constraint_sequence : str
            size: L
            A 1-letter amino acid string matching the length of the energy table
        selected_res: list
            A list of residueGroup objects. The union of these defines the set of residues that will be allowed to vary during MCMC.
        """
        assert self.residue_info != None, "Must set residue info before applying constraints"
        assert len(constraint_sequence) == len(self.residue_info), f"Number of residues does not match: {len(constraint_sequence)} vs {self.residue_info}"
        self.constraint_sequence = self.seq2array(constraint_sequence)
        variable_pos = []
        for i,(chain_id,res_num) in enumerate(self.residue_info):
            for res_group in selected_var_res:
                if res_group.is_res_selected(chain_id,res_num):
                    variable_pos.append(i)
                    continue
        self.variable_pos = np.array(variable_pos)
        print(f"self.variable_pos: {self.variable_pos}")
        if len(self.variable_pos) == 0:
            raise ValueError("No variable positions were selected")
        
    def unsetConstraints(self):
        self.constraint_sequence = None
        self.variable_pos = np.array([i for i in range(len(self.residue_info))])

    def seq2array(self,seq:str):
        return np.array([AA_to_int[aa1] for aa1 in seq])

    def array2seq(self,seq_array:np.array):
        return ''.join([int_to_AA[aa_i] for aa_i in seq_array])

    def initSeq(self) -> np.array:
        if self.constraint_sequence is None:
            return np.random.randint(20,size=len(self.residue_info))
        else:
            # initialize whole sequence to constraint sequence
            sequence = self.constraint_sequence 
            # set variable positions to random residues
            sequence[self.variable_pos] = np.random.randint(20,size=len(self.variable_pos))
            return sequence
    
    def getTempSchedule(self,n_it,Ti,Tf,type='linear'):
        assert Tf > 0, "Final temp cannot be zero"
        # linear for now
        return np.linspace(Ti,Tf,n_it)

    # def proposeMutation(self,sequence) -> tuple:
    #     # select a variable position
    #     res_idx = self.variable_pos[np.random.randint(len(self.variable_pos))]
    #     # select a *new* residue
    #     # (this was the most efficient way I could think of to randomly select a new aa (excluding the current aa))
    #     aa_idx = (sequence[res_idx] + np.random.randint(self.aa_size-1)) % 20
    #     return (res_idx,aa_idx)

    def acceptByMetropolisCriteria(self,deltaE,T,rand):
        if deltaE < 0:
            return True
        else:
            # a = -(deltaE)/T
            # P_accept = ne.evaluate('exp(a)')
            P_accept = np.exp(-(deltaE)/T)
            if rand < P_accept:
                return True
            else:
                False

    def getTotalEnergy(self,sequence: np.array, complexity_weight):
        if complexity_weight < 0:
            return self.etab.getEnergy(sequence)
        else:
            return self.etab.getEnergy(sequence)-complexity_weight*sequenceComplexity(sequence[self.variable_pos])

    def getMutationEnergy(self,sequence: np.array, position, mutation, complexity_weight):
        if complexity_weight < 0:
            return self.etab.getResidueEnergy(sequence, position, mutation)
        else:
            return self.etab.getResidueEnergy(sequence, position, mutation)-complexity_weight*sequenceComplexity(sequence[self.variable_pos])

    def getPerResidueSequenceEnergy(self,sequence, complexity_weight):
        seq_E = np.zeros(sequence.shape)
        for i, s in enumerate(sequence):
            seq_E[i] = self.getMutationEnergy(sequence, i, s, complexity_weight)
        return seq_E

    def sample(self,params=None, early_stopping=-1, verbose=False):
        ''' Sample a sequence by MCMC
        Parameters
        ----------
        params : mcmcParams
            A simple class for storing the parameters for MCMC
        Note: there are five separate arrays stored in memory
        -init_seq
        -current_seq
        -propose_seq
        -best_seq
        -overall_best_seq
        Care is taken to ensure that these are copied as infrequently as possible (instead, they are modified in place) 
        '''
        assert self.etab != None, "Must set energy table before sampling"
        # print("Early stopping: ", early_stopping, flush=True)
        try:
            # reset
            self.cycle_store = []
            self.init_seq_store = []
            self.best_seq_store = []
            self.init_seq_e_store = []
            self.best_seq_e_store = []
            self.best_seq_complexity_store = []
            self.top_100_seq_store = []
            self.top_100_E_store = []
            self.top_100_marker = -1
            if params == None:
                if verbose:
                    print("Using default parameters")
                params = mcmcParams()
            if verbose:
                params.printParams()
            if params.seed is not None and params.seed > 0:
                np.random.seed(params.seed)
            start = time.time()
            T_schedule = self.getTempSchedule(params.n_it,params.Ti,params.Tf)
            allcycles_best_seq = None
            repeat_count = 0
            for cycle in range(params.n_cycles):
                rands = np.random.rand(params.n_it)
                seq_update_idx = np.random.randint(low=1, high=self.aa_size, size=params.n_it)
                seq_update_pos = np.random.randint(low=0, high=len(self.variable_pos), size=params.n_it)

                self.cycle_store.append(cycle)

                # initialize sequence
                init_seq = self.initSeq()
                self.init_seq_store.append(self.array2seq(init_seq))
                current_seq = init_seq.copy()
                best_seq = init_seq.copy()
                # current_E = self.getTotalEnergy(current_seq,params.complexity_weight)
                # self.seq_E = self.getPerResidueSequenceEnergy(current_seq,params.complexity_weight)
                best_E = self.getTotalEnergy(current_seq,params.complexity_weight)
                self.init_seq_e_store.append(best_E.item())
                diff_E_best = 0
                for it in range(params.n_it):
                    if it + 1 % 10000 == 0:
                        diff_E_best = self.getTotalEnergy(current_seq, params.complexity_weight) -  best_E

                    # propose random mutation
                    # 1) select a variable position
                    r_idx = self.variable_pos[seq_update_pos[it]]
                    # 2) select a *new* residue identity
                    aa_idx = (current_seq[r_idx] + seq_update_idx[it]) % 20
                    # propose_seq[r_idx] = aa_idx
                    delta_E = self.etab.getResidueEnergy(current_seq, r_idx, aa_idx)
                    accept = self.acceptByMetropolisCriteria(delta_E, T_schedule[it], rands[it])

                    if accept:
                        # update current sequence
                        current_seq[r_idx] = aa_idx
                        diff_E_best += delta_E
                        
                        if diff_E_best < 0:
                            best_seq = current_seq.copy()
                            best_E += diff_E_best
                            diff_E_best = 0


                best_E = self.getTotalEnergy(best_seq,params.complexity_weight)
                if allcycles_best_seq is None or best_E < self.getTotalEnergy(allcycles_best_seq,params.complexity_weight):
                    allcycles_best_seq = best_seq.copy()
                    # print("resetting best seq")
                    repeat_count = 0
                elif cycle > 0:
                    # if the best sequence from this cycle is the same as the last, then we count this as convergence
                    # print("increasing repeat count")
                    repeat_count += np.array_equal(best_seq,allcycles_best_seq)
                    # print(best_seq)
                    # print(overall_best_seq)
                    # print("new repeat: ", repeat_count)
                    # repeat_count += 1

                self.best_seq_store.append(self.array2seq(best_seq))
                self.best_seq_e_store.append(best_E.item())
                self.best_seq_complexity_store.append(sequenceComplexity(best_seq[self.variable_pos]))
                
                # Report results
                if verbose:
                    print(f"Cycle: {cycle}")
                    print(f"Initial sequence:      {self.array2seq(init_seq)}, E = {self.etab.getEnergy(init_seq)}, complexity = {sequenceComplexity(init_seq[self.variable_pos])}")
                    print(f"Current sequence:      {self.array2seq(current_seq)}, E = {self.etab.getEnergy(current_seq)}, complexity = {sequenceComplexity(current_seq[self.variable_pos])}")
                    print(f"Best sequence:         {self.array2seq(best_seq)}, E = {self.etab.getEnergy(best_seq)}, complexity = {sequenceComplexity(best_seq[self.variable_pos])}")
                    print("repeat count: ", repeat_count, " | early stopping: ", early_stopping)
                if early_stopping > 0 and repeat_count >= early_stopping - 1:
                    print(f"The lowest energy sequence has now been observed {early_stopping} times, indicating convergence. Terminate search early.")
                    break
            if verbose:
                print(f"Overall best sequence: {self.array2seq(allcycles_best_seq)}, E = {self.etab.getEnergy(allcycles_best_seq)}, complexity = {sequenceComplexity(allcycles_best_seq[self.variable_pos])}")
            stop = time.time()
            if verbose:
                print(f"Elapsed: {stop - start}")
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)
            print(e)

    def getDataFrame(self, native_seq, nsr):
        data_dict = {'cycle':self.cycle_store,
              'init_seq':self.init_seq_store,
              'best_seq':self.best_seq_store,
              'native_seq': native_seq,
              'best_NSR': nsr,
              'init_seq_e':self.init_seq_e_store,
              'best_seq_e':self.best_seq_e_store,
              'best_seq_complexity':self.best_seq_complexity_store}

        return pd.DataFrame(data_dict)