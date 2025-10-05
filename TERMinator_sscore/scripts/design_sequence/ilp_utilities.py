import regex as re
import numpy as np
import pulp

import torch

from terminator.utils.common import int_to_AA,int_to_3lt_AA,AA_to_int,aa_three_to_one

'''
Code for loading a Potts model and optimizing over sequence space.

The Potts model may either be in a KNN-sparse tensor format (output by the network) or
as a long-form "energy table" where each line describes a self or pair energy between two residues.

The optimization process itself is handled by an external solver (e.g. CPLEX). The code here 
functions as a wrapper to set up the problem and call the solver via the API.
'''

class PottsModel:
    """A class for storing and accessing self and pair energies for a given protein structure"""

    def __init__(self):
        # list of tuples: the (chain id,resnum) of each residue (position in this list defines the index of each residue)
        self.residue_info = []

        # list of lists: (L,a) where L is the number of residues in the structure and a is the number of amino acid types
        self.self_energies = []

        # dictionary of lists: (P,a^2) where P is the number of unique pairs. 
        self.pair_energies = dict()

        self.aa_size = 20

        print('aa_size',self.aa_size)

    def loadFromKNNSparse(self, etab: np.ndarray, E_idx: np.ndarray, _residue_info=None):
        """ Converts the Potts Model output by the network into a format that is convenient for setting up the LP problem
        
        Parameters
        ----------
        etab : np.ndarray
            The self and pair energies
            size: L x k x 20 x 20
        E_idx : np.ndarray
            The k-NN graph providing the indices of the residues 
            size: L x k
        _residue_info : list
            A list of tuples providing the chain ID and number of each residue in the structure (e.g. [('A',0),('A',1),...,('A',100)]
            size: L
        """
        print('etab',type(etab),etab.shape)
        print('E_idx',type(E_idx),E_idx)
        assert etab.shape[0] == E_idx.shape[0]
        if _residue_info == None:
            # create dummy residue info
            self.residue_info = [('A',resnum) for resnum in range(etab.shape[0])]
        else:
            self.residue_info = _residue_info
        print('residue_info',type(self.residue_info),len(self.residue_info))
        assert len(self.residue_info) == etab.shape[0]

        # get self energies
        selfE = etab[:,0:1].squeeze(1) # L x 20 x 20
        selfE = np.diagonal(selfE,offset=0,axis1=-2,axis2=-1) # L x 20
        for res_idx in range(selfE.shape[0]):
            res_self_e = [selfE[res_idx][aa_idx] for aa_idx in range(self.aa_size)]
            self.self_energies.append(res_self_e)
        
        # get pair energies
        pairE = etab[:,1:] # L x k-1 x 20 x 20
        pair_E_idx = E_idx[:,1:] # L x k-1
        for res_i_idx in range(pairE.shape[0]):
            for j,res_j_idx in enumerate(pair_E_idx[res_i_idx]):
                if res_i_idx > res_j_idx:
                    continue
                self.pair_energies[(res_i_idx,res_j_idx)] = [pairE[res_i_idx][j][aa_i_idx][aa_j_idx] for aa_i_idx in range(self.aa_size) for aa_j_idx in range(self.aa_size)]

    def readEtabLine(self,line):
        self_e_pattern = "(\w),(\d+)\s(\D{3})\s([+-]?([0-9]*[.])?[0-9]+)"
        pair_e_pattern = "(\w),(\d+)\s(\w),(\d+)\s(\D{3})\s(\D{3})\s([+-]?([0-9]*[.])?[0-9]+)"
        # try to parse as a self energy line
        match = re.match(self_e_pattern,line)
        if match:
            return {'chainid':match.group(1),'resnum':int(match.group(2)),
                    'aa':aa_three_to_one(match.group(3)),'energy':float(match.group(4))}

        # try to parse as a pair energy line
        match = re.match(pair_e_pattern,line)
        if match:
            return {'chainid_i':match.group(1),'resnum_i':int(match.group(2)),
                    'chainid_j':match.group(3),'resnum_j':int(match.group(4)),
                    'aa_i':aa_three_to_one(match.group(5)),'aa_j':aa_three_to_one(match.group(6)),
                    'energy':float(match.group(7))}

        raise ValueError("Line formatting is not consistent with self or pair energy:",line)
    
    def readEtabLineSingleA(self,line):
            self_e_pattern = "(\w),(\d+)\s(\D)\s([+-]?([0-9]*[.])?[0-9]+)"
            pair_e_pattern = "(\w),(\d+)\s(\w),(\d+)\s(\D)\s(\D)\s([+-]?([0-9]*[.])?[0-9]+)"
            # try to parse as a self energy line
            match = re.match(self_e_pattern,line)
            if match:
                return {'chainid':match.group(1),'resnum':int(match.group(2)),
                        'aa':match.group(3),'energy':float(match.group(4))}

            # try to parse as a pair energy line
            match = re.match(pair_e_pattern,line)
            if match:
                return {'chainid_i':match.group(1),'resnum_i':int(match.group(2)),
                        'chainid_j':match.group(3),'resnum_j':int(match.group(4)),
                        'aa_i':match.group(5),'aa_j':match.group(6),
                        'energy':float(match.group(7))}

            raise ValueError("Line formatting is not consistent with self or pair energy:",line)

    def loadFromEtabFile(self,path):
        """ Converts the Potts Model output by the dTERMen into a format that is convenient for setting up the LP problem
        
        Parameters
        ----------
        path : str
            The path to a dTERMen energy table file
        """
        # temp storage for the parameters
        self_table = dict()
        pair_table = dict()
        with open(path,'r') as file:
            lines = [line.rstrip() for line in file]

        # get the residue info
        try:
            line_data = [self.readEtabLine(line) for line in lines]
        except ValueError:
            line_data = [self.readEtabLineSingleA(line) for line in lines]
        unique_res = set()
        for data in line_data:
            if len(data) == 4:
                # self energy
                unique_res.add((data['chainid'],data['resnum']))
        self.residue_info = [res for res in sorted(list(unique_res))]
        res2idx = {res:idx for idx,res in enumerate(self.residue_info)} # key: (chainid,resnum), val: res_idx

        # get the self and pair energies
        for data in line_data:
            if len(data) == 4: 
                #self energy term
                res_idx = res2idx[(data['chainid'],data['resnum'])]
                aa_idx = AA_to_int[data['aa']]
                if res_idx not in self_table:
                    self_table[res_idx] = dict()
                self_table[res_idx][aa_idx] = data['energy']
            elif len(data) == 7: 
                #pair energy term
                res_i_idx = res2idx[(data['chainid_i'],data['resnum_i'])]
                res_j_idx = res2idx[(data['chainid_j'],data['resnum_j'])]
                aa_i_idx = AA_to_int[data['aa_i']]
                aa_j_idx = AA_to_int[data['aa_j']]

                pair_index = (res_i_idx,res_j_idx) if (res_i_idx < res_j_idx) else (res_j_idx,res_i_idx)
                pair_amino_acid = (aa_i_idx,aa_j_idx) if (res_i_idx < res_j_idx) else (aa_j_idx,aa_i_idx)

                if pair_index not in pair_table:
                    pair_table[pair_index] = dict()

                pair_table[pair_index][pair_amino_acid] = data['energy']
            else:
                print("line does not match the expected .etab formatting")

        # store in class attributes
        for res_idx in self_table:
            res_self_e = [self_table[res_idx][aa_idx] for aa_idx in range(self.aa_size)]
            self.self_energies.append(res_self_e)

        for pair in pair_table:
            res_i_idx = pair[0]
            res_j_idx = pair[1]
            pair_energies = []
            for aa_i_idx in range(self.aa_size):
                for aa_j_idx in range(self.aa_size):
                    pair_e = pair_table[pair][(aa_i_idx,aa_j_idx)] if pair in pair_table and (aa_i_idx,aa_j_idx) in pair_table[pair] else 0.0
                    pair_energies.append(pair_e)
            self.pair_energies[(res_i_idx,res_j_idx)] = pair_energies

    # def writeEtab(self,name):
    #     with open(name+'.etab','w') as file:
    #         # self energies
    #         residueIterator = self.residuesIterator()
    #         for res_idx in residueIterator:
    #             for aa_idx in range(self.aa_size):
    #                 file.write(f"{self.residue_info[res_idx][0]},{self.residue_info[res_idx][1]} {int_to_3lt_AA[aa_idx]} {self.self_energies[res_idx][aa_idx]}")
                
    #         # pair energies

    def getEnergy(self,sequence):
        assert len(sequence) == len(self.residue_info)
        total_energy = 0.0
        residueIterator = self.residuesIterator()
        for res_idx in residueIterator:
            # self energy
            aa_idx = AA_to_int[sequence[res_idx]]
            total_energy += self.getSelfEnergy(res_idx,aa_idx)

            # pair energy
            pairs = self.getPairsContainingResidue(res_idx)
            for res_i_idx,res_j_idx in pairs:
                # avoid double counting pair energies
                if (res_i_idx > res_j_idx):
                    continue
                aa_i_idx = AA_to_int[sequence[res_i_idx]]
                aa_j_idx = AA_to_int[sequence[res_j_idx]]
                total_energy += self.getPairEnergy(res_i_idx,res_j_idx,aa_i_idx,aa_j_idx)
        return total_energy
        

    def getSelfEnergy(self,position_idx,aa_idx):
        return self.self_energies[position_idx][aa_idx]

    def getPairEnergy(self,position_i_idx,position_j_idx,aa_i_idx,aa_j_idx):
        if position_i_idx > position_j_idx:
            position_i_idx, position_j_idx = position_j_idx, position_i_idx
        return self.pair_energies[(position_i_idx,position_j_idx)][self.getAAPairIdx(aa_i_idx,aa_j_idx)]

    def getAAPairIdx(self,aa_i_idx,aa_j_idx):
        # (2,2) -> 4
        #    1  2  3
        # 1 -> -> ->
        # 2 ->  X
        return aa_i_idx*self.aa_size + aa_j_idx

    def getPairsContainingResidue(self,res_idx):
        return [tuple(sorted([res_i_idx,res_j_idx])) for (res_i_idx,res_j_idx) in self.pair_energies if (res_idx == res_i_idx or res_idx == res_j_idx)]

    def residuesIterator(self):
        for idx,res in enumerate(self.self_energies):
            yield idx

    def selfEnergiesIterator(self):
        for idx,res in enumerate(self.self_energies):
            for aa_idx in range(self.aa_size):
                yield (idx,aa_idx,self.getSelfEnergy(idx,aa_idx))

    def pairEnergiesIterator(self):
        for pair in self.pair_energies:
            for aa_i_idx in range(self.aa_size):
                for aa_j_idx in range(self.aa_size):
                    yield (pair[0],pair[1],aa_i_idx,aa_j_idx,self.getPairEnergy(pair[0],pair[1],aa_i_idx,aa_j_idx))

class LPDesign:
    """
    A class for defining a LP optimization problem given a Potts Model and constraints.

    The general framework is borrowed from https://github.com/KeatingLab/sortcery_design/blob/master/notebooks/08_NOTEBOOK_design_06%20-%20specificity%20-%20poly-%202%20vs%201.ipynb
    Each self/pair energy is defined as a distinct binary variable with the energy assigned as 
    the weight. Constraints are applied to ensure consistency in the final solution (e.g. each
    position may have only one amino acid selected).    
    ...

    Attributes
    ----------
    model : PottsModel
        A potts model describing a given protein structure
    problem : pulp.LpProblem
        The optimization problem
    solver : pulp.solver
        The solver that will perform the optimization
    
    Methods
    -------
    method_name(var=Type)
        description
    """

    def __init__(self, _model: PottsModel, _solverName='CPLEX_PY'):
        """
        Parameters
        ----------
        _model : PottsModel
            A potts model describing the structure for which sequence is to be designed

        _solver : str
            The solver that should be used (default is the python API for CPLEX).

        """
        self.model = _model
        self.problem = pulp.LpProblem("problem", pulp.LpMinimize)
        self.solver = pulp.getSolver(_solverName)
        self.verbose = True
        self.n_variable_res = None
        # see https://coin-or.github.io/pulp/guides/how_to_configure_solvers.html for help with installing/configuring the solver

        # dictionary[variable_name] = lp.LpVariable
        # self.all_variables = {}
        self.self_variables = {}
        self.pair_variables = {}

        self.defineVariables()
        self.setConstraints()
        self.setObjective()

    def defineVariables(self):
        """Define the variables corresponding to self and pair energies at each position in the energy table"""
        self.n_variable_res = len(self.model.residue_info)

        # Define self energies
        if self.verbose:
            print("Define self energies")
        selfe_iterator = self.model.selfEnergiesIterator()
        for idx,aa_idx,_ in selfe_iterator:
            selfe_name = self.getSelfEName(idx,aa_idx)
            var = pulp.LpVariable(selfe_name,0,1,pulp.LpBinary)
            # self.variable_names[selfe_name] = var
            self.self_variables[selfe_name] = var
            
        # Define pair energies
        if self.verbose:
            print("Define pair energies")
        paire_iterator = self.model.pairEnergiesIterator()
        for res_i_idx,res_j_idx,aa_i_idx,aa_j_idx,_ in paire_iterator:
            paire_name = self.getPairEName(res_i_idx,res_j_idx,aa_i_idx,aa_j_idx)
            var = pulp.LpVariable(paire_name,0,1,pulp.LpBinary)
            # self.variable_names[paire_name] = var
            self.pair_variables[paire_name] = var

    def setConstraints(self):
        """
        Apply constraints to ensure consistency in the final solution
        
        Three forms of consistency:
        1) self energy consistency: a single self energy is selected (i.e. non-zero) per position
        2) pair energy consistency: a single pair energy is selected per pair of positions
        3) self-pair energy consistency: the pair energy must match the amino acids selected at each position
        """
        
        ## self energy consistency
        if self.verbose:
            print("Apply self energy constraints")
        residue_iterator = self.model.residuesIterator()
        for res_idx in residue_iterator:
            # ensures that only one energy can have a value of 1, per position
            self.problem += pulp.lpSum(self.self_variables[self.getSelfEName(res_idx,aa_idx)] for aa_idx in range(self.model.aa_size)) == 1

        if self.verbose:
            print("Apply pair and self-pair energy constraints")
        residue_iterator = self.model.residuesIterator()
        for res_idx in residue_iterator:
            # get pairs by their residue indices, (i,j)
            pairs = self.model.getPairsContainingResidue(res_idx)
            # only consider i < j to avoid double counting
            pairs = [(i,j) for (i,j) in pairs if i < j]

            for (res_i_idx,res_j_idx) in pairs:
                ## pair energy consistency
                self.problem += pulp.lpSum(self.pair_variables[self.getPairEName(res_i_idx,res_j_idx,aa_i_idx,aa_j_idx)] for aa_i_idx in range(self.model.aa_size) for aa_j_idx in range(self.model.aa_size)) == 1

                ## self-pair energy consistency
                # ensure that res_i is consistent with the pair energy
                for aa_i_idx in range(self.model.aa_size):
                    self_pair_constraints = [(self.self_variables[self.getSelfEName(res_i_idx,aa_i_idx)],-1)]
                    self_pair_constraints += [(self.pair_variables[self.getPairEName(res_i_idx,res_j_idx,aa_i_idx,aa_j_idx)],1) for aa_j_idx in range(self.model.aa_size)]
                    self.problem += pulp.LpAffineExpression(var for var in self_pair_constraints) == 0

                # ensure that res_j is consistent with the pair energy
                for aa_j_idx in range(self.model.aa_size):
                    self_pair_constraints = [(self.self_variables[self.getSelfEName(res_j_idx,aa_j_idx)],-1)]
                    self_pair_constraints += [(self.pair_variables[self.getPairEName(res_i_idx,res_j_idx,aa_i_idx,aa_j_idx)],1) for aa_i_idx in range(self.model.aa_size)]
                    self.problem += pulp.LpAffineExpression(var for var in self_pair_constraints) == 0

    def setSequenceConstraint(self,constraint_sequence):
        '''Constraint certain positions to a single sequence'''
        assert len(constraint_sequence) == len(self.model.self_energies)

        counter = 0
        for res_idx,aa in enumerate(constraint_sequence):
            if aa == 'X':
                # no constraint at this position
                counter+=1
            else:
                # ensures that the self energy for the select amino acid must be used in the optimal sequence
                self.problem += self.self_variables[self.getSelfEName(res_idx,AA_to_int[aa])] == 1
        self.n_variable_res = counter
        print(f"{self.n_variable_res} of {len(self.model.residue_info)} residues are variable ")

    def setObjective(self):
        """Define the LP optimization objective"""
        if self.verbose:
            print("Set objective function")
        selfe_iterator = self.model.selfEnergiesIterator()
        all_energies = [(self.self_variables[self.getSelfEName(idx,aa_idx)],E) for idx,aa_idx,E in selfe_iterator]
        paire_iterator = self.model.pairEnergiesIterator()
        all_energies += [(self.pair_variables[self.getPairEName(res_i_idx,res_j_idx,aa_i_idx,aa_j_idx)],E)
                            for res_i_idx,res_j_idx,aa_i_idx,aa_j_idx,E in paire_iterator]
        self.problem += pulp.LpAffineExpression(var for var in all_energies)

    def findSolution(self):
        """Call the external solver to find a solution to the problem"""
        if self.verbose:
            print("Try to find solution")
        status = self.problem.solve(self.solver)
        print('CPLEX complete')
        optSeq = ''
        if pulp.LpStatus[status]=='Optimal':
            optSeq = self.getSequenceFromVariables()
            print('Optimal sequence:',optSeq)

            # outFile.write(topX[-1])
            # outFile.write("\n")
            # outFile.flush()
        else:
            print(pulp.LpStatus[status])
        return optSeq

    def getSequenceFromVariables(self):
        """Get sequence and verify that self constraints are satisfied"""
        sequence = ''
        residue_iterator = self.model.residuesIterator()
        for res_idx in residue_iterator:
            res_aa = ''
            for aa_idx in range(self.model.aa_size):
                var = self.self_variables[self.getSelfEName(res_idx,aa_idx)]
                if pulp.value(var) > 10**-6:
                    if res_aa == '':
                        res_aa = int_to_AA[aa_idx]
                    else:
                        raise ValueError('Only a single amino acid type should be selected per residue')
            sequence+=res_aa
        return sequence
            
    def getSelfEName(self,idx,aa_idx):
        return f"{idx}_{int_to_AA[aa_idx]}"

    def getPairEName(self,res_i_idx,res_j_idx,aa_i_idx,aa_j_idx):
        return f"{res_i_idx}_{res_j_idx}_{int_to_AA[aa_i_idx]}_{int_to_AA[aa_j_idx]}"

    def addCurrentSolutionToConstraints(self, hamming_distance):
        """Add the current solution as a constraint such that future solutions have at least hamming_distance mutations relative to this solution"""
        print("Adding current solution as a constraint prior to finding another solution")
        residue_iterator = self.model.residuesIterator()
        solution_self_e_list = []
        for res_idx in residue_iterator:
            res_aa = ''
            for aa_idx in range(self.model.aa_size):
                var = self.self_variables[self.getSelfEName(res_idx,aa_idx)]
                if pulp.value(var) > 10**-6:
                    solution_self_e_list.append(var)
        self.problem += pulp.lpSum(var for var in solution_self_e_list) < len(self.model.residue_info) - hamming_distance

def NSR(seq1,seq2):
    assert len(seq1) == len(seq2)
    return np.sum(np.array([x==y for x,y in zip(seq1,seq2)]))/len(seq1)

# # pylint: disable=broad-except
# def _knn_etab_to_sequence_wrapper(etab, E_idx, res_info, out_file, native_seq):
#     """Wrapper for _etab_to_sequence_wrapper that does error handling"""
#     try:
#         return knn_etab_to_sequence(etab, E_idx, res_info, out_file, native_seq)
#     except Exception:
#         eprint(out_file.name)
#         eprint(idx_dict)
#         traceback.print_exc()
#         return False, out_file.name

class optimizeEnergyTableViaLP:
    def __init__(self,process_number=1):
        """
        """
        # Open file to write output to
        # self.out_file = open(file_prefix+'_optimalsequences.csv','w')
        self.out_file = open(f"optimalsequences_{process_number}.csv",'w')
        self.out_file.write('name,constraint,designed_sequence,native_sequence,nsr,designed_sequence_energy,native_sequence_energy\n')

    def knn_etab_to_sequence(self, name, etab, E_idx, res_info, constraint = "", constraint_seq = None, native_seq = None):
        """Find the optimal sequence given the etab

        If a native sequence is provided, will compute NSR

        Args
        ====
        etab : np.ndarray
            Etab outputted by TERMinator
        E_idx : np.ndarray
            Indexing matrix associated with :code:`etab_matrix`
        res_info : dict
            Chain ID and resnum for all residues in the structure
        native_seq : str
            The sequence on the PDB file that was used for design

        Returns
        =======
        bool
            Whether or not the parsing occured without errors
        out_path : str
            The output path for the provided file handle
        """

        potts_model = PottsModel()
        potts_model.loadFromKNNSparse(etab,E_idx,res_info)

        lp_optimizer = LPDesign(potts_model)
        lp_optimizer.setSequenceConstraint(constraint_seq)
        optimal_sequence = lp_optimizer.findSolution()
        optimal_sequence_e = potts_model.getEnergy(optimal_sequence)

        if native_seq == None:
            # dummy native sequence
            native_seq = ''.join(['A']*len(optimal_sequence))
        native_sequence_e = potts_model.getEnergy(native_seq)
        
        # write to file
        self.out_file.write(f"{name},{constraint},{optimal_sequence},{native_seq},{NSR(optimal_sequence,native_seq)},{optimal_sequence_e},{native_sequence_e}"+'\n')
        return True

    def dtermen_etab_to_sequence(self, etab_path, name = None, native_seq = None):
        if name == None:
            name = etab_path.split('/')[-1][:-5]
        potts_model = PottsModel()
        potts_model.loadFromEtabFile(etab_path)

        lp_optimizer = LPDesign(potts_model)
        optimal_sequence = lp_optimizer.findSolution()

        optimal_sequence_e = potts_model.getEnergy(optimal_sequence)
        if native_seq == None:
            # dummy native sequence
            native_seq = ['A']*len(optimal_sequence)
        native_sequence_e = potts_model.getEnergy(native_seq)
        
        # write to file
        self.out_file.write(f"{name},{''},{optimal_sequence},{native_seq},{NSR(optimal_sequence,native_seq)},{optimal_sequence_e},{native_sequence_e}"+'\n')

        return True