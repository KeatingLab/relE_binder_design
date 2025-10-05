# zero is used as padding
AA_to_int = {
    'A': 1,
    'ALA': 1,
    'C': 2,
    'CYS': 2,
    'D': 3,
    'ASP': 3,
    'E': 4,
    'GLU': 4,
    'F': 5,
    'PHE': 5,
    'G': 6,
    'GLY': 6,
    'H': 7,
    'HIS': 7,
    'I': 8,
    'ILE': 8,
    'K': 9,
    'LYS': 9,
    'L': 10,
    'LEU': 10,
    'M': 11,
    'MET': 11,
    'N': 12,
    'ASN': 12,
    'P': 13,
    'PRO': 13,
    'Q': 14,
    'GLN': 14,
    'R': 15,
    'ARG': 15,
    'S': 16,
    'SER': 16,
    'T': 17,
    'THR': 17,
    'V': 18,
    'VAL': 18,
    'W': 19,
    'TRP': 19,
    'Y': 20,
    'TYR': 20,
    'X': 21
}

AA_to_int = {key: val - 1 for key, val in AA_to_int.items()}

int_to_AA = {y: x for x, y in AA_to_int.items() if len(x) == 1}

int_to_3lt_AA = {y: x for x, y in AA_to_int.items() if len(x) == 3}

def seq_to_ints(sequence):
    """
    Given a string of one-letter encoded AAs, return its corresponding integer encoding
    """
    return [AA_to_int[residue] for residue in sequence]

def ints_to_seq(int_list):
    return [int_to_AA[i] for i in int_list]

def aa_three_to_one(residue):
    return int_to_AA[AA_to_int[residue]]

def aa_to_similar(residue):
    if residue == 'MSE':  # convert MSE (seleno-met) to MET
        residue = 'MET'
    elif residue == 'FME':  # convert MSE (n-formylmethionine) to MET
        residue = 'MET'
    elif residue == 'SEP':  # convert SEP (phospho-ser) to SER
        residue = 'SER'
    elif residue == 'SAC':  # convert SAC (n-acetyl-ser) to SER
        residue = 'SER'
    elif residue == 'OAS':  # convert OAS (o-acetyl-ser) to SER
        residue = 'SER'
    elif residue == 'TPO':  # convert TPO (phospho-thr) to THR
        residue = 'THR'
    elif residue == 'IYT':  # convert IYT (n-alpha-acetyl-3,5-diiodotyrosyl-d-threonine) to THR
        residue = 'THR'
    elif residue == 'PTR':  # convert PTR (phospho-tyr) to TYR
        residue = 'TYR'
    elif residue == 'TYS':  # convert TYS (o-sulfo-l-tyr) to TYR
        residue = 'TYR'
    elif residue == 'CSO':  # convert CSO (hydroxy-cys) to CYS
        residue = 'CYS'
    elif residue == 'SEC':  # convert SEC (seleno-cys) to CYS
        residue = 'CYS'
    elif residue == 'CSS':  # convert CSS (s-mercaptocysteine) to CYS
        residue = 'CYS'
    elif residue == 'CAS':  # convert CAS (s-(dimethylarsenic)cysteine) to CYS
        residue = 'CYS'
    elif residue == 'CAF':  # convert CAF (s-dimethylarsinoyl-cysteine) to CYS
        residue = 'CYS'
    elif residue == 'OCS':  # convert OCS (cysteine sulfonic acid) to CYS
        residue = 'CYS'
    elif residue == 'CSD':  # convert CSD (3-sulfinoalanine) to CYS
        residue = 'CYS'
    elif residue == 'CME':  # convert CME (s,s-(2-hydroxyethyl)thiocysteine) to CYS
        residue = 'CYS'
    elif residue == 'YCM':  # convert YCM (s-(2-amino-2-oxoethyl)-l-cysteine) to CYS
        residue = 'CYS'
    elif residue == 'SAH':  # convert SAH (s-adenosyl-l-homocysteine) to CYS
        residue = 'CYS'
    elif residue == 'HYP':  # convert HYP (4-hydroxyproline) to PRO
        residue = 'PRO'
    elif residue == 'M3L':  # convert M3L (n-trimethyllysine) to LYS
        residue = 'LYS'
    elif residue == 'LLP':  # convert LLP (n'-pyridoxyl-lysine-5'-monophosphate) to LYS
        residue = 'LYS'
    elif residue == 'KPI':  # convert KPI ((2s)-2-amino-6-[(1-hydroxy-1-oxo-propan-2-ylidene)amino]hexanoic acid) to LYS
        residue = 'LYS'
    elif residue == 'KPX':  # convert KPX (lysine nz-corboxylic acid) to LYS
        residue = 'LYS'
    elif residue == 'MLY':  # convert MLY (n-dimethyl-lysine) to LYS
        residue = 'LYS'
    elif residue == 'KCX':  # convert KCX (lysine nz-carboxylic acid) to LYS
        residue = 'LYS'
    elif residue == 'PCA':  # convert PCA (pyroglutamic acid) to GLN
        residue = 'GLN'
    elif residue == 'DGL':  # convert DGL (d-glutamic acid) to GLU
        residue = 'GLU'
    elif residue == 'BHD':  # convert BHD (beta-hydroxyaspartic acid) to ASP
        residue = 'ASP'
    elif residue == 'IAS':  # convert IAS (beta-l-aspartic acid) to ASP
        residue = 'ASP'
    elif residue == 'ABA':  # convert ABA (alpha-aminobutyric acid) to ALA
        residue = 'ALA'
    elif residue == '0A9':  # convert 0A9 (methyl l-phenylalaninate) to PHE
        residue = 'PHE'
    elif residue == 'KYN':  # convert KYN (l-kynurenine) to TRP
        residue = 'TRP'
    elif residue == '0Q4':
        residue = '0Q4'
    elif residue == 'FOL':
        residue = 'FOL'
    return residue
        
