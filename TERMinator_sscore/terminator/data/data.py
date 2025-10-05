"""Datasets and dataloaders for loading TERMs.

This file contains dataset and dataloader classes
to be used when interacting with TERMs.
"""
import glob
import math
import multiprocessing as mp
import os
import pdb
import pickle
import random
import sys
from tkinter import N
import numpy as np
import torch
import torch.nn.functional as F
import torch_cluster
import torch_geometric
import torch.distributed as dist
import torch.linalg as LA
# from nrgten.encom import ENCoM #swans
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import Dataset, Sampler , IterableDataset
from torch.utils.data.distributed import DistributedSampler
from tqdm import tqdm
import copy

from scripts.data.preprocessing.cleanStructs import extractBackbone
from scripts.data.preprocessing.packageTensors import concatFeatureDicts, dumpCoordsTensors

# pylint: disable=no-member, not-callable


# Jing featurization functions


def _normalize(tensor, dim=-1):
    '''Normalizes a `torch.Tensor` along dimension `dim` without `nan`s.'''
    return torch.nan_to_num(torch.div(tensor, torch.norm(tensor, dim=dim, keepdim=True)))


def _rbf(D, D_min=0., D_max=20., D_count=16, device='cpu'):
    '''Returns an RBF embedding of `torch.Tensor` `D` along a new axis=-1.

    That is, if `D` has shape [...dims], then the returned tensor will have
    shape [...dims, D_count].

    From https://github.com/jingraham/neurips19-graph-protein-design
    '''
    D_mu = torch.linspace(D_min, D_max, D_count, device=device)
    D_mu = D_mu.view([1, -1])
    D_sigma = (D_max - D_min) / D_count
    D_expand = torch.unsqueeze(D, -1)

    rbf = torch.exp(-((D_expand - D_mu) / D_sigma)**2)
    return rbf


def _dihedrals(X, eps=1e-7):
    """ Compute dihedral angles between residues given atomic backbone coordinates

    Args
    ----
    X : torch.FloatTensor
        Tensor specifying atomic backbone coordinates
        Shape: num_res x 4 x 3

    Returns
    -------
    D_features : torch.FloatTensor
        Dihedral angles, lifted to the 3-torus
        Shape: num_res x 7
    """
    # From https://github.com/jingraham/neurips19-graph-protein-design

    X = torch.reshape(X[:, :3], [3 * X.shape[0], 3])
    dX = X[1:] - X[:-1]
    U = _normalize(dX, dim=-1)
    u_2 = U[:-2]
    u_1 = U[1:-1]
    u_0 = U[2:]

    # Backbone normals
    n_2 = _normalize(torch.cross(u_2, u_1), dim=-1)
    n_1 = _normalize(torch.cross(u_1, u_0), dim=-1)

    # Angle between normals
    cosD = torch.sum(n_2 * n_1, -1)
    cosD = torch.clamp(cosD, -1 + eps, 1 - eps)
    D = torch.sign(torch.sum(u_2 * n_1, -1)) * torch.acos(cosD)

    # This scheme will remove phi[0], psi[-1], omega[-1]
    D = F.pad(D, [1, 2])
    D = torch.reshape(D, [-1, 3])
    # Lift angle representations to the circle
    D_features = torch.cat([torch.cos(D), torch.sin(D)], 1)
    return D_features


def _positional_embeddings(edge_index, num_embeddings=16, dev='cpu'):
    """ Sinusoidally encode sequence distances for edges.

    Args
    ----
    edge_index : torch.LongTensor
        Edge indices for sparse representation of protein graph
        Shape: 2 x num_edges
    num_embeddings : int or None, default=128
        Dimensionality of sinusoidal embedding.

    Returns
    -------
    E : torch.FloatTensor
        Sinusoidal encoding of sequence distances
        Shape: num_edges x num_embeddings

    """
    # From https://github.com/jingraham/neurips19-graph-protein-design
    d = edge_index[0] - edge_index[1]

    frequency = torch.exp(
        torch.arange(0, num_embeddings, 2, dtype=torch.float32, device=dev) * -(np.log(10000.0) / num_embeddings))
    angles = d.unsqueeze(-1) * frequency
    E = torch.cat((torch.cos(angles), torch.sin(angles)), -1)
    return E


def _orientations(X_ca):
    """ Compute forward and backward vectors per residue.

    Args
    ----
    X_ca : torch.FloatTensor
        Tensor specifying atomic backbone coordinates for CA atoms.
        Shape: num_res x 3

    Returns
    -------
    torch.FloatTensor
        Pairs of forward, backward vectors per residue.
        Shape: num_res x 2 x 3
    """
    # From https://github.com/drorlab/gvp-pytorch
    forward = _normalize(X_ca[1:] - X_ca[:-1])
    backward = _normalize(X_ca[:-1] - X_ca[1:])
    forward = F.pad(forward, [0, 0, 0, 1])
    backward = F.pad(backward, [0, 0, 1, 0])
    return torch.cat([forward.unsqueeze(-2), backward.unsqueeze(-2)], -2)


def _sidechains(X):
    """ Compute vectors pointing in the approximate direction of the sidechain.

    Args
    ----
    X : torch.FloatTensor
        Tensor specifying atomic backbone coordinates.
        Shape: num_res x 4 x 3

    Returns
    -------
    vec : torch.FloatTensor
        Sidechain vectors.
        Shape: num_res x 3
    """
    # From https://github.com/drorlab/gvp-pytorch
    n, origin, c = X[:, 0], X[:, 1], X[:, 2]
    c, n = _normalize(c - origin), _normalize(n - origin)
    bisector = _normalize(c + n)
    perp = _normalize(torch.cross(c, n))
    vec = -bisector * math.sqrt(1 / 3) - perp * math.sqrt(2 / 3)
    return vec


def _jing_featurize(protein, dev='cpu'):
    """ Featurize individual proteins for use in torch_geometric Data objects,
    as done in https://github.com/drorlab/gvp-pytorch

    Args
    ----
    protein : dict
        Dictionary of protein features

        - :code:`name` - PDB ID of the protein
        - :code:`coords` - list of dicts specifying backbone atom coordinates
        in the format of that outputted by :code:`parseCoords.py`
        - :code:`seq` - protein sequence
        - :code:`chain_idx` - an integer per residue such that each unique integer represents a unique chain

    Returns
    -------
    torch_geometric.data.Data
        Data object containing
        - :code:`x` - CA atomic coordinates
        - :code:`seq` - sequence of protein
        - :code:`name` - PDB ID of protein
        - :code:`node_s` - Node scalar features
        - :code:`node_v` - Node vector features
        - :code:`edge_s` - Edge scalar features
        - :code:`edge_v` - Edge vector features
        - :code:`edge_index` - Sparse representation of edge
        - :code:`mask` - Residue mask specifying residues with incomplete coordinate sets
    """
    name = protein['name']
    with torch.no_grad():
        coords = torch.as_tensor(protein['coords'], device=dev, dtype=torch.float32)
        seq = torch.as_tensor(protein['seq'], device=dev, dtype=torch.long)

        mask = torch.isfinite(coords.sum(dim=(1, 2)))
        coords[~mask] = np.inf

        X_ca = coords[:, 1]
        edge_index = torch_cluster.knn_graph(X_ca, k=30, loop=True)  # TODO: make param

        pos_embeddings = _positional_embeddings(edge_index)
        # generate mask for interchain interactions
        pos_chain = (protein['chain_idx'][edge_index.view(-1)]).view(2, -1)
        pos_mask = (pos_chain[0] != pos_chain[1])
        # zero out all interchain positional embeddings
        pos_embeddings = pos_mask.unsqueeze(-1) * pos_embeddings

        E_vectors = X_ca[edge_index[0]] - X_ca[edge_index[1]]
        rbf = _rbf(E_vectors.norm(dim=-1), D_count=16, device=dev)  # TODO: make param

        dihedrals = _dihedrals(coords)
        orientations = _orientations(X_ca)
        sidechains = _sidechains(coords)

        node_s = dihedrals
        node_v = torch.cat([orientations, sidechains.unsqueeze(-2)], dim=-2)
        edge_s = torch.cat([rbf, pos_embeddings], dim=-1)
        edge_v = _normalize(E_vectors).unsqueeze(-2)

        node_s, node_v, edge_s, edge_v = map(torch.nan_to_num, (node_s, node_v, edge_s, edge_v))

    data = torch_geometric.data.Data(x=X_ca,
                                     seq=seq,
                                     name=name,
                                     node_s=node_s,
                                     node_v=node_v,
                                     edge_s=edge_s,
                                     edge_v=edge_v,
                                     edge_index=edge_index,
                                     mask=mask)
    return data


# Ingraham featurization functions


def _ingraham_featurize(batch, device="cpu", data_type="coord"):
    """ Pack and pad data in batch into torch tensors
    as done in https://github.com/jingraham/neurips19-graph-protein-design

    Args
    ----
    batch : list of dict
        list of protein backbone coordinate dictionaries,
        in the format of that outputted by :code:`parseCoords.py`
    device : str
        device to place torch tensors
    type : str
        type of data to featurize

    Returns
    -------
    X : torch.Tensor
        Batched coordinates tensor
    mask : torch.Tensor
        Mask for X
    lengths : np.ndarray
        Array of lengths of batched proteins
    """
    B = len(batch)
    lengths = np.array([b.shape[0] for b in batch], dtype=np.int32)
    l_max = max(lengths)
    if data_type == "coord":
        X = np.zeros([B, l_max, 4, 3])
    elif data_type.find("dynamic_signatures") > -1 or data_type.find("conformational_pos_variation") > -1:
        X = np.zeros([B, l_max, 1])
    elif data_type.find("conformational_ensembles") > -1:
        num_ensembles = np.array([b.shape[3] for b in batch], dtype=np.int32)
        if not np.all(num_ensembles == num_ensembles[0]):
            raise ValueError("All conformational ensembles must have the same number of ensembles.")
        X = np.zeros([B, l_max, 4, 3, num_ensembles[0]])
    else:
        raise ValueError("Type for ingraham featurization must be \'coord\' or \'dynamic_signatures\' or \'conformational_ensembles\' or \'conformational_pos_variation\'.")

    # Build the batch
    for i, x in enumerate(batch):
        l = x.shape[0]
        if data_type == "coord":
            x_pad = np.pad(x, [[0, l_max - l], [0, 0], [0, 0]], 'constant', constant_values=(np.nan, ))
            X[i, :, :, :] = x_pad
        elif data_type.find("dynamic_signatures") > -1 or data_type.find("conformational_pos_variation") > -1:
            print("padding")
            print(x.shape)
            x_pad = np.pad(x, [[0, l_max - l], [0, 0]], 'constant', constant_values=(np.nan, ))
            X[i, :, :] = x_pad
        elif data_type.find("conformational_ensembles") > -1:
            x_pad = np.pad(x, [[0, l_max - l], [0, 0], [0, 0], [0, 0]], 'constant', constant_values=(np.nan, ))
            X[i, :, :, :, :] = x_pad
        

    # Mask
    isnan = np.isnan(X)
    if data_type == "coord":
        mask = np.isfinite(np.sum(X, (2, 3))).astype(np.float32)
    elif data_type.find("dynamic_signatures") > -1 or data_type.find("conformational_pos_variation") > -1:
        mask = np.isfinite(np.sum(X, 2)).astype(np.float32)
    elif data_type.find("conformational_ensembles") > -1:
        mask = np.isfinite(np.sum(X, (2, 3, 4))).astype(np.float32)
    X[isnan] = 0.

    # Conversion
    X = torch.from_numpy(X).to(dtype=torch.float32, device=device)
    mask = torch.from_numpy(mask).to(dtype=torch.float32, device=device)
    return X, mask, lengths

def _quaternions(R, eps=1e-10):
    """ Convert a batch of 3D rotations [R] to quaternions [Q]
        R [...,3,3]
        Q [...,4]
    """
    def _R(i, j):
        return R[..., i, j]

    # Simple Wikipedia version
    # en.wikipedia.org/wiki/Rotation_matrix#Quaternion
    # For other options see math.stackexchange.com/questions/2074316/calculating-rotation-axis-from-rotation-matrix
    diag = torch.diagonal(R, dim1=-2, dim2=-1)
    Rxx, Ryy, Rzz = diag.unbind(-1)
    magnitudes = 0.5 * torch.sqrt(
        torch.abs(1 + torch.stack([Rxx - Ryy - Rzz, -Rxx + Ryy - Rzz, -Rxx - Ryy + Rzz], -1)) + eps)
    signs = torch.sign(torch.stack([_R(2, 1) - _R(1, 2), _R(0, 2) - _R(2, 0), _R(1, 0) - _R(0, 1)], -1))
    xyz = signs * magnitudes
    # The relu enforces a non-negative trace
    w = torch.sqrt(F.relu(1 + diag.sum(-1, keepdim=True))) / 2.
    Q = torch.cat((xyz, w), -1)
    Q = F.normalize(Q, dim=-1)

    return Q


def _orientations_coarse(X, edge_index, eps=1e-6):
    # Pair features

    # Shifted slices of unit vectors
    dX = X[1:, :] - X[:-1, :]
    U = F.normalize(dX, dim=-1)
    u_2 = U[:-2, :]
    u_1 = U[1:-1, :]
    u_0 = U[2:, :]
    # Backbone normals
    n_2 = F.normalize(torch.cross(u_2, u_1), dim=-1)
    n_1 = F.normalize(torch.cross(u_1, u_0), dim=-1)

    # Bond angle calculation
    cosA = -(u_1 * u_0).sum(-1)
    cosA = torch.clamp(cosA, -1 + eps, 1 - eps)
    A = torch.acos(cosA)
    # Angle between normals
    cosD = (n_2 * n_1).sum(-1)
    cosD = torch.clamp(cosD, -1 + eps, 1 - eps)
    D = torch.sign((u_2 * n_1).sum(-1)) * torch.acos(cosD)
    # Backbone features
    AD_features = torch.stack((torch.cos(A), torch.sin(A) * torch.cos(D), torch.sin(A) * torch.sin(D)), -1)
    AD_features = F.pad(AD_features, (0, 0, 1, 2), 'constant', 0)

    # Build relative orientations
    o_1 = F.normalize(u_2 - u_1, dim=-1)
    O = torch.stack((o_1, n_2, torch.cross(o_1, n_2)), 2)
    O = O.view(list(O.shape[:1]) + [9])
    O = F.pad(O, (0, 0, 1, 2), 'constant', 0)

    # DEBUG: Viz [dense] pairwise orientations
    # O = O.view(list(O.shape[:2]) + [3,3])
    # dX = X.unsqueeze(2) - X.unsqueeze(1)
    # dU = torch.matmul(O.unsqueeze(2), dX.unsqueeze(-1)).squeeze(-1)
    # dU = dU / torch.norm(dU, dim=-1, keepdim=True)
    # dU = (dU + 1.) / 2.
    # plt.imshow(dU.data.numpy()[0])
    # plt.show()
    # print(dX.size(), O.size(), dU.size())
    # exit(0)

    O_pairs = O[edge_index]
    X_pairs = X[edge_index]

    # Re-view as rotation matrices
    O_pairs = O_pairs.view(list(O_pairs.shape[:-1]) + [3,3])

    # Rotate into local reference frames
    dX = X_pairs[0] - X_pairs[1]
    dU = torch.matmul(O_pairs[1], dX.unsqueeze(-1)).squeeze(-1)
    dU = F.normalize(dU, dim=-1)
    R = torch.matmul(O_pairs[1].transpose(-1, -2), O_pairs[0])
    Q = _quaternions(R)

    # Orientation features
    O_features = torch.cat((dU, Q), dim=-1)

    # DEBUG: Viz pairwise orientations
    # IMG = Q[:,:,:,:3]
    # # IMG = dU
    # dU_full = torch.zeros(X.shape[0], X.shape[1], X.shape[1], 3).scatter(
    #     2, E_idx.unsqueeze(-1).expand(-1,-1,-1,3), IMG
    # )
    # print(dU_full)
    # dU_full = (dU_full + 1.) / 2.
    # plt.imshow(dU_full.data.numpy()[0])
    # plt.show()
    # exit(0)
    # print(Q.sum(), dU.sum(), R.sum())
    return AD_features, O_features

def _ingraham_geometric_featurize(protein, dev='cpu'):
    """ Featurize individual proteins for use in torch_geometric Data objects,
    as done in https://github.com/drorlab/gvp-pytorch

    Args
    ----
    protein : dict
        Dictionary of protein features

        - :code:`name` - PDB ID of the protein
        - :code:`coords` - list of dicts specifying backbone atom coordinates
        in the format of that outputted by :code:`parseCoords.py`
        - :code:`seq` - protein sequence
        - :code:`chain_idx` - an integer per residue such that each unique integer represents a unique chain

    Returns
    -------
    torch_geometric.data.Data
        Data object containing
        - :code:`x` - CA atomic coordinates
        - :code:`seq` - sequence of protein
        - :code:`name` - PDB ID of protein
        - :code:`node_s` - Node scalar features
        - :code:`node_v` - Node vector features
        - :code:`edge_s` - Edge scalar features
        - :code:`edge_v` - Edge vector features
        - :code:`edge_index` - Sparse representation of edge
        - :code:`mask` - Residue mask specifying residues with incomplete coordinate sets
    """
    name = protein['name']
    with torch.no_grad():
        coords = torch.as_tensor(protein['coords'], device=dev, dtype=torch.float32)
        seq = torch.as_tensor(protein['seq'], device=dev, dtype=torch.long)

        mask = torch.isfinite(coords.sum(dim=(1, 2)))
        coords[~mask] = np.inf

        X_ca = coords[:, 1]
        edge_index = torch_cluster.knn_graph(X_ca, k=30, loop=True)  # TODO: make param

        pos_embeddings = _positional_embeddings(edge_index)
        # generate mask for interchain interactions
        pos_chain = (protein['chain_idx'][edge_index.view(-1)]).view(2, -1)
        pos_mask = (pos_chain[0] != pos_chain[1])
        # zero out all interchain positional embeddings
        pos_embeddings = pos_mask.unsqueeze(-1) * pos_embeddings

        E_vectors = X_ca[edge_index[0]] - X_ca[edge_index[1]]
        rbf = _rbf(E_vectors.norm(dim=-1), D_count=16, device=dev)  # TODO: make param

        dihedrals = _dihedrals(coords)
        _, orientations = _orientations_coarse(X_ca, edge_index)

        node_features = dihedrals
        edge_features = torch.cat([pos_embeddings, rbf, orientations], dim=-1)

        node_features, edge_features, = map(torch.nan_to_num, (node_features, edge_features))

    data = torch_geometric.data.Data(x=X_ca,
                                     seq=seq,
                                     name=name,
                                     node_features=node_features,
                                     edge_features=edge_features,
                                     edge_index=edge_index,
                                     mask=mask)
    return data


# Batching functions


def convert(tensor):
    """Converts given tensor from numpy to pytorch."""
    return torch.from_numpy(tensor)


def _package(b_idx, flex_type="", noise_level=0.0, bond_length_noise_level=0.0):
    """Package the given datapoints into tensors based on provided indices.

    Tensors are extracted from the data and padded. Coordinates are featurized
    and the length of TERMs and chain IDs are added to the data.

    residue_idx code adapted from https://github.com/dauparas/ProteinMPNN

    Args
    ----
    b_idx : list of tuples (dicts, int)
        The feature dictionaries, as well as an int for the sum of the lengths of all TERMs,
        for each datapoint to package.
    flex_type : str
        Methodology to calculate flex data.
    noise_level : float
        Std of noise to add.
    bond_length_noise_level : float
        Std of noise to add to bond lengths.

    Returns
    -------
    dict
        Collection of batched features required for running TERMinator. This contains:

        - :code:`msas` - the sequences for each TERM match to the target structure

        - :code:`features` - the :math:`\\phi, \\psi, \\omega`, and environment values of the TERM matches

        - :code:`ppoe` - the :math:`\\phi, \\psi, \\omega`, and environment values of the target structure

        - :code:`seq_lens` - lengths of the target sequences

        - :code:`focuses` - the corresponding target structure residue index for each TERM residue

        - :code:`contact_idxs` - contact indices for each TERM residue

        - :code:`src_key_mask` - mask for TERM residue padding

        - :code:`X` - coordinates

        - :code:`x_mask` - mask for the target structure

        - :code:`seqs` - the target sequences

        - :code:`ids` - the PDB ids

        - :code:`chain_idx` - the chain IDs
    """
    # wrap up all the tensors with proper padding and masks
    batch = [data[0] for data in b_idx]
    focus_lens = [data[1] for data in b_idx]
    features, msas, focuses, seq_lens, coords, flexes = [], [], [], [], [], []
    term_lens = []
    seqs = []
    ids = []
    chain_lens = []
    ppoe = []
    contact_idxs = []
    gvp_data = []
    # geometric_data = []

    sortcery_seqs = []
    sortcery_nrgs = []
    sscore_vals = []
    
    for _, data in enumerate(batch):
        # have to transpose these two because then we can use pad_sequence for padding
        features.append(convert(data['features']).transpose(0, 1))
        msas.append(convert(data['msas']).transpose(0, 1))

        ppoe.append(convert(data['ppoe']))
        focuses.append(convert(data['focuses']))
        contact_idxs.append(convert(data['contact_idxs']))
        seq_lens.append(data['seq_len'])
        term_lens.append(data['term_lens'].tolist())
        if flex_type and flex_type.find("batch") > -1:
            noise = generate_noise(flex_type=flex_type, noise_level=noise_level, size=data['coords'].shape, X=data['coords'], bond_length_noise_level=bond_length_noise_level, chain_lens=data['chain_lens'])
            data['coords'] =  np.add(data['coords'], noise)
        coords.append(data['coords'])
        seqs.append(convert(data['sequence']))
        ids.append(data['pdb'])
        chain_lens.append(data['chain_lens'])
        if (flex_type is not None) and (flex_type != ""):
            flexes.append(data['flex'])
        if flex_type and flex_type.find("random") == -1 and data['flex'] is None:
            print(f"flex is none for {data['pdb']}.")

        if 'sortcery_seqs' in data:
            assert len(batch) == 1, "batch_size for SORTCERY fine-tuning should be set to 1"
            sortcery_seqs = convert(data['sortcery_seqs']).unsqueeze(0)
        if 'sortcery_nrgs' in data:
            sortcery_nrgs = convert(data['sortcery_nrgs']).unsqueeze(0)
        if 'sscore' in data:
            sscore_vals.append(convert(np.nan_to_num(data['sscore'],0).astype(dtype="float32")))

        chain_idx = []
        for i, c_len in enumerate(data['chain_lens']):
            chain_idx.append(torch.ones(c_len) * i)
        chain_idx = torch.cat(chain_idx, dim=0)
        gvp_data.append(
            _jing_featurize({
                'name': data['pdb'],
                'coords': data['coords'],
                'seq': data['sequence'],
                'chain_idx': chain_idx
            }))
        # geometric_data.append(
        #     _ingraham_geometric_featurize({
        #         'name': data['pdb'],
        #         'coords': data['coords'],
        #         'seq': data['sequence'],
        #         'chain_idx': chain_idx
        #     }))

    # transpose back after padding
    features = pad_sequence(features, batch_first=True).transpose(1, 2)
    msas = pad_sequence(msas, batch_first=True).transpose(1, 2).long()

    # we can pad these using standard pad_sequence
    ppoe = pad_sequence(ppoe, batch_first=True)
    focuses = pad_sequence(focuses, batch_first=True)
    contact_idxs = pad_sequence(contact_idxs, batch_first=True)
    src_key_mask = pad_sequence([torch.zeros(l) for l in focus_lens], batch_first=True, padding_value=1).bool()
    seqs = pad_sequence(seqs, batch_first=True)
    if 'sscore' in data:
        sscore_vals = pad_sequence(sscore_vals, batch_first=True)

    # we do some padding so that tensor reshaping during batchifyTERM works
    # TODO(alex): explain this since I have no idea what's going on
    max_aa = focuses.size(-1)
    for lens in term_lens:
        max_term_len = max(lens)
        diff = max_aa - sum(lens)
        lens += [max_term_len] * (diff // max_term_len)
        lens.append(diff % max_term_len)

    # featurize coordinates same way as ingraham et al
    X, x_mask, _ = _ingraham_featurize(coords)

    # featurize flex data same way as ingraham et al
    if flex_type and (flex_type.find("dynamic_signatures") > -1 or flex_type.find("conformational_ensembles") > -1 or flex_type.find("conformational_pos_variation") > -1):
        flex, flex_mask, _ = _ingraham_featurize(flexes, data_type=flex_type)
    else:
        flex, flex_mask = None, None
    if flex_type and flex_type.find("random") == -1 and flex is None:
        print(flex_type, "flex is none!")

    # pad with -1 so we can store term_lens in a tensor
    seq_lens = torch.tensor(seq_lens)
    max_all_term_lens = max([len(term) for term in term_lens])
    for i, _ in enumerate(term_lens):
        term_lens[i] += [-1] * (max_all_term_lens - len(term_lens[i]))
    term_lens = torch.tensor(term_lens)

    # generate chain_idx from chain_lens
    chain_idx = []
    for c_lens in chain_lens:
        arrs = []
        for i, chain_len in enumerate(c_lens):
            arrs.append(torch.ones(chain_len) * i)
        chain_idx.append(torch.cat(arrs, dim=-1))
    chain_idx = pad_sequence(chain_idx, batch_first=True)

    return {
        'msas': msas,
        'features': features.float(),
        'ppoe': ppoe.float(),
        'seq_lens': seq_lens,
        'focuses': focuses,
        'contact_idxs': contact_idxs,
        'src_key_mask': src_key_mask,
        'term_lens': term_lens,
        'X': X,
        'x_mask': x_mask,
        'seqs': seqs,
        'ids': ids,
        'chain_idx': chain_idx,
        'gvp_data': gvp_data,
        'sortcery_seqs': sortcery_seqs,
        'sortcery_nrgs': sortcery_nrgs,
        'sscore': sscore_vals,
        # 'geometric_data': geometric_data,
        'flex': flex,
        'flex_mask': flex_mask
    }

### Source for dihedral calculations: https://www.cgl.ucsf.edu/Outreach/pc204/lecture_notes/phipsi/structured/phipsi.py"""

def Vector(x, y, z):
        return (x, y, z)

def subtract(u, v):
    "Return difference between two vectors."
    x = u[0] - v[0]
    y = u[1] - v[1]
    z = u[2] - v[2]
    return Vector(x, y, z)

def length(v):
    "Return length of a vector."
    sum = 0.0
    for c in v:
        sum += c * c
    return math.sqrt(sum)

def cross(u, v):
    "Return the cross product of two vectors."
    x = u[1] * v[2] - u[2] * v[1]
    y = u[2] * v[0] - u[0] * v[2]
    z = u[0] * v[1] - u[1] * v[0]
    return Vector(x, y, z)

def angle(v0, v1):
    "Return angle [0..pi] between two vectors."
    cosa = dot(v0, v1) / length(v0) / length(v1)
    cosa = round(cosa, 10)
    if cosa > 1:
        cosa = 1
    if cosa < -1:
        cosa = -1
    return math.acos(cosa)

def dot(u, v):
    "Return dot product of two vectors."
    sum = 0.0
    for cu, cv in zip(u, v):
        sum += cu * cv
    return sum

def calc_dihedral(p0, p1, p2, p3):
    """Return angle [0..2*pi] formed by vertices p0-p1-p2-p3."""

    v01 = subtract(p0, p1)
    v32 = subtract(p3, p2)
    v12 = subtract(p1, p2)
    v0 = cross(v12, v01)
    v3 = cross(v12, v32)
    # The cross product vectors are both normal to the axis
    # vector v12, so the angle between them is the dihedral
    # angle that we are looking for.  However, since "angle"
    # only returns values between 0 and pi, we need to make
    # sure we get the right sign relative to the rotation axis
    a = angle(v0, v3)
    if dot(cross(v0, v3), v12) > 0:
        a = -a
    return a

def update_pos(p0, p1, p2, dihedral, bond_length):
    """Return vertex p0 consistent with old p0, current p1, p2, p3, dihedral angle, and bond length."""
    
    # define vectors and angle to use in updating pos
    v01 = subtract(p0, p1)
    v10 = subtract(p1, p0)
    v12 = subtract(p1, p2)
    v0 = cross(v12, v01)
    t1 = cross(v12, v0)
    t2 = cross(v12, t1)
    a1 = angle(v01, v12)

    # define translation and rotation matrices
    v12_norm = np.divide(v12, length(v12))
    t1_norm = np.divide(t1, length(t1))
    t2_norm = np.divide(t2, length(t2))
    matrix_transform = np.transpose(np.array([v12_norm, t1_norm, t2_norm]))
    matrix_transform_inv = np.linalg.inv(matrix_transform)
    R = [[np.cos(dihedral), -1*np.sin(dihedral)], [np.sin(dihedral), np.cos(dihedral)]]

    # transform to local reference frame and update position, then transform back
    p0_local = np.matmul(matrix_transform_inv, p0)
    p1_local = np.matmul(matrix_transform_inv, p1)
    center_offset = np.cos(a1)*bond_length
    center_local = p1_local + [center_offset, 0, 0]
    p0_local_centered = p0_local - center_local
    p0_local_centered[1] = p0_local_centered[1]
    p0_local_centered[1:] = np.matmul(R, p0_local_centered[1:])
    p0_local_updated = p0_local_centered + center_local
    p0_updated = np.matmul(matrix_transform, p0_local_updated)
    return p0_updated

def calc_dihedrals_bond_lengths(X=None,chain_lens=None, expected_bond_lengths=None, expected_dihedrals=None):
    """Calculate noise

    Args
    ----
    X : torch.tensor
        Position of backbone residues
        size: num_residues x 4 x 3
    chain_lens : list
        list of chain lengths

    Returns
    -------
    dihedrals : pairs of (phi, psi) for each residue
    bond_lengths : tuples of bond length distances for each residue
    """

    X = X.reshape((X.shape[0]*X.shape[1], X.shape[2]))
    start_idx = 0
    dihedrals = []
    bond_lengths = []
    for c in chain_lens:
        c *= 4
        cur_dihedrals = [np.nan, np.nan]
        cur_bond_lengths = [np.nan, np.nan, np.nan]
        for atom_id in range(start_idx, start_idx + c - 4):
            p0 = X[atom_id,:]
            p1 = X[atom_id+1,:]
            p2 = X[atom_id+2,:]
            p3 = X[atom_id+3,:]
            if atom_id % 4 == 0:
                p3 = X[atom_id+4,:]
            elif atom_id % 4 == 1:
                p2 = X[atom_id+3,:]
                p3 = X[atom_id+4,:]
            elif atom_id % 4 == 2:
                p1 = X[atom_id+2,:]
                p2 = X[atom_id+3,:]
                p3 = X[atom_id+4,:]
            if atom_id % 4 == 3:
                continue
            dihedral = calc_dihedral(p0, p1, p2, p3)
            bond_length = length(subtract(p0,p1))
            if atom_id % 4 == 0:
                cur_dihedrals[1] = dihedral
                cur_bond_lengths[0] = bond_length
            elif atom_id % 4 == 1:
                cur_bond_lengths[1] = bond_length
            elif atom_id % 4 == 2:
                cur_dihedrals[0] = dihedral
                cur_bond_lengths[2] = bond_length
            if not np.isnan(cur_dihedrals[0]) and not np.isnan(cur_dihedrals[1]):
                dihedrals.append(cur_dihedrals)
                cur_dihedrals = [np.nan, np.nan]
            if not np.isnan(cur_bond_lengths[0]) and not np.isnan(cur_bond_lengths[1]) and not np.isnan(cur_bond_lengths[2]):
                bond_lengths.append(cur_bond_lengths)
                cur_bond_lengths = [np.nan, np.nan, np.nan]
        start_idx += c
    dihedrals = np.swapaxes(np.array(dihedrals), 0, 1)
    bond_lengths = np.swapaxes(np.array(bond_lengths), 0, 1)
    # if expected_bond_lengths is not None and expected_dihedrals is not None:
    #     for atom_idx in range(bond_lengths.shape[1]):
    #         print(f"expected dihedrals: {expected_dihedrals[:, atom_idx]}")
    #         print(f"dihedrals: {dihedrals[:, atom_idx]}")
    #         print(f"dihedrals diff: {dihedrals[:, atom_idx] - expected_dihedrals[:, atom_idx]}")
    return dihedrals.round(5), bond_lengths.round(5)

def generate_noise(flex_type, noise_level, size, X=None, bond_length_noise_level=0.0, chain_lens=None, expected_dihedrals=None):
    """Calculate noise

    Args
    ----
    flex_type : str
        methodology to calculate flex data
    noise_level : float
        std of noise to add
    size : tuple
        shape of coordinate entry in TERM data
    X : torch.tensor
        Position of backbone residues
        size: num_residues x 4 x 3
    bond_length_noise_level : float
        std of noise to add to bond lengths
    chain_lens : list
        list of chain lengths

    Returns
    -------
    noise : noise to be added to backbone atoms
    """
    print(f"noise level: {noise_level}, size={size}")
    if flex_type.find("fixed") == -1 and flex_type.find("torsion") == -1:
        noise = np.random.normal(loc=0, scale=noise_level, size=size)
    elif flex_type.find("fixed") > -1:
        flex_size = (size[0], size[2])
        noise = np.random.normal(loc=0, scale=noise_level, size=flex_size)
        noise = np.repeat(noise, 4)
        noise = np.reshape(noise, size)
    elif flex_type.find("torsion") > -1:
        X = X.reshape((X.shape[0]*X.shape[1], X.shape[2]))
        noise = np.zeros(X.shape)
        start_idx = 0
        for c in chain_lens:
            c *= 4
            noise[start_idx,:] = [0, 0, 0]# np.random.normal(loc=0, scale=0.02, size=(1,3))
            # X[start_idx, :] += noise[start_idx, :]
            for atom_id in range(start_idx + 1, start_idx + c - 4):
                p0 = X[atom_id,:]
                p1 = X[atom_id+1,:]
                p2 = X[atom_id+2,:]
                p3 = X[atom_id+3,:]
                if atom_id % 4 == 0:
                    p3 = X[atom_id+4,:]
                elif atom_id % 4 == 1:
                    p2 = X[atom_id+3,:]
                    p3 = X[atom_id+4,:]
                    continue
                elif atom_id % 4 == 2:
                    p1 = X[atom_id+2,:]
                    p2 = X[atom_id+3,:]
                    p3 = X[atom_id+4,:]
                if atom_id % 4 == 3:
                    continue
                dihedral_update = np.random.normal(loc=0, scale=noise_level, size=1)[0]
                if atom_id % 4 == 0 or atom_id % 4 == 2:
                    dihedral_update = 0.1 #dihedral_update #0.1
                else:
                    dihedral_update = 0
                bond_length = length(subtract(p0, p1))
                new_pos = update_pos(p0, p1, p2, dihedral_update, bond_length)
                new_v01 = subtract(new_pos, p1)
                bond_length = length(new_v01)
                bond_length_updated = bond_length + np.random.normal(loc=0, scale=bond_length_noise_level, size=1)[0]
                bond_length_multiplier = bond_length_updated / bond_length
                new_v01 = np.multiply(new_v01, bond_length_multiplier)
                new_pos = np.add(new_v01, p1)
                noise[atom_id, :] = new_pos - X[atom_id, :]
                for prev_atom_id in range(start_idx, atom_id):
                    prev_p0 = X[prev_atom_id, :] + noise[prev_atom_id, :]
                    prev_bond_length = length(subtract(prev_p0, p1))
                    prev_new_pos = update_pos(prev_p0, p1, p2, dihedral_update, prev_bond_length)
                    noise[prev_atom_id, :] = (prev_new_pos - X[prev_atom_id, :])
            noise[start_idx + c - 3,:] = [0, 0, 0] #np.random.normal(loc=0, scale=0.02, size=(1,3))
            # X[start_idx + c - 3,:] += noise[start_idx + c - 3,:]
            noise[start_idx + c - 2,:] = [0, 0, 0]# np.random.normal(loc=0, scale=0.02, size=(1,3))
            # X[start_idx + c - 1,:] += noise[start_idx + c - 2,:]
            noise[start_idx + c - 1,:] = [0, 0, 0]# np.random.normal(loc=0, scale=0.02, size=(1,3))
            # X[start_idx + c - 2,:] += noise[start_idx + c - 1,:]
            start_idx += c
        X = X.reshape((int(X.shape[0]/4), 4, 3))
        noise = noise.reshape(X.shape)

    return noise


def get_flex(pdb_id, flex_type, size, nrgten_folder=None, noise_level=0.0, bond_length_noise_level=0.0):
    """Calculate and save or load flex data to add for each pdb input

    Args
    ----
    pdb_id : str
        PDB ID to load.
    flex_type : str
        methodology to calculate flex data
    size : tuple
        shape of coordinate entry in TERM data
    nrgten_folder : str
        folder to find nrgten files
    noise_level : float
        std of noise to add
    bond_length_noise_level : float
        std of noise to add to bond lengths 
    Returns
    -------
    flex : list
        flex or dynamical signature for all backbone atoms 
    """
    if flex_type.find("random") > -1 and flex_type.find("batch") == -1:
        flex = generate_noise(flex_type, noise_level, size, bond_length_noise_level=bond_length_noise_level)
    elif flex_type == "dynamic_signature":
        try:
            nrgten_path = f"{nrgten_folder}/{pdb_id}_dynamic_signatures.txt"
            with open(nrgten_path, 'r') as f:
                flex_data = f.readlines()
            flex = np.array([int(sig[:-1]) for sig in flex_data])
            flex = flex[:, None]
        except:
            print("Must run nrgten on " + pdb_id + " before training.")
            return None
    elif flex_type.find("conformational_ensembles") > -1 or flex_type.find("conformational_pos_variation") > -1:
        try:
            nrgten_path = f"{nrgten_folder}/{pdb_id}_conformational_ensembles.features"
            with open(nrgten_path, 'rb') as fp:
                data = pickle.load(fp)
                flex = data['coords']
            if flex_type.find("conformational_pos_variation") > -1:
                flex = np.std(flex, axis=-1)
                flex = np.mean(flex, axis=(1, 2))
                flex = flex.flatten()
                flex = flex[:, None]
        except:
            print("Must run nrgten on " + pdb_id + " before training.")
            return None
    else:
        raise Exception(f"flex type '{flex_type}' not recognized")
    return flex

# Non-lazy data loading functions

def load_file(in_folder, pdb_folder, flex_folder, pdb_id, flex_type="", noise_level=0.0, bond_length_noise_level=0, min_protein_len=30, max_protein_len=0, num_ensembles=1):
    """Load the data specified in the proper .features file and return them.
    If the read sequence length is less than :code:`min_protein_len`, instead return None.

    Args
    ----
    in_folder : str
        folder to find TERM file.
    pdb_folder : str
        folder to find raw pdb files.
    flex_folder : str
        folder to find flex files.
    pdb_id : str
        PDB ID to load.
    flex_type : str
        methodology to calculate flex data
    noise_level : float
        std of noise to add
    bond_length_noise_level : float
        std of noise to add to bond lengths 
    min_protein_len : int
        minimum cutoff for loading TERM file.
    max_protein_len : int
        maximumum cutoff for loading protein to fit on GPU.
    num_ensembles : int
        number of conformational ensembles for each protein.

    Returns
    -------
    data : dict
        Data from TERM file (as dict)
    total_term_len : int
        Sum of lengths of all TERMs
    seq_len : int
        Length of protein sequence
    """
    path = f"{in_folder}/{pdb_id}/{pdb_id}.features"
    with open(path, 'rb') as fp:
        data = pickle.load(fp)
        if flex_type:
            if not pdb_folder and flex_type == "dynamic_signature":
                raise Exception("pdb folder must be input if using dynamic signature flex input")
            flex = get_flex(pdb_id, flex_type, data['coords'].shape, flex_folder, noise_level, bond_length_noise_level=bond_length_noise_level)
            if flex is None:
                return None
            if flex_type.find("random") > -1 and flex_type.find("batch") == -1:
                data['coords'] = np.add(data['coords'], flex)
        else:
            flex = None
        seq_len = data['seq_len']
        total_term_length = data['term_lens'].sum()
        if seq_len < min_protein_len or (max_protein_len and seq_len * num_ensembles > max_protein_len):
            print(f"{pdb_id} has length problems: seq_len {seq_len}, min_protein_len: {min_protein_len}, max_protein_len: {max_protein_len}, num_ensembes: {num_ensembles}!")
            return None
    return data, total_term_length, seq_len, flex


class TERMDataset(Dataset):
    """TERM Dataset that loads all feature files into a Pytorch Dataset-like structure.

    Attributes
    ----
    dataset : list
        list of tuples containing features, TERM length, and sequence length
    shuffle_idx : list
        array of indices for the dataset, for shuffling
    """
    def __init__(self, in_folder, pdb_folder, flex_folder=None, flex_type="", noise_level=0.0, bond_length_noise_level=0.0, pdb_ids=None, min_protein_len=30, max_protein_len=None, num_ensembles=1, num_processes=32):
        """
        Initializes current TERM dataset by reading in feature files.

        Reads in all feature files from the given directory, using multiprocessing
        with the provided number of processes. Stores the features, the TERM length,
        and the sequence length as a tuple representing the data. Can read from PDB ids or
        file paths directly. Uses the given protein length as a cutoff.

        Args
        ----
        in_folder : str
            path to directory containing feature files generated by :code:`scripts/data/preprocessing/generateDataset.py`
        pdb_folder : str
            path to directory containing pdb files
        flex_folder : str
            path to directory containing flex files
        flex_type : str
            methodology to calculate flex data
        noise_level : float
            std of noise to add
        bond_length_noise_level : float
            std of noise to add to bond lengths 
        pdb_ids: list, optional
            list of pdbs from `in_folder` to include in the dataset
        min_protein_len: int, default=30
            minimum length of a protein in the dataset
        max_protein_len : int
            maximumum cutoff for loading protein to fit on GPU.
        num_ensembles : int
            number of conformational ensembles for each protein.
        num_processes: int, default=32
            number of processes to use during dataloading
        """
        self.dataset = []

        with mp.Pool(num_processes) as pool:

            if pdb_ids:      
                print("Loading feature files from pdb_ids subset list")
                progress = tqdm(total=len(pdb_ids))
            else:
                print("Loading feature file paths from features folder")

                filelist = list(glob.glob(f'{in_folder}/*/*.features'))
                print(filelist)
                progress = tqdm(total=len(filelist))
                # get pdb_ids
                pdb_ids = [os.path.basename(path)[:-len(".features")] for path in filelist]
                print(pdb_ids)

            def update_progress(res):
                del res
                if torch.distributed.is_initialized():
                    if torch.distributed.get_rank() == 0:
                        progress.update(1)
                else:
                    progress.update(1)

            res_list = [
                pool.apply_async(load_file, (in_folder, pdb_folder, flex_folder, id, flex_type, noise_level, bond_length_noise_level),
                                    kwds={"min_protein_len": min_protein_len, "max_protein_len": max_protein_len, "num_ensembles": num_ensembles},
                                    callback=update_progress) for id in pdb_ids
            ]
            pool.close()
            pool.join()
            progress.close()
            for i, res in enumerate(res_list):
                data = res.get()
                if data is not None:
                    features, total_term_length, seq_len, flex = data
                    if flex is not None and flex.shape[0] != features['coords'].shape[0]:
                        print(f"Skipping {pdb_ids[i]} because sequence lengths of flex and coord input are not consistent ({flex.shape}, {features['coords'].shape})")
                        continue
                    features['flex'] = flex
                    self.dataset.append((features, total_term_length, seq_len))
                else:
                    print(f"Error processing {pdb_ids[i]}. Skipping file.")
            self.shuffle_idx = np.arange(len(self.dataset))
        print(f"Done with init, have {len(self.dataset)} examples to use.")

    def shuffle(self):
        """Shuffle the current dataset."""
        np.random.shuffle(self.shuffle_idx)

    def __len__(self):
        """Returns length of the given dataset.

        Returns
        -------
        int
            length of dataset
        """
        return len(self.dataset)

    def __getitem__(self, idx):
        """Extract a given item with provided index.

        Args
        ----
        idx : int
            Index of item to return.
        Returns
        ----
        data : dict
            Data from TERM file (as dict)
        total_term_len : int
            Sum of lengths of all TERMs
        seq_len : int
            Length of protein sequence
        """
        data_idx = self.shuffle_idx[idx]
        if isinstance(data_idx, list):
            return np.array([self.dataset[i] for i in data_idx])
        return self.dataset[data_idx]              

def TERMBatchSampler_init(self, ddp,
                dataset,
                dev,
                batch_size=4,
                sort_data=False,
                shuffle=True,
                semi_shuffle=False,
                semi_shuffle_cluster_size=500,
                batch_shuffle=True,
                drop_last=False,
                max_term_res=55000,
                max_seq_tokens=None,
                flex_type=None,
                noise_level=0.0,
                bond_length_noise_level=0.0,
                num_ensembles=1):
    """
    Reads in and processes a given dataset.

    Given the provided dataset, load all the data. Then cluster the data using
    the provided method, either shuffled or sorted and then shuffled.

    Args
    ----
    ddp : bool
        Indicator for usage of distributed data parallel in creating TERMBatchSampler.
    dataset : TERMDataset
        Dataset to batch.
    dev : str
        Computation device to use
    batch_size : int or None, default=4
        Size of batches created. If variable sized batches are desired, set to None.
    sort_data : bool, default=False
        Create deterministic batches by sorting the data according to the
        specified length metric and creating batches from the sorted data.
        Incompatible with :code:`shuffle=True` and :code:`semi_shuffle=True`.
    shuffle : bool, default=True
        Shuffle the data completely before creating batches.
        Incompatible with :code:`sort_data=True` and :code:`semi_shuffle=True`.
    semi_shuffle : bool, default=False
        Sort the data according to the specified length metric,
        then partition the data into :code:`semi_shuffle_cluster_size`-sized partitions.
        Within each partition perform a complete shuffle. The upside is that
        batching with similar lengths reduces padding making for more efficient computation,
        but the downside is that it does a less complete shuffle.
    semi_shuffle_cluster_size : int, default=500
        Size of partition to use when :code:`semi_shuffle=True`.
    batch_shuffle : bool, default=True
        If set to :code:`True`, shuffle samples within a batch.
    drop_last : bool, default=False
        If set to :code:`True`, drop the last samples if they don't form a complete batch.
    max_term_res : int or None, default=55000
        When :code:`batch_size=None, max_term_res>0, max_seq_tokens=None`,
        batch by fitting as many datapoints as possible with the total number of
        TERM residues included below `max_term_res`.
        Calibrated using :code:`nn.DataParallel` on two V100 GPUs.
    max_seq_tokens : int or None, default=None
        When :code:`batch_size=None, max_term_res=None, max_seq_tokens>0`,
        batch by fitting as many datapoints as possible with the total number of
        sequence residues included below `max_seq_tokens`. Exactly one of :code:`max_term_res`
        and :code:`max_seq_tokens` must be None.
    flex_type : str
        methodology to calculate flex data
    noise_level : float
        std of noise to add
    bond_length_noise_level : float
        std of noise to add to bond lengths
    num_ensembles : int
        number of conformational ensembles for each training example
    """
    if ddp:
        DistributedSampler.__init__(self, dataset)
    else:
        Sampler.__init__(self, dataset)
    self.size = len(dataset)
    self.dataset, self.total_term_lengths, self.seq_lengths = zip(*dataset)
    assert not (max_term_res is None
                and max_seq_tokens is None), "Exactly one of max_term_res and max_seq_tokens must be None"
    if max_term_res is None and max_seq_tokens > 0:
        self.lengths = self.seq_lengths
    elif max_term_res > 0 and max_seq_tokens is None:
        self.lengths = self.total_term_lengths
    else:
        raise ValueError("Exactly one of max_term_res and max_seq_tokens must be None")
    self.shuffle = shuffle
    self.sort_data = sort_data
    self.batch_shuffle = batch_shuffle
    self.batch_size = batch_size
    self.drop_last = drop_last
    self.max_term_res = max_term_res
    self.max_seq_tokens = max_seq_tokens
    self.semi_shuffle = semi_shuffle
    self.semi_shuffle_cluster_size = semi_shuffle_cluster_size
    self.ddp = ddp
    self.dev = dev
    self.flex_type = flex_type
    self.noise_level = noise_level
    self.bond_length_noise_level = bond_length_noise_level
    self.num_ensembles = num_ensembles

    assert not (shuffle and semi_shuffle), "Lazy Dataloader shuffle and semi shuffle cannot both be set"
    # initialize clusters
    self._cluster()
    print(f"Done with clustering. Have {len(self.clusters)} clusters.")


def TERMBatchSampler_cluster(self):
    """ Shuffle data and make clusters of indices corresponding to batches of data.

    This method speeds up training by sorting data points with similar TERM lengths
    together, if :code:`sort_data` or :code:`semi_shuffle` are on. Under `sort_data`,
    the data is sorted by length. Under `semi_shuffle`, the data is broken up
    into clusters based on length and shuffled within the clusters. Otherwise,
    it is randomly shuffled. Data is then loaded into batches based on the number
    of proteins that will fit into the GPU without overloading it, based on
    :code:`max_term_res` or :code:`max_seq_tokens`.
    """
    if self.ddp and dist.get_rank() > 0:
        dims = [0, 0]
        dist.broadcast_object_list(dims, src=0, device=torch.device(self.dev))
        clusters = [list(-1*np.ones(dims[1], dtype=int))]*dims[0]
        # clusters = [torch.tensor(range(len(self.dataset)+1)).to(self.dev)]*len(self.dataset)
        dist.broadcast_object_list(clusters, src=0, device=torch.device(self.dev))
        # while True:
        #     batch = torch.tensor(range(len(self.dataset)), dtype=torch.uint8).to(self.dev)
        #     dist.recv(batch, src=0)
        #     batch = list((batch.cpu()).numpy())
        #     if not batch:
        #         break
        #     clusters.append(batch)
        final_clusters = []
        for batch in clusters:
            if batch[0] > 0:
                final_clusters.append(batch[1:batch[0]+1])
        self.clusters = final_clusters
    else:
        # if we sort data, use sorted indexes instead
        if self.sort_data:
            idx_list = np.argsort(self.lengths)
        elif self.semi_shuffle:
            # trying to speed up training
            # by shuffling points with similar term res together
            idx_list = np.argsort(self.lengths)
            shuffle_borders = []

            # break up datapoints into large clusters
            border = 0
            while border < len(self.lengths):
                shuffle_borders.append(border)
                border += self.semi_shuffle_cluster_size

            # shuffle datapoints within clusters
            last_cluster_idx = len(shuffle_borders) - 1
            for cluster_idx in range(last_cluster_idx + 1):
                start = shuffle_borders[cluster_idx]
                if cluster_idx < last_cluster_idx:
                    end = shuffle_borders[cluster_idx + 1]
                    np.random.shuffle(idx_list[start:end])
                else:
                    np.random.shuffle(idx_list[start:])

        else:
            idx_list = list(range(len(self.dataset)))
            np.random.shuffle(idx_list)

        # Cluster into batches of similar sizes
        clusters, batch = [], []

        # if batch_size is None, fit as many proteins we can into a batch
        # without overloading the GPU
        if self.batch_size is None:
            if self.max_term_res is None and self.max_seq_tokens > 0:
                cap_len = self.max_seq_tokens
            elif self.max_term_res > 0 and self.max_seq_tokens is None:
                cap_len = self.max_term_res
            cap_len /= self.num_ensembles

            current_batch_lens = []
            total_data_len = 0
            for count, idx in enumerate(idx_list):
                current_batch_lens.append(self.lengths[idx])
                total_data_len = max(current_batch_lens) * len(current_batch_lens) * self.num_ensembles
                if count != 0 and total_data_len > cap_len:
                    clusters.append(batch)
                    batch = [idx]
                    current_batch_lens = [self.lengths[idx]]
                else:
                    batch.append(idx)

        else:  # used fixed batch size
            for count, idx in enumerate(idx_list):
                if count != 0 and count % self.batch_size == 0:
                    clusters.append(batch)
                    batch = [idx]
                else:
                    batch.append(idx)

        if len(batch) > 0 and not self.drop_last:
            clusters.append(batch)
        self.clusters = clusters
        if self.ddp:
            send_clusters = []
            max_len = 0
            for batch in clusters:
                batch = [len(batch)] + batch
                send_clusters.append(batch)
                max_len = max(max_len, len(batch))
            for i, batch in enumerate(send_clusters):
                send_clusters[i] = batch + list(-1*np.ones(max_len - len(batch), dtype=int))
            dims = [len(send_clusters), len(send_clusters[0])]
            dist.broadcast_object_list(dims, src=0, device=torch.device(self.dev))
            dist.broadcast_object_list(send_clusters, src=0, device=torch.device(self.dev))
            # _send_batch(send_clusters, self.num_replicas, self.dev)

def TERMBatchSampler_package(self, b_idx):
    """Package the given datapoints into tensors based on provided indices.

    Tensors are extracted from the data and padded. Coordinates are featurized
    and the length of TERMs and chain IDs are added to the data.

    Args
    ----
    b_idx : list of tuples (dicts, int, int)
        The feature dictionaries, the sum of the lengths of all TERMs, and the sum of all sequence lengths
        for each datapoint to package.

    Returns
    -------
    dict
        Collection of batched features required for running TERMinator. This contains:

        - :code:`msas` - the multiple sequence alignment

        - :code:`features` - the TERM features

        - :code:`ppoe` - the :math:`\\phi, \\psi, \\omega`, and environment values of the target structure

        - :code:`seq_lens` - the lengths of the sequences

        - :code:`focuses` - the corresponding target structure residue index for each TERM residue

        - :code:`contact_idxs` - contact indices

        - :code:`src_key_mask` - mask for TERM residue padding

        - :code:`X` - coordinates

        - :code:`x_mask` - mask for the target structure

        - :code:`seqs` - the sequences

        - :code:`ids` - the PDB ids

        - :code:`chain_idx` - the chain IDs
    """
    return _package([b[0:2] for b in b_idx], self.flex_type, self.noise_level, self.bond_length_noise_level)


def TERMBatchSampler_len(self):
    """Returns length of dataset, i.e. number of batches.

    Returns
    -------
    int
        length of dataset.
    """
    if not self.ddp:
        return len(self.clusters)
    if len(self.clusters) % self.num_replicas != 0:
        if not self.drop_last:
            return math.floor((len(self.clusters) + self.num_replicas) / self.num_replicas)
        else:
            return math.floor((len(self.clusters)) / self.num_replicas)
    return math.ceil(len(self.clusters) / self.num_replicas)


def TERMBatchSampler_iter(self):
    """Allows iteration over dataset."""
    if not self.ddp:
        if self.shuffle or self.semi_shuffle:
            self._cluster()
            np.random.shuffle(self.clusters)
        for batch in self.clusters:
            yield batch
    else:
        if self.shuffle or self.semi_shuffle:
            self._cluster()
            g = torch.Generator()
            g.manual_seed(self.seed + self.epoch)
            indices = torch.randperm(len(self.clusters), generator=g).tolist()
        else:
            indices = list(range(len(self.clusters)))
        
        self.num_samples = self.__len__()
        self.total_size = self.num_samples * self.num_replicas
        if not self.drop_last:
            padding_size = self.total_size - len(indices)
            if padding_size <= len(indices):
                indices += indices[:padding_size]
            else:
                indices += (indices * math.ceil(padding_size / len(indices)))[:padding_size]
        else:
            indices = indices[:self.total_size]
        assert len(indices) == self.total_size
        indices = indices[self.rank:self.total_size:self.num_replicas]
        assert len(indices) == self.num_samples
        for index in indices:
            yield self.clusters[index]
        

class TERMBatchSamplerWrapper(object):
    """BatchSampler/Dataloader helper class for TERM data using TERMDataset.

    Attributes in wrapped TERMBatchSampler class
    ----
    size: int
        Length of the dataset
    dataset: List
        List of features from TERM dataset
    total_term_lengths: List
        List of TERM lengths from the given dataset
    seq_lengths: List
        List of sequence lengths from the given dataset
    lengths: List
        TERM lengths or sequence lengths, depending on
        whether :code:`max_term_res` or :code:`max_seq_tokens` is set.
    batch_size : int or None, default=4
        Size of batches created. If variable sized batches are desired, set to None.
    sort_data : bool, default=False
        Create deterministic batches by sorting the data according to the
        specified length metric and creating batches from the sorted data.
        Incompatible with :code:`shuffle=True` and :code:`semi_shuffle=True`.
    shuffle : bool, default=True
        Shuffle the data completely before creating batches.
        Incompatible with :code:`sort_data=True` and :code:`semi_shuffle=True`.
    semi_shuffle : bool, default=False
        Sort the data according to the specified length metric,
        then partition the data into :code:`semi_shuffle_cluster_size`-sized partitions.
        Within each partition perform a complete shuffle. The upside is that
        batching with similar lengths reduces padding making for more efficient computation,
        but the downside is that it does a less complete shuffle.
    semi_shuffle_cluster_size : int, default=500
        Size of partition to use when :code:`semi_shuffle=True`.
    batch_shuffle : bool, default=True
        If set to :code:`True`, shuffle samples within a batch.
    drop_last : bool, default=False
        If set to :code:`True`, drop the last samples if they don't form a complete batch.
    max_term_res : int or None, default=55000
        When :code:`batch_size=None, max_term_res>0, max_seq_tokens=None`,
        batch by fitting as many datapoints as possible with the total number of
        TERM residues included below `max_term_res`.
        Calibrated using :code:`nn.DataParallel` on two V100 GPUs.
    max_seq_tokens : int or None, default=None
        When :code:`batch_size=None, max_term_res=None, max_seq_tokens>0`,
        batch by fitting as many datapoints as possible with the total number of
        sequence residues included below `max_seq_tokens`.
    """
    def __init__(self,
                 ddp):
        self.ddp = ddp
        if ddp:
            class TERMBatchSampler(DistributedSampler):
                pass
        else:
            class TERMBatchSampler(Sampler):
                pass
        TERMBatchSampler.__init__ = TERMBatchSampler_init
        TERMBatchSampler._cluster = TERMBatchSampler_cluster
        TERMBatchSampler.package = TERMBatchSampler_package
        TERMBatchSampler.__len__ = TERMBatchSampler_len
        TERMBatchSampler.__iter__ = TERMBatchSampler_iter

        self.sampler = TERMBatchSampler

# needs to be outside of object for pickling reasons (?)
def read_lens(in_folder, pdb_id, min_protein_len=30):
    """ Reads the lengths specified in the proper .length file and return them.

    If the read sequence length is less than :code:`min_protein_len`, instead return None.

    Args
    ----
    in_folder : str
        folder to find TERM file.
    pdb_id : str
        PDB ID to load.
    min_protein_len : int
        minimum cutoff for loading TERM file.
    Returns
    -------
    pdb_id : str
        PDB ID that was loaded
    total_term_length : int
        number of TERMS in file
    seq_len : int
        sequence length of file, or None if sequence length is less than :code:`min_protein_len`
    """
    path = f"{in_folder}/{pdb_id}/{pdb_id}.length"
    # pylint: disable=unspecified-encoding
    with open(path, 'rt') as fp:
        total_term_length = int(fp.readline().strip())
        seq_len = int(fp.readline().strip())
        if seq_len < min_protein_len:
            return None
    return pdb_id, total_term_length, seq_len


class TERMLazyDataset(Dataset):
    """TERM Dataset that loads all feature files into a Pytorch Dataset-like structure.

    Unlike TERMDataset, this just loads feature filenames, not actual features.

    Attributes
    ----
    dataset : list
        list of tuples containing feature filenames, TERM length, and sequence length
    shuffle_idx : list
        array of indices for the dataset, for shuffling
    """
    def __init__(self, in_folder, pdb_ids=None, min_protein_len=30, num_processes=32):
        """
        Initializes current TERM dataset by reading in feature files.

        Reads in all feature files from the given directory, using multiprocessing
        with the provided number of processes. Stores the feature filenames, the TERM length,
        and the sequence length as a tuple representing the data. Can read from PDB ids or
        file paths directly. Uses the given protein length as a cutoff.

        Args
        ----
        in_folder : str
            path to directory containing feature files generated by :code:`scripts/data/preprocessing/generateDataset.py`
        pdb_ids: list, optional
            list of pdbs from `in_folder` to include in the dataset
        min_protein_len: int, default=30
            minimum length of a protein in the dataset
        num_processes: int, default=32
            number of processes to use during dataloading
        """
        self.dataset = []

        with mp.Pool(num_processes) as pool:

            if pdb_ids:
                print("Loading feature file paths")
                progress = tqdm(total=len(pdb_ids))

                def update_progress(res):
                    del res
                    if torch.distributed.is_initialized():
                        if torch.distributed.get_rank() == 0:
                            progress.update(1)
                    else:
                        progress.update(1)

                res_list = [
                    pool.apply_async(read_lens, (in_folder, pdb_id),
                                     kwds={"min_protein_len": min_protein_len},
                                     callback=update_progress) for pdb_id in pdb_ids
                ]
                pool.close()
                pool.join()
                progress.close()
                for res in res_list:
                    data = res.get()
                    if data is not None:
                        pdb_id, total_term_length, seq_len = data
                        filename = f"{in_folder}/{pdb_id}/{pdb_id}.features"
                        self.dataset.append((os.path.abspath(filename), total_term_length, seq_len))
            else:
                print("Loading feature file paths")

                filelist = list(glob.glob(f'{in_folder}/*/*.features'))
                progress = tqdm(total=len(filelist))

                def update_progress(res):
                    del res
                    if torch.distributed.is_initialized():
                        if torch.distributed.get_rank() == 0:
                            progress.update(1)
                    else:
                        progress.update(1)

                # get pdb_ids
                pdb_ids = [os.path.basename(path)[:-len(".features")] for path in filelist]

                res_list = [
                    pool.apply_async(read_lens, (in_folder, pdb_id),
                                     kwds={"min_protein_len": min_protein_len},
                                     callback=update_progress) for pdb_id in pdb_ids
                ]
                pool.close()
                pool.join()
                progress.close()
                for res in res_list:
                    data = res.get()
                    if data is not None:
                        pdb_id, total_term_length, seq_len = data
                        filename = f"{in_folder}/{pdb_id}/{pdb_id}.features"
                        self.dataset.append((os.path.abspath(filename), total_term_length, seq_len))

        self.shuffle_idx = np.arange(len(self.dataset))

    def shuffle(self):
        """Shuffle the dataset"""
        np.random.shuffle(self.shuffle_idx)

    def __len__(self):
        """Returns length of the given dataset.

        Returns
        -------
        int
            length of dataset
        """
        return len(self.dataset)

    def __getitem__(self, idx):
        """Extract a given item with provided index.

        Args
        ----
        idx : int
            Index of item to return.
        Returns
        ----
        data : dict
            Data from TERM file (as dict)
        total_term_len : int
            Sum of lengths of all TERMs
        seq_len : int
            Length of protein sequence
        """
        data_idx = self.shuffle_idx[idx]
        if isinstance(data_idx, list):
            return np.array([self.dataset[i] for i in data_idx])
        return self.dataset[data_idx]


def TERMLazyBatchSampler_init(self, ddp,
                dataset,
                batch_size=4,
                sort_data=False,
                shuffle=True,
                semi_shuffle=False,
                semi_shuffle_cluster_size=500,
                batch_shuffle=True,
                drop_last=False,
                max_term_res=55000,
                max_seq_tokens=None,
                term_matches_cutoff=None,
                term_dropout=None,
                noise=0.0):
    
    """
    Reads in and processes a given dataset.

    Given the provided dataset, load all the data. Then cluster the data using
    the provided method, either shuffled or sorted and then shuffled.

    Args
    ----
    ddp : bool
        Indicator for usage of distributed data parallel in creating TERMBatchSampler.
    dataset : TERMLazyDataset
        Dataset to batch.
    batch_size : int or None, default=4
        Size of batches created. If variable sized batches are desired, set to None.
    sort_data : bool, default=False
        Create deterministic batches by sorting the data according to the
        specified length metric and creating batches from the sorted data.
        Incompatible with :code:`shuffle=True` and :code:`semi_shuffle=True`.
    shuffle : bool, default=True
        Shuffle the data completely before creating batches.
        Incompatible with :code:`sort_data=True` and :code:`semi_shuffle=True`.
    semi_shuffle : bool, default=False
        Sort the data according to the specified length metric,
        then partition the data into :code:`semi_shuffle_cluster_size`-sized partitions.
        Within each partition perform a complete shuffle. The upside is that
        batching with similar lengths reduces padding making for more efficient computation,
        but the downside is that it does a less complete shuffle.
    semi_shuffle_cluster_size : int, default=500
        Size of partition to use when :code:`semi_shuffle=True`.
    batch_shuffle : bool, default=True
        If set to :code:`True`, shuffle samples within a batch.
    drop_last : bool, default=False
        If set to :code:`True`, drop the last samples if they don't form a complete batch.
    max_term_res : int or None, default=55000
        When :code:`batch_size=None, max_term_res>0, max_seq_tokens=None`,
        batch by fitting as many datapoints as possible with the total number of
        TERM residues included below `max_term_res`.
        Calibrated using :code:`nn.DataParallel` on two V100 GPUs.
    max_seq_tokens : int or None, default=None
        When :code:`batch_size=None, max_term_res=None, max_seq_tokens>0`,
        batch by fitting as many datapoints as possible with the total number of
        sequence residues included below `max_seq_tokens`.
    term_matches_cutoff : int or None, default=None
        Use the top :code:`term_matches_cutoff` TERM matches for featurization.
        If :code:`None`, apply no cutoff.
    term_dropout : str or None, default=None
        Let `t` be the number of TERM matches in the given datapoint.
        Select a random int `n` from 1 to `t`, and take a random subset `n`
        of the given TERM matches to keep. If :code:`term_dropout='keep_first'`,
        keep the first match and choose `n-1` from the rest.
        If :code:`term_dropout='all'`, choose `n` matches from all matches.
    noise : float
        std of noise to add to all backbone atoms.
    """
    if ddp:
        DistributedSampler.__init__(self, dataset)
    else:
        Sampler.__init__(self, dataset)
    self.dataset = dataset
    self.size = len(dataset)
    self.filepaths, self.total_term_lengths, self.seq_lengths = zip(*dataset)
    assert not (max_term_res is None
                and max_seq_tokens is None), "Exactly one of max_term_res and max_seq_tokens must be None"
    if max_term_res is None and max_seq_tokens > 0:
        self.lengths = self.seq_lengths
    elif max_term_res > 0 and max_seq_tokens is None:
        self.lengths = self.total_term_lengths
    else:
        raise Exception("Exactly one of max_term_res and max_seq_tokens must be None")
    self.shuffle = shuffle
    self.sort_data = sort_data
    self.batch_shuffle = batch_shuffle
    self.batch_size = batch_size
    self.drop_last = drop_last
    self.max_term_res = max_term_res
    self.max_seq_tokens = max_seq_tokens
    self.semi_shuffle = semi_shuffle
    self.semi_shuffle_cluster_size = semi_shuffle_cluster_size
    self.term_matches_cutoff = term_matches_cutoff
    assert term_dropout in ["keep_first", "all", None], f"term_dropout={term_dropout} is not a valid argument"
    self.term_dropout = term_dropout
    self.noise = noise
    assert not (shuffle and semi_shuffle), "Lazy Dataloader shuffle and semi shuffle cannot both be set"

    # initialize clusters
    self._cluster()

def TERMLazyBatchSampler_cluster(self):
    """ Shuffle data and make clusters of indices corresponding to batches of data.

    This method speeds up training by sorting data points with similar TERM lengths
    together, if :code:`sort_data` or :code:`semi_shuffle` are on. Under `sort_data`,
    the data is sorted by length. Under `semi_shuffle`, the data is broken up
    into clusters based on length and shuffled within the clusters. Otherwise,
    it is randomly shuffled. Data is then loaded into batches based on the number
    of proteins that will fit into the GPU without overloading it, based on
    :code:`max_term_res` or :code:`max_seq_tokens`.
    """

    # if we sort data, use sorted indexes instead
    if self.sort_data:
        idx_list = np.argsort(self.lengths)
    elif self.semi_shuffle:
        # trying to speed up training
        # by shuffling points with similar term res together
        idx_list = np.argsort(self.lengths)
        shuffle_borders = []

        # break up datapoints into large clusters
        border = 0
        while border < len(self.lengths):
            shuffle_borders.append(border)
            border += self.semi_shuffle_cluster_size

        # shuffle datapoints within clusters
        last_cluster_idx = len(shuffle_borders) - 1
        for cluster_idx in range(last_cluster_idx + 1):
            start = shuffle_borders[cluster_idx]
            if cluster_idx < last_cluster_idx:
                end = shuffle_borders[cluster_idx + 1]
                np.random.shuffle(idx_list[start:end])
            else:
                np.random.shuffle(idx_list[start:])

    else:
        idx_list = list(range(len(self.dataset)))
        np.random.shuffle(idx_list)

    # Cluster into batches of similar sizes
    clusters, batch = [], []

    # if batch_size is None, fit as many proteins we can into a batch
    # without overloading the GPU
    if self.batch_size is None:
        if self.max_term_res is None and self.max_seq_tokens > 0:
            cap_len = self.max_seq_tokens
        elif self.max_term_res > 0 and self.max_seq_tokens is None:
            cap_len = self.max_term_res

        current_batch_lens = []
        total_data_len = 0
        for count, idx in enumerate(idx_list):
            current_batch_lens.append(self.lengths[idx])
            total_data_len = max(current_batch_lens) * len(current_batch_lens)
            if count != 0 and total_data_len > cap_len:
                clusters.append(batch)
                batch = [idx]
                current_batch_lens = [self.lengths[idx]]
            else:
                batch.append(idx)

    else:  # used fixed batch size
        for count, idx in enumerate(idx_list):
            if count != 0 and count % self.batch_size == 0:
                clusters.append(batch)
                batch = [idx]
            else:
                batch.append(idx)

    if len(batch) > 0 and not self.drop_last:
        clusters.append(batch)
    self.clusters = clusters

def TERMLazyBatchSampler_package(self, b_idx):
    """Package the given datapoints into tensors based on provided indices.

    Tensors are extracted from the data and padded. Coordinates are featurized
    and the length of TERMs and chain IDs are added to the data.

    Args
    ----
    b_idx : list of (str, int, int)
        The path to the feature files, the sum of the lengths of all TERMs, and the sum of all sequence lengths
        for each datapoint to package.

    Returns
    -------
    dict
        Collection of batched features required for running TERMinator. This contains:

        - :code:`msas` - the multiple sequence alignment

        - :code:`features` - the TERM features

        - :code:`ppoe` - the :math:`\\phi, \\psi, \\omega`, and environment values of the target structure

        - :code:`seq_lens` - the lengths of the sequences

        - :code:`focuses` - the corresponding target structure residue index for each TERM residue

        - :code:`contact_idxs` - contact indices

        - :code:`src_key_mask` - mask for TERM residue padding

        - :code:`X` - coordinates

        - :code:`x_mask` - mask for the target structure

        - :code:`seqs` - the sequences

        - :code:`ids` - the PDB ids

        - :code:`chain_idx` - the chain IDs
    """
    if self.batch_shuffle:
        b_idx_copy = b_idx[:]
        random.shuffle(b_idx_copy)
        b_idx = b_idx_copy

    # load the files specified by filepaths
    batch = []
    for data in b_idx:
        filepath = data[0]
        with open(filepath, 'rb') as fp:
            protein_data = pickle.load(fp)
            protein_data['coords'] = protein_data['coords'] + np.random.normal(loc=0, scale=self.noise, size=protein_data['coords'].shape)
            batch.append((protein_data, data[1]))
            if 'ppoe' not in batch[-1][0].keys():
                print(filepath)
        fp.close()

    # package batch
    packaged_batch = _package(batch)

    features = packaged_batch["features"]
    msas = packaged_batch["msas"]
    # apply TERM matches cutoff
    if self.term_matches_cutoff:
        features = features[:, :self.term_matches_cutoff]
        msas = msas[:, :self.term_matches_cutoff]
    # apply TERM matches dropout
    if self.term_dropout:
        # sample a random number of alignments to keep
        n_batch, n_align, n_terms, n_features = features.shape
        if self.term_dropout == 'keep_first':
            n_keep = torch.randint(0, n_align, [1]).item()
        elif self.term_dropout == 'all':
            n_keep = torch.randint(1, n_align, [1]).item()
        # sample from a multinomial distribution
        weights = torch.ones([1, 1]).expand([n_batch * n_terms, n_keep])
        if n_keep == 0:
            sample_idx = torch.ones(1)
        else:
            sample_idx = torch.multinomial(weights, n_keep)
            sample_idx = sample_idx.view([n_batch, n_terms, n_keep]).transpose(-1, -2)
            sample_idx_features = sample_idx.unsqueeze(-1).expand([n_batch, n_keep, n_terms, n_features])
            sample_idx_msas = sample_idx

        if self.term_dropout == 'keep_first':
            if n_keep == 0:
                features = features[:, 0:1]
                msas = msas[:, 0:1]
            else:
                sample_features = torch.gather(features[:, 1:], 1, sample_idx_features)
                sample_msas = torch.gather(msas[:, 1:], 1, sample_idx_msas)
                features = torch.cat([features[:, 0:1], sample_features], dim=1)
                msas = torch.cat([msas[:, 0:1], sample_msas], dim=1)
        elif self.term_dropout == 'all':
            features = torch.gather(features, 1, sample_idx_features)
            msas = torch.gather(msas, 1, sample_idx_msas)

    packaged_batch["features"] = features
    packaged_batch["msas"] = msas

    return packaged_batch

def TERMLazyBatchSampler_len(self):
    """Returns length of dataset, i.e. number of batches.

    Returns
    -------
    int
        length of dataset.
    """
    if not self.ddp:
        return len(self.clusters)
    if len(self.clusters) % self.num_replicas != 0:
        if not self.drop_last:
            return math.floor((len(self.clusters) + self.num_replicas) / self.num_replicas)
        else:
            return math.floor((len(self.clusters)) / self.num_replicas)
    return math.ceil(len(self.clusters) / self.num_replicas)

def separateTargetAndBinder(target_data,binder_data):
    # Moves the binder 1000 angstroms away from the target and then combines their data
    # lazy,
    binder_data['coords'] = binder_data['coords'] + 1000
    return concatFeatureDicts(target_data,binder_data)

def separateTargetAndBinder(target_data,binder_data):
    # Moves the binder 1000 angstroms away from the target and then combines their data
    # lazy,
    binder_data['coords'] = binder_data['coords'] + 1000
    return concatFeatureDicts(target_data,binder_data)

class BinderScoringIterableDataset(IterableDataset):
    def __init__(self, filename, target_pdb_path='', max_res_num = 27500, binder_subset=None, mode="binder_and_complex", skip_package=False):
        # the max number of residues per batch (assuming single V100)
        self.max_res_num = max_res_num

        # the subset of binders that should be loaded
        self.binder_subset = binder_subset
        # print(f"BinderScoringIterableDataset subset length: {len(self.binder_subset)}")

        if mode == "binder_and_complex":
            self.load_complex = True
            self.load_binder = True
        elif mode == "complex_only":
            self.load_complex = True
            self.load_binder = False
        elif mode == "binder_only":
            self.load_complex = False
            self.load_binder = True
        else:
            raise ValueError(f"'{mode}' is not a valid modes")
        
        self.skip_package = skip_package
        
        #Store the filename
        self.filename = filename

        with open(self.filename,'r') as file:
            start_pos = file.tell()
            data = self.readPDBFromOpenFile(file,'',True)
            if data['pdb'][-7:] != '_target':
                if target_pdb_path != '':
                    print('No target structure found at the beginning of the file, instead will be reading from: ',target_pdb_path)
                    with open(target_pdb_path) as target_file:
                        target_name = target_pdb_path.split('/')[-1][:-len('.pdb')]
                        self.target_data = self.readPDBFromOpenFile(target_file,target_name,True)
                    self.target_packaged_data = self.prepareTargetFeatureDict()
                    self.binder_start_pos = start_pos # binders start at the beginning of the file
                else:
                    raise ValueError('The first PDB in the file should be a target protein structure or provided as a separate file')
            else:
                self.target_data = data
                self.target_packaged_data = self.prepareTargetFeatureDict()
                self.binder_start_pos = file.tell() # binders start after target

        self.current_complex_data = None
        self.current_complex_packaged_data = None

        self.current_binder_data = None
        self.current_binder_packaged_data = None

    def readPDBFromOpenFile(self,file,name='',target=False):
        """Read an individual structure from a multi-entry PDB file"""
        data = dict()

        # First line should be the name of the file
        line = file.readline()
        # If line is empty, we know we've reached the end of the file
        if not line:
            return None

        if name == '':
            if line[0:10] != 'HEADER    ':
                raise AssertionError('Each section of the multi-entry PDB file should begin with "HEADER    "')
            name = line[10:].rstrip()
        if not target and self.binder_subset is not None:
            if name not in self.binder_subset:
                # print(f"{name} not in binder subset")
                while line[:3] != "END":
                    line = file.readline()
                return False
            else:
                print(f"Loading {name} from binder subset")

        # Read the following lines until 'END', extracting lines corresponding to backbone atoms
        pdb_lines = extractBackbone("",None,False,'',file)

        # Parse the PDB lines; get input data
        data = dumpCoordsTensors(None,None,False,pdb_lines)
        data['pdb'] = name # Manually set the name
        # print('coords shape:',data['coords'].shape)

        # Add 
        # print('Loaded ',data['pdb']) #ss debug
        return data

    def prepareTargetFeatureDict(self):
        packaged_target_data = _package([(self.target_data,1)])

        # Augment with features that are important for scoring
        packaged_target_data['chain_lens'] = [self.target_data['chain_lens']]
        packaged_target_data['res_info'] = [self.target_data['res_info']]
        packaged_target_data['binder_chain_id'] = ''
        return packaged_target_data

    def createBatch(self,file):
        # load binder/complex data until the max_res_num is reached
        self.batch = []
        batch_feature_dicts = []
        n_complexes_in_batch = 0
        total_data_len = 0
        max_num_res_per_structure = 0
        n_structure_in_batch = 0
        self.current_binder_data = None

        # read through file, loading the binder/complex data and adding to batch 
        while total_data_len < self.max_res_num - max_num_res_per_structure:
            # # record start position in case we can't fit this binder/complex in the batch
            # self.binder_start_pos = file.tell()
            self.current_binder_data = self.readPDBFromOpenFile(file)

            # if at the end of the file return the batch, if there is one
            if self.current_binder_data == False:
                # skip section, not in binder list
                print('skip')
                continue
            elif self.current_binder_data == None:
                # end of the binder file
                print('none')
                if self.batch == []:
                    # nothing loaded
                    print('batch empty')
                    return None
                else:
                    # structures in batch, break loop and package
                    # file.seek(self.binder_start_pos) # ? why reset position
                    print('done loading, pack batch')
                    break

            # concatenate the binder features with the target to get the complex
            self.current_complex_data = concatFeatureDicts(self.target_data,self.current_binder_data)
            if self.current_complex_data['seq_len'] > max_num_res_per_structure:
                max_num_res_per_structure = self.current_complex_data['seq_len']
            if self.load_complex and self.load_binder:
                n_structure_in_batch += 2
                total_data_len = max_num_res_per_structure * n_structure_in_batch
                n_complexes_in_batch += 1
            else:
                n_structure_in_batch += 1
                total_data_len = max_num_res_per_structure * n_structure_in_batch
                n_complexes_in_batch += 1                

            # add binder/complex to batch
            if self.load_binder:
                self.batch.append(self.current_binder_data)
            if self.load_complex:
                self.batch.append(self.current_complex_data)

            if self.load_binder:
                batch_feature_dicts.append((self.current_binder_data,1))
            if self.load_complex:
                batch_feature_dicts.append((self.current_complex_data,1))

        if self.skip_package:
            return batch_feature_dicts
        
        # Package together the feature dictionaries and manually add the extra features necessary for binder scoring
        self.packaged_batch = _package(batch_feature_dicts)
        self.packaged_batch['chain_lens'] = [structure['chain_lens'] for structure in self.batch]
        self.packaged_batch['res_info'] = [structure['res_info'] for structure in self.batch]
        self.packaged_batch['binder_chain_id'] = [structure['chain_ids'][-1] for structure in self.batch]

        print(f"Packaged batch with {n_complexes_in_batch} unique binders. Batch has {total_data_len} residues when including padding")
        return self.packaged_batch

        
    def loadBatch(self,file):
        packaged_batch = self.createBatch(file)
        while packaged_batch != None:
            yield packaged_batch
            packaged_batch = self.createBatch(file)
            
        print('No more binder data')
        file.close()

    def __iter__(self):
        # Create an iterator
        file = open(self.filename)

        # Go to the line where the binder PDBs begin
        file.seek(self.binder_start_pos)

        # Generator iterates over structures in the multi-entry PDB file
        getFiles = self.loadBatch(file)

        return getFiles


class ComplexScoringDataset(Dataset):
    def __init__(self, pdbListPath = "", pdbList = [], targetChainID = "", binderChainID = "", mode = "binder_and_complex"):
        
        if mode == "binder_and_complex":
            self.load_complex = True
            self.load_binder = True
        elif mode == "complex_only":
            self.load_complex = True
            self.load_binder = False
        elif mode == "binder_only":
            self.load_complex = False
            self.load_binder = True
        else:
            raise ValueError(f"'{mode}' is not a valid modes")

        #Store the list of paths to PDB files in object's memory
        if pdbListPath != "":
            with open(pdbListPath,'r') as file:
                self.pdb_paths = [line.rstrip() for line in file]
        elif pdbList != []:
            self.pdb_paths = [line.rstrip() for line in pdbList]
        else:
            raise ValueError
        
        self.pdb_names = [path.split('/')[-1][:-4] for path in self.pdb_paths]
        # print(self.pdb_names)
        if targetChainID == "":
            self.target_chains = [self.get_target_chains(pdb) for pdb in self.pdb_names]
        else:
            self.target_chains = [targetChainID]*len(self.pdb_names)
        if binderChainID == "":
            self.binder_chains = [self.get_binder_chains(pdb) for pdb in self.pdb_names]
        else:
            self.binder_chains = [binderChainID]*len(self.pdb_names)
        
        # print(self.binder_chains)

        self.target_data = None
        self.target_packaged_data = None
        self.binder_data = None
        self.binder_packaged_data = None
        self.complex_data = None
        self.complex_packaged_data = None

    def get_binder_chains(self, pdb):
        split = pdb.split('_')
        assert len(split) >= 3
        # target_chain_ids = list(split[1])
        binder_chain_id = split[2]
        assert len(binder_chain_id) == 1
        return binder_chain_id

    def get_target_chains(self, pdb):
        split = pdb.split('_')
        assert len(split) >= 3
        # target_chain_ids = list(split[1])
        target_chain_ids = split[1]
        return target_chain_ids

    def __len__(self):
        """Returns length of the given dataset.
        Returns
        -------
        int
            length of dataset
        """
        return len(self.pdb_paths)

    def __getitem__(self, idx):
        """Extract a given item with provided index.
        Args
        ----
        idx : int
            Index of item to return.
        Returns
        ----
        data : dict
            Data from TERM file (as dict)
        """
        # Load the PDB file section describing the target
        pdb_lines = extractBackbone(self.pdb_paths[idx],None,False,self.binder_chains[idx],None)
        self.target_data = dumpCoordsTensors(None,None,False,pdb_lines)
        self.target_data['pdb'] = self.pdb_names[idx]+'_target' # Manually set the name
        print('Loaded target only',self.target_data['pdb'])

        # Package the data
        self.packaged_target_data = _package([(self.target_data,1)])

        # Augment with features that are important for scoring
        self.packaged_target_data['chain_lens'] = [self.target_data['chain_lens']]
        self.packaged_target_data['res_info'] = [self.target_data['res_info']]
        self.packaged_target_data['binder_chain_id'] = ''

        # Load the PDB file section describing the binder
        pdb_lines = extractBackbone(self.pdb_paths[idx],None,False,self.target_chains[idx],None)
        self.binder_data = dumpCoordsTensors(None,None,False,pdb_lines)
        self.binder_data['pdb'] = self.pdb_names[idx]+'_binder' # Manually set the name
        print('Loaded binder only',self.binder_data['pdb'])

        # Now create the complex by concatenating the feature dictionaries
        self.complex_data = concatFeatureDicts(self.target_data,self.binder_data)
        self.complex_data['pdb'] = self.pdb_names[idx]
        print('Loaded complex',self.complex_data['pdb'])

        # For consistency with BinderScoringDataset class, we will combine the binder/complex data before packaging
        batch_feature_dicts = []
        if self.load_binder:
            batch_feature_dicts.append((self.binder_data,1))
        if self.load_complex:
            batch_feature_dicts.append((self.complex_data,1))
      
        self.packaged_binder_complex_batch = _package(batch_feature_dicts)
        self.packaged_binder_complex_batch['chain_lens'] = [structure[0]['chain_lens'] for structure in batch_feature_dicts]
        self.packaged_binder_complex_batch['res_info'] = [structure[0]['res_info'] for structure in batch_feature_dicts]
        self.packaged_binder_complex_batch['binder_chain_id'] = [structure[0]['chain_ids'][-1] for structure in batch_feature_dicts]

        return self.packaged_target_data,self.packaged_binder_complex_batch

class SingleChainScoringDataset(Dataset):
    def __init__(self, pdbList):
        #Store the list of paths to PDB files in object's memory
        with open(pdbList,'r') as file:
            self.pdb_paths = [line.rstrip() for line in file]
            self.pdb_names = [path.split('/')[-1][:-4] for path in self.pdb_paths]
        
        self.structure_data = None
        self.structure_packaged_data = None

    def __len__(self):
        """Returns length of the given dataset.
        Returns
        -------
        int
            length of dataset
        """
        return len(self.pdb_paths)

    def __getitem__(self, idx):
        """Extract a given item with provided index.
        Args
        ----
        idx : int
            Index of item to return.
        Returns
        ----
        data : dict
            Data from TERM file (as dict)
        """
        # Load the PDB file section describing the target
        pdb_lines = extractBackbone(self.pdb_paths[idx],None,False,'',None)
        self.structure_data = dumpCoordsTensors(None,None,False,pdb_lines)
        self.structure_data['pdb'] = self.pdb_names[idx]
        print('Loaded structure:',self.structure_data['pdb'])

        # Package the data
        self.structure_packaged_data = _package([(self.structure_data,1)])
        self.structure_packaged_data['res_info'] = [self.structure_data['res_info']]

        return self.structure_packaged_data

def TERMLazyBatchSampler_iter(self):
    """Allows iteration over dataset."""
    if not self.ddp:
        if self.shuffle or self.semi_shuffle:
            self._cluster()
            np.random.shuffle(self.clusters)
        for batch in self.clusters:
            yield batch
    else:
        if self.shuffle or self.semi_shuffle:
            self._cluster()
            g = torch.Generator()
            g.manual_seed(self.seed + self.epoch)
            indices = torch.randperm(len(self.clusters), generator=g).tolist()
        else:
            indices = list(range(len(self.clusters)))
        
        self.num_samples = self.__len__()
        self.total_size = self.num_samples * self.num_replicas
        if not self.drop_last:
            padding_size = self.total_size - len(indices)
            if padding_size <= len(indices):
                indices += indices[:padding_size]
            else:
                indices += (indices * math.ceil(padding_size / len(indices)))[:padding_size]
        else:
            indices = indices[:self.total_size]
        assert len(indices) == self.total_size
        indices = indices[self.rank:self.total_size:self.num_replicas]
        assert len(indices) == self.num_samples

        for index in indices:
            yield self.clusters[index]
        
    

class TERMLazyBatchSamplerWrapper(object):
    """BatchSampler/Dataloader helper class for TERM data using TERMLazyDataset.

    Attributes in wrapped TERMLazyBatchSampler class
    ----------
    dataset : TERMLazyDataset
        Dataset to batch.
    size : int
        Length of dataset
    batch_size : int or None, default=4
        Size of batches created. If variable sized batches are desired, set to None.
    sort_data : bool, default=False
        Create deterministic batches by sorting the data according to the
        specified length metric and creating batches from the sorted data.
        Incompatible with :code:`shuffle=True` and :code:`semi_shuffle=True`.
    shuffle : bool, default=True
        Shuffle the data completely before creating batches.
        Incompatible with :code:`sort_data=True` and :code:`semi_shuffle=True`.
    semi_shuffle : bool, default=False
        Sort the data according to the specified length metric,
        then partition the data into :code:`semi_shuffle_cluster_size`-sized partitions.
        Within each partition perform a complete shuffle. The upside is that
        batching with similar lengths reduces padding making for more efficient computation,
        but the downside is that it does a less complete shuffle.
    semi_shuffle_cluster_size : int, default=500
        Size of partition to use when :code:`semi_shuffle=True`.
    batch_shuffle : bool, default=True
        If set to :code:`True`, shuffle samples within a batch.
    drop_last : bool, default=False
        If set to :code:`True`, drop the last samples if they don't form a complete batch.
    max_term_res : int or None, default=55000
        When :code:`batch_size=None, max_term_res>0, max_seq_tokens=None`,
        batch by fitting as many datapoints as possible with the total number of
        TERM residues included below `max_term_res`.
        Calibrated using :code:`nn.DataParallel` on two V100 GPUs.
    max_seq_tokens : int or None, default=None
        When :code:`batch_size=None, max_term_res=None, max_seq_tokens>0`,
        batch by fitting as many datapoints as possible with the total number of
        sequence residues included below `max_seq_tokens`.
    term_matches_cutoff : int or None, default=None
        Use the top :code:`term_matches_cutoff` TERM matches for featurization.
        If :code:`None`, apply no cutoff.
    term_dropout : str or None, default=None
        Let `t` be the number of TERM matches in the given datapoint.
        Select a random int `n` from 1 to `t`, and take a random subset `n`
        of the given TERM matches to keep. If :code:`term_dropout='keep_first'`,
        keep the first match and choose `n-1` from the rest.
        If :code:`term_dropout='all'`, choose `n` matches from all matches.
    """

    def __init__(self,
                 ddp):
        self.ddp = ddp
        if ddp:
            class TERMLazyBatchSampler(DistributedSampler):
                pass
        else:
            class TERMLazyBatchSampler(Sampler):
                pass
        TERMLazyBatchSampler.__init__ = TERMLazyBatchSampler_init
        TERMLazyBatchSampler._cluster = TERMLazyBatchSampler_cluster
        TERMLazyBatchSampler.package = TERMLazyBatchSampler_package
        TERMLazyBatchSampler.__len__ = TERMLazyBatchSampler_len
        TERMLazyBatchSampler.__iter__ = TERMLazyBatchSampler_iter

        self.sampler = TERMLazyBatchSampler
