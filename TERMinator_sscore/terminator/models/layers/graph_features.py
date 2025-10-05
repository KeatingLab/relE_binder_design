""" Backbone featurization modules

This file contains modules which featurize the protein backbone graph via
its backbone coordinates. Adapted from https://github.com/jingraham/neurips19-graph-protein-design and https://github.com/dauparas/ProteinMPNN
"""
from audioop import bias
from traceback import print_tb
import numpy as np
import torch
from torch import nn
import torch.nn.functional as F

from .utils import gather_edges, gather_nodes

# pylint: disable=no-member


class PositionalEncodings(nn.Module):
    """ Module to generate differential positional encodings for protein graph edges """
    def __init__(self, num_embeddings):
        super().__init__()
        self.num_embeddings = num_embeddings

    def forward(self, E_idx):
        """ Generate directional differential positional encodings for edges

        Args
        ----
        E_idx : torch.LongTensor
            Protein kNN edge indices
            Shape: n_batches x seq_len x k

        Returns
        -------
        E : torch.Tensor
            Directional Diffential positional encodings for edges
            Shape: n_batches x seq_len x k x num_embeddings
        """
        dev = E_idx.device
        # i-j
        N_nodes = E_idx.size(1)
        ii = torch.arange(N_nodes, dtype=torch.float32).view((1, -1, 1)).to(dev)
        d = (E_idx.float() - ii).unsqueeze(-1)
        # Original Transformer frequencies
        frequency = torch.exp(
            torch.arange(0, self.num_embeddings, 2, dtype=torch.float32) *
            -(np.log(10000.0) / self.num_embeddings)).to(dev)
        # Grid-aligned
        # frequency = 2. * np.pi * torch.exp(
        #     -torch.linspace(
        #         np.log(self.period_range[0]),
        #         np.log(self.period_range[1]),
        #         self.num_embeddings / 2
        #     )
        # )
        angles = d * frequency.view((1, 1, 1, -1))
        E = torch.cat((torch.cos(angles), torch.sin(angles)), -1)
        return E

class PositionalChainEncodings(nn.Module):
    """ Module to generate differential positional encodings for chain info for graph edges
        Adapted from https://github.com/dauparas/ProteinMPNN
    """
    def __init__(self, num_embeddings, max_relative_feature=32):
        super(PositionalChainEncodings, self).__init__()
        self.num_embeddings = num_embeddings
        self.max_relative_feature = max_relative_feature
        self.linear = nn.Linear(2*max_relative_feature+1+1, num_embeddings)

    def forward(self, offset, mask):
        d = torch.clip(offset + self.max_relative_feature, 0, 2*self.max_relative_feature)*mask + (1-mask)*(2*self.max_relative_feature+1)
        d_onehot = torch.nn.functional.one_hot(d, 2*self.max_relative_feature+1+1)
        E = self.linear(d_onehot.float())
        return E


class ProteinFeatures(nn.Module):
    """ Protein backbone featurization based on Ingraham et al NeurIPS

    Attributes
    ----------
    embeddings : PositionalEncodings
        Module to generate differential positional embeddings for edges
    dropout : nn.Dropout
        Dropout module
    node_embeddings, edge_embeddings : nn.Linear
        Embedding layers for nodes and edges
    norm_nodes, norm_edges : nn.LayerNorm
        Normalization layers for node and edge features
    """
    def __init__(self,
                 edge_features,
                 node_features,
                 num_positional_embeddings=16,
                 num_rbf=16,
                 top_k=30,
                 features_type='full',
                 flex_type='',
                 num_ensembles=1,
                 augment_eps=0.,
                 dropout=0.1,):
        """ Extract protein features """
        super().__init__()
        self.edge_features = edge_features
        self.node_features = node_features
        self.top_k = top_k
        self.augment_eps = augment_eps
        self.num_rbf = num_rbf
        self.num_positional_embeddings = num_positional_embeddings

        # Feature types
        self.features_type = features_type
        self.feature_dimensions = {
            'coarse': (3, num_positional_embeddings + num_rbf + 7),
            'full': (6, num_positional_embeddings + num_rbf + 7),
            'dist': (6, num_positional_embeddings + num_rbf),
            'hbonds': (3, 2 * num_positional_embeddings),
            'pairwise_dist': (6, 2*num_positional_embeddings + 7 + num_rbf*25)
        }

        # Positional encodings
        self.chain_embeddings = PositionalChainEncodings(num_positional_embeddings)
        self.embeddings = PositionalEncodings(num_positional_embeddings)
        self.dropout = nn.Dropout(dropout)

        # Normalization and embedding
        node_in, edge_in = self.feature_dimensions[features_type]
        if flex_type.find("conformational_ensembles") > -1:
            if flex_type.find("raw") > -1:
                if flex_type.find("stack") > -1:
                    self.conf_embedding = nn.Linear(num_ensembles, 1, bias=True)
            elif flex_type.find("embedded") > -1:
                if flex_type.find("stack") > -1:
                    self.conf_node_embedding = nn.Linear(num_ensembles, 1, bias=True)
                    self.conf_edge_embedding = nn.Linear(num_ensembles, 1, bias=True)
                else:
                    node_in = node_in * num_ensembles
                    edge_in = edge_in * num_ensembles
        
        self.flex_type = flex_type
        self.node_embedding = nn.Linear(node_in, node_features, bias=True)
        self.edge_embedding = nn.Linear(edge_in, edge_features, bias=True)
        self.norm_nodes = nn.LayerNorm(node_features)  # Normalize(node_features)
        self.norm_edges = nn.LayerNorm(edge_features)  # Normalize(edge_features)


    def _dist(self, X, mask, eps=1E-6):
        """ Pairwise euclidean distances """
        # Convolutional network on NCHW
        mask_2D = torch.unsqueeze(mask, 1) * torch.unsqueeze(mask, 2)
        dX = torch.unsqueeze(X, 1) - torch.unsqueeze(X, 2)
        D = mask_2D * torch.sqrt(torch.sum(dX**2, 3) + eps)

        # Identify k nearest neighbors (including self)
        D_max, _ = torch.max(D, -1, keepdim=True)
        D_adjust = D + (1. - mask_2D) * (D_max + 1.0)
        if (X.shape[1] >= self.top_k):
            D_neighbors, E_idx = torch.topk(D_adjust, self.top_k, dim=-1, largest=False)
        else:
            D_neighbors, E_idx = torch.topk(D_adjust, X.shape[1], dim=-1, largest=False)
        mask_neighbors = gather_edges(mask_2D.unsqueeze(-1), E_idx)

        # Debug plot KNN
        # print(E_idx[:10,:10])
        # D_simple = mask_2D * torch.zeros(D.size()).scatter(-1, E_idx, torch.ones_like(knn_D))
        # print(D_simple)
        # fig = plt.figure(figsize=(4,4))
        # ax = fig.add_subplot(111)
        # D_simple = D.data.numpy()[0,:,:]
        # plt.imshow(D_simple, aspect='equal')
        # plt.axis('off')
        # plt.tight_layout()
        # plt.savefig('D_knn.pdf')
        # exit(0)
        return D_neighbors, E_idx, mask_neighbors

    def _rbf(self, D):
        dev = D.device
        # Distance radial basis function
        D_min, D_max, D_count = 0., 20., self.num_rbf
        D_mu = torch.linspace(D_min, D_max, D_count).to(dev)
        D_mu = D_mu.view([1, 1, 1, -1])
        D_sigma = (D_max - D_min) / D_count
        D_expand = torch.unsqueeze(D, -1)
        RBF = torch.exp(-((D_expand - D_mu) / D_sigma)**2)

        # for i in range(D_count):
        #     fig = plt.figure(figsize=(4,4))
        #     ax = fig.add_subplot(111)
        #     rbf_i = RBF.data.numpy()[0,i,:,:]
        #     # rbf_i = D.data.numpy()[0,0,:,:]
        #     plt.imshow(rbf_i, aspect='equal')
        #     plt.axis('off')
        #     plt.tight_layout()
        #     plt.savefig('rbf{}.pdf'.format(i))
        #     print(np.min(rbf_i), np.max(rbf_i), np.mean(rbf_i))
        # exit(0)
        return RBF

    def _get_rbf(self, A, B, E_idx):
        """Adapted from https://github.com/dauparas/ProteinMPNN
        """
        D_A_B = torch.sqrt(torch.sum((A[:,:,None,:] - B[:,None,:,:])**2,-1) + 1e-6) #[B, L, L]
        D_A_B_neighbors = gather_edges(D_A_B[:,:,:,None], E_idx)[:,:,:,0] #[B,L,K]
        RBF_A_B = self._rbf(D_A_B_neighbors)
        return RBF_A_B

    def _quaternions(self, R, eps=1e-10):
        """ Convert a batch of 3D rotations [R] to quaternions [Q]
            R [...,3,3]
            Q [...,4]
        """
        def _R(i, j):
            return R[:, :, :, i, j]

        # Simple Wikipedia version
        # en.wikipedia.org/wiki/Rotation_matrix#Quaternion
        # For other options see math.stackexchange.com/questions/2074316/calculating-rotation-axis-from-rotation-matrix
        diag = torch.diagonal(R, dim1=-2, dim2=-1)
        Rxx, Ryy, Rzz = diag.unbind(-1)
        magnitudes = 0.5 * torch.sqrt(
            torch.abs(1 + torch.stack([Rxx - Ryy - Rzz, -Rxx + Ryy - Rzz, -Rxx - Ryy + Rzz], -1) + eps))
        signs = torch.sign(torch.stack([_R(2, 1) - _R(1, 2), _R(0, 2) - _R(2, 0), _R(1, 0) - _R(0, 1)], -1))
        xyz = signs * magnitudes
        # The relu enforces a non-negative trace
        w = torch.sqrt(F.relu(1 + diag.sum(-1, keepdim=True))) / 2.
        Q = torch.cat((xyz, w), -1)
        Q = F.normalize(Q, dim=-1)

        return Q

    def _contacts(self, D_neighbors, mask_neighbors, cutoff=8):
        """ Contacts """
        D_neighbors = D_neighbors.unsqueeze(-1)
        neighbor_C = mask_neighbors * (D_neighbors < cutoff).type(torch.float32)
        return neighbor_C

    def _hbonds(self, X, E_idx, mask_neighbors, eps=1E-3):
        """ Hydrogen bonds and contact map
        """
        X_atoms = dict(zip(['N', 'CA', 'C', 'O'], torch.unbind(X, 2)))

        # Virtual hydrogens
        X_atoms['C_prev'] = F.pad(X_atoms['C'][:, 1:, :], (0, 0, 0, 1), 'constant', 0)
        X_atoms['H'] = X_atoms['N'] + F.normalize(
            F.normalize(X_atoms['N'] - X_atoms['C_prev'], -1) + F.normalize(X_atoms['N'] - X_atoms['CA'], -1), -1)

        def _distance(X_a, X_b):
            return torch.norm(X_a[:, None, :, :] - X_b[:, :, None, :], dim=-1)

        def _inv_distance(X_a, X_b):
            return 1. / (_distance(X_a, X_b) + eps)

        # DSSP vacuum electrostatics model
        U = (0.084 * 332) * (_inv_distance(X_atoms['O'], X_atoms['N']) + _inv_distance(X_atoms['C'], X_atoms['H']) -
                             _inv_distance(X_atoms['O'], X_atoms['H']) - _inv_distance(X_atoms['C'], X_atoms['N']))

        HB = (U < -0.5).type(torch.float32)
        neighbor_HB = mask_neighbors * gather_edges(HB.unsqueeze(-1), E_idx)
        # print(HB)
        # HB = F.sigmoid(U)
        # U_np = U.cpu().data.numpy()
        # # plt.matshow(np.mean(U_np < -0.5, axis=0))
        # plt.matshow(HB[0,:,:])
        # plt.colorbar()
        # plt.show()
        # D_CA = _distance(X_atoms['CA'], X_atoms['CA'])
        # D_CA = D_CA.cpu().data.numpy()
        # plt.matshow(D_CA[0,:,:] < contact_D)
        # # plt.colorbar()
        # plt.show()
        # exit(0)
        return neighbor_HB

    # def _pairwise_dist(self, X, E_idx, mask_neighbors, chain_idx, eps):
        # ## Calculate pairwise distances
        # ## Calculation of virtual beta-C distance adapted from Daupras et al., 2022
        # X_ca = X[:,:,1,:]
        # X_n = X[:,:,0,:]
        # X_c = X[:,:,2,:]
        # b = X_ca - X_n
        # c = X_c - X_ca
        # a = torch.cross(b, c)
        # X_b = -0.58273431*a + 0.56802827*b - 0.54067466*c + X_ca
        # X_b = torch.unsqueeze(X_b, -1)
        # X = torch.swapaxes(X, 2, 3)
        # X = torch.cat((X, X_b), dim=-1)
        # X = torch.swapaxes(X, 2, 3)
        # X_shift1 = torch.roll(X, shifts=1, dims=2)
        # X_shift2 = torch.roll(X, shifts=2, dims=2)
        # X_shift3 = torch.roll(X, shifts=3, dims=2)
        # X_shift4 = torch.roll(X, shifts=4, dims=2)
        # X = torch.cat((X, X_shift1, X_shift2, X_shift3, X_shift4), dim=2)
        # X = torch.swapaxes(X, 2, 3)
        # dX = torch.unsqueeze(X, 1) - torch.unsqueeze(X, 2)
        # D = torch.sqrt(torch.sum(dX**2, 3) + eps)
        # D_rba_list = []
        # rba_params = np.linspace(0, 20, num=16)
        # for rba_param in rba_params:
        #     D_rba_list.append(torch.exp(-1*(D - rba_param)**2))
        # pairwise_dist = torch.cat(D_rba_list, -1)
        # pairwise_dist = mask_neighbors * gather_edges(pairwise_dist.unsqueeze(-1), E_idx)

        # # Calculate chain identity
        # dev = E_idx.device
        # N_batch = E_idx.size(0)
        # N_terms = E_idx.size(1)
        # N_neighbors = E_idx.size(3)
        # chain_idx_expand = chain_idx.view(N_batch, 1, -1, 1).expand((-1, N_terms, -1, N_neighbors))
        # E_chain_idx = torch.gather(chain_idx_expand.to(dev), 2, E_idx)
        # same_chain = (E_chain_idx == E_chain_idx[:, :, :, 0:1]).to(dev)

        # pairwise_dist = torch.cat((pairwise_dist, same_chain), dim=-1)

        # print(f"Pairwise dist shape: {pairwise_dist.shape}.")
        # return pairwise_dist

    def _pairwise_dist(self, X, mask, residue_idx, chain_labels, D_neighbors, E_idx):
        """Adapted from https://github.com/dauparas/ProteinMPNN
        """
        b = X[:,:,1,:] - X[:,:,0,:]
        c = X[:,:,2,:] - X[:,:,1,:]
        a = torch.cross(b, c, dim=-1)
        Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + X[:,:,1,:]
        Ca = X[:,:,1,:]
        N = X[:,:,0,:]
        C = X[:,:,2,:]
        O = X[:,:,3,:]

        RBF_all = []
        RBF_all.append(self._rbf(D_neighbors)) #Ca-Ca
        RBF_all.append(self._get_rbf(N, N, E_idx)) #N-N
        RBF_all.append(self._get_rbf(C, C, E_idx)) #C-C
        RBF_all.append(self._get_rbf(O, O, E_idx)) #O-O
        RBF_all.append(self._get_rbf(Cb, Cb, E_idx)) #Cb-Cb
        RBF_all.append(self._get_rbf(Ca, N, E_idx)) #Ca-N
        RBF_all.append(self._get_rbf(Ca, C, E_idx)) #Ca-C
        RBF_all.append(self._get_rbf(Ca, O, E_idx)) #Ca-O
        RBF_all.append(self._get_rbf(Ca, Cb, E_idx)) #Ca-Cb
        RBF_all.append(self._get_rbf(N, C, E_idx)) #N-C
        RBF_all.append(self._get_rbf(N, O, E_idx)) #N-O
        RBF_all.append(self._get_rbf(N, Cb, E_idx)) #N-Cb
        RBF_all.append(self._get_rbf(Cb, C, E_idx)) #Cb-C
        RBF_all.append(self._get_rbf(Cb, O, E_idx)) #Cb-O
        RBF_all.append(self._get_rbf(O, C, E_idx)) #O-C
        RBF_all.append(self._get_rbf(N, Ca, E_idx)) #N-Ca
        RBF_all.append(self._get_rbf(C, Ca, E_idx)) #C-Ca
        RBF_all.append(self._get_rbf(O, Ca, E_idx)) #O-Ca
        RBF_all.append(self._get_rbf(Cb, Ca, E_idx)) #Cb-Ca
        RBF_all.append(self._get_rbf(C, N, E_idx)) #C-N
        RBF_all.append(self._get_rbf(O, N, E_idx)) #O-N
        RBF_all.append(self._get_rbf(Cb, N, E_idx)) #Cb-N
        RBF_all.append(self._get_rbf(C, Cb, E_idx)) #C-Cb
        RBF_all.append(self._get_rbf(O, Cb, E_idx)) #O-Cb
        RBF_all.append(self._get_rbf(C, O, E_idx)) #C-O
        RBF_all = torch.cat(tuple(RBF_all), dim=-1)

        offset = residue_idx[:,:,None]-residue_idx[:,None,:]
        offset = gather_edges(offset[:,:,:,None], E_idx)[:,:,:,0] #[B, L, K]

        d_chains = ((chain_labels[:, :, None] - chain_labels[:,None,:])==0).long() #find self vs non-self interaction
        E_chains = gather_edges(d_chains[:,:,:,None], E_idx)[:,:,:,0]
        E_positional = self.chain_embeddings(offset.long(), E_chains)
        E = torch.cat((E_positional, RBF_all), -1)
        return E 

    def _orientations_coarse(self, X, E_idx, eps=1e-6):
        # Pair features

        # Shifted slices of unit vectors
        dX = X[:, 1:, :] - X[:, :-1, :]
        U = F.normalize(dX, dim=-1)
        u_2 = U[:, :-2, :]
        u_1 = U[:, 1:-1, :]
        u_0 = U[:, 2:, :]
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
        AD_features = torch.stack((torch.cos(A), torch.sin(A) * torch.cos(D), torch.sin(A) * torch.sin(D)), 2)
        AD_features = F.pad(AD_features, (0, 0, 1, 2), 'constant', 0)

        # Build relative orientations
        o_1 = F.normalize(u_2 - u_1, dim=-1)
        O = torch.stack((o_1, n_2, torch.cross(o_1, n_2)), 2)
        O = O.view(list(O.shape[:2]) + [9])
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

        O_neighbors = gather_nodes(O, E_idx)
        X_neighbors = gather_nodes(X, E_idx)

        # Re-view as rotation matrices
        O = O.view(list(O.shape[:2]) + [3, 3])
        O_neighbors = O_neighbors.view(list(O_neighbors.shape[:3]) + [3, 3])

        # Rotate into local reference frames
        dX = X_neighbors - X.unsqueeze(-2)
        dU = torch.matmul(O.unsqueeze(2), dX.unsqueeze(-1)).squeeze(-1)
        dU = F.normalize(dU, dim=-1)
        R = torch.matmul(O.unsqueeze(2).transpose(-1, -2), O_neighbors)
        Q = self._quaternions(R)

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

    def _dihedrals(self, X, eps=1e-7):
        # First 3 coordinates are N, CA, C
        X = X[:, :, :3, :].reshape(X.shape[0], 3 * X.shape[1], 3)

        # Shifted slices of unit vectors
        dX = X[:, 1:, :] - X[:, :-1, :]
        U = F.normalize(dX, dim=-1)
        u_2 = U[:, :-2, :]
        u_1 = U[:, 1:-1, :]
        u_0 = U[:, 2:, :]
        # Backbone normals
        n_2 = F.normalize(torch.cross(u_2, u_1), dim=-1)
        n_1 = F.normalize(torch.cross(u_1, u_0), dim=-1)

        # Angle between normals
        cosD = (n_2 * n_1).sum(-1)
        cosD = torch.clamp(cosD, -1 + eps, 1 - eps)
        D = torch.sign((u_2 * n_1).sum(-1)) * torch.acos(cosD)

        # This scheme will remove phi[0], psi[-1], omega[-1]
        D = F.pad(D, (1, 2), 'constant', 0)
        D = D.view((D.size(0), int(D.size(1) / 3), 3))

        # print(cosD.cpu().data.numpy().flatten())
        # print(omega.sum().cpu().data.numpy().flatten())

        # Bond angle calculation
        # A = torch.acos(-(u_1 * u_0).sum(-1))

        # DEBUG: Ramachandran plot
        # x = phi.cpu().data.numpy().flatten()
        # y = psi.cpu().data.numpy().flatten()
        # plt.scatter(x * 180 / np.pi, y * 180 / np.pi, s=1, marker='.')
        # plt.xlabel('phi')
        # plt.ylabel('psi')
        # plt.axis('square')
        # plt.grid()
        # plt.axis([-180,180,-180,180])
        # plt.show()

        # Lift angle representations to the circle
        D_features = torch.cat((torch.cos(D), torch.sin(D)), 2)
        return D_features

    def forward(self, X, mask):
        """ Featurize coordinates as an attributed graph

        Args
        ----
        X : torch.Tensor
            Backbone coordinates
            Shape: n_batch x seq_len x 4 x 3
        mask : torch.ByteTensor
            Mask for residues
            Shape: n_batch x seq_len

        Returns
        -------
        V : torch.Tensor
            Node embeddings
            Shape: n_batches x seq_len x n_hidden
        E : torch.Tensor
            Edge embeddings in kNN dense form
            Shape: n_batches x seq_len x k x n_hidden
        E_idx : torch.LongTensor
            Edge indices
            Shape: n_batches x seq_len x k x n_hidden
        """
        # Data augmentation
        if self.training and self.augment_eps > 0:
            X = X + self.augment_eps * torch.randn_like(X)

        # Build k-Nearest Neighbors graph
        X_ca = X[:, :, 1, :]
        D_neighbors, E_idx, mask_neighbors = self._dist(X_ca, mask)

        # Pairwise features
        AD_features, O_features = self._orientations_coarse(X_ca, E_idx)
        RBF = self._rbf(D_neighbors)

        # Pairwise embeddings
        E_positional = self.embeddings(E_idx)

        if self.features_type == 'coarse':
            # Coarse backbone features
            V = AD_features
            E = torch.cat((E_positional, RBF, O_features), -1)
        elif self.features_type == 'hbonds':
            # Hydrogen bonds and contacts
            neighbor_HB = self._hbonds(X, E_idx, mask_neighbors)
            neighbor_C = self._contacts(D_neighbors, mask_neighbors)
            # Dropout
            neighbor_C = self.dropout(neighbor_C)
            neighbor_HB = self.dropout(neighbor_HB)
            # Pack
            V = mask.unsqueeze(-1) * torch.ones_like(AD_features)
            neighbor_C = neighbor_C.expand(-1, -1, -1, int(self.num_positional_embeddings / 2))
            neighbor_HB = neighbor_HB.expand(-1, -1, -1, int(self.num_positional_embeddings / 2))
            E = torch.cat((E_positional, neighbor_C, neighbor_HB), -1)
        elif self.features_type == 'full':
            # Full backbone angles
            V = self._dihedrals(X)
            E = torch.cat((E_positional, RBF, O_features), -1)
        elif self.features_type == 'dist':
            # Full backbone angles
            V = self._dihedrals(X)
            E = torch.cat((E_positional, RBF), -1)
        elif self.features_type == 'pairwise_dist':
            # Full backbone angles + pairwise distances
            E_pairwise = self._pairwise_dist(X, mask, _, _, D_neighbors, E_idx)
            V = self._dihedrals(X)
            E = torch.cat((E_positional, RBF), -1)

        # Embed the nodes
        V = self.node_embedding(V)
        V = self.norm_nodes(V)
        E = self.edge_embedding(E)
        E = self.norm_edges(E)

        # DEBUG
        # U = (np.nan * torch.zeros(X.size(0),X.size(1),X.size(1),3)).scatter(2, E_idx.unsqueeze(-1).expand(-1,-1,-1,3), E[:,:,:,:3])
        # plt.imshow(U.data.numpy()[0,:,:,0])
        # plt.show()
        # exit(0)
        return V, E, E_idx


class IndexDiffEncoding(nn.Module):
    """ Module to generate differential positional encodings for multichain protein graph edges

    Similar to ProteinFeatures, but zeros out features between interchain interactions """
    def __init__(self, num_embeddings):
        super().__init__()
        self.num_embeddings = num_embeddings

    def forward(self, E_idx, chain_idx):
        """ Generate directional differential positional encodings for edges

        Args
        ----
        E_idx : torch.LongTensor
            Protein kNN edge indices
            Shape: n_batches x seq_len x k
        chain_idx : torch.LongTensor
            Indices for residues such that each chain is assigned a unique integer
            and each residue in that chain is assigned that integer
            Shape: n_batches x seq_len

        Returns
        -------
        E : torch.Tensor
            Directional Diffential positional encodings for edges
            Shape: n_batches x seq_len x k x num_embeddings
        """
        dev = E_idx.device
        # i-j
        N_batch = E_idx.size(0)
        N_terms = E_idx.size(1)
        N_nodes = E_idx.size(2)
        N_neighbors = E_idx.size(3)
        ii = torch.arange(N_nodes, dtype=torch.float32).view((1, -1, 1)).to(dev)
        d = (E_idx.float() - ii).unsqueeze(-1)

        # Original Transformer frequencies
        frequency = torch.exp(
            torch.arange(0, self.num_embeddings, 2, dtype=torch.float32) *
            -(np.log(10000.0) / self.num_embeddings)).to(dev)
        # Grid-aligned
        # frequency = 2. * np.pi * torch.exp(
        #     -torch.linspace(
        #         np.log(self.period_range[0]),
        #         np.log(self.period_range[1]),
        #         self.num_embeddings / 2
        #     )
        # )
        angles = d * frequency.view((1, 1, 1, -1))
        E = torch.cat((torch.cos(angles), torch.sin(angles)), -1)

        # we zero out positional frequencies from inter-chain edges
        # the idea is, the concept of "sequence distance"
        # between two residues in different chains doesn't
        # make sense :P
        chain_idx_expand = chain_idx.view(N_batch, 1, -1, 1).expand((-1, N_terms, -1, N_neighbors))
        E_chain_idx = torch.gather(chain_idx_expand.to(dev), 2, E_idx)
        same_chain = (E_chain_idx == E_chain_idx[:, :, :, 0:1]).to(dev)

        E *= same_chain.unsqueeze(-1)
        return E


class MultiChainProteinFeatures(ProteinFeatures):
    """ Protein backbone featurization which accounts for differences
    between inter-chain and intra-chain interactions.

    Attributes
    ----------
    embeddings : IndexDiffEncoding
        Module to generate differential positional embeddings for edges
    dropout : nn.Dropout
        Dropout module
    node_embeddings, edge_embeddings : nn.Linear
        Embedding layers for nodes and edges
    norm_nodes, norm_edges : nn.LayerNorm
        Normalization layers for node and edge features
    """
    def __init__(self,
                 edge_features,
                 node_features,
                 num_positional_embeddings=16,
                 num_rbf=16,
                 top_k=30,
                 features_type='full',
                 flex_type=None,
                 num_ensembles=1,
                 augment_eps=0.,
                 dropout=0.1):
        """ Extract protein features """
        super().__init__(edge_features,
                         node_features,
                         num_positional_embeddings=num_positional_embeddings,
                         num_rbf=num_rbf,
                         top_k=top_k,
                         features_type=features_type,
                         flex_type=flex_type,
                         num_ensembles=num_ensembles,
                         augment_eps=augment_eps,
                         dropout=dropout)

        # so uh this is designed to work on the batched TERMS
        # but if we just treat the whole sequence as one big TERM
        # the math is the same so i'm not gonna code a new module lol
        self.embeddings = IndexDiffEncoding(num_positional_embeddings)

    # pylint: disable=arguments-differ
    def forward(self, X, chain_idx, mask):
        """ Featurize coordinates as an attributed graph

        Args
        ----
        X : torch.Tensor
            Backbone coordinates
            Shape: n_batch x seq_len x 4 x 3
        chain_idx : torch.LongTensor
            Indices for residues such that each chain is assigned a unique integer
            and each residue in that chain is assigned that integer
            Shape: n_batches x seq_len
        mask : torch.ByteTensor
            Mask for residues
            Shape: n_batch x seq_len

        Returns
        -------
        V : torch.Tensor
            Node embeddings
            Shape: n_batches x seq_len x n_hidden
        E : torch.Tensor
            Edge embeddings in kNN dense form
            Shape: n_batches x seq_len x k x n_hidden
        E_idx : torch.LongTensor
            Edge indices
            Shape: n_batches x seq_len x k x n_hidden
        """

        # Data augmentation
        if self.training and self.augment_eps > 0:
            X = X + self.augment_eps * torch.randn_like(X)
        # Setup input to handle flex info
        X_list = [X]
        if self.flex_type and self.flex_type.find("conformational_ensembles") > -1:
            if self.flex_type.find("raw") > -1:
                old_shape = X.shape
                X = self.conf_embedding(X)
                X = torch.squeeze(X, dim=-1)
                print(f"Transformed coord input from {old_shape} to {X.shape}.")
                X_list = [X]
            elif self.flex_type.find("embedded") > -1:
                X_list = torch.split(X, 1, 4)
                X_list = [X[:,:,:,:,0] for X in X_list]
                print(f"Spliting coord input into {len(X_list)} tensors for separate embedding.")
            elif self.flex_type.find("graph") > -1:
                X_list = [X]
        V_list = []
        E_list = []
        for X in X_list:
            # Build k-Nearest Neighbors graph
            X_ca = X[:, :, 1, :]
            D_neighbors, E_idx, mask_neighbors = self._dist(X_ca, mask)

            # Pairwise features
            AD_features, O_features = self._orientations_coarse(X_ca, E_idx)
            RBF = self._rbf(D_neighbors)

            # Pairwise embeddings
            # we unsqueeze to generate "1 TERM" per sequence,
            # then squeeze it back to get rid of it
            E_positional = self.embeddings(E_idx.unsqueeze(1), chain_idx).squeeze(1)

            if self.features_type == 'coarse':
                # Coarse backbone features
                V = AD_features
                E = torch.cat((E_positional, RBF, O_features), -1)
            elif self.features_type == 'hbonds':
                # Hydrogen bonds and contacts
                neighbor_HB = self._hbonds(X, E_idx, mask_neighbors)
                neighbor_C = self._contacts(D_neighbors, mask_neighbors)
                # Dropout
                neighbor_C = self.dropout(neighbor_C)
                neighbor_HB = self.dropout(neighbor_HB)
                # Pack
                V = mask.unsqueeze(-1) * torch.ones_like(AD_features)
                neighbor_C = neighbor_C.expand(-1, -1, -1, int(self.num_positional_embeddings / 2))
                neighbor_HB = neighbor_HB.expand(-1, -1, -1, int(self.num_positional_embeddings / 2))
                E = torch.cat((E_positional, neighbor_C, neighbor_HB), -1)
            elif self.features_type == 'full':
                # Full backbone angles
                V = self._dihedrals(X)
                E = torch.cat((E_positional, RBF, O_features), -1)
            elif self.features_type == 'dist':
                # Full backbone angles
                V = self._dihedrals(X)
                E = torch.cat((E_positional, RBF), -1)
            elif self.features_type == 'pairwise_dist':
                # Full backbone angles
                residue_idx = -100*torch.ones((chain_idx.shape[0], chain_idx.shape[1]), dtype=torch.int32)
                dev = "cpu"
                if X.is_cuda:
                    dev = "cuda:" + str(X.get_device())
                residue_idx = residue_idx.to(dtype=torch.long,device=dev)
                for i_batch, chain in enumerate(torch.split(chain_idx, 1, 0)):
                    l0 = 0
                    l1 = 0
                    list_idx, counts_idx = torch.unique_consecutive(chain, return_counts=True)
                    for c, idx in enumerate(list_idx):
                        if idx == 0 and c > 0:
                            break 
                        l1 += counts_idx[c]
                        arange_list = torch.arange(l0, l1).to(device=dev)
                        residue_idx[i_batch, l0:l1] = 100*(idx)+arange_list
                        l0 += counts_idx[c]
                
                pairwise_dist = self._pairwise_dist(X, mask, residue_idx, chain_idx, D_neighbors, E_idx)
                print(f"\t\tUsing pairwise distances; pairwise_dist_embedding shape: {pairwise_dist.shape}")
                V = self._dihedrals(X)
                E = torch.cat((E_positional, O_features, pairwise_dist), -1)
            V_list.append(V)
            E_list.append(E)
        if self.flex_type and self.flex_type.find("embedded") > -1 and self.flex_type.find("stack") > -1:
            V = torch.squeeze(self.conf_node_embedding(torch.stack(V_list, -1)), dim=-1)
            E = torch.squeeze(self.conf_node_embedding(torch.stack(E_list, -1)), dim=-1)
        else:
            V = torch.cat(V_list, -1)
            E = torch.cat(E_list, -1)
        if self.flex_type and self.flex_type.find("conformational_ensembles") > -1 and self.flex_type.find("embedded") > -1:
            print(f"Combined {len(X_list)} tensors to form node and edge tensors of shape {V.shape}")
        # V = torch.squeeze(torch.stack(V_list))
        # E = torch.squeeze(torch.stack(E_list))

        # Embed the nodes
        V = self.node_embedding(V)
        V = self.norm_nodes(V)
        E = self.edge_embedding(E)
        E = self.norm_edges(E)

        # DEBUG
        # U = (np.nan * torch.zeros(X.size(0),X.size(1),X.size(1),3)).scatter(2, E_idx.unsqueeze(-1).expand(-1,-1,-1,3), E[:,:,:,:3])
        # plt.imshow(U.data.numpy()[0,:,:,0])
        # plt.show()
        # exit(0)
        return V, E, E_idx
