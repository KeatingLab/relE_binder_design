""" GNN Potts Model Encoder modules

This file contains the GNN Potts Model Encoder, as well as an ablated version of
itself. """
from __future__ import print_function

import torch
from torch import nn

from terminator.models.layers.graph_features import MultiChainProteinFeatures
from terminator.models.layers.s2s_modules import (EdgeMPNNLayer, EdgeTransformerLayer, NodeMPNNLayer,
                                                  NodeTransformerLayer)
from terminator.models.layers.utils import (cat_edge_endpoints, cat_neighbors_nodes, gather_edges, gather_nodes,
                                            merge_duplicate_pairE)

# pylint: disable=no-member, not-callable


class AblatedPairEnergies(nn.Module):
    """Ablated GNN Potts Model Encoder

    Attributes
    ----------
    dev: str
        Device representing where the model is held
    hparams: dict
        Dictionary of parameter settings (see :code:`terminator/utils/model/default_hparams.py`)
    features : MultiChainProteinFeatures
        Module that featurizes a protein backbone (including multimeric proteins)
    W : nn.Linear
        Output layer that projects edge embeddings to proper output dimensionality
    """
    def __init__(self, hparams):
        """ Graph labeling network """
        super().__init__()
        hdim = hparams['energies_hidden_dim']
        self.hparams = hparams

        # Featurization layers
        self.features = MultiChainProteinFeatures(node_features=hdim,
                                                  edge_features=hdim,
                                                  top_k=hparams['k_neighbors'],
                                                  features_type=hparams['energies_protein_features'],
                                                  augment_eps=hparams['energies_augment_eps'],
                                                  dropout=hparams['energies_dropout'])

        self.W = nn.Linear(hparams['energies_input_dim'] * 3, hparams['energies_output_dim'])

    def forward(self, V_embed, E_embed, X, x_mask, chain_idx):
        """ Create kNN etab from TERM features, then project to proper output dimensionality.

        Args
        ----
        V_embed : torch.Tensor
            TERM node embeddings
            Shape: n_batch x n_res x n_hidden
        E_embed : torch.Tensor
            TERM edge embeddings
            Shape : n_batch x n_res x n_res x n_hidden
        X : torch.Tensor
            Backbone coordinates
            Shape: n_batch x n_res x 4 x 3
        x_mask : torch.ByteTensor
            Mask for X.
            Shape: n_batch x n_res
        chain_idx : torch.LongTensor
            Indices such that each chain is assigned a unique integer and each residue in that chain
            is assigned that integer.
            Shape: n_batch x n_res

        Returns
        -------
        etab : torch.Tensor
            Energy table in kNN dense form
            Shape: n_batch x n_res x k x n_hidden
        E_idx : torch.LongTensor
            Edge index for `etab`
            Shape: n_batch x n_res x k
        """
        # compute the kNN etab
        _, _, E_idx = self.features(X, chain_idx, x_mask)  # notably, we throw away the backbone features
        E_embed_neighbors = gather_edges(E_embed, E_idx)
        h_E = cat_edge_endpoints(E_embed_neighbors, V_embed, E_idx)
        etab = self.W(h_E)

        # merge duplicate pairEs
        n_batch, n_res, k, out_dim = etab.shape
        # ensure output etab is masked properly
        etab = etab * x_mask.view(n_batch, n_res, 1, 1)
        etab = etab.unsqueeze(-1).view(n_batch, n_res, k, 20, 20)
        etab[:, :, 0] = etab[:, :, 0] * torch.eye(20).to(etab.device) # zero off-diagonal energies
        etab = merge_duplicate_pairE(etab, E_idx)
        etab = etab.view(n_batch, n_res, k, out_dim)

        return etab, E_idx


class PairEnergies(nn.Module):
    """GNN Potts Model Encoder

    Attributes
    ----------
    dev: str
        Device representing where the model is held
    hparams: dict
        Dictionary of parameter settings (see :code:`terminator/utils/model/default_hparams.py`)
    features : MultiChainProteinFeatures
        Module that featurizes a protein backbone (including multimeric proteins)
    W_v : nn.Linear
        Embedding layer for incoming TERM node embeddings
    W_e : nn.Linear
        Embedding layer for incoming TERM edge embeddings
    edge_encoder : nn.ModuleList of EdgeTransformerLayer or EdgeMPNNLayer
        Edge graph update layers
    node_encoder : nn.ModuleList of NodeTransformerLayer or NodeMPNNLayer
        Node graph update layers
    W_out : nn.Linear
        Output layer that projects edge embeddings to proper output dimensionality
    W_proj : nn.Linear (optional)
        Output layer that projects node embeddings to proper output dimensionality.
        Enabled when :code:`hparams["node_self_sub"]=True`
    """
    def __init__(self, hparams):
        """ Graph labeling network """
        super().__init__()
        self.hparams = hparams

        hdim = hparams['energies_hidden_dim']

        # Hyperparameters
        self.node_features = hdim
        self.edge_features = hdim
        self.input_dim = hdim
        hidden_dim = hdim
        output_dim = hparams['energies_output_dim']
        dropout = hparams['energies_dropout']
        num_encoder_layers = hparams['energies_encoder_layers']

        # Featurization layers
        self.features = MultiChainProteinFeatures(node_features=hdim,
                                                  edge_features=hdim,
                                                  top_k=hparams['k_neighbors'],
                                                  features_type=hparams['energies_protein_features'],
                                                  flex_type=hparams['flex_type'],
                                                  num_ensembles=hparams['num_ensembles'],
                                                  augment_eps=hparams['energies_augment_eps'],
                                                  dropout=hparams['energies_dropout'])

        # Embedding layers
        self.W_v = nn.Linear(hdim + hparams['energies_input_dim'] + hparams['flex_input_dim'], hdim, bias=True)
        self.W_e = nn.Linear(hdim + hparams['energies_input_dim'], hdim, bias=True)
        edge_layer = EdgeTransformerLayer if not hparams['energies_use_mpnn'] else EdgeMPNNLayer
        node_layer = NodeTransformerLayer if not hparams['energies_use_mpnn'] else NodeMPNNLayer


        self.num_encoder_ensemble_layers = 1
        self.hidden_size_multiplier = 1
        # if enabled, condense conformational ensemble graph embeddings
        if hparams['use_flex'] and hparams['flex_type'].find("conformational_ensembles") > -1 and hparams['flex_type'].find("graph") > -1:
            if hparams['flex_type'].find("stack") > -1:
                self.condense_graphs = nn.Linear(hparams['num_ensembles'], 1)
            else:
                self.condense_graphs = nn.Linear(hidden_dim*hparams['num_ensembles'], hidden_dim)
            if hparams['flex_type'].find('unique') > -1:
                self.num_encoder_ensemble_layers = hparams['num_ensembles']
            if hparams['flex_type'].find('shared') > -1:
                self.hidden_size_multiplier = 1 # hparams['num_ensembles']
        # Encoder layers 
        if hparams['use_flex'] and hparams['flex_type'].find("conformational_ensembles") > -1 and hparams['flex_type'].find("graph") > -1 and ( hparams['flex_type'].find("unique") > -1 or  hparams['flex_type'].find("shared") > -1):
            self.edge_encoder = nn.ModuleList([nn.ModuleList(
                [edge_layer(hidden_dim, hidden_dim * 3 * self.hidden_size_multiplier, dropout=dropout, merge_pair=hparams['merge_pair_embeddings']) for _ in range(num_encoder_layers)]) 
                                for _ in range(self.num_encoder_ensemble_layers)]) 
            self.node_encoder = nn.ModuleList([nn.ModuleList(
                [node_layer(hidden_dim, hidden_dim * 2 * self.hidden_size_multiplier, dropout=dropout, num_ensembles=hparams['num_ensembles']) for _ in range(num_encoder_layers)])
                                for _ in range(self.num_encoder_ensemble_layers)]) 
        else:
            self.edge_encoder = nn.ModuleList(
                [edge_layer(hidden_dim, hidden_dim * 3 * self.hidden_size_multiplier, dropout=dropout,  merge_pair=hparams['merge_pair_embeddings']) for _ in range(num_encoder_layers)])
            self.node_encoder = nn.ModuleList(
                [node_layer(hidden_dim, hidden_dim * 2 * self.hidden_size_multiplier, dropout=dropout, num_ensembles=hparams['num_ensembles']) for _ in range(num_encoder_layers)])

        # self.edge_encoder = nn.ModuleList(
        #     [edge_layer(hidden_dim, hidden_dim * 3, dropout=dropout) for _ in range(num_encoder_layers)]) 
                            
        # self.node_encoder = nn.ModuleList(
        #     [node_layer(hidden_dim, hidden_dim * 2, dropout=dropout) for _ in range(num_encoder_layers)])

        # if enabled, generate self energies in etab from node embeddings
        if "node_self_sub" in hparams.keys() and hparams["node_self_sub"] is True:
            self.W_proj = nn.Linear(hidden_dim, 20)

        # project edges to proper output dimensionality
        self.W_out = nn.Linear(hidden_dim, output_dim, bias=True)
        self.neighborhood_aggregator = None
        if hparams["sscore_module"] == "linear":
            self.S_out = nn.Linear(hidden_dim, 1, bias=True)
        elif hparams["sscore_module"] == "mlp":
            mlp_n_hidden = hparams["sscore_module_nhidden"]
            self.S_out = nn.Sequential(
                nn.Linear(hidden_dim, mlp_n_hidden),
                nn.ReLU(),
                nn.Linear(mlp_n_hidden, 1, bias=True),
            )
        elif hparams["sscore_module"] == "aggnode":
            self.S_out = nn.ModuleList(
                [node_layer(hidden_dim, hidden_dim * 2 * self.hidden_size_multiplier, dropout=dropout, num_ensembles=hparams['num_ensembles']),
                nn.Linear(hidden_dim, 1, bias=True)]
            )
        else:
            raise ValueError(f"sscore_module: {self.hparams['sscore_module']} not recognized")     

        # Initialization
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(self, V_embed, E_embed, X, x_mask, chain_idx, pos_variation=None):
        """ Create kNN etab from backbone and TERM features, then project to proper output dimensionality.

        Args
        ----
        V_embed : torch.Tensor or None
            TERM node embeddings. None only accepted if :code:`hparams['energies_input_dim']=0`.
            Shape: n_batch x n_res x n_hidden
        E_embed : torch.Tensor or None
            TERM edge embeddings. None only accepted if :code:`hparams['energies_input_dim']=0`.
            Shape : n_batch x n_res x n_res x n_hidden
        X : torch.Tensor
            Backbone coordinates
            Shape: n_batch x n_res x 4 x 3
        x_mask : torch.ByteTensor
            Mask for X.
            Shape: n_batch x n_res
        chain_idx : torch.LongTensor
            Indices such that each chain is assigned a unique integer and each residue in that chain
            is assigned that integer.
            Shape: n_batch x n_res
        pos_variation : torch.LongTensor
            Positional variation for each residue calculated by NRGTEN
            Shape : n_batch * n_res * 1

        Returns
        -------
        etab : torch.Tensor
            Energy table in kNN dense form
            Shape: n_batch x n_res x k x n_hidden
        E_idx : torch.LongTensor
            Edge index for `etab`
            Shape: n_batch x n_res x k
        sscore : torch.Tensor
            Structure scores per residue
            Shape: n_batch x n_res
        """
        # Prepare node and edge embeddings
        if self.hparams['use_flex'] and self.hparams['flex_type'].find("conformational_ensembles") > -1 and self.hparams['flex_type'].find("graph") > -1:
            X_list = torch.split(X, 1, 4)
            X_list = [X[:,:,:,:,0] for X in X_list]
            print(f"Spliting coord input into {len(X_list)} tensors for separate graph generation.")
        else:
            X_list = [X]
        h_V_list_encode = []
        h_E_list_encode = []
        E_idx_list_encode = []
        multiple_ensembles = False
        for X in X_list:
            
            V, E, E_idx = self.features(X, chain_idx, x_mask)
            if self.hparams['energies_input_dim'] != 0:
                if not self.hparams['use_coords']:  # this is hacky/inefficient but i am lazy
                    V = torch.zeros_like(V)
                    E = torch.zeros_like(E)
                # fuse backbone and TERM embeddings
                h_V = self.W_v(torch.cat([V, V_embed], dim=-1))
                E_embed_neighbors = gather_edges(E_embed, E_idx)
                h_E = self.W_e(torch.cat([E, E_embed_neighbors], dim=-1))
            elif self.hparams['use_flex'] and self.hparams['flex_type'].find('dynamic_signatures') > -1:
                h_V = self.W_v(torch.cat([V, V_embed], dim=-1))
                h_E = self.W_e(E)
            elif self.hparams['use_flex'] and self.hparams['flex_type'].find('conformational_pos_variation') > -1:
                h_V = self.W_v(torch.cat([V, pos_variation], dim=-1))
                h_E = self.W_e(E)
            else:
                # just use backbone features
                h_V = self.W_v(V)
                h_E = self.W_e(E)
            h_V_list_encode.append(h_V)
            h_E_list_encode.append(h_E)
            E_idx_list_encode.append(E_idx)
        
        h_V_list = []
        h_E_list = []
        if self.hparams['use_flex'] and self.hparams['flex_type'].find("conformational_ensembles") > -1 and self.hparams['flex_type'].find("graph") > -1:
            E_idx = E_idx_list_encode[round(self.hparams['num_ensembles'] / 2)]
        else:
            E_idx = E_idx_list_encode[0]
        k = E_idx.shape[2]
        shared_graph_calculation = False

        for i, (h_V, h_E) in enumerate(zip(h_V_list_encode, h_E_list_encode)):
            encoder_ensemble_layer_id = i % self.num_encoder_ensemble_layers
            if self.hparams['use_flex'] and self.hparams['flex_type'].find("conformational_ensembles") > -1 and self.hparams['flex_type'].find("graph") > -1 and self.hparams['flex_type'].find('unique') > -1:
                print(f"\tUsing layer id {encoder_ensemble_layer_id}.")
            # Graph updates
            mask_attend = gather_nodes(x_mask.unsqueeze(-1), E_idx).squeeze(-1)
            mask_attend = x_mask.unsqueeze(-1) * mask_attend
            if self.hparams['use_flex'] and self.hparams['flex_type'].find("conformational_ensembles") > -1 and self.hparams['flex_type'].find("graph") > -1 and ( self.hparams['flex_type'].find("unique") > -1 or  self.hparams['flex_type'].find("shared") > -1):
                edge_iter = self.edge_encoder[encoder_ensemble_layer_id]
                node_iter =  self.node_encoder[encoder_ensemble_layer_id]
            else:
                edge_iter = self.edge_encoder
                node_iter = self.node_encoder
            for edge_layer, node_layer in zip(edge_iter,node_iter):
                h_EV_edges = cat_edge_endpoints(h_E, h_V, E_idx)
                if self.hparams['use_flex'] and self.hparams['flex_type'].find("conformational_ensembles") > -1 and self.hparams['flex_type'].find("graph") > -1 and self.hparams['flex_type'].find('shared') > -1:
                    h_EV_edges_list = [h_EV_edges]
                    for j, (h_V_other, h_E_other) in enumerate(zip(h_V_list_encode, h_E_list_encode)):
                        if j == i:
                            continue
                        h_EV_edges_list.append(cat_edge_endpoints(h_E_other, h_V_other, E_idx))
                    h_EV_edges = torch.cat(h_EV_edges_list, -2)
                    h_E = torch.cat(h_E_list_encode, -2)
                    h_E = edge_layer(h_E, h_EV_edges, E_idx, mask_E=x_mask, mask_attend=mask_attend, multiple_ensembles=True, k_neighbors=k)
                    h_E_list_encode = torch.split(h_E, k, -2)
                    h_E = h_E_list_encode[i]
                else:
                    h_E = edge_layer(h_E, h_EV_edges, E_idx, mask_E=x_mask, mask_attend=mask_attend)
                h_EV_nodes = cat_neighbors_nodes(h_V, h_E, E_idx)
                h_V_expand = h_V.unsqueeze(-2).expand(-1, -1, h_EV_nodes.size(-2), -1)
                if self.hparams['use_flex'] and self.hparams['flex_type'].find("conformational_ensembles") > -1 and self.hparams['flex_type'].find("graph") > -1 and self.hparams['flex_type'].find('shared') > -1:
                    h_EV_nodes_list = [h_EV_nodes]
                    h_V_expand_list = [h_V_expand]
                    for j, (h_V_other, h_E_other) in enumerate(zip(h_V_list_encode, h_E_list_encode)):
                        if j == i:
                            continue
                        h_EV_nodes_list.append(cat_neighbors_nodes(h_V_other, h_E_other, E_idx))
                        h_V_expand_list.append(h_V_other.unsqueeze(-2).expand(-1, -1, h_EV_nodes_list[-1].size(-2), -1))
                    h_EV_nodes = torch.cat(h_EV_nodes_list, -2)
                    h_V_expand = torch.cat(h_V_expand_list, -2)
                    h_V = node_layer(h_V_expand, h_EV_nodes, mask_V=x_mask, mask_attend=mask_attend, multiple_ensembles=True, k_neighbors=k)
                    h_V_list_encode = torch.split(h_V, 1, -2)
                    h_V_list_encode = [torch.squeeze(h_V_tensor, dim=2) for h_V_tensor in h_V_list_encode]
                    h_V = h_V_list_encode[i]
                    shared_graph_calculation = True
                else:
                    h_V = node_layer(h_V_expand, h_EV_nodes, mask_V=x_mask, mask_attend=mask_attend, multiple_ensembles=multiple_ensembles, multiple_ensembles_idx=i)
            if shared_graph_calculation:
                h_V_list = h_V_list_encode
                h_E_list = h_E_list_encode
                break
            h_V_list.append(h_V)
            h_E_list.append(h_E)
        # condense conformational ensemble graph embeddings
        if self.hparams['use_flex'] and self.hparams['flex_type'].find("conformational_ensembles") > -1 and self.hparams['flex_type'].find("graph") > -1:
            if self.hparams['flex_type'].find("stack") > -1:
                h_V = torch.stack(h_V_list, -1)
                h_E = torch.stack(h_E_list, -1)
            else:
                h_V = torch.cat(h_V_list, -1)
                h_E = torch.cat(h_E_list, -1)
            orig_V_shape = h_V.shape
            h_V = self.condense_graphs(h_V)
            h_E = self.condense_graphs(h_E)
            print(f"Concattenated {len(h_V_list)} graph tensors to form node and edge graph tensors of shape {h_V.shape} from tensors of shape {orig_V_shape}")

        h_V = torch.squeeze(h_V, dim=-1)
        h_E = torch.squeeze(h_E, dim=-1)

        if self.hparams['sscore_from_embedding'] == "node":
            if self.hparams["sscore_module"] == "aggnode":
                h_EV_nodes = cat_neighbors_nodes(h_V, h_E, E_idx)
                h_V_expand = h_V.unsqueeze(-2).expand(-1, -1, h_EV_nodes.size(-2), -1)
                h_V = self.S_out[0](h_V_expand, h_EV_nodes, mask_V=x_mask, mask_attend=mask_attend)
                sscore = self.S_out[1](h_V).squeeze(-1)
            else:
                # project node embedding to structure score output
                sscore = self.S_out(h_V).squeeze(-1) # b x n
        elif self.hparams['sscore_from_embedding'] == "final_edge":
            # project last edge embedding to structure score output (shouldn't work too well)
            h_E_final = h_E[:,:,-1,:] # b x n x h
            sscore = self.S_out(h_E_final).squeeze(-1) # b x n
        elif self.hparams['sscore_from_embedding'] == "self_edge":
            # project self edge embedding to structure score output (default)
            h_E_self = h_E[:,:,0,:] # b x n x h
            sscore = self.S_out(h_E_self).squeeze(-1) # b x n
        elif self.hparams['sscore_from_embedding'] == 'edge_average':
            # average embeddings around each residue and then project to structure score output
            h_E_mean = torch.mean(h_E,-2) # b x n x h
            sscore = self.S_out(h_E_mean).squeeze(-1) # b x n
        else:
            raise ValueError(f"sscore_from_embedding: {self.hparams['sscore_from_embedding']} not recognized")

        # project to output and merge duplicate pairEs
        h_E = self.W_out(h_E)
        n_batch, n_res, k, out_dim = h_E.shape
        h_E = h_E * x_mask.view(n_batch, n_res, 1, 1) # ensure output etab is masked properly
        h_E = h_E.unsqueeze(-1).view(n_batch, n_res, k, 20, 20)
        h_E[:, :, 0] = h_E[:, :, 0] * torch.eye(20).to(h_E.device) # zero off-diagonal energies
        h_E = merge_duplicate_pairE(h_E, E_idx)

        # if specified, use generate self energies from node embeddings
        if "node_self_sub" in self.hparams.keys() and self.hparams["node_self_sub"] is True:
            h_V = self.W_proj(h_V)
            h_E[..., 0, :, :] = torch.diag_embed(h_V, dim1=-2, dim2=-1)

        # reshape to fit kNN output format
        h_E = h_E.view(n_batch, n_res, k, out_dim)
        return h_E, E_idx, sscore
