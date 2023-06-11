import math
from re import X
import numpy as np

import torch
import torch.nn as nn
from numpy import linalg as LA

class residueBackbone:
    def __init__(self, coords: np.array, res_idx: int, target: bool, prev_res = None):
        # res_coords: 4 x 3
        assert coords.shape == (4,3)
        self.res_idx = res_idx
        self.target = target
        # print("res_idx",self.res_idx)

        # get the backbone atom coordinates (x,y,z)
        self.N = coords[0]
        self.Ca = coords[1]
        self.C = coords[2]
        self.O = coords[3] 

        # in case coords consist of multiple chains, need to check if residue atoms are close enough to be bonded
        self.prev_res = prev_res if not None and self.resBonded(prev_res) else None

        # define import bond vectors
        self.placeAmideHydrogen()
        self.getCarbonyl()

    def resBonded(self, prev_res):
        if prev_res is None:
            return False
        peptideBondDistance = LA.norm(prev_res.C - self.N)
        return True if peptideBondDistance < 2.0 else False # 1.3 is the optimal distance, but we can safely use a higher value in case there are odd structures

    def placeAmideHydrogen(self):
        if self.prev_res == None:
            self.H = None
            self.N_to_H = None
        else:
            # find the negative bisector of normalized N_i -> C_i-1 and  N_i -> Ca_i
            N_to_prevC = self.prev_res.C - self.N
            N_to_prevC = N_to_prevC / LA.norm(N_to_prevC)
            N_to_Ca = self.Ca - self.N
            N_to_Ca = N_to_Ca / LA.norm(N_to_Ca)

            nBisector = -(N_to_prevC + N_to_Ca)

            self.N_to_H = nBisector/LA.norm(nBisector) * 1.0 #1.0 Å is the N-H bond length and can be adjusted if need be
            self.H = self.N + self.N_to_H

    def getCarbonyl(self):
        return self.O - self.C

    def getLonePairPlaneNormal(self):
        # print("Ca",self.Ca)
        # print("C",self.C)
        # print("O",self.O)
        Ca_to_C = self.C - self.Ca
        C_to_O = self.O - self.C
        return np.cross(Ca_to_C,C_to_O)

    def distance(self,other_R):
        return LA.norm(self.Ca - other_R.Ca)

    def __hash__(self):
        return hash(self.res_idx)

class boundingBox:
    def __init__(self, pad):
        self.pad = pad
        self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max = math.inf, -math.inf, math.inf, -math.inf, math.inf, -math.inf

    def addPoint(self, x, y, z):
        if x < self.x_min:
            self.x_min = x - self.pad
        if x > self.x_max:
            self.x_max = x + self.pad
        if y < self.y_min:
            self.y_min = y - self.pad
        if y > self.y_max:
            self.y_max = y + self.pad
        if z < self.z_min:
            self.z_min = z - self.pad
        if z > self.z_max:
            self.z_max = z + self.pad

    def addResidue(self, R):
        self.addPoint(R.Ca[0],R.Ca[1],R.Ca[2])

    def addResidues(self, R_list):
        for R in R_list:
            self.addResidue(R)

    def getXWidth(self):
        return self.x_max - self.x_min

    def getYWidth(self):
        return self.y_max - self.y_min
    
    def getZWidth(self):
        return self.z_max - self.z_min

    def reportBoundaries(self):
        print(f"Bounding box with boundaries: x_min = {self.x_min}, x_max = {self.x_max}, y_min = {self.y_min}, y_max = {self.y_max}, z_min = {self.z_min}, z_max = {self.z_max}")


class CartesianLookupTable:
    def __init__(self, R_list, pad, bin_width):
        self.bin_width = bin_width
        self.bb = boundingBox(pad)
        self.bb.addResidues(R_list)
        self.addResidues(R_list)
        self.bb.reportBoundaries()
        return

    def addResidues(self, R_list):
        # lets try storing the residues with a dictionary (average lookup time is O(1), which should be fast enough)
        self.residue_map = dict()
        for R in R_list:
            R_hash = self.hashResidue(R)
            if R_hash not in self.residue_map:
                self.residue_map[R_hash] = [R]
            else:
                self.residue_map[R_hash].append(R)

    def getBin(self, val, min_val):
        return int((val - min_val) // self.bin_width)

    def hashResidue(self, R):
        x, y, z = R.Ca[0], R.Ca[1], R.Ca[2]
        x_bin = self.getBin(x, self.bb.x_min)
        y_bin = self.getBin(y, self.bb.y_min)
        z_bin = self.getBin(z, self.bb.z_min)
        # print(f"x = {x}, y = {y}, z = {z}. hash = {x_bin},{y_bin},{z_bin}")
        return (x_bin,y_bin,z_bin)

    def getNeighborhood(self, query_R, distance_cutoff):
        query_x, query_y, query_z  = self.hashResidue(query_R)
        n_neighbor_bins = int(distance_cutoff // self.bin_width)
        neighbour_hash_list = [(x,y,z) 
                                for x in range(query_x-n_neighbor_bins,query_x+n_neighbor_bins+1)
                                for y in range(query_y-n_neighbor_bins,query_y+n_neighbor_bins+1)
                                for z in range(query_z-n_neighbor_bins,query_z+n_neighbor_bins+1)]
        return neighbour_hash_list

    def getNearbyRes(self, query_R, distance_cutoff):
        nearby_res_list = []
        neighbour_hashes = self.getNeighborhood(query_R,distance_cutoff)
        for hash in neighbour_hashes:
            if hash in self.residue_map:
                for R in self.residue_map[hash]:
                    if query_R.distance(R) < distance_cutoff:
                        nearby_res_list.append(R)
        return nearby_res_list


class bbHBond:
    def __init__(self, target_coords: np.array, verbose: bool = False):
        # target_coords: n_res x 4 x 3
        self.target_residue_list = []
        self.binder_residue_list = []
        self.min_dist = 7.5
        self.verbose = verbose

        self.hbonds_file = None
        self.hbonds_energy_file = None
        self.exposed_groups_file = None

        self.hbond_info_dict = None # dict[(donor_res,acceptor_res)] = E
        self.res_hbond_satisfied = None # dict[res] = [Bool,Bool] whether the NH or CO is already participating in a h-bond
        self.exposed_donor_res = None # dict[res] = E
        self.exposed_acceptor_res = None # dict[res] = E

        self.setTargetResidues(target_coords)

    def setTargetResidues(self, target_coords: np.array):
        print(f"setting {target_coords.shape[0]} target residues")
        self.target_residue_list = []
        prev_res = None
        for i,res_coords in enumerate(target_coords):
            self.target_residue_list.append(residueBackbone(res_coords,i,True,prev_res))
            prev_res = self.target_residue_list[-1]
        self.n_target_res = len(self.target_residue_list)
            
        # create 3D hash table and place residues by Ca coordinates
        self.res_hash_table = CartesianLookupTable(self.target_residue_list,1.0,3.0)

    def setBinderResidues(self, binder_coords: np.array):
        # binder_coords: n_res x 4 x 3
        print(f"setting {binder_coords.shape[0]} binder residues")
        self.binder_residue_list = []
        prev_res = None
        for i,res_coords in enumerate(binder_coords):
            self.binder_residue_list.append(residueBackbone(res_coords,i+self.n_target_res,False,prev_res))
            prev_res = self.binder_residue_list[-1]
        self.n_binder_res = len(self.binder_residue_list)

    def findHBonds(self):
        self.hbond_info_dict = dict()
        nearby_res = dict() #dict[resA] = [resB,....,resC] residues with Ca within min_dist
        self.res_hbond_satisfied = {res:[False,False] for res in self.binder_residue_list}
        self.exposed_donor_res = dict() 
        self.exposed_acceptor_res = dict()
        # find internal hydrogen bonds
        for i,rA in enumerate(self.binder_residue_list):
            for j,rB in enumerate(self.binder_residue_list):
                if i >= j:
                    # we only want unique pairs of residues
                    continue
                if LA.norm(rA.Ca - rB.Ca) > self.min_dist:
                    # too far to interact
                    continue
                if rA not in nearby_res:
                    nearby_res[rA] = [rB]
                else:
                    nearby_res[rA].append(rB)

                if rB not in nearby_res:
                    nearby_res[rB] = [rA]
                else:
                    nearby_res[rB].append(rA)
                self.checkForHBond(rA,rB)

        # find interface hydrogen bonds
        for rBinder in self.binder_residue_list:
            nearby_target_res = self.res_hash_table.getNearbyRes(rBinder,self.min_dist)
            nearby_res[rBinder].extend(nearby_target_res)
            for rTarget in nearby_target_res:
                self.checkForHBond(rBinder,rTarget)

        # now check for binder residues that are free to hydrogen bond with water
        for i,rBinder in enumerate(self.binder_residue_list):
            is_donor,is_acceptor = self.res_hbond_satisfied[rBinder]
            if not is_donor and i != 0:
                self.checkIfNHExposed(rBinder,nearby_res[rBinder])
            if not is_acceptor and i != self.n_binder_res-1:
                self.checkIfCOExposed(rBinder,nearby_res[rBinder])

    def reportHBonds(self):
        print(f"There are {len(self.hbond_info_dict)} hydrogen bonds")
        for hbond_res in self.hbond_info_dict:
            donorR = hbond_res[0]
            acceptorR = hbond_res[1]
            print(f"donor residue (target? {donorR.target}): {donorR.res_idx}, acceptor residue (target? {acceptorR.target}): {acceptorR.res_idx}, energy: {self.hbond_info_dict[hbond_res]}")

    def writeHBonds(self, name, target_res_info, binder_res_info):
        if self.hbonds_file == None:
            self.hbonds_file = open('all_hbonds.csv','w')
            self.hbonds_file.write(f"name,donor_chain,donor_resnum,acceptor_chain,acceptor_resnum,energy"+'\n')
        for hbond_res in self.hbond_info_dict:
            donorR = hbond_res[0]
            acceptorR = hbond_res[1]
            donor_chain,donor_resnum = target_res_info[donorR.res_idx] if donorR.target else binder_res_info[donorR.res_idx - self.n_target_res]
            acceptor_chain,acceptor_resnum = target_res_info[acceptorR.res_idx] if acceptorR.target else binder_res_info[acceptorR.res_idx - self.n_target_res]
            self.hbonds_file.write(f"{name},{donor_chain},{donor_resnum},{acceptor_chain},{acceptor_resnum},{self.hbond_info_dict[hbond_res]}"+'\n')

    def writeExposed(self, name, binder_res_info):
        if self.exposed_groups_file == None:
            self.exposed_groups_file = open('all_exposed_groups.csv','w')
            self.exposed_groups_file.write(f"name,chain,resnum,group"+'\n')
        for R in self.exposed_donor_res:
            chain,resnum = binder_res_info[R.res_idx - self.n_target_res]
            self.exposed_groups_file.write(f"{name},{chain},{resnum},NH"+'\n')
        for R in self.exposed_acceptor_res:
            chain,resnum = binder_res_info[R.res_idx - self.n_target_res]
            self.exposed_groups_file.write(f"{name},{chain},{resnum},CO"+'\n')

    def writeHBondsTotalE(self, name):
        if self.hbonds_energy_file == None:
            self.hbonds_energy_file = open('seed_hbond_energy.csv','w')
            self.hbonds_energy_file.write(f"name,n_hbonds,n_exposed_nh,n_exposed_co,energy"+'\n')
        total_E = 0
        for hbond_res in self.hbond_info_dict:
            total_E+=self.hbond_info_dict[hbond_res]
        for res in self.exposed_donor_res:
            total_E+=self.exposed_donor_res[res]
        for res in self.exposed_acceptor_res:
            total_E+=self.exposed_acceptor_res[res]
        self.hbonds_energy_file.write(f"{name},{len(self.hbond_info_dict)},{len(self.exposed_donor_res)},{len(self.exposed_acceptor_res)},{total_E}"+'\n')

    def close_files(self):
        self.hbonds_file.close()
        self.hbonds_energy_file.close()
        self.exposed_groups_file.close()

    def checkForHBond(self, R1, R2):
        # print(f"checking {R1.res_idx} and {R2.res_idx}")
        # only consider interacting residues with energies of < -0.25 kcal/mol (this is taken from the paper)
        # R1 is donor, R2 is acceptor
        E1 = self.hBondEnergy(R1,R2)
        if E1 < -0.25:
            self.hbond_info_dict[(R1,R2)] = E1
            if not R1.target:
                self.res_hbond_satisfied[R1][0] = True
            if not R2.target:
                self.res_hbond_satisfied[R2][1] = True

        # R2 is donor, R1 is acceptor
        E2 = self.hBondEnergy(R2,R1)
        if E2 < -0.25:
            self.hbond_info_dict[(R2,R1)] = E2
            if not R2.target:
                self.res_hbond_satisfied[R2][0] = True
            if not R1.target:
                self.res_hbond_satisfied[R1][1] = True

    def hBondEnergy(self, donor: residueBackbone, acceptor: residueBackbone):
        # approach from https://pubmed.ncbi.nlm.nih.gov/2709375/ (note that STRIDE paper seems to use slightly different notation)
        # the energy depends on the distance between the donor and acceptor atoms and the orientation of the hydrogen bond
        if donor.prev_res == None:
            # ignore if the donor is at the N-terminus (unclear how to handle the sp3 hybridized NH3+)
            return 0

        r,t,t_o,t_1 = self.getGeometricParameters(donor, acceptor)

        E_m = -2.8 
        r_m = 3
        C = -3*E_m*r_m**8 # kcal Å^8 / mol
        D = -4*E_m*r_m**6 # kcal Å^6 / mol
        E_r = (C/r**8) - (D/r**6)
        if E_r > -0.25:
            if self.verbose:
                print("Atoms are at a distance such that E_r > -0.25, skip the rest of calculation")
            return 0.0

        E_t = math.cos(t)**2 if t < (math.pi/2) else 0

        degrees_to_radians_110 = 110 * (math.pi/180)
        if t_1 < (math.pi/2):
            E_p = (0.9 + 0.1 * math.sin(2*t_1)) * math.cos(t_o)
        elif t_1 < degrees_to_radians_110:
            K_1 = 0.9/(math.cos(degrees_to_radians_110)**6)
            K_2 = math.cos(degrees_to_radians_110)**2
            E_p = K_1*(K_2-math.cos(t_1)**2)**3 * math.cos(t_o)
        else:
            E_p = 0

        E = E_r * E_t * E_p
        
        if self.verbose:
            print(f"donor: target = {donor.target}, res_idx = {donor.res_idx}, acceptor: target = {acceptor.target}, res_idx = {acceptor.res_idx}")
            print(f"r = {r}, t = {t}, t_o = {t_o}, t_1 = {t_1}")
            print(f"E = {E}, E_r = {E_r}, E_t = {E_t}, E_p = {E_p}")
        return E

    def getGeometricParameters(self, donor: residueBackbone, acceptor: residueBackbone):
        # r is the distance between donor and acceptor
        r = LA.norm(donor.N - acceptor.O) 

        # t is the angle between donor N -> donor H and donor N -> acceptor O vectors (ideal t is 0)
        t = self.angleBetweenVectors(donor.N_to_H,acceptor.O-donor.H)

        # t_o is the angle between the acceptor O -> donor H and component of acceptor O -> donor H in the plane of the donor O lone pairs (ideal t_o is 0)
        # find using vector normal to the plane
        acceptorO_to_donorH = donor.H-acceptor.O
        lonePairNormal = acceptor.getLonePairPlaneNormal()
        # print('acceptorO_to_donorH',acceptorO_to_donorH)
        # print('lonePairNormal',lonePairNormal)
        alpha = self.angleBetweenVectors(acceptorO_to_donorH,lonePairNormal)
        # print('alpha',alpha)
        t_o = min(abs((math.pi/2)-alpha),abs(alpha-(math.pi/2)))

        # t_1 is the angle between acceptor O -> donor H in the plane of the donor O lone pairs and acceptor C -> and acceptor O (ideal t_1 is ~pi/3 or 60 degrees)
        # 1) project acceptor O -> donor H onto the normal vector
        acceptorO_to_donorH_lonePairNormal_projection = self.vectorProjection(acceptorO_to_donorH,lonePairNormal)
        # 2) find the acceptor O -> donor H in the plane of the donor O lone pairs
        acceptorO_to_donorH_inLonePairPlane = acceptorO_to_donorH - acceptorO_to_donorH_lonePairNormal_projection
        # 3) find t_1 between the O -> H in the plane and the carbonyl bond
        t_1 = self.angleBetweenVectors(acceptorO_to_donorH_inLonePairPlane,donor.getCarbonyl())

        return r,t,t_o,t_1

    def checkIfNHExposed(self,R,nearby_R):
        # place water optimally
        W = self.placeWaterAroundNH(R)

        # check if any atoms have significant overlap
        if not self.hasOverlap(W,nearby_R):
            self.exposed_donor_res[R] = -2.8
            if self.verbose:
                print(f"NH in {R.res_idx} is exposed")

    def checkIfCOExposed(self,R,nearby_R):
        # place waters optimally
        W1,W2 = self.placeWatersAroundCarbonyl(R)

        # check if any atoms have significant overlap
        if not self.hasOverlap(W1,nearby_R):
            self.exposed_acceptor_res[R] = -4.0
            if self.verbose:
                print(f"CO in {R.res_idx} is exposed")
        elif not self.hasOverlap(W2,nearby_R):
            self.exposed_acceptor_res[R] = -4.0
            if self.verbose:
                print(f"CO in {R.res_idx} is exposed")

    def hasOverlap(self,W,other_R):
        r_w = 1.4
        r_n = 1.5
        r_o = 1.5
        r_c = 1.7
        water_to_N_cutoff = 1.4+1.5
        water_to_O_cutoff = 1.4+1.5
        water_to_C_cutoff = 1.4+1.7
        for R in other_R:
            if (LA.norm(W - R.N) < water_to_N_cutoff):
                return True
            if (LA.norm(W - R.Ca) < water_to_C_cutoff):
                return True
            if (LA.norm(W - R.C) < water_to_C_cutoff):
                return True
            if (LA.norm(W - R.O) < water_to_O_cutoff):
                return True
        return False

    def placeWaterAroundNH(self,res):
        # place a single water along the N->H (optimal distance between N and O is 3.0 Å)
        NH = res.H-res.N
        return res.N+3.0*(NH/LA.norm(NH))

    def placeWatersAroundCarbonyl(self,res):
        # place two waters, one near each lone pair (optimal distance between O is 2.8 Å)
        w1_global = self.sphericalToCartesian(2.8,math.pi/2,math.radians(30))
        w2_global = self.sphericalToCartesian(2.8,math.pi/2,math.radians(120))
        x_2,y_2,z_2 = self.constructReferenceFrameAroundNCCa(res)
        w1_local = self.globalFrameToLocal(w1_global,x_2,y_2,z_2)
        w2_local = self.globalFrameToLocal(w2_global,x_2,y_2,z_2)
        return (res.O+w1_local,res.O+w2_local)

    def angleBetweenVectors(self, v1: np.array, v2: np.array):
        return np.arccos(np.dot(v1,v2)/(LA.norm(v1)*LA.norm(v2)))

    def vectorProjection(self, v1, v2):
        v2_hat = v2/LA.norm(v2)
        v2_component_of_v1 = np.dot(v1,v2)/LA.norm(v2)
        return v2_component_of_v1*v2_hat

    def sphericalToCartesian(self,r,theta,psi):
        x = r*math.sin(theta)*math.cos(psi)
        y = r*math.sin(theta)*math.sin(psi)
        z = r*math.cos(theta)
        return np.array([x,y,z])

    def constructReferenceFrameAroundNCCa(self,res):
        '''
        Local reference frame is defined for easy conversion from spherical coordinates
        x_2 = y_2 x z_2 
        y_2 = C->O (normalized)
        z_2 = N->C x C->Ca (normalized)
        '''
        z_2 = res.getLonePairPlaneNormal()
        z_2 = z_2/LA.norm(z_2)
        y_2 = res.O - res.C
        y_2 = y_2/LA.norm(y_2)
        x_2 = np.cross(y_2,z_2)
        return x_2,y_2,z_2

    def globalFrameToLocal(self,v1,x,y,z):
        '''transform v1 (which is in the global frame) to a frame defined by M = [x,y,z]'''
        M = np.stack([x,y,z],axis=1)
        return np.matmul(M,v1)
