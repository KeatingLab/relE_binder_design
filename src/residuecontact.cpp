#include "residuecontact.h"

/* --- --- --- --- --- checkVDWRadii --- --- --- --- --- */

mstreal checkVDWRadii::maxRadius() {
    return max_radius;
}

mstreal checkVDWRadii::maxSumRadii() {
    return checkVDWRadii::maxRadius() * 2;
}

mstreal checkVDWRadii::getRadius(const Atom& a, bool strict) {
    return getRadius(a.getParent()->getName(),a.getName(),strict);
}

mstreal checkVDWRadii::getRadius(const string& resName, const string& atomName, bool strict) {
    bool found = radii.count(resName) && radii[resName].count(atomName);
    if (!found) {
        if (strict) MstUtils::error("Atom name ("+atomName+") not found for residue ("+resName+")","checkVDWRadii::getRadius"); 
        else return 0.0; 
    }
    return radii[resName][atomName];
}

mstreal checkVDWRadii::sumRadii(const Atom& a1, const Atom& a2) {
    return (getRadius(a1)+getRadius(a2));
}

bool checkVDWRadii::clash(const Atom& a1, const Atom& a2, double lb) {
    return (a1.distance(a2) < checkVDWRadii::sumRadii(a1, a2) * lb);
}

bool checkVDWRadii::contact(const Atom& a1, const Atom& a2, double lb, double ub) {
    double dist = a1.distance(a2);
    double s = checkVDWRadii::sumRadii(a1, a2);
    return (dist >= s * lb && dist < s * ub);
}

bool checkVDWRadii::independent(const Atom& a1, const Atom& a2, double ub) {
    return (a1.distance(a2) >= checkVDWRadii::sumRadii(a1, a2) * ub);
}

atomInteractionType checkVDWRadii::interactionType(const Atom& a1, const Atom& a2, double lb, double ub) {
    double dist = a1.distance(a2);
    double s = checkVDWRadii::sumRadii(a1, a2);
    if (dist < s * lb) {
        return ATOMCLASH;
    } else if (dist < s * ub) {
        return ATOMCONTACT;
    } else {
        return NOINTERACTION;
    }
}

// bool checkVDWRadii::initConstants() {
//     // https://pubs.rsc.org/en/content/articlelanding/2013/DT/c3dt50599e
//     radii = {{"N",  1.60}, {"NT", 1.60}, {"CA", 2.365},
//              {"C", 2.10},  {"CT", 2.10}, {"O", 1.60}, 
//              {"CB", 2.2350}, {"S", 1.80}};
//     ignore_atoms = {"H"};
//     max_radius = 0;
//     for (auto it : radii) if (it.second > max_radius) max_radius = it.second;
//     // cout << "max radius: " << max_radius << endl;
//     return true;
//     naccess vdw.radii
// }

string checkVDWRadii::getName(string atomName, bool strict) {
    while (atomName.length() >= 1) {
        if (radii.count(atomName) > 0) {
            return atomName;
        } else {
            if (strict) break;
            atomName = atomName.substr(0, atomName.size() - 1);
        }
    }
    if (strict) MstUtils::error("VDW radius not defined for atom name"+atomName,"checkVDWRadii::getName");
    return "";
}

bool checkVDWRadii::initConstants() {
    // copied from naccess vdw.radii using script
    radii["ALA"]["N"] = 1.65;
    radii["ALA"]["CA"] = 1.87;
    radii["ALA"]["C"] = 1.76;
    radii["ALA"]["O"] = 1.40;
    radii["ALA"]["CB"] = 1.87;
    radii["ARG"]["N"] = 1.65;
    radii["ARG"]["CA"] = 1.87;
    radii["ARG"]["C"] = 1.76;
    radii["ARG"]["O"] = 1.40;
    radii["ARG"]["CB"] = 1.87;
    radii["ARG"]["CG"] = 1.87;
    radii["ARG"]["CD"] = 1.87;
    radii["ARG"]["NE"] = 1.65;
    radii["ARG"]["CZ"] = 1.76;
    radii["ARG"]["NH1"] = 1.65;
    radii["ARG"]["NH2"] = 1.65;
    radii["ASP"]["N"] = 1.65;
    radii["ASP"]["CA"] = 1.87;
    radii["ASP"]["C"] = 1.76;
    radii["ASP"]["O"] = 1.40;
    radii["ASP"]["CB"] = 1.87;
    radii["ASP"]["CG"] = 1.76;
    radii["ASP"]["OD1"] = 1.40;
    radii["ASP"]["OD2"] = 1.40;
    radii["ASN"]["N"] = 1.65;
    radii["ASN"]["CA"] = 1.87;
    radii["ASN"]["C"] = 1.76;
    radii["ASN"]["O"] = 1.40;
    radii["ASN"]["CB"] = 1.87;
    radii["ASN"]["CG"] = 1.76;
    radii["ASN"]["OD1"] = 1.40;
    radii["ASN"]["ND2"] = 1.65;
    radii["CYS"]["N"] = 1.65;
    radii["CYS"]["CA"] = 1.87;
    radii["CYS"]["C"] = 1.76;
    radii["CYS"]["O"] = 1.40;
    radii["CYS"]["CB"] = 1.87;
    radii["CYS"]["SG"] = 1.85;
    radii["GLU"]["N"] = 1.65;
    radii["GLU"]["CA"] = 1.87;
    radii["GLU"]["C"] = 1.76;
    radii["GLU"]["O"] = 1.40;
    radii["GLU"]["CB"] = 1.87;
    radii["GLU"]["CG"] = 1.87;
    radii["GLU"]["CD"] = 1.76;
    radii["GLU"]["OE1"] = 1.40;
    radii["GLU"]["OE2"] = 1.40;
    radii["GLN"]["N"] = 1.65;
    radii["GLN"]["CA"] = 1.87;
    radii["GLN"]["C"] = 1.76;
    radii["GLN"]["O"] = 1.40;
    radii["GLN"]["CB"] = 1.87;
    radii["GLN"]["CG"] = 1.87;
    radii["GLN"]["CD"] = 1.76;
    radii["GLN"]["OE1"] = 1.40;
    radii["GLN"]["NE2"] = 1.65;
    radii["GLY"]["N"] = 1.65;
    radii["GLY"]["CA"] = 1.87;
    radii["GLY"]["C"] = 1.76;
    radii["GLY"]["O"] = 1.40;
    radii["HIS"]["N"] = 1.65;
    radii["HIS"]["CA"] = 1.87;
    radii["HIS"]["C"] = 1.76;
    radii["HIS"]["O"] = 1.40;
    radii["HIS"]["CB"] = 1.87;
    radii["HIS"]["CG"] = 1.76;
    radii["HIS"]["ND1"] = 1.65;
    radii["HIS"]["CD2"] = 1.76;
    radii["HIS"]["CE1"] = 1.76;
    radii["HIS"]["NE2"] = 1.65;
    radii["ILE"]["N"] = 1.65;
    radii["ILE"]["CA"] = 1.87;
    radii["ILE"]["C"] = 1.76;
    radii["ILE"]["O"] = 1.40;
    radii["ILE"]["CB"] = 1.87;
    radii["ILE"]["CG1"] = 1.87;
    radii["ILE"]["CG2"] = 1.87;
    radii["ILE"]["CD1"] = 1.87;
    radii["LEU"]["N"] = 1.65;
    radii["LEU"]["CA"] = 1.87;
    radii["LEU"]["C"] = 1.76;
    radii["LEU"]["O"] = 1.40;
    radii["LEU"]["CB"] = 1.87;
    radii["LEU"]["CG"] = 1.87;
    radii["LEU"]["CD1"] = 1.87;
    radii["LEU"]["CD2"] = 1.87;
    radii["LYS"]["N"] = 1.65;
    radii["LYS"]["CA"] = 1.87;
    radii["LYS"]["C"] = 1.76;
    radii["LYS"]["O"] = 1.40;
    radii["LYS"]["CB"] = 1.87;
    radii["LYS"]["CG"] = 1.87;
    radii["LYS"]["CD"] = 1.87;
    radii["LYS"]["CE"] = 1.87;
    radii["LYS"]["NZ"] = 1.50;
    radii["MET"]["N"] = 1.65;
    radii["MET"]["CA"] = 1.87;
    radii["MET"]["C"] = 1.76;
    radii["MET"]["O"] = 1.40;
    radii["MET"]["CB"] = 1.87;
    radii["MET"]["CG"] = 1.87;
    radii["MET"]["SD"] = 1.85;
    radii["MET"]["CE"] = 1.87;
    radii["PHE"]["N"] = 1.65;
    radii["PHE"]["CA"] = 1.87;
    radii["PHE"]["C"] = 1.76;
    radii["PHE"]["O"] = 1.40;
    radii["PHE"]["CB"] = 1.87;
    radii["PHE"]["CG"] = 1.76;
    radii["PHE"]["CD1"] = 1.76;
    radii["PHE"]["CD2"] = 1.76;
    radii["PHE"]["CE1"] = 1.76;
    radii["PHE"]["CE2"] = 1.76;
    radii["PHE"]["CZ"] = 1.76;
    radii["PRO"]["N"] = 1.65;
    radii["PRO"]["CA"] = 1.87;
    radii["PRO"]["C"] = 1.76;
    radii["PRO"]["O"] = 1.40;
    radii["PRO"]["CB"] = 1.87;
    radii["PRO"]["CG"] = 1.87;
    radii["PRO"]["CD"] = 1.87;
    radii["SER"]["N"] = 1.65;
    radii["SER"]["CA"] = 1.87;
    radii["SER"]["C"] = 1.76;
    radii["SER"]["O"] = 1.40;
    radii["SER"]["CB"] = 1.87;
    radii["SER"]["OG"] = 1.40;
    radii["THR"]["N"] = 1.65;
    radii["THR"]["CA"] = 1.87;
    radii["THR"]["C"] = 1.76;
    radii["THR"]["O"] = 1.40;
    radii["THR"]["CB"] = 1.87;
    radii["THR"]["OG1"] = 1.40;
    radii["THR"]["CG2"] = 1.87;
    radii["TRP"]["N"] = 1.65;
    radii["TRP"]["CA"] = 1.87;
    radii["TRP"]["C"] = 1.76;
    radii["TRP"]["O"] = 1.40;
    radii["TRP"]["CB"] = 1.87;
    radii["TRP"]["CG"] = 1.76;
    radii["TRP"]["CD1"] = 1.76;
    radii["TRP"]["CD2"] = 1.76;
    radii["TRP"]["NE1"] = 1.65;
    radii["TRP"]["CE2"] = 1.76;
    radii["TRP"]["CE3"] = 1.76;
    radii["TRP"]["CZ2"] = 1.76;
    radii["TRP"]["CZ3"] = 1.76;
    radii["TRP"]["CH2"] = 1.76;
    radii["TYR"]["N"] = 1.65;
    radii["TYR"]["CA"] = 1.87;
    radii["TYR"]["C"] = 1.76;
    radii["TYR"]["O"] = 1.40;
    radii["TYR"]["CB"] = 1.87;
    radii["TYR"]["CG"] = 1.76;
    radii["TYR"]["CD1"] = 1.76;
    radii["TYR"]["CD2"] = 1.76;
    radii["TYR"]["CE1"] = 1.76;
    radii["TYR"]["CE2"] = 1.76;
    radii["TYR"]["CZ"] = 1.76;
    radii["TYR"]["OH"] = 1.40;
    radii["VAL"]["N"] = 1.65;
    radii["VAL"]["CA"] = 1.87;
    radii["VAL"]["C"] = 1.76;
    radii["VAL"]["O"] = 1.40;
    radii["VAL"]["CB"] = 1.87;
    radii["VAL"]["CG1"] = 1.87;
    radii["VAL"]["CG2"] = 1.87;
    radii["ASX"]["N"] = 1.65;
    radii["ASX"]["CA"] = 1.87;
    radii["ASX"]["C"] = 1.76;
    radii["ASX"]["O"] = 1.40;
    radii["ASX"]["CB"] = 1.87;
    radii["ASX"]["CG"] = 1.76;
    radii["ASX"]["AD1"] = 1.50;
    radii["ASX"]["AD2"] = 1.50;
    radii["GLX"]["N"] = 1.65;
    radii["GLX"]["CA"] = 1.87;
    radii["GLX"]["C"] = 1.76;
    radii["GLX"]["O"] = 1.40;
    radii["GLX"]["CB"] = 1.87;
    radii["GLX"]["CG"] = 1.76;
    radii["GLX"]["CD"] = 1.87;
    radii["GLX"]["AE1"] = 1.50;
    radii["GLX"]["AE2"] = 1.50;
    radii["ACE"]["C"] = 1.76;
    radii["ACE"]["O"] = 1.40;
    radii["ACE"]["CH3"] = 1.87;
    radii["__A"]["P"] = 1.90;
    radii["__A"]["O1P"] = 1.40;
    radii["__A"]["O2P"] = 1.40;
    radii["__A"]["O5*"] = 1.40;
    radii["__A"]["C5*"] = 1.80;
    radii["__A"]["C4*"] = 1.80;
    radii["__A"]["O4*"] = 1.40;
    radii["__A"]["C3*"] = 1.80;
    radii["__A"]["O3*"] = 1.40;
    radii["__A"]["C2*"] = 1.80;
    radii["__A"]["C1*"] = 1.80;
    radii["__A"]["N9"] = 1.60;
    radii["__A"]["C8"] = 1.80;
    radii["__A"]["N7"] = 1.60;
    radii["__A"]["C5"] = 1.80;
    radii["__A"]["C6"] = 1.80;
    radii["__A"]["N6"] = 1.60;
    radii["__A"]["N1"] = 1.60;
    radii["__A"]["C2"] = 1.80;
    radii["__A"]["N3"] = 1.60;
    radii["__A"]["C4"] = 1.80;
    radii["__C"]["P"] = 1.90;
    radii["__C"]["O1P"] = 1.40;
    radii["__C"]["O2P"] = 1.40;
    radii["__C"]["O5*"] = 1.40;
    radii["__C"]["C5*"] = 1.80;
    radii["__C"]["C4*"] = 1.80;
    radii["__C"]["O4*"] = 1.40;
    radii["__C"]["C3*"] = 1.80;
    radii["__C"]["O3*"] = 1.40;
    radii["__C"]["C2*"] = 1.80;
    radii["__C"]["C1*"] = 1.80;
    radii["__C"]["N1"] = 1.60;
    radii["__C"]["C2"] = 1.80;
    radii["__C"]["O2"] = 1.40;
    radii["__C"]["N3"] = 1.60;
    radii["__C"]["C4"] = 1.80;
    radii["__C"]["N4"] = 1.60;
    radii["__C"]["C5"] = 1.80;
    radii["__C"]["C6"] = 1.80;
    radii["__G"]["P"] = 1.90;
    radii["__G"]["O1P"] = 1.40;
    radii["__G"]["O2P"] = 1.40;
    radii["__G"]["O5*"] = 1.40;
    radii["__G"]["C5*"] = 1.80;
    radii["__G"]["C4*"] = 1.80;
    radii["__G"]["O4*"] = 1.40;
    radii["__G"]["C3*"] = 1.80;
    radii["__G"]["O3*"] = 1.40;
    radii["__G"]["C2*"] = 1.80;
    radii["__G"]["C1*"] = 1.80;
    radii["__G"]["N9"] = 1.60;
    radii["__G"]["C8"] = 1.80;
    radii["__G"]["N7"] = 1.60;
    radii["__G"]["C5"] = 1.80;
    radii["__G"]["C6"] = 1.80;
    radii["__G"]["O6"] = 1.40;
    radii["__G"]["N1"] = 1.60;
    radii["__G"]["C2"] = 1.80;
    radii["__G"]["N2"] = 1.60;
    radii["__G"]["N3"] = 1.60;
    radii["__G"]["C4"] = 1.80;
    radii["__T"]["P"] = 1.90;
    radii["__T"]["O1P"] = 1.40;
    radii["__T"]["O2P"] = 1.40;
    radii["__T"]["O5*"] = 1.40;
    radii["__T"]["C5*"] = 1.80;
    radii["__T"]["C4*"] = 1.80;
    radii["__T"]["O4*"] = 1.40;
    radii["__T"]["C3*"] = 1.80;
    radii["__T"]["O3*"] = 1.40;
    radii["__T"]["C2*"] = 1.80;
    radii["__T"]["C1*"] = 1.80;
    radii["__T"]["N1"] = 1.60;
    radii["__T"]["C2"] = 1.80;
    radii["__T"]["O2"] = 1.40;
    radii["__T"]["N3"] = 1.60;
    radii["__T"]["C4"] = 1.80;
    radii["__T"]["O4"] = 1.40;
    radii["__T"]["C5"] = 1.80;
    radii["__T"]["C5M"] = 1.80;
    radii["__T"]["C6"] = 1.80;
    radii["__U"]["P"] = 1.90;
    radii["__U"]["O1P"] = 1.40;
    radii["__U"]["O2P"] = 1.40;
    radii["__U"]["O5*"] = 1.40;
    radii["__U"]["C5*"] = 1.80;
    radii["__U"]["C4*"] = 1.80;
    radii["__U"]["O4*"] = 1.40;
    radii["__U"]["C3*"] = 1.80;
    radii["__U"]["O3*"] = 1.40;
    radii["__U"]["C2*"] = 1.80;
    radii["__U"]["O2*"] = 1.40;
    radii["__U"]["C1*"] = 1.80;
    radii["__U"]["N1"] = 1.60;
    radii["__U"]["C2"] = 1.80;
    radii["__U"]["O2"] = 1.40;
    radii["__U"]["N3"] = 1.60;
    radii["__U"]["C4"] = 1.80;
    radii["__U"]["O4"] = 1.40;
    radii["__U"]["C5"] = 1.80;
    radii["__U"]["C6"] = 1.80;
    radii["2MG"]["P"] = 1.90;
    radii["2MG"]["O1P"] = 1.40;
    radii["2MG"]["O2P"] = 1.40;
    radii["2MG"]["O5*"] = 1.40;
    radii["2MG"]["C5*"] = 1.80;
    radii["2MG"]["C4*"] = 1.80;
    radii["2MG"]["O4*"] = 1.40;
    radii["2MG"]["C3*"] = 1.80;
    radii["2MG"]["O3*"] = 1.40;
    radii["2MG"]["C2*"] = 1.80;
    radii["2MG"]["O2*"] = 1.40;
    radii["2MG"]["C1*"] = 1.80;
    radii["2MG"]["N9"] = 1.60;
    radii["2MG"]["C8"] = 1.80;
    radii["2MG"]["N7"] = 1.60;
    radii["2MG"]["C5"] = 1.80;
    radii["2MG"]["C6"] = 1.80;
    radii["2MG"]["O6"] = 1.40;
    radii["2MG"]["N1"] = 1.60;
    radii["2MG"]["C2"] = 1.80;
    radii["2MG"]["N2"] = 1.60;
    radii["2MG"]["C2A"] = 1.80;
    radii["2MG"]["N3"] = 1.60;
    radii["2MG"]["C4"] = 1.80;
    radii["H2U"]["P"] = 1.90;
    radii["H2U"]["O1P"] = 1.40;
    radii["H2U"]["O2P"] = 1.40;
    radii["H2U"]["O5*"] = 1.40;
    radii["H2U"]["C5*"] = 1.80;
    radii["H2U"]["C4*"] = 1.80;
    radii["H2U"]["O4*"] = 1.40;
    radii["H2U"]["C3*"] = 1.80;
    radii["H2U"]["O3*"] = 1.40;
    radii["H2U"]["C2*"] = 1.80;
    radii["H2U"]["O2*"] = 1.40;
    radii["H2U"]["C1*"] = 1.80;
    radii["H2U"]["N1"] = 1.60;
    radii["H2U"]["C2"] = 1.80;
    radii["H2U"]["O2"] = 1.40;
    radii["H2U"]["N3"] = 1.60;
    radii["H2U"]["C4"] = 1.80;
    radii["H2U"]["O4"] = 1.40;
    radii["H2U"]["C5"] = 1.80;
    radii["H2U"]["C6"] = 1.80;
    radii["HEM"]["FE"] = 1.47;
    radii["HEM"]["CHA"] = 2.00;
    radii["HEM"]["CHB"] = 2.00;
    radii["HEM"]["CHC"] = 2.00;
    radii["HEM"]["CHD"] = 2.00;
    radii["HEM"]["N A"] = 1.55;
    radii["HEM"]["C1A"] = 1.78;
    radii["HEM"]["C2A"] = 1.78;
    radii["HEM"]["C3A"] = 1.78;
    radii["HEM"]["C4A"] = 1.78;
    radii["HEM"]["CMA"] = 1.90;
    radii["HEM"]["CAA"] = 1.90;
    radii["HEM"]["CBA"] = 1.90;
    radii["HEM"]["CGA"] = 1.90;
    radii["HEM"]["N B"] = 1.55;
    radii["HEM"]["C1B"] = 1.78;
    radii["HEM"]["C2B"] = 1.78;
    radii["HEM"]["C3B"] = 1.78;
    radii["HEM"]["C4B"] = 1.78;
    radii["HEM"]["CMB"] = 1.90;
    radii["HEM"]["CAB"] = 1.90;
    radii["HEM"]["CBB"] = 1.90;
    radii["HEM"]["N C"] = 1.55;
    radii["HEM"]["C1C"] = 1.78;
    radii["HEM"]["C2C"] = 1.78;
    radii["HEM"]["C3C"] = 1.78;
    radii["HEM"]["C4C"] = 1.78;
    radii["HEM"]["CMC"] = 1.90;
    radii["HEM"]["CAC"] = 1.90;
    radii["HEM"]["CBC"] = 1.90;
    radii["HEM"]["N D"] = 1.55;
    radii["HEM"]["C1D"] = 1.78;
    radii["HEM"]["C2D"] = 1.78;
    radii["HEM"]["C3D"] = 1.78;
    radii["HEM"]["C4D"] = 1.78;
    radii["HEM"]["CMD"] = 1.90;
    radii["HEM"]["CAD"] = 1.90;
    radii["HEM"]["CBD"] = 1.90;
    radii["HEM"]["CGD"] = 1.90;
    radii["HEM"]["O1A"] = 1.35;
    radii["HEM"]["O2A"] = 1.35;
    radii["HEM"]["O1D"] = 1.35;
    radii["HEM"]["O2D"] = 1.35;
    radii["HOH"]["O"] = 1.40;

    for (auto it : radii) for (auto it2 : it.second) if (it2.second > max_radius) max_radius = it2.second;
    return true;
}

/* --- --- --- --- --- vdwContacts --- --- --- --- --- */

vdwContacts::vdwContacts(vector<Residue*> S_res) {
    // By default set all residues to be checked
    setResidues(S_res,S_res);

    // Add each residue in the structure to the proximity search
    preparePS(S_res);
}

vdwContacts::vdwContacts(vector<Chain*> resIChains, vector<Chain*> resJChains) {
    vector<Residue*> resIVec, resJVec;
    for (Chain* C : resIChains) {
        vector<Residue*> cRes = C->getResidues();
        resIVec.insert(resIVec.end(),cRes.begin(),cRes.end());
    }
    for (Chain* C : resJChains) {
        vector<Residue*> cRes = C->getResidues();
        resJVec.insert(resJVec.end(),cRes.begin(),cRes.end());
    }
    vector<Residue*> checkForOverlap = MstUtils::setintersect(resIVec,resJVec);
    if (!checkForOverlap.empty()) MstUtils::error("When providing two sets of residues, both sets must not have overlapping members","vdwContacts::vdwContacts");
    setResidues(resIVec,resJVec);
    preparePS();
}

vdwContacts::vdwContacts(vector<Residue*> resIVec, vector<Residue*> resJVec) {
    vector<Residue*> checkForOverlap = MstUtils::setintersect(resIVec,resJVec);
    if (!checkForOverlap.empty()) MstUtils::error("When providing two sets of residues, both sets must not have overlapping members","vdwContacts::vdwContacts");
    setResidues(resIVec,resJVec);
    preparePS();
}

void vdwContacts::setResidues(vector<Residue*> resIVec, vector<Residue*> resJVec) {
    resISet = set<Residue*>(resIVec.begin(),resIVec.end());
    resJSet = set<Residue*>(resJVec.begin(),resJVec.end());
    cout << "set I has " << resISet.size() << " and set J has " << resJSet.size() << endl;
}

void vdwContacts::preparePS(vector<Residue*> toAdd) {
    if (!toAdd.empty()) {
        for (Residue* R : toAdd) for (Atom* A : R->getAtoms()) allAtoms.push_back(A); 
        allAtomsPS = ProximitySearch(allAtoms,1.5);
    } else {
        // for (Residue* R : resISet) for (Atom* A : R->getAtoms()) allAtoms.push_back(A); 
        for (Residue* R : resJSet) for (Atom* A : R->getAtoms()) allAtoms.push_back(A); 
        allAtomsPS = ProximitySearch(allAtoms,1.5);
    }
}

set<Residue*> vdwContacts::getInteractingRes(Residue* Ri) {
    vector<Atom*> resAtoms = Ri->getAtoms();
    set<Residue*> contacts;
    for (Atom* Ai : resAtoms) {
        mstreal A_maxRadius_sum = vdwR.getRadius(Ai,strictAtomName) + vdwR.maxRadius();
        vector<int> atomsToCheck = allAtomsPS.getPointsWithin(Ai->getCoor(), 0, A_maxRadius_sum);
        for (int atomIdx : atomsToCheck) {
            Atom* Aj = allAtoms[atomIdx];
            Residue* Rj = Aj->getResidue();
            
            if (!isContact(Ri,Rj)) continue;
            
            // Check if contacting, given the Aj radius
            atomInteractionType intType = vdwR.interactionType(Ai,Aj);
            if (intType == ATOMCONTACT) {
                contacts.insert(Rj);
            }
        }
    }
    return contacts;
}

vector<pair<Residue*,Residue*>> vdwContacts::getInteractingRes() {
    vector<pair<Residue*,Residue*>> result;
    int numContacts = 0;
    for (Residue* Ri : resISet) {
        set<Residue*> interactingRes;
        for (Residue* Rj : getInteractingRes(Ri)) {
            interactingRes.insert(Rj);
            numContacts++;
        }
        for (Residue* Rj : interactingRes) {
            result.emplace_back(pair<Residue*,Residue*>(Ri,Rj));
        }
    }
    // if (verbose) cout << "In total, selected residues have " << numContacts << " VDW contacts" << endl;
    return result;
}

map<int,set<int> > vdwContacts::getAllInteractingRes(bool verbose) {
    map<int,set<int> > allInteracting;
    int numContacts = 0;
    for (Residue* Ri : resISet) {
        // Check if map is already populated
        if (interacting.count(Ri) == 1) continue;
        // If not, compute the contacts anew
        set<int> interactingResIdx;
        for (Residue* Rj : getInteractingRes(Ri)) {
            interactingResIdx.insert(Rj->getResidueIndex());
            numContacts++;
        }
        allInteracting[Ri->getResidueIndex()] = interactingResIdx;
    }
    if (verbose) cout << "In total, structure has " << numContacts << " VDW contacts" << endl;
    return allInteracting;
}

bool vdwContacts::isContact(Residue* Ri, Residue* Rj) {
    // Check if should ignore, since in the same residue
    if (Ri == Rj) return false;

    // Check that both residues are in the sets under consideration
    if (resISet.find(Ri) == resISet.end()) return false;
    if (resJSet.find(Rj) == resJSet.end()) return false;

    // Check that residues are not too close in the chain
    if (Ri->getParent() == Rj->getParent()) {
        int Ri_pos = Ri->getResidueIndexInChain();
        int Rj_pos = Rj->getResidueIndexInChain();
        bool tooClose = (abs(Ri_pos-Rj_pos) <= ignoreDistance);
        if (tooClose) return false;
    }
    return true;
}

/* --- --- --- --- --- potentialContacts --- --- --- --- --- */

potentialContacts::potentialContacts(vector<Residue*> _targetResidues, vector<Residue*> _binderResidues, bool _sameChain, bool strict) {
    targetResidues = _targetResidues;
    binderResidues = _binderResidues;
    sameChain = _sameChain;

    for (Residue* R : targetResidues) {
        targetCA.push_back(R->findAtom("CA"));
        string aa3 = R->getName();
        if (!SeqToolsExtension::AAinSet(aa3)) {
            aa3 = SeqToolsExtension::findEquivalentResidueInAlphabet(aa3,strict);
            if (aa3 == "") continue;
        }
        res_t aaIdx = SeqTools::aaToIdx(aa3);
        targetResidueAAIdentity[R] = aaIdx;
    }
    for (Residue* R : binderResidues) binderCA.push_back(R->findAtom("CA",strict));
}

vector<pair<Residue*,Residue*>> potentialContacts::getContacts(bool simpleDistanceCheck) {
    vector<pair<Residue*,Residue*>> contacts;
    for (Atom* targetA : targetCA) {
        CartesianPoint targetCoor = targetA->getCoor();
        Residue* targetRes = targetA->getResidue();
        for (Atom* binderA : binderCA) {
            CartesianPoint binderCoor = binderA->getCoor();
            mstreal distance = targetCoor.distance(binderCoor);
            Residue* binderRes = binderA->getResidue();
            if (distance > defaultDistanceToCheck) continue;
            if (sameChain && ignoreContact(targetRes,binderRes)) continue; 
            else if (simpleDistanceCheck) contacts.push_back(pair<Residue*,Residue*>(targetRes,binderRes));
            else {
                cout << "Not defined!" << endl;
            }
        }
    }
    return contacts;
}

CartesianPoint potentialContacts::getCbFromRes(Residue* R) {
    // Place Cb given spherical coordinates
    CartesianPoint Cb;
    Cb.setPositionBySphericalCoordinates(radius,polarAngle*(M_PI/180),azimuthalAngle*(M_PI/180));

    // Transform point to local coordinate frame
    residueFrame rF(R);
    Transform globalToLocal = TransformFactory::switchFrames(rF,Frame());
    globalToLocal.apply(Cb);

    // Return vector corresponding to Cb position
    return Cb;
}

mstreal potentialContacts::getNormalizedCbDistance(Residue* Ri, Residue* Rj) {
    // Get all vectors
    CartesianPoint RiCa = Ri->findAtom("CA")->getCoor();
    CartesianPoint RjCa = Rj->findAtom("CA")->getCoor();
    CartesianPoint RiCb = getCbFromRes(Ri);
    CartesianPoint RjCb = getCbFromRes(Rj);

    // Compute distance between Ca and Cb atoms
    mstreal CaDistance = RiCa.distance(RjCa);
    mstreal CbDistance = RiCb.distance(RjCb);

    // Compute the bounds (min and max Cb distance possible)
    // note: does not assume CaCb vector length
    mstreal RiCaCbR = RiCb.distance(RiCa);
    mstreal RjCaCbR = RjCb.distance(RjCa);
    mstreal minDistance = max(CaDistance - RiCaCbR - RjCaCbR, 0.0);
    mstreal maxDistance = CaDistance + RiCaCbR + RjCaCbR;

    return (CbDistance - minDistance) / (maxDistance - minDistance);
}

// void potentialContacts::writeBackboneAndCb(string name) {
//     vector
//     for ()
// }

bool potentialContacts::ignoreContact(Residue* Ri, Residue* Rj) {
    // Check if should ignore, because is the same residue
    if (Ri == Rj) return true;

    // Check that residues are not too close in the chain
    if (Ri->getParent() == Rj->getParent()) {
        int Ri_pos = Ri->getResidueIndexInChain();
        int Rj_pos = Rj->getResidueIndexInChain();
        bool tooClose = (abs(Ri_pos-Rj_pos) <= ignoreDistance);
        if (tooClose) return true;
    }
    return false;
}