#include "alignframes.h"

void alignFrames::findInteractingRes()
{
    if (aa == SeqTools::unknownIdx())
        MstUtils::error("Must set amino acid before finding interacting residues", "alignFrames::findInteractingRes");
    allInteractingRes.clear();
    for (int target_idx = 0; target_idx < targets.size(); target_idx++)
    {
        augmentedStructure *aS = targets[target_idx];
        for (int res_i = 0; res_i < aS->residueSize(); res_i++)
        {
            res_t res_i_aa = SeqTools::aaToIdx(aS->getResidue(res_i).getName());
            if (res_i_aa != aa)
            {
                continue;
            }
            const set<int> &interacting_res = vdwContacts[target_idx][res_i];
            for (int res_j : interacting_res)
                allInteractingRes.emplace_back(target_idx, res_i, res_j);
        }
    }
}

Structure *alignFrames::getAlignedInteractingRes(int i)
{
    if ((i < 0) || (i >= allInteractingRes.size()))
        MstUtils::error("Provided value " + MstUtils::toString(i) + " is out of range: (0," + MstUtils::toString(allInteractingRes.size() - 1), "alignFrames::getAlignedInteractingRes");
    const interactingRes &interaction = allInteractingRes[i];

    // get the residues from the target structure
    Residue *rI = &targets[interaction.target]->getResidue(interaction.res_i);
    Residue *rJ = &targets[interaction.target]->getResidue(interaction.res_j);

    Structure *interactingResiduePair = constructStructureFromResiduePair(rI, rJ, interaction);

    // get the transformation that aligns residue i to the reference frame
    residueFrame *rFi = targets[interaction.target]->getResidueFrame(interaction.res_i);
    residueFrame *rFf = targets[interaction.target]->getResidueFrame(interaction.res_j);
    Transform rFi_to_ref = TransformFactory::switchFrames(rFrame, *rFi);
    // Transform rFi_to_ref = TransformFactory::switchFrames(*rFi,rFrame);

    // apply the transformation to the selected residues
    rFi_to_ref.apply(interactingResiduePair);

    return interactingResiduePair;
};

void alignFrames::writeAlignedInteractingResToPDB(string pdbPath)
{
    // open file
    fstream pdb_out;
    MstUtils::openFile(pdb_out, pdbPath, fstream::out, "alignFrames::writeAlignedInteractingResToPDB");

    // transform interacting residues
    for (int i = 0; i < getNumInteracting(); i++)
    {
        Structure *interactingResiduePair = getAlignedInteractingRes(i);

        // write to file
        pdb_out << "HEADER    " << interactingResiduePair->getName() << endl;
        interactingResiduePair->writePDB(pdb_out);

        delete interactingResiduePair;
    }
    pdb_out.close();
}

void alignFrames::addTarget(augmentedStructure *S)
{
    targets.push_back(S);
}

void alignFrames::addTarget(const augmentedStructure &S)
{
    augmentedStructure *newTarget = new augmentedStructure(S);
    targets.push_back(newTarget);
}

void alignFrames::writeDBFile(string dbPath)
{
    fstream ofs;
    MstUtils::openFile(ofs, dbPath, fstream::out | fstream::binary, "alignFrames::writeDBFile");
    MstUtils::writeBin(ofs, 'V');
    MstUtils::writeBin(ofs, (int)1); // format version
    for (int ti = 0; ti < targets.size(); ti++)
    {
        if (targets[ti] == NULL) MstUtils::error("cannot write a database, in which full structures are not populated", "alignFrames::writeDBFile");
        MstUtils::writeBin(ofs, 'S'); // marks the start of a structure section
        targets[ti]->writeData(ofs);

        if (vdwContacts.find(ti) != vdwContacts.end())
        {
            map<int, set<int>> &vals = vdwContacts[ti];
            MstUtils::writeBin(ofs, 'I'); // marks the start of a residue pair interaction property section
            MstUtils::writeBin(ofs, (string)"vdw");
            MstUtils::writeBin(ofs, (int)vals.size());
            for (auto i = vals.begin(); i != vals.end(); ++i)
            {
                MstUtils::writeBin(ofs, (int)i->first);
                MstUtils::writeBin(ofs, (int)(i->second).size());
                for (int j : i->second) 
                {
                    MstUtils::writeBin(ofs, (int)j);
                }
            }
        }
    }
    ofs.close();
}

void alignFrames::readDBFile(string dbPath)
{
    fstream ifs;
    MstUtils::openFile(ifs, dbPath, fstream::in | fstream::binary, "alignFrames::readDBFile");
    char sect;
    string name;
    int ver = 0;
    int ti = numTargets();
    MstUtils::readBin(ifs, sect);
    if (sect == 'V')
    {
        MstUtils::readBin(ifs, ver);
        MstUtils::readBin(ifs, sect);
    }
    if (sect != 'S')
        MstUtils::error("first section must be a structure one, while reading database file " + dbPath, "alignFrames::readDBFile");
    while (ifs.peek() != EOF)
    {
        Structure targetStruct;
        streampos loc = ifs.tellg();
        targetStruct.readData(ifs);
        augmentedStructure *targetAugS = new augmentedStructure(targetStruct);
        // targetSource.push_back(targetInfo(dbFile, targetFileType::BINDATABASE, loc, memSave));
        int L = targetStruct.residueSize();
        addTarget(targetAugS);
        while (ifs.peek() != EOF)
        {
            MstUtils::readBin(ifs, sect);
            if (sect == 'I')
            {
                MstUtils::readBin(ifs, name); // ignore name for now, but could use in the future
                map<int, set<int>> &vals = vdwContacts[ti];
                int ri, rj, N, n;
                mstreal cd;
                MstUtils::readBin(ifs, N);
                for (int i = 0; i < N; i++)
                {
                    MstUtils::readBin(ifs, ri);
                    MstUtils::readBin(ifs, n);
                    for (int j = 0; j < n; j++)
                    {
                        MstUtils::readBin(ifs, rj);
                        vals[ri].insert(rj);
                    }
                }
            }
            else if (sect == 'S')
            {
                break;
            }
            else
            {
                MstUtils::error("unknown section type" + MstUtils::toString(sect) + ", while reading database file " + dbPath, "alignFrames::readDBFile");
            }
        }
        ti++;
    }
    ifs.close();
    cout << "Loaded database with " << targets.size() << " structure" << endl;
}

Structure *alignFrames::constructStructureFromResiduePair(Residue *Ri, Residue *Rj, const interactingRes &data)
{
    Structure *interactingResiduePair = new Structure;

    // add each residue
    Chain *chainA = new Chain("A", "");
    chainA->insertResidueCopy(Ri);
    bool allow_rename = false;
    interactingResiduePair->appendChain(chainA, allow_rename);

    Chain *chainB = new Chain("B", "");
    chainB->insertResidueCopy(Rj);
    interactingResiduePair->appendChain(chainB, allow_rename);

    string name = aaName + "_" + MstUtils::toString(data.target) + "_" + MstUtils::toString(data.res_i) + "_" + MstUtils::toString(data.res_j);
    interactingResiduePair->setName(name);
    return interactingResiduePair;
}