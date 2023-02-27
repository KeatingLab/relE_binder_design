#include "msttypes.h"
#include "mstoptions.h"

#include "alignframes.h"
#include "generateseeds.h"
#include "scoreinterface.h"
#include "residuecontact.h"
#include "residueframe.h"

#include <chrono>

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Finds potential interface contacts");
    op.addOption("frameDB","Path to mobile frame database",true);
    op.addOption("contactData","Path to a JSON file describing probability density of contacts",true);
    op.addOption("complexPDB","A PDB file defining the complex between the protein target and designed binder. The target must have a defined sequence",true);
    op.addOption("binderChains","The binder chain ID(s) delimited by '_', e.g. '0_1'. If in --targetPDB mode, will be removed from the structure before scoring",true);
    op.addOption("targetChains","--complexPDB mode only. The protein chain(s) delimited by '_', e.g. 'A_B'",false);
    op.setOptions(argc,argv);

    if (op.isGiven("complexPDB")) {
        if (!op.isGiven("targetChains") || !op.isGiven("binderChains")) MstUtils::error("If scoring a complex, must provide '--targetChains' and '--binderChains'");
    }
    string mobileFrameDB = op.getString("frameDB");
    string contactData = op.getString("contactData");
    string complexPDB = op.getString("complexPDB","");
    string targetChainIDsString = op.getString("targetChains");
    string binderChainIDsString = op.getString("binderChains");

    binderScorer* scorer = nullptr;

    binderScorerParams params;
    params.frameDBPath = mobileFrameDB;
    params.potentialContactsJSONPath = contactData;
    params.posCut = 0.5;
    params.oriCut = 10;

    augmentedStructure complex(complexPDB);
    scorer = new binderScorer(params,complex,binderChainIDsString,targetChainIDsString);

    scorer->defineTargetBindingSiteResiduesByrSASA();
    vector<Residue*> targetResidues = scorer->getTargetResidues();
    vector<Residue*> binderResidues = scorer->getBinderResidues();

    potentialContacts pConts(targetResidues,binderResidues);
    pConts.load2DProbabilityDensities(contactData);

    // Keep track of contacts with a few different files
    fstream out;
    MstUtils::openFile(out,"allPotentialContacts.csv",fstream::out);
    out << "resIchainID,resInum,resIname,resJchainID,resJnum,resJname,CaDistance,normalizedCbdistance,RiCbOrientation,RjCbOrientation,numBBHeavyAtomsBetweenResidues,potentialContactSS,potentialContactSB,potentialContactBS,potentialContactBB" << endl;

    fstream ss_out;
    MstUtils::openFile(ss_out,"potentialContactsSS.tsv",fstream::out);

    fstream sb_out;
    MstUtils::openFile(sb_out,"potentialContactsSB.tsv",fstream::out);

    fstream bs_out;
    MstUtils::openFile(bs_out,"potentialContactsBS.tsv",fstream::out);

    fstream bb_out;
    MstUtils::openFile(bb_out,"potentialContactsBB.tsv",fstream::out);

    int nSS = 0, nSB = 0, nBS = 0, nBB = 0, nTotal = 0;
    cout << targetResidues.size() << " target residues and " << binderResidues.size() << " binder residues" << endl;
    for (Residue* targetR : targetResidues) {
        for (Residue* binderR : binderResidues) {
            mstreal CaDistance = pConts.getCaDistance(targetR,binderR);
            mstreal getNormalizedCbDistance = pConts.getNormalizedCbDistance(targetR,binderR);
            mstreal RiCbOrientation = pConts.getCaCbtoRiCaRjCaAngle(targetR,binderR);
            mstreal RjCbOrientation = pConts.getCaCbtoRiCaRjCaAngle(binderR,targetR);
            int numBBHeavyAtomsBetweenResidues = pConts.bbHeavyAtomsBetweenResidues(targetR,binderR);

            bool SS = pConts.isPotentialSSContact(targetR,binderR,0.05);
            bool SB = pConts.isPotentialSBContact(targetR,binderR,0.05);
            bool BS = pConts.isPotentialSBContact(binderR,targetR,0.05);
            bool BB = pConts.isPotentialBBContact(targetR,binderR);

            out << targetR->getChainID() << "," << targetR->getNum() << "," << targetR->getName() << ",";
            out << binderR->getChainID() << "," << binderR->getNum() << "," << binderR->getName() << ",";
            out << CaDistance << "," << getNormalizedCbDistance << "," << RiCbOrientation << "," << RjCbOrientation << ",";
            out << numBBHeavyAtomsBetweenResidues << ",";
            out << SS << "," << SB << "," << BS << "," << BB << endl;

            if (SS) {
                ss_out << targetR->getChainID() << targetR->getNum() << "\t" << binderR->getChainID() << binderR->getNum() << "\t" << 1.0 << endl;
                nSS++;
            } if (SB) {
                sb_out << targetR->getChainID() << targetR->getNum() << "\t" << binderR->getChainID() << binderR->getNum() << "\t" << 1.0 << endl;
                nSB++;
            } if (BS) {
                bs_out << targetR->getChainID() << targetR->getNum() << "\t" << binderR->getChainID() << binderR->getNum() << "\t" << 1.0 << endl;
                nBS++;
            } if (BB) {
                bb_out << targetR->getChainID() << targetR->getNum() << "\t" << binderR->getChainID() << binderR->getNum() << "\t" << 1.0 << endl;
                nBB++;
            } if (SS | SB | BS | BB) nTotal++;
        }
    }
    cout << "SS = " << nSS << endl;
    cout << "SB = " << nSB << endl;
    cout << "BS = " << nBS << endl;
    cout << "BB = " << nBB << endl;
    cout << "Total (union) = " << nTotal << endl;

    cout << "Done" << endl;

    return 0;
}