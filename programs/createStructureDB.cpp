#include "mstrotlib.h"
#include "msttypes.h"
#include "mstoptions.h"
#include "mstfasst.h"

#include "alignframes.h"
#include "residuecontact.h"
#include "residueframe.h"

// #include <chrono>

int main(int argc, char *argv[])
{
    MstOptions op;
    op.setTitle("Creates a structure database from input PDB files.");
    op.addOption("pL", "a file with a list of PDB files.");
    op.addOption("db", "a previously-written database.");
    op.addOption("dL", "a file with a list of databases (will consolidate into one).");
    op.addOption("vdw","If provided, will compute the van der Waals contacts between residues");
    op.addOption("o", "output database file name.", true);
    op.addOption("c", "clean up PDB files, so that only protein residues with enough of a backbone to support rotamer building survive.");
    op.addOption("s", "split final PDB files into chains by connectivity. Among other things, this avoids \"gaps\" within chains (where missing residues would go)");
    op.addOption("batch", "an integer. If specified, instead of building the database will spread all the work across this many "
                          "jobs, writing corresponding batch files for submission to the cluster. Will also produce a file called "
                          "<out>.fin.sh (where <out> is the base of the name specified in --o), which is to be run after all jobs "
                          "finish to complete the database building process.");
    op.addOption("p", "the name of the partition to submit slurm batch scripts to. If not provided, will use (default = defq)");

    op.setOptions(argc, argv);

    if (!op.isGiven("batch"))
    {
        proteinFrameDB frameDB;
        cout << "Reading structures..." << endl;
        if (op.isGiven("pL"))
        {
            vector<string> pdbFiles = MstUtils::fileToArray(op.getString("pL"));
            for (int i = 0; i < pdbFiles.size(); i++)
            {
                Structure P(pdbFiles[i]);
                if (op.isGiven("c"))
                {
                    Structure C;
                    RotamerLibrary::extractProtein(C, P);
                    if (P.residueSize() != C.residueSize())
                    {
                        cout << pdbFiles[i] << ", had " << P.residueSize() << " residues, and " << C.residueSize() << " residues after cleaning..." << endl;
                    }
                    C.setName(P.getName());
                    P = C;
                }
                if (op.isGiven("s")) {
                    P = P.reassignChainsByConnectivity();
                    P.deleteShortChains();
                }
                if (P.residueSize() != 0)
                    frameDB.addTarget(P);
                else
                    cout << "skipping " << pdbFiles[i] << " as it ends up having no residues..." << endl;
            }
        }
        if (op.isGiven("db"))
        {
            frameDB.readDBFile(op.getString("db"));
            // If adding to an existing DB, print what properties already exist within it
            cout << "structures DB with " << frameDB.numTargets() << " frames" << endl;
        }
        if (op.isGiven("dL"))
        {
            vector<string> dbFiles = MstUtils::fileToArray(op.getString("dL"));
            for (int i = 0; i < dbFiles.size(); i++)
            {
                frameDB.readDBFile(dbFiles[i]);
            }
        }
        if (op.isGiven("vdw")) {
            cout << "Computing vdW contacts..." << endl;
            // compute and add some properties
            for (int ti = 0; ti < frameDB.numTargets(); ti++)
            {
                cout << "\ttarget " << ti + 1 << "/" << frameDB.numTargets() << "..." << endl;
                const augmentedStructure& P = frameDB.getTarget(ti);
                vdwContacts vdwC(P.getResidues());
                bool verbose = true;
                frameDB.setVDWContacts(ti,vdwC.getAllInteractingRes(verbose));
            }
        }
        frameDB.writeDBFile(op.getString("o")+".db");
    }
    else
    {
        if (!op.isGiven("pL"))
            MstUtils::error("--pL must be given with --batch");
        if (!op.isInt("batch") || (op.getInt("batch") <= 0))
            MstUtils::error("--batch must be a positive integer!");
        int nJ = op.getInt("batch");
        vector<string> pdbFiles = MstUtils::fileToArray(op.getString("pL"));
        srand(time(NULL) + (int)getpid());
        MstUtils::shuffle(pdbFiles);
        vector<pair<int, int>> tasks = MstUtils::splitTasks(pdbFiles.size(), nJ);
        fstream outf;
        string dbListFile = MstSys::pathBase(op.getString("o")) + ".dL";
        fstream dblf;
        MstUtils::openFile(dblf, dbListFile, ios::out);
        vector<string> toClean;
        toClean.push_back(dbListFile);
        for (int i = 0; i < nJ; i++)
        {
            string base = MstSys::pathBase(op.getString("o")) + "." + MstUtils::toString(i);

            // dump subset of PDB files this job will work on
            string listFile = base + ".list";
            MstUtils::openFile(outf, listFile, ios::out);
            for (int k = tasks[i].first; k <= tasks[i].second; k++)
            {
                outf << pdbFiles[k] << endl;
            }
            outf.close();

            // create a job script file that will work on this subset
            string batchFile = base + ".sh";
            string dbFile = base;
            MstUtils::openFile(outf, batchFile, ios::out);
            // time per structure depends on what properties need to be calculated
            int time_per_structure = 1;
            vector<string> allOpts = op.getAllGivenOptions();
            cout << "Given options, calculating that " << time_per_structure << " min are needed per structure (on average)" << endl;
            int hrs = (int)ceil(time_per_structure * (tasks[i].second - tasks[i].first + 1) / 60.0); // fifteen minutes per structure should be plenty
            string partition_name = op.getString("p", "defq");
            outf << "#!/bin/bash\n";
            outf << "#SBATCH -J fasstDB.%A\n"
                 << "#SBATCH -o fasstDB.%A.log\n";
            outf << "#SBATCH -p " << partition_name << "\n"
                 << "#SBATCH -n 1\n"
                 << "#SBATCH --mem=2G\n";
            outf << "#SBATCH -t 0-" << hrs << ":00:00\n";
            outf << op.getExecName() << " --pL " << listFile << " --o " << dbFile;
            // keep all other options from the call to self
            for (int j = 0; j < allOpts.size(); j++)
            {
                if ((allOpts[j].compare("batch") == 0) || (allOpts[j].compare("pL") == 0) || (allOpts[j].compare("o") == 0))
                    continue;
                outf << " --" << allOpts[j] << " " << op.getString(allOpts[j]);
            }
            outf << endl;
            outf.close();
            dblf << dbFile << ".db" << endl;
            toClean.push_back(base + ".sh*");
            toClean.push_back(base);
            toClean.push_back(base + ".list");
        }
        dblf.close();
        fstream fin;
        MstUtils::openFile(fin, "fin." + MstSys::pathBase(op.getString("o")) + ".sh", ios::out);
        fin << op.getExecName() << " --dL " << dbListFile << " --o " << op.getString("o");
        fin << endl;
        fin << "if [ $? -eq 0 ]; then # only clean up if database creation was successful" << endl;
        for (int k = 0; k < toClean.size(); k++)
        {
            fin << "  rm " << toClean[k] << endl;
        }
        fin << "fi" << endl;
        fin.close();
    }
    cout << "Done" << endl; // this makes it easy to check if all the jobs were completed
    return 0;
}