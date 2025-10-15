# RelE peptide binder design scripts

This repository contains all of the code used to design extensions of RelB peptide. The majority of the scripts are written in C++, but the scoring and sequence design scripts are written in python. The examples directory shows how the different steps of design can be run.

---

## Installation

### Python Programs

#### Environments

Most of the dependency management is handled with [Poetry](https://python-poetry.org/). To isolate the environment and control the python version, we have used [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main). Note that this environment was created and tested on linux (not on mac).

First create and activate your conda enviroment with python 3.9
```bash
conda create env --name rele_binder_design python==3.9 -y && conda activate rele_binder_design
```

Next install the dependencies using poetry. 
```bash
cd $REPO
poetry install # installs the exact environment as defined by the poetry.lock file
```

Note that this environment contains CPU-only pytorch. A legacy conda environment yaml `rele_binder_design.yml` with GPU-enabled torch is also provided, but may be difficult to install on newer GPUs as the pytorch/pytorch geometric versions are old.

### C++ Programs

#### Dependencies

All C++ specific source code can be found in `/fragment_tools_cpp`. Before building these programs, first set up the following dependencies
- "Mosaist" or [MST](https://github.com/swanss/Mosaist), a library for working with protein structures and sequences.

- [FreeSASA](https://github.com/mittinatten/freesasa/releases/tag/2.0.3), for calculating the solvent accessible surface area of biomolecules. Configure with the following command `./configure --disable-json --disable-xml --disable-threads`. 

- [JSON](https://github.com/nlohmann/json), for working with JSON formatted files.

You will need to build Mosaist and FreeSASA according to the instructions provided with the repo. JSON does not need to be compiled.

Edit the `makefile` to provide the paths to the respective installation directories. For example, if each of these has been cloned/installed in the same parent directory containing `/rele_binder_design`, then the DIR variables would be set as follows:

```makefile
MSTDIR = ../Mosaist
SASADIR = ../freesasa-2.0.3
JSONDIR = ../json
```

#### Building the programs

The following compilers have been used to build these programs.

- Apple clang version 11.0.0 (clang-1100.0.33.17)

- gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-36)

Use the following commands to build

- `make all` - builds all programs
- `make test` - builds only programs in the `tests` directory
- `make bin/[executable name]` - builds the specific executable with its dependencies
- `make clean` - removes build intermediates and products

---

# Known bugs
- FreeSASA currently crashes when presented with structures containing hydrogens (e.g. rosetta models)
