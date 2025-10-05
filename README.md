# RelE peptide binder design scripts

This repository contains all of the code used to design extensions of RelB peptide. The majority of the scripts are written in C++, but the scoring and sequence design scripts are written in python. The examples directory shows how the different steps of design can be run.

---

## Installation

### Python Programs

#### Environments

Most of the dependency management is handled with [Poetry](https://python-poetry.org/). To isolate the environment and control the python version, we have used [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main). Note that this environment was created and tested on linux (not on mac)

First create and activate your conda enviroment with python 3.9
```bash
conda create env --name rele_binder_design python==3.9 -y && conda activate rele_binder_design
```

Next install the dependencies using poetry. 
```bash
cd $REPO
poetry install # installs the exact environment as defined by the poetry.lock file

# it is easier to install the torch-geometric dependencies directly with pip
poetry run pip install \
  --find-links https://data.pyg.org/whl/torch-2.4.0%2Bcu124.html \
  torch-scatter torch-sparse torch-cluster pyg-lib
```

#### Running Programs

See examples in `peptide_binder_design/examples/python_programs`

### C++ Programs

#### Dependencies

Before building `interfaceGenerator` programs, you must download the following dependencies

- "Mosaist" or [MST](https://github.com/Grigoryanlab/Mosaist), a library for working with protein structures and sequences.

- [FreeSASA](https://github.com/mittinatten/freesasa), for calculating the solvent accessible surface area of biomolecules. Configure with the following command `./configure --disable-json --disable-xml --disable-threads`. 

- [JSON](https://github.com/nlohmann/json), for working with JSON formatted files.

You will need to build MST according to the instructions provided with the repo. FreeSASA and JSON do not need to be compiled separately.

Edit the `makefile` to provide the paths to their respective installation directories. For example, if each of these has been cloned/installed in the same parent directory containing interfaceGenerator, then the DIR variables would be set as follows:

```makefile
MSTDIR = ../MST
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

## Main programs

### `generateSeeds`

Generate short segments of protein backbone, i.e. 'interface seeds', around the target protein.

### `scoreBinder`

Scores one, or many, binding structures relative to a target protein. The score is computed per interface contact and is equal to the log of the probability of the binder residue, given the target residue.

# Known bugs to fix
- FreeSASA currently crashes when presented with structures containing hydrogens (e.g. rosetta models)
