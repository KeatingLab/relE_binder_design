# interfaceGenerator

A suite of C++ programs for *de novo* design of protein-binding peptides or mini-proteins. 

---

## Installation

### Python Programs

#### Environments

Using conda or [mamba](https://github.com/mamba-org/mamba), set up the environment for [TERMinator_sscore](https://github.com/swanss/TERMinator_sscore), following the directions in the README.md

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
