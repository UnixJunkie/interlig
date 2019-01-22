# README #

InterComp is a software for the sequence-order alignment and superposition of biological molecules.

Currently InterComp supports two input formats:

* PDB
* mol2

## Installation ##

To compile InterComp just clone this repo and run the command:

```
make
```

A binary named InterComp will be generated in the cloned folder.

## Usage ##
To run the binary in the folder where it is installed:

    ./InterfaceComparison [-PDB, -list, -mol2] target1 target2 <options>

    ### Main options: ###
    
    -PDB   : Compare two PDB files
             (./InterfaceComparison -list target1.pdb target2.pdb)
             
    -list  : Compare a PDB file with a list of PDB files (one path per line)
             (./InterfaceComparison -list target.pdb pdb_paths.list)
             
    -mol2  : Compare a mol2 file with a mol2 database (multiple molecules in the database)
             (./InterfaceComparison -mol2 target.mol2 database.mol2)

    ### Other options: ###
    
    -d0           : d0 parameter in the Levitt-Gerstein score (default: 0.5)
    -dW           : weight of sequence similarity (e.g. BLOSUM62) in the optimization procedure (default: 0.5)
    -eps          : weight of structural similarity in the optimization procedure (default: 0.5)
    -nullP        : percentage of ignored atoms in the smaller molecule (default: 0)
    -v            : verbose output
    -anneal <int> : number of annealing rounds (default: 1)
    -seed <int>   : seed for the stochastic process (random number generation)
    -ch <int>     : Markov chain length multiplier (default: 10)\n
    -matrix <path>: optional path to the sequence similarity matrix (e.g. BLOSUM62)
    -super        : calculate optimal superposition between two aligned PDB files (only -PDB mode)
                    (output to molA.pdb and molB.pdb files)
    -apply        : apply the same superposition as in -super to a third molecule (only -PDB mode)
                    (output to molC.pdb file)
    -h            : print this help

## Usage examples ##

When aligning two PDB files with default optimization parameters, InterComp is called with:

```
./InterComp -PDB examples/target1.pdb examples/target2.pdb
```

The order of the two files following the -PDB flag is irrelevant.

In this case, InterComp only outputs the structural and sequence score of the alignment.

If a superposition of the two molecules is needed, use:

```
./InterComp -PDB examples/target1.pdb examples/target2.pdb -super
```

This will generate two PDB files:

* molA.pdb, which contains the CA atoms from the smallest of the two input molecules
* molB.pdb, which contains only those CA atoms from the largest of the input molecules that superimpose to the atoms in molA.pdb. These are also rotated and translated accordingly onto the atoms in molA.pdb.

If we want to run on mol2 files with d0 = 0.1 and eps = 1:

```
./InterfaceComparison -mol2 examples/target.mol2 examples/database.mol2 -d0 0.1 -eps 1
```

The database file here is *always* the second file following the mol2 flag.
