# Quaternary Symmetry

QuatSymm is a tool to determine the stoichiometry and symmetry of an assembly. It is able to efficiently handle large complexes, such as icosohedral virus particles.

QuatSymm was initially developed by Peter Rose for the RCSB. It was extended by Aleix Lafita and used during the Critical Assessment of Protein Structure 2016 (CASP 12). If you use QuatSymm, please cite:

> Lafita A, Bliven S, Kryshtafovych A, et al. Assessment of protein assembly prediction in CASP12. Proteins. (2018) 86:247–256. https://doi.org/10.1002/prot.25408

Algorithmic details are given in Aleix Lafita's thesis:

> Lafita A. Assessment of protein assembly prediction in CASP12 & Conformational dynamics of integrin alpha-l domains. ETH Zurich. (2017)


## Installation

QuatSymm requires Java 8 or newer to run.

The latest release can be downloaded from [Github](https://github.com/rcsb/symmetry/releases) and contains the quatsymm jar file and the `runQuatSymm.sh` wrapper script. The script may be called directly or added to your PATH.

## Usage


The easiest way to run QuatSymm is via the included wrapper script.

```
usage:  runQuatSymm.sh [OPTIONS] [structures...]
```

QuatSymm can also run directly through java without using the wrapper script.
From the command line, it can be run as

```bash
java -Xmx2G -jar quatsymm-*.jar [OPTIONS] [structures...]
```

## Options

Options are specified in gnu style. Boolean options can be negated by prefixing with "no".
Common options have short forms.

| Short Option | Long Option | Description
| :----------: | :---------- | :----------
| -h | --help                           | Print usage information
|    | --version                        | Print CE-Symm version
| -i | --input <file>                   | File listing whitespace-delimited query structures
| -q | --noverbose                      | Disable verbose logging information, as well as the default (--simple) output.
| -v | --verbose                        | Output verbose logging information.
| -o | --stats <file>                   | Output a tsv file with detailed symmetry information (default)
| -j | --show3d                         | Force jMol display for each structure [default for <10 structures when specified on command line]
| -J | --noshow3d                       | Disable jMol display [default with --input or for >=10 structures]
|    | --pdbfilepath <dir>              | Download directory for new structures [default tmp folder]. Can also be set with the PDB_DIR environmental variable.
|    | --threads <arg>                  | Number of threads [default cores-1] --minSeqLen <int>                The minimum subunit length to be considered for clustering and symmetry analysis (default: 20)
|    | --minSeqId <float>               | Sequence identity threshold to consider for the sequence subunit clustering. Two subunits with sequence identity equal or higher than the threshold will be clustered together (range: [0,1], default: 0.95)
|    | --minSequenceCoverage <float>    | The minimum coverage of the sequence alignment between two subunits to be clustered together (range: [0,1], default: 0.9)
|    | --minStructureCoverage <float>   | The minimum coverage of the structure alignment between two subunits to be clustered together (range: [0,1], default: 0.9)
|    | --maxClustRmsd <float>           | Structure similarity threshold (measured with RMSD) to consider for the structural subunit clustering.  (default: 3.0 A)
|    | --clustMethod <str>              | Method to cluster the subunits: SEQUENCE, STRUCTURE, or SEQUENCE_STRUCTURE(default: STRUCTURE)
|    | --maxSymmRmsd <float>            | Structure similarity threshold (measured with RMSD) in the symmetry detection (default: 7.0 A)


## Examples

Analyze an assembly and display the results using Jmol:

```bash
runQuatSymm.sh 1vym
```

To run QuatSymm on many structures, specify input (`-i`) and output (`-o`) files. The input
should contain a list of whitespace separated structure identifiers. Lines
beginning with '#' are ignored. Note that using <tt>--input</tt> disables Jmol
display by default. Use <tt>--jmol</tt> to override this behavior.

```bash
runQuatSymm.sh -i queries.txt -o output.tsv
```
An important feature in the batch mode is the multithreading implementation.
By *default*, all the allowed CPUs are used to run the symmetry analyzes in 
parallel. The option `-threads N` can be used to set the number N of threads
to use.

Other examples:
```
# Pseudosymmetry - hemoglobin alpha+beta
runQuatSymm.sh 4hhb
# Local symmetry - ABC transporter BtuCD+F
runQuatSymm.sh 4fi3
```

## Structure Names

QuatSymm accepts a wide variety of ways to specify structures. Some examples:
- PDB Code: `1TGH`
- Single chain: `1ITB.A`
- Residue Range: `2NWX.A:255-416`
  - Multiple ranges can be concatenated by commas if needed
- SCOP/CATH/ECOD/PDP domains: `d1u6dx_`, `e2j5aA1`
- Files: `file.cif`, `file:///path/to/file.pdb?chainId=A`
- URLs: `http://files.rcsb.org/download/1TIM.pdb`
- Biological Assemblies: `BIO:3HDP:2`
  - Biological Assembly support is currently experimental. Some visualization features may not behave properly.

## Output

The output is tab-delimited with the following fields:

* __Name__: Name of the structure
* __Size__: Total number of subunits
* __Subunits__: List of subunits included in the symmetry
* __Stoichiometry__: Stoichiometry of the assembly. Letters do not correspond to chain names, but are assigned to entities starting with the most common entity.
* __Pseudostoichiometry__: Whether some non-identical entities where considered as equivalent
* __Symmetry__: Point group (Cn, Dn, T, O, or I)
* __Local__: Indicates that some subunits were from the reported symmetry 
* __Method__: Distinguishes which of the internal methods was used: NO_ROTATION (C1),
	C2_ROTATION (C2), ROTATION (Cn), or ROTO_TRANSLATION (helical).
* __SymmRMSD__: RMSD over all symmetric transformations
* __SymmTMscore__: TM-score over all symmetric transformations

## License & Availability

QuatSymm is licensed under LGPL 2.1 (see LICENSE, which should have come bundled
with the executable).

The source code is available at https://github.com/rcsb/symmetry and as part
of the [BioJava library](https://github.com/biojava/biojava).


--------------------
[Home](../README.md)
