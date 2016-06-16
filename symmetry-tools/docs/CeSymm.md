# CE-Symm

## Usage


The easiest way to run CE-Symm is via the included wrapper script.

```
usage:  runCESymm.sh [OPTIONS] [structures...]
```

CE-Symm can also run directly through java without using the wrapper script.
On many systems, double clicking the JAR file should open the program in
interactive mode. From the command line, it can be run as

```bash
java -Xmx500M -jar cesymm-*.jar [OPTIONS] [structures...]
```

## Options

Options are specified in gnu style. Boolean options can be negated by prefixing with "no".
All common options have short forms.

Short Option | Long Option | Description
:----------: | :---------- | :----------
-h  | --help            | Print usage information
    | --version         | Print CE-Symm version
-i  | --input=file      | File listing whitespace-delimited query structures
-v  | --verbose         | Output verbose logging information.
-q  | --noverbose       | Disable verbose logging information, as well as the default (--simple) output.
-o  | --simple=file     | Output result in a simple format (default)
    | --stats=file      | Output a tsv file with detailed symmetry information
    | --tsv=file        | Output alignment as a tsv-formated list of aligned residues.
    | --xml=file        | Output alignment as XML
    | --fatcat=file     | Output alignment as FATCAT output
    | --fasta=file      | Output alignment as FASTA alignment output
    | --axes=file       | Output each axis of symmetry found in the structure as two points and a rotation angle
-j  | --show3d          | Force Jmol display for each structure [default for <10 structures when specified on command line]
-J  | --noshow3d        | Disable Jmol display [default with --input or for >=10 structures]
    | --ordermethod=Class   | Order detection method: SEQUENCE_FUNCTION (default), GRAPH_COMPONENT, ANGLE, or USER_INPUT
    | --order <int>     | Force a particular order. If positive, implies --ordermethod=USER_INPUT.
    | --refinemethod=Class  | Refiner method: SEQUENCE_FUNCTION (default), NOT_REFINED, or GRAPH_COMPONENT
    | --symmtype=Class      | Restrict symmetry to: CLOSED, OPEN, or AUTO (default)
    | --pdbfilepath=dir | Download directory for new structures [default tmp folder]. Can also be set with the PDB_DIR environmental variable.
    | --threads=int     | Number of threads
    | --maxgapsize=float| This parameter configures the maximum gap size G, that is applied during the AFP extension. The larger the value, the longer the calculation time can become, Default value is 30. Set to 0 for no limit.
    | --scoringstrategy=str |   Which scoring function to use: CA_SCORING, SIDE_CHAIN_SCORING, SIDE_CHAIN_ANGLE_SCORING, CA_AND_SIDE_CHAIN_ANGLE_SCORING, or SEQUENCE_CONSERVATION
    | --winsize=int     | This configures the fragment size m of Aligned Fragment Pairs (AFPs).
    | --maxrmsd=float   | The maximum RMSD at which to stop alignment optimization. (default: unlimited=99)
    | --nointernalgaps  | Force alignment to include a residue from all repeats. (By default only 50% of repeats must be aligned in each column.)
    | --gapopen=float   | Gap opening penalty during alignment optimization [default: 5.0].
    | --gapextension=float  | Gap extension penalty during alignment optimization [default: 0.5].
    | --symmlevels=int  | Run iteratively the algorithm to find multiple symmetry levels. The parameter controls the maximum symmetry levels allowed. 0 means unbounded. [default: 0].
    | --noopt           | Disable optimization of the resulting symmetry alignment.
    | --unrefinedscorethreshold=float | The TM-score threshold before refinement of the self-alignment. TM-scores above this value will be considered significant results [default: 0.4, interval [0.0,1.0]].
    | --refinedscorethreshold=float | The TM-score threshold of the multiple alignment of repeats after refinement. TM-scores above this value will be considered significant results [default: 0.36, interval [0.0,1.0]].
    | --ssethrehold=int | The minimum number of secondary structure elements (SSE) for each symmetric subunit, for the result to be singificant [default: 2].
    | --maxorder=int    | The maximum number of symmetric subunits [default: 8].
    | --rndseed=int     | The random seed used in optimization, for reproducibility of the results [default: 0].
    | --minlen=int      | The minimum length, expressed in number of core aligned residues, of a symmetric subunit [default: 15].
    | --dcutoff=float   | The maximum distance, in A, allowed between any two aligned residue positions [default: 7.0].
    | --scopversion=version | Version of SCOP or SCOPe to use when resolving SCOP identifiers [defaults to latest SCOPe]

## Interactive mode

By default, running CE-Symm will enter interactive mode. It will prompt you
for an input structure and display the result in a Jmol window.

```bash
runCESymm.sh
```

You can also pass one or more structures to the program to bypass the input
dialog. A number of common ways to specify protein structures are supported,
including PDB IDs, SCOP domain identifiers, PDP domains, and filenames. See
BioJava's [StructureID.getStructure()](
http://www.biojava.org/docs/api/org/biojava3/structure/StructureIO.html#getStructure%28java.lang.String%29)
method for the complete syntax.

```bash
runCESymm.sh 1HIV
runCESymm.sh 1GEN.A d1tl2a_ PDP:1RI6Aa
```

## Batch mode

To run CE-Symm on many structures, specify input and output files. The input
should contain a list of whitespace separated structure identifiers. Lines
beginning with '#' are ignored. Note that using <tt>--input</tt> disables Jmol
display by default. Use <tt>--jmol</tt> to override this behavior.

```bash
runCESymm.sh --input=queries.txt --xml=output.xml
```
An important feature in the batch mode is the multithreading implementation.
By *default*, all the allowed CPUs are used to run the symmetry analyses in 
parallel. The option `-threads N` can be used to set the number N of threads
to use. Note that the scaling efficiency of the multithreading is not perfect,
and the **thread overhead** becomes significant with more than 8 threads.

## Specifying Structures

CE-Symm accepts a wide variety of ways to specify structures. Some examples:
- PDB Code: `1TGH`
- Single chain: `1ITB.A`
- Residue Range: `2NWX.A:255-416`
  - Multiple ranges can be concatenated by commas if needed
- SCOP/CATH/ECOD domains: `d1u6dx_`, `e2j5aA1`
- Files: `file.cif`, `file:///path/to/file.pdb?chainId=A`
- URLs: `http://files.rcsb.org/download/1TIM.pdb`
- Biological Assemblies: `BIO:3HDP:2`
  - Biological Assembly support is currently experimental. Some visualization features may not behave properly.

## Output

CE-Symm can output results in a number of formats. Format options may be followed
by a filename (`--stats=out.txt`). If the filename is empty or '-', output will
be sent to the terminal. To prevent structure names from being interpreted as
filenames, it is recommended to use the hyphen explicitly (e.g. `--stats=-`,
although `--stats` may work alone in some contexts).

If no format is specified, CE-Symm will default to printing the simple format to
standard out, although this can be surpressed with the `-q` option.

The following formats are supported. Most formats represent the symmetry as an
alignment from the structure to itself.

* __Simple__: Intended to be a simple human-readable summary of results. Lists
  one line per query, giving the symmetry determination and information about
  how that conclusion was reached.
* __Stats__: Provides more detailed statistics about each result. See below
  for a description of the fields.
* __Axes__: Provides detailed infomation about the axes of symmetry found. 
* See below for a description of the fields.
* __TSV__: A list of aligned residues for each structure,  with a line containing 
  only '//' to separate records.
* __FATCAT__: A single file containing all alignments in FATCAT's traditional output
  format, with a line containing only '//' to separate records.
* __XML__: All the alignments in a custom XML format suitable for machine parsing.


### Statistics Output

The `--stats` option outputs a tab-delimited file, where each row corresponds to an input structure, with the following columns:

- __Name__ Name of the structure
- __NumRepeats__ Total number of repeats detected by CE-Symm, including multiple
  levels of symmetry if detected
- __SymmGroup__ Symmetry Group of the top level of symmetry. This includes
  point group symmetry (Cn or Dn), helical symmetry (H), and translational
  repeats (R).
- __Refined__ 'true' or 'false', indicating whether refinement was successful
- __SymmLevels__ Number of symmetry levels detected
- __SymmType__ OPEN or CLOSED for symmetry or NONE for asymmetry. Multiple levels of symmetry are comma separated, in order from the first level to the last
- __RotationAngle__ Angle of rotation at the principal axis (degrees). Closed symmetry may deviate from ideal values due to the superposition procedure. Multiple axes of symmetry are comma separated, in order from the first symmetry level to the last
- __ScrewTranslation__ Translation parallel the principal axis (Å). Multiple axes of symmetry are comma separated, in order from the first symmetry level to the last
- __UnrefinedTMscore__ TM-Score of the self-alignment prior to refinement
- __UnrefinedRMSD__ RMSD of the self-alignment prior to refinement
- __SymmTMscore__ Average pairwise TM-Score of all repeats in the alignment
  after refinement
- __SymmRMSD__ Average RMSD of all repeats in the alignment after refinement
- __RepeatLength__ Total number of alignment positions in the multiple repeats alignment
- __CoreLength__ Number of _ungapped_ columns in the multiple alignment. Equivalent to (RepeatLength*NumRepeats)-Gaps
- __Length__ Total length of the protein
- __Coverage__ Fraction of the protein aligned
- __Repeats__ Substructure identifiers of each repeat range, semicolon separated. The numbers correspond to residue numbers in the original structure

### Axes Output

Information about the rotation axes can be extracted from the `--axes` output format. This option also outputs a tab-delimited file, where each row corresponds to one axis of symmetry. In the case of simple cyclic or open symmetry, only one axis is output. For more complex cases (dihedral symmetry or axes that don't form a point group), a hierarchy of axes is constructed that align smaller and smaller portions of the structure. Axes with the same `SymmLevel` are related to each other by combinations of lower-level axes. Note that the axes generated by CE-Symm (displayed in yellow in Jmol) are optimized independently, so for dihedral cases they may differ slightly from the blue point-group axes.

- __Name__ Name of the structure
- __SymmLevel__ Level of symmetry of the axis
- __SymmType__ CLOSED or OPEN depending on whether applying the axis to the last repeat superimposes with the first repeat
- __SymmOrder__ Number of repeats in the structure to which the axis applies. For CLOSED symmetry axes, this is equivalent to 360/RotationAngle.
- __RotationAngle__ Angle of rotation at the principal axis (degrees). Closed symmetry may deviate from ideal values due to the superposition procedure.
- __ScrewTranslation__ Translation parallel the principal axis (Å)
- __Point1__ One end of the axis
- __Point2__ Other end of the axis
- __AlignedRepeats__ Defines which repeats are aligned. Each parenthetical group gives a set of aligned repeats. This is similar to the cycle notation used in group theory. For example, `(A;C;E)(B;D;F)` indicates a three-fold rotation superimposing `A` and `B` onto `C` and `D`, `C` and `D` onto `E` and `F`, and (if closed) `E` and `F` onto `A` and `B`.

### XML output

The XML output is primarily suitable for loading results into BioJava. Other formats are generally preferable for common uses. However, a few clarifications are in order for advanced users:

- The `<Structures>` element gives residue ranges for each repeat of the multiple alignment. This information is also available in the `--simple` format.
- The `<eqrX>` elements give the indices of aligned residues for each repeat. Use of the `--tsv` output is strongly recommended rather than parsing these, since tsv uses standard PDB residue numbers.
- The `<Matrix4d>` elements give the 4D transformation matrices required to superimpose each repeat onto the first repeat (thus the first entry is always the identity matrix). It is possible to express these matrices in terms of the rotation axes from the `--axes` output, but doing so requires careful tracking of which repeats are aligned by each axis.

## License & Availability

CE-Symm is licensed under LGPL 2.1 (see LICENSE, which should have come bundled
with the executable).

The source code is available at https://github.com/rcsb/symmetry and as part
of the [BioJava library](https://github.com/biojava/biojava).
A webserver is also provided at http://source.rcsb.org/jfatcatserver/symmetry.jsp

If you use CE-Symm in your research, please cite:

Douglas Myers-Turnbull, Spencer E Bliven, Peter W Rose, Zaid K Aziz, Philippe
  Youkharibache, Philip E Bourne, and Andreas Prlic. Systematic detection of
  internal symmetry in proteins using CE-Symm. 2014. Awaiting publication.

## Examples

```
# Get help
runCESymm.sh -h

# Full interactive mode. Prompts for an input structure & runs CE-Symm,
# displaying the result using Jmol
runCESymm.sh

# Run CE-Symm on a single structure or a short list of structures
# and display the results using Jmol
runCESymm.sh 1HIV
runCESymm.sh 1GEN.A d1tl2a_ PDP:1RI6Aa
runCESymm.sh ./myprotein.pdb

# Print result summary to standard out, repressing Jmol display
runCESymm.sh --noshow3d 1GEN.A 1TL2.A 1RI6.A

# Print detailed alignment information
runCESymm.sh --stats=- 1GEN.A

# Output XML file with full alignments and detailed results
runCESymm.sh --xml=output.xml -J 1GEN.A 1TL2.A 1RI6.A

# Run in batch mode.
# The input file should list one structure per line.
# Accepts PDB IDs ("4hhb"), SCOP identifiers ("d4hhba_"), or ranges ("4hhb.A","4hhb.A_1-141")
# Lines beginning with '#' are ignored.
runCESymm.sh --input=queries.txt --xml=output.xml
```


--------------------
[Home](../README.md)
