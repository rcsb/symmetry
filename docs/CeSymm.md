CE-Symm
=======

Usage
-----

The easiest way to run CE-Symm is via the included wrapper script.

```
usage:  runCESymm.sh [OPTIONS] [structures...]
```

CE-Symm can also run through java. The main class is <tt>demo.CeSymmMain</tt>, and it requires classpath entries for BioJava and the RCSB symmetry package.

```bash
java -Xmx500M -cp "jars/*" demo.CeSymmMain [OPTIONS] [structures...]
```

Options
-------

Options are specified in gnu style. Boolean options can be negated by prefixing with "no".
All common options have short forms.

Short Option | Long Option | Description
:----------: | :---------- | :----------
-h  | --help            | Print usage information
    | --version         | Print CE-Symm version
-i  | --input=file      | File listing whitespace-delimited query structures
-o  | --xml=file        | Output alignment as XML
    | --fatcat=file     | Output alignment as FATCAT output
    | --ce=file         | Output alignment as CE output
    | --pdb=file        | Output alignment as two-model PDB file
    | --tsv=file        | Output alignment as tab-separated file
-v  | --verbose         | Print detailed output (equivalent to "--tsv=-")
-j  | --show3d          | Force jMol display for each structure [default for <10 structures when specified on command line]
-J  | --noshow3d        | Disable jMol display [default with --input or for >=10 structures]
-t  | --order           | Use TM-Score with order for deciding significance. [default]
-T  | --noorder         | Use TM-Score alone for deciding significance.
    | --ordermethod=Class   | Order detection method. Can be a full class name or a short class name from the org.biojava3.structure.align.symm.order package. [default SequenceFunctionOrderDetector]
    | --pdbfilepath=dir | Download directory for new structures [default temp folder]
    | --pdbdirsplit     | Indicates that --pdbfilepath is split into multiple subdirs, like the ftp site. [default]
    | --nopdbdirsplit   | Indicates that --pdbfilepath should be a single directory.
    | --threads         | Number of threads
    | --scopversion=version | Version of SCOP or SCOPe to use when resolving SCOP identifiers [defaults to latest SCOPe]

Interactive mode
----------------

By default, running CE-Symm will enter interactive mode. It will prompt you
for an input structure and display the result in a jMol window.

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

Batch mode
----------

To run CE-Symm on many structures, specify input and output files. The input
should contain a list of whitespace separated structure identifiers. Lines
beginning with '#' are ignored. Note that using <tt>--input</tt> disables jMol
display by default. Use <tt>--jmol</tt> to override this behavior.

```bash
runCESymm.sh --input=queries.txt --xml=output.xml
```


Examples
--------

```
# Get help
runCESymm.sh -h

# Full interactive mode. Prompts for an input structure & runs CE-Symm,
# displaying the result using jMol
runCESymm.sh

# Run CE-Symm on a single structure or a short list of structures
# and display the results using jMol
runCESymm.sh 1HIV
runCESymm.sh 1GEN.A 1TL2.A 1RI6.A

# Print result summary to standard out, repressing jMol display
runCESymm.sh --nojmol 1GEN.A 1TL2.A 1RI6.A

# Print detailed alignment information
runCESymm.sh --nojmol -v 1GEN.A 1TL2.A 1RI6.A

# Output XML file with full alignments and detailed results
runCESymm.sh --xml=output.xml 1GEN.A 1TL2.A 1RI6.A

# Run in batch mode.
# The input file should list one structure per line.
# Accepts PDB IDs ("4hhb"), SCOP identifiers ("d4hhba_"), or ranges ("4hhb.A","4hhb.A_1-141")
# Lines beginning with '#' are ignored.
runCESymm.sh --input=queries.txt --xml=output.xml

```


--------------------
[Home](../README.md)