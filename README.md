[![Build Status](https://travis-ci.org/rcsb/symmetry.png)](https://travis-ci.org/rcsb/symmetry)

# Symmetry in Biomolecular Structures

This project collects tools to detect, analyze, and visualize **protein symmetry**. This includes the CE-Symm tool for **internal symmetry**, a tool for **quaternary symmetry**, and other experiments relating to symmetry. 

This project serves dual purposes. First, it acts as an incubator for new methods and contains RCSB-specific or other non-production ready code. Stable methods have been moved to the [biojava](http://github.com/biojava/biojava) library, where they are accessible for programatic use. Second, user interfaces are built here following each BioJava release. See [releases](https://github.com/rcsb/symmetry/releases) for the latest CE-Symm executable.

## Tools

User interfaces are available within the [symmetry-tools](https://github.com/rcsb/symmetry/blob/master/symmetry-tools) module.

### CE-Symm

<img src="docu/img/1u6d_symmetry.png" align="left" width="150" alt="C6 internal symmetry in PDB:1U6D" title="PDB:1U6D" />

CE-Symm is a tool for detecting internal symmetry in protein structures. CE-Symm 2.0 is able to detect both open and closed symmetry and provide a multiple alignment of all repeats. 

See [CE-Symm documentation](symmetry-tools/docs/CeSymm.md) for more details.

When using CE-Symm, please cite:

> Douglas Myers-Turnbull, Spencer E Bliven, Peter W Rose, Zaid K Aziz, Philippe
> Youkharibache, Philip E Bourne, Andreas Prlic (2013) Systematic detection of
> internal symmetry in proteins using CE-Symm. *J Mol Biol*, 426(11), 2255–2268.
> [PMID 24681267]

or the following preprint about CE-Symm 2.0:

> Spencer E Bliven, Aleix Lafita,Peter W Rose, Guido Capitani, Andreas  Prlić,
> Philip E Bourne. Analyzing the symmetrical arrangement of structural repeats
> in proteins with CE-Symm. Preprint on *BioRxiv*: 
> https://doi.org/10.1101/297960

### Quaternary Symmetry

<img src="docu/img/1G63.jpg" alt="Tetrahedral symmetry of PDB:1G63" title="PDB:1G63" style="float: left; width:150px" align="left" width="150"/>

The QuatSymm tool is used for identifying quaternary symmetry in protein complexes. It reports the symmetry and stoichiometry for a complex. It is also able to tolerate mutations and thereby detect pseudosymmetry in complexes with closely related subunits.

See [QuatSymm documentation](symmetry-tools/docs/QuatSymm.md) for more details.

The QuatSymm algorithms are used by the RCSB Protein Data Bank's [quaternary symmetry analysis](http://www.rcsb.org/pdb/browse/stoichiometry.do). The tool was also utilized for biological assembly assessment during the Critical Assessment of Protein Structure 2016 (CASP 12).


## Dependencies

This project requires Java 1.8 or newer. Release versions of all tools can be run with no other dependencies.

To build from source, the library requires [Apache Maven](http://maven.apache.org/).

Additional dependencies ([BioJava](http://www.biojava.org), [Jmol](http://www.jmol.org), etc) should be automatically resolved by maven.


## Building

To build, run

    mvn package
    
This will build the CE-Symm release packages (e.g. `symmetry-tools/target/cesymm-X.X.X.jar`) and jar files containing the other code.
