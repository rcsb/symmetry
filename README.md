# Symmetry in Biomolecular Structures

This project collects tools to detect, analyze, and visualize **protein symmetry**. This includes the CE-Symm tool for **internal symmetry**, a tool for **quaternary symmetry**, and other experiments relating to symmetry. 

This project serves dual purposes. First, it acts as an incubator for new methods and contains RCSB-specific or other non-production ready code. Stable methods have been moved to the [biojava](http://github.com/biojava/biojava) library, where they are accessible for programatic use. Second, user interfaces are built here following each BioJava release. See [releases](https://github.com/rcsb/symmetry/releases) for the latest CE-Symm executable.

## Tools

User interfaces are available within the [symmetry-tools](https://github.com/rcsb/symmetry/blob/master/symmetry-tools) module.

### Internal Symmetry

<img src="docu/img/1u6d_symmetry.png" align="left" width="150" alt="C6 internal symmetry in PDB:1U6D" title="PDB:1U6D" />

The stable code to analyze and visualize the internal symmetry of biomolecular structures was ported to the [biojava repository](http://github.com/biojava/biojava). However, this repository contains code in development for the analysis of internal symmetry (symmetry-core module) and command line binaries of the algorithms (symmetry-tools module).

See [CE-Symm usage](symmetry-tools/docs/CeSymm.md) for more details.

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

Source code for detecting quaternary symmetry was incubated here until [Jan 15, 2015](https://github.com/rcsb/symmetry/releases/tag/quaternary), when the code to analyze and visualize the symmetry in biological assemblies was ported to the [biojava repository](http://github.com/biojava/biojava).


This code was used for the RCSB PDB's [quaternary symmetry analysis](http://www.rcsb.org/pdb/browse/stoichiometry.do).

## Dependencies

This project requires Java 1.8. The library is configured to build using [Apache Maven](http://maven.apache.org/).

The following dependencies should be automatically resolved by maven:

- [BioJava](http://www.biojava.org)
- [Jmol](http://www.jmol.org)

## Building

To build, run

    mvn package
    
This will build the CE-Symm release package (`symmetry-tools/target/CeSymm-X.X.X.jar`) and jar files containing the other code.

### Build Status
[![Build Status](https://travis-ci.org/rcsb/symmetry.png)](https://travis-ci.org/rcsb/symmetry)
