# Symmetry in Biomolecular Structures

Detect, analyze, and visualize **protein symmetry** in biological assemblies (**quaternary symmetry**) and inside single chains or domains (**internal symmetry**). The results can be visualized in [Jmol](http://www.jmol.org).

### Quaternary Symmetry

The code to analyze and visualize the symmetry in biological assemblies was ported to the [biojava repository](http://github.com/biojava/biojava).

![PDB ID 1G63](https://raw.github.com/rcsb/symmetry/master/docu/img/1G63.jpg)

View the results when using this code base to [analyze all of the PDB](http://www.rcsb.org/pdb/browse/stoichiometry.do).

### Internal Symmetry

The stable code to analyze and visualize the internal symmetry of biomolecular structures was also ported to the [biojava repository](http://github.com/biojava/biojava). However, this repository contains code in development for the analysis of internal symmetry (symmetry-core module) and command line binaries of the algorithms (symmetry-tools module).

![1u6d](https://raw.github.com/rcsb/symmetry/master/docu/img/1u6d_symmetry.png)

## Dependencies

This project requires Java 1.6. The library is configured to build using [Apache Maven](http://maven.apache.org/).

The following dependencies should be automatically resolved by maven:

- [BioJava](http://www.biojava.org)
- [Jmol](http://www.jmol.org)

### Build Status
[![Build Status](https://travis-ci.org/rcsb/symmetry.png)](https://travis-ci.org/rcsb/symmetry)
