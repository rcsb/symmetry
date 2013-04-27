### Build Status
![](https://travis-ci.org/rcsb/symmetry.png)


# Symmetry in Protein Structures

> Detect, analyze, and visualize **protein symmetry**
 - Detect symmetry in **biological assemblies**
 - Detect **internal pseudo-symmetry** in protein chains
 - Visualize results in [Jmol](http://www.jmol.org)

![PDB ID 1G63](https://raw.github.com/rcsb/symmetry/master/docu/img/1G63.jpg)
 
## Dependencies

This library depends on 

- [BioJava](http://www.biojava.org)
- [Jmol](http://www.jmol.org)

## Analysis of whole PDB

- View the results when using this code base to [analyze all of the PDB](http://www.rcsb.org/pdb/browse/stoichiometry.do)

## Getting Started

 - Symmetry in Biological Assemblies
  View the symmetry of the biological assembly of ![PDB ID 1a0s](https://raw.github.com/rcsb/symmetry/master/docu/img/BioAssemblySymmetryScreenshot1a0s.png)

View [DemoOrientBioAssembly.java](https://github.com/rcsb/symmetry/blob/master/src/main/java/demo/DemoOrientBioAssembly.java)
 
 - Intra-chain symmetry (algorithm name: CE-SYMM)
   View the internal pseudo-symmetry for ![SCOP ID d1jlya1](https://raw.github.com/rcsb/symmetry/master/docu/img/CeSymmScreenshotd1jlya1.png)

View [DemoCeSymm.java](https://github.com/rcsb/symmetry/blob/master/src/main/java/demo/DemoCeSymm.java)



