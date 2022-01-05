RCSB Symmetry
-------------

## CE-Symm 2.2.1 (symmetry-2.2.1)

Released 5 January 2022

Bug Fixes:

- Upgrade to BioJava 6.0.4

## CE-Symm 2.2.0 (symmetry-2.2.0)

Released 4 January 2022.

This version requires Java 8 or newer, and is tested up to Java 17.

Bug Fixes:

- Fix Log4j vulnerabilities (CVE-2021-44228 and related)
- Fix errors when using RCSB services (e.g. fetching domains)
- Fix issue with Jmol not displaying correctly in european locales
- Fix issues with older Java versions due to various dependencies (#104)
- Upgrade to BioJava 6.0.3

Removed Features:

- Remove ScanSymmetry class

## CE-Symm 2.1.0 (symmetry-2.1.0)

Released 9 January 2019.

Includes both CE-Symm and QuatSymm.

This version requires Java 8 or newer.

New features:

- Use BioJava 5.2.0, which brings lots of improvements:
  - Support Java 8-11 (#97)
  - Improved structure parsing (e.g. better symmetry operator support, mmCIF & MMTF parsing, better distinction between asymm_id and auth_id for identifying chains, etc)
  - Newer Jmol version
- Better documentation (#96)
- Separate QuatSymm package


## CE-Symm 2.0.0 (symmetry-2.0.0)

Released 3 July 2018 after a long period of beta release candidates.

This version requires Java 7-8.

New features:

- Multiple alignment between all repeats of the protein
- Detection of multiple axes of symmetry (both dihedral point groups and hierarchical symmetry)
- Detection and visualization of point groups
- Better symmetry detection
- Monte-Carlo alignment optimization
- GUI/visualization improvements
- Major command line option changes
- Based on BioJava 4.2.12
  - Diverse input formats (pdb, mmcif, local & remote files, chain & residue selections)
  - Numerous bug fixes & improvements


## CE-Symm 1.0.0 (symmetry-tools-1.0.0)

Released 18 August 2014

This is the initial release of the CE-Symm tool. It corresponds to the version
used in Myers-Turnbull (2014).

This version runs on Java 6-8.

Note that the git history at this time is convoluted due to the merging of
several projects. The similarly-named symmetry-1.0.0 tag contains only
quaternary symmetry code, rather than CE-Symm. Further reorganizations followed
the migration of stable quaternary symmetry and ce-symm code into BioJava.