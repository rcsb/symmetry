#!/bin/bash

# example:

# show help:
#  bash runCESymm.sh -h 

# run interactively & display jmol
#  bash runCESymm.sh

# Print alignment
#  bash runCESymm.sh --nojmol --alignment 1hiv

# Run on several structures (allows PDB ID, SCOP ID, ranges, etc)
#  bash runCESymm.sh --nojmol 1hiv d1ijqa1 1MER.A

# send the arguments to the java app
java -Xmx500M -jar "${project.build.finalName}.jar" "$@"