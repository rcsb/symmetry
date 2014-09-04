#!/bin/bash
# CE-Symm version ${project.version}
#
# Source:   https://github.com/rcsb/symmetry-tools
# Server:   http://source.rcsb.org/jfatcatserver/symmetry.jsp

### EXAMPLES ###

# Get help
#runCESymm.sh -h

# Full interactive mode. Prompts for an input structure & runs CE-Symm,
# displaying the result using jMol
#runCESymm.sh

# Run CE-Symm on a single structure or a short list of structures
# and display the results using jMol
#runCESymm.sh 1HIV
#runCESymm.sh 1GEN.A d1tl2a_ PDP:1RI6Aa
#runCESymm.sh ./myprotein.pdb

# Print result summary to standard out, repressing jMol display
#runCESymm.sh --noshow3d 1GEN.A 1TL2.A 1RI6.A

# Print detailed alignment information
#runCESymm.sh -v 1GEN.A

# Output XML file with full alignments and detailed results
#runCESymm.sh --xml=output.xml -J 1GEN.A 1TL2.A 1RI6.A

# Run in batch mode.
# The input file should list one structure per line.
# Accepts PDB IDs ("4hhb"), SCOP identifiers ("d4hhba_"), or ranges ("4hhb.A","4hhb.A_1-141")
# Lines beginning with '#' are ignored.
#runCESymm.sh --input=queries.txt --xml=output.xml

# Alignments can be output as two-model PDB files.
# The --pdb option can take a directory
#runCESymm.sh --input=queries.txt --pdb=.
# It can also take a format string, where "%s" gets substituted for the structure name.
#runCESymm.sh --input=queries.txt --pdb=%s.cesymm.pdb


### Execute jar ###


# Get the base directory of the argument.
# Can resolve single symlinks if readlink is installed
function scriptdir {
    cd "$(dirname "$1")"
    cd "$(dirname "$(readlink "$1" 2>/dev/null || basename "$1" )")"
    pwd
}
DIR="$(scriptdir "$0" )"
# send the arguments to the java app
java -Xmx500M -jar "$DIR/${project.build.finalName}.jar" "$@"
