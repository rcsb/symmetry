package writers;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.cluster.SubunitCluster;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;

/**
 * Writes the QuatSymm multiple alignment in FASTA format in a single
 * file. Different entries are split by the characters '//'.
 * 
 * @author Aleix Lafita
 *
 */
public class QuatSymmFastaWriter extends QuatSymmWriter {
	public QuatSymmFastaWriter(String filename) throws IOException {
		super(filename);
	}

	@Override
	public synchronized void writeResult(String identifier,
			QuatSymmetryResults result) throws StructureException {
		if (result != null ) {
			for (SubunitCluster cluster:result.getSubunitClusters()){
				// There is bug because Structure Identifiers are null - quick fix here
				List<StructureIdentifier> structident = cluster.getSubunits().stream()
						.map(n -> new StructureName(identifier + "_" + n.getName()))
						.collect(Collectors.toList());
				
				MultipleAlignment alignment = cluster.getMultipleAlignment();
				alignment.getEnsemble().setStructureIdentifiers(structident);
				if(alignment != null) {
					writer.write(MultipleAlignmentWriter.toFASTA(alignment));
				}
				writer.println("//");
				writer.flush();
			}
		}
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		// No header for Fasta files
	}
	
	

}