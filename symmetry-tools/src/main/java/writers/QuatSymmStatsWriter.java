package writers;

import java.io.IOException;

import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;

/**
 * Writes a stats summary of a Quaternary Symmetry result in TSV format
 * 
 * @author Aleix Lafita
 * 
 */
public class QuatSymmStatsWriter extends QuatSymmWriter {

	public QuatSymmStatsWriter(String filename) throws IOException {
		super(filename);
	}

	/**
	 * Writes a line to the file with the Combined Symmetry results of an entry.
	 * 
	 * @param result
	 * @throws IOException
	 */
	public synchronized void writeResult(QuatSymmetryResults result)
			throws IOException {

		writer.println("identifier" + "\t" + result.getSubunits().getSubunitCount() + "\t"
				+ result.getSubunits().getStoichiometry() + "\t" + result.getSymmetry()
				+ "\t" + result.getMethod() + "\t" + result.getScores().getRmsd() + "\t"
				+ result.getScores().getTm());
		writer.flush();
	}
	
	public synchronized void writeBlank(StructureIdentifier id) {
		writer.println(id.getIdentifier() + "\t0\tNA\tNA\tNA\t0\t0");
		writer.flush();
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		writer.println("Name\t" + "Subunits\t" + "Stoichiometry\t" + "Symmetry\t"
				+ "Type\t" + "RMSD" + "TMscore");
		writer.flush();
	}

}