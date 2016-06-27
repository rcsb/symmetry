package writers;

import java.io.IOException;

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

	@Override
	public synchronized void writeResult(String identifier, QuatSymmetryResults result)
			throws IOException {

		writer.println(identifier + "\t" + result.getSubunits().getSubunitCount() + "\t"
				+ result.getSubunits().getStoichiometry() + "\t" + result.getSymmetry()
				+ "\t" + result.getMethod() + "\t" + result.getScores().getRmsd() + "\t"
				+ result.getScores().getTm());
		writer.flush();
	}

	@Override
	public synchronized void writeHeader() throws IOException {
		writer.println("Name\t" + "Subunits\t" + "Stoichiometry\t" + "Symmetry\t"
				+ "Type\t" + "RMSD\t" + "TMscore");
		writer.flush();
	}

}