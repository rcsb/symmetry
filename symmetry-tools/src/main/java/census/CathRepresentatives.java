package census;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import org.biojava.nbio.structure.cath.CathCategory;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathFactory;

/**
 * This script loops through all the CATH domains in the database and extracts a
 * representative of each level (class, architecture, topology or homology).
 * 
 * @author Aleix Lafita
 *
 */
public class CathRepresentatives {

	private final static String output = "resources/cath_H_repr.list";

	public static void main(String[] args) throws IOException {

		// Get an instance of the CATH database
		CathDatabase cath = CathFactory.getCathDatabase();
		List<CathDomain> domains = cath.getByCategory(CathCategory.Homology);

		StringBuffer buf = new StringBuffer();

		for (CathDomain domain : domains)
			buf.append(domain.getDomainName() + "\n");

		File file = new File(output);
		FileWriter writer = new FileWriter(file);
		writer.write(buf.toString());
		writer.close();

	}

}
