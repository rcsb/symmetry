package census;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.ecod.EcodDatabase;
import org.biojava.nbio.structure.ecod.EcodDomain;
import org.biojava.nbio.structure.ecod.EcodFactory;

/**
 * This script loops through all the ECOD domains in the latest database and
 * extracts a representative of a specified level (X, H, T or F groups).
 * <p>
 * Representatives are chosen randomly from the set of Manual representatives,
 * if present, or from the set of automated classifications, otherwise.
 * 
 * @author Aleix Lafita
 *
 */
public class EcodRepresentatives {

	private final static String output = "resources/ecod_F_repr.tsv";

	public static void main(String[] args) throws IOException {

		// Get an instance of the ECOD database
		EcodDatabase ecod = EcodFactory.getEcodDatabase();
		List<EcodDomain> domains = ecod.getAllDomains();

		Map<String, String> repr = new HashMap<String, String>();
		Map<String, Boolean> famSeen = new HashMap<String, Boolean>();

		for (EcodDomain domain : domains) {

			// Check if the family already has a representative
			if (famSeen.get(getFamId(domain)) != null) {
				if (famSeen.get(getFamId(domain)))
					continue; // The family has a manual representative
				else if (!domain.getManual())
					continue; // The domain is not manual, do not override
			}

			repr.put(getFamId(domain), domain.getDomainId());
			famSeen.put(getFamId(domain), domain.getManual());
		}

		StringBuffer buf = new StringBuffer();
		buf.append("Family\tDomain\n");
		for (String key : repr.keySet())
			buf.append(key + "\t" + repr.get(key) + "\n");

		File file = new File(output);
		FileWriter writer = new FileWriter(file);
		writer.write(buf.toString());
		writer.close();

	}

	private static String getFamId(EcodDomain d) {
		return d.getXGroup() + "." + d.getHGroup() + "." + d.getTGroup() + "."
				+ d.getFGroup();
		//return d.getXGroup().toString(); //+ "." + d.getHGroup() + "." + d.getTGroup();
	}

}
