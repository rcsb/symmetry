package org.biojava.nbio.structure.align.symm.ecodcensus;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.biojava.nbio.structure.ecod.EcodDatabase;
import org.biojava.nbio.structure.ecod.EcodDomain;
import org.biojava.nbio.structure.ecod.EcodFactory;

public class EcodRepresentatives {

	private final EcodDatabase ecod;
	private static final int DEFAULT_LEVEL = 3; // T-Group
	
	public EcodRepresentatives() {
		this(EcodFactory.getEcodDatabase());
	}
	public EcodRepresentatives(EcodDatabase ecod) {
		super();
		this.ecod = ecod;
	}
	
	public List<EcodDomain> getDomains(int level) throws IOException {
		Map<String,EcodDomain> reps = getRepresentatives(ecod, level);
		return new ArrayList<EcodDomain>(reps.values());
	}
	public List<EcodDomain> getDomains() throws IOException {
		return getDomains(DEFAULT_LEVEL);
	}
	
	/**
	 * Generates a representative for each Ecod topology group
	 * @param ecod
	 * @param level Ecod hierarchy level, from 1 (X-group) through 4 (F-group)
	 * @return
	 * @throws IOException
	 */
	public static Map<String,EcodDomain> getRepresentatives(EcodDatabase ecod,int level) throws IOException {
		SortedMap<String,EcodDomain> reps = new TreeMap<String,EcodDomain>();
		for(EcodDomain d : ecod.getAllDomains()) {
			String classification = getHierarchicalName(d, level);
			if( reps.containsKey(classification)) {
				// Use first manual representative encountered
				Boolean isRep = d.getManual();
				if( isRep != null && isRep) {
					EcodDomain currentRep = reps.get(classification);
					Boolean currentIsRep = currentRep.getManual();
					if(currentIsRep == null || !currentIsRep) {
						// only replace a representative if it was not manual and this is
						reps.put(classification, d);
					}
				}

			} else {
				reps.put(classification, d);
			}
		}
		return reps;
	}
	
	/**
	 * Joins the hierarchical classification of an ecod domain with periods.
	 * 
	 * Levels are: X,H,T,F, and domain
	 * Example:
	 * getHierarchicalName(e2gs2A1,4) -> "206.1.1.25"
	 * getHierarchicalName(e2gs2A1,1) -> "206"
	 * 
	 * @param d Domain
	 * @param level Number of hierarchical levels to include
	 * @return dot-separated string of first levels of the group
	 */
	public static String getHierarchicalName(EcodDomain d,int level) {
		final String sep = ".";
		if(level < 1 || level > 5) {
			throw new IllegalArgumentException("Unrecognized ecod hierarchy level: "+level);
		}
		StringBuilder str = new StringBuilder();
		if(level >= 1) {
			str.append( d.getXGroup() );
		}
		if(level >= 2) {
			str.append(sep).append( d.getHGroup() );
		}
		if(level >= 3) {
			str.append(sep).append( d.getTGroup() );
		}
		if(level >= 4) {
			str.append(sep).append( d.getFGroup() );
		}
		if(level >= 5) {
			str.append(sep).append(d.getUid());
		}

		return str.toString();
	}
	
	public static void main(String[] args) throws IOException {
		EcodRepresentatives ecodreps = new EcodRepresentatives();
		List<EcodDomain> reps;
//		reps = ecodreps.getDomains(1);
//		System.out.println("Found "+reps.size()+" X-group representatives");
//		reps = ecodreps.getDomains(2);
//		System.out.println("Found "+reps.size()+" H-group representatives");
//		reps = ecodreps.getDomains(3);
//		System.out.println("Found "+reps.size()+" T-group representatives");
//		reps = ecodreps.getDomains(4);
//		System.out.println("Found "+reps.size()+" F-group representatives");
//		reps = ecodreps.getDomains(5);
//		System.out.println("Found "+reps.size()+" domains");

		reps = ecodreps.getDomains(3);
		for(EcodDomain d: reps) {
			String rangeStr = String.format("%s.%s",d.getPdbId(),d.getRange());
			System.out.println(rangeStr);
		}
	}


}
