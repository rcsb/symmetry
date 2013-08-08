package org.biojava3.structure.align.symm.census;

import java.util.concurrent.Callable;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava3.changeux.IdentifyAllSymmetries;
import org.biojava3.structure.align.symm.CeSymm;
import org.biojava3.structure.align.symm.census2.CensusJob;
import org.biojava3.structure.align.symm.protodomain.Protodomain;
import org.biojava3.structure.utils.SymmetryTools;

/**
 * @deprecated Replaced by {@link CensusJob}
 */
@Deprecated
public class SymmetryCalculationJob implements Callable<CensusResult> {


	public String name;
	AtomCache cache;
	Integer count;
	ScopDescription scopDescription;


	public String getName() {
		return name;
	}



	public void setName(String name) {
		this.name = name;
	}

	public AtomCache getCache() {
		return cache;
	}



	public void setCache(AtomCache cache) {
		this.cache = cache;
	}

	public Integer getCount() {
		return count;
	}



	public void setCount(Integer count) {
		this.count = count;
	}


	public ScopDescription getScopDescription() {
		return scopDescription;
	}

	public void setScopDescription(ScopDescription scopDescription) {
		this.scopDescription = scopDescription;
	}

	public CensusResult call() throws Exception {

		System.out.println("# " + count + " " + name);
		try {

			String name1 = name;
			String name2 = name;

			Atom[] ca1 = cache.getAtoms(name1);

			if ( ca1.length < 20) {
				System.out.println(" ... ignoring " + name);
				return null;
			}

			Atom[] ca2 = cache.getAtoms(name2);

			boolean isSymmetric = false;

			StructureAlignment ceSymm = ScanSCOPForSymmetry.getCeSymm();
			AFPChain afpChain = ceSymm.align(ca1,ca2);

			double angle = -1;

			if ( afpChain == null) { 
				return null;
			}
			angle = SymmetryTools.getAngle(afpChain,ca1,ca2);

			if ( IdentifyAllSymmetries.isSignificant(afpChain)) {

				if ( angle > 20) {


					isSymmetric = true;
				}
			}

			int order = CeSymm.getSymmetryOrder(afpChain);
			
			String protodomain = Protodomain.fromSymmetryAlignment(afpChain, ca1, 1, cache).toString();
			
			/*
			String protodomain = "";
			if ( afpChain.getOptLength() > 0 ) {
				


				int[]optLen = afpChain.getOptLen();
				int[][][] optAln = afpChain.getOptAln();

				int p1 = optAln[0][0][0];
				int lenBlock1 = optLen[0] - 1;
				int p2 = optAln[0][0][lenBlock1];

				Atom a1 = ca1[p1];
				Atom a2 = ca1[p2];

				String chainId  = a1.getGroup().getChain().getChainID();
				Group g1 = a1.getGroup();
				Group g2 = a2.getGroup();
				PdbChainKey key = PdbChainKey.fromName(name);
				protodomain = key.getPdbId();
				protodomain += "." + chainId + "_";
				protodomain+= g1.getResidueNumber().toString()  + "-" + g2.getResidueNumber().toString();
			}
			*/
			
			//System.out.println(isSymmetric + " : " + name + " " +  domain + " : "  );

			//StringBuffer str = printTabbedResult(afpChain, isSymmetric, superfamily,name1, count);
			//StringBuffer str = printHTMLResult(afpChain, isSymmetric, superfamily,name1, count, angle);

			CensusResult result = convertResult(afpChain, isSymmetric, scopDescription, name1, count, angle, order, protodomain);

			//System.out.println(result);
			return result;
		} catch (Exception e){
			e.printStackTrace();
			System.err.println("ERROR processing " + name + " " + e.getMessage());


		}
		System.err.println("returning null");
		return null;

	}

	private CensusResult convertResult(AFPChain afpChain, boolean isSymmetric, ScopDescription superfamily, 
			String name, int count, double angle , int order, String protodomain){

		String description = "";
		Character scopClass = '?';
		if ( superfamily != null){
			description  = superfamily.getDescription();
			scopClass = superfamily.getClassificationId().charAt(0);
		}
		CensusResult r = new CensusResult();

		r.setRank(count);
		r.setIsSignificant(isSymmetric);
		r.setName(name);
		if ( superfamily != null)
			r.setClassificationId(superfamily.getClassificationId());
		else 
			r.setClassificationId("-");
		r.setzScore(afpChain.getProbability());
		r.setRmsd(afpChain.getTotalRmsdOpt());
		r.setTmScore(afpChain.getTMScore());
		r.setAlignScore(afpChain.getAlignScore());
		r.setIdentity((float)afpChain.getIdentity());
		r.setSimilarity((float)afpChain.getSimilarity());
		r.setLength1(afpChain.getCa1Length());
		r.setAligLength(afpChain.getOptLength());
		r.setAngle((float)angle);
		r.setDescription(description);
		r.setScopClass(scopClass);
		r.setOrder(order);
		r.setProtoDomain(protodomain);

		return r;
	}

}
