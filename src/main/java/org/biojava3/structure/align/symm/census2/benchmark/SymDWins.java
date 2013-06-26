package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava3.structure.align.symm.census2.Result;
import org.biojava3.structure.align.symm.census2.Significance;
import org.biojava3.structure.align.symm.census2.SignificanceFactory;
import org.biojava3.structure.align.symm.protodomain.Protodomain;

public class SymDWins {

	public static List<String> findWins(Sample ceSymm, Sample symD, Significance ceSymmSig, Significance symDSig) {
		List<String> scopIds = new ArrayList<String>();
		HashSet<String> all = new HashSet<String>();
		HashSet<String> a = new HashSet<String>();
		HashSet<String> b = new HashSet<String>();
		for (Case c : ceSymm.getData()) {
			boolean sig = ceSymmSig.isSignificant(c.getResult());
			if (sig) a.add(c.getScopId());
			all.add(c.getScopId());
		}
		for (Case c : symD.getData()) {
			boolean sig = symDSig.isSignificant(c.getResult());
			if (sig) b.add(c.getScopId());
			if (!all.contains(c.getScopId())) throw new IllegalArgumentException("SCOP Id " + c.getScopId() + " is missing");
		}
		for (Case c : ceSymm.getData()) {
			if (!c.getKnownInfo().hasRotationalSymmetry()) continue;
			String s = c.getScopId();
			if (b.contains(s) && !a.contains(s)) {
				scopIds.add(s);
			}
		}
		return scopIds;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if (args.length != 2) {
			System.err.println("Usage: " + SymDWins.class.getSimpleName() + " cesymm-file symd-file");
			return;
		}
		Significance symDSig = new Significance() {
			private static final double cutoff = 9.2;
			@Override
			public boolean isPossiblySignificant(AFPChain afpChain) {
				return true; // whatever
			}

			@Override
			public boolean isSignificant(Protodomain protodomain, int order, double angle, AFPChain afpChain) {
				return true; // whatever
			}

			@Override
			public boolean isSignificant(Result result) {
				if (result.getAlignment() == null) return false;
				return result.getAlignment().getzScore() >= cutoff;
			}
		};
		Sample ceSymm = Sample.fromXML(new File(args[0]));
//		for (Case c : ceSymm.getData()) {
//			if (c.getAlignment().getTmScore() >= 0.38 && c.getAlignment().getTmScore() <= 0.4) System.out.println(c.getScopId());
//		}
		ScopDatabase scop = ScopFactory.getSCOP(ScopFactory.VERSION_1_75A);
		List<String> wins = findWins(ceSymm, Sample.fromXML(new File(args[1])), SignificanceFactory.symmetric(), symDSig);
		for (String s : wins) {
			ScopDomain d = scop.getDomainByScopID(s);
			int sf = d.getFoldId();
			String sfDescription = scop.getScopDescriptionBySunid(sf).getDescription();
			System.out.println(s + "\t\t" + d.getClassificationId() + "\t\t" + sfDescription);
		}
	}

}
