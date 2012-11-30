/**
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on Nov 9, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package demo;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopFactory;
import org.biojava.bio.structure.scop.ScopInstallation;

import org.rcsb.ava.dbsearch.domainassign.DomainAssigner;
import org.rcsb.fatcat.server.PdbChainKey;

public class FetchPDPDomain {
	public static void main(String[] args){
		//String pdpId = "PDP:3ANUAa"; // a multi fragment one
		//String pdpId = "PDP:2VDCGa"; // a large one
		String pdpId = "PDP:2UV8Ag" ; // shuffled order of ranges
		String range = DomainAssigner.getDomainRangeForPDPId(pdpId);
		System.out.println("got range: " + pdpId + " : " + range);
		//A:345-343,A:983-1092,A:1182-1452,A:1505-1701
		//2UV8Ag	A:345-343,A:983-1092,A:1182-1452,A:1505-1701
		PdbChainKey key  = PdbChainKey.fromName(pdpId);
		try {
			System.out.println(DomainAssigner.getDomainRangesForName(key));
		} catch (Exception e){
			e.printStackTrace();
		}

		AtomCache cache = new AtomCache();
		//ScopInstallation scop = new ScopInstallation();
		//ScopFactory.setScopDatabase(scop);

		//		ScopDatabase scop = ScopFactory.getSCOP();
		//		
		//		System.out.println(scop.getDomainsForPDB("1bxz"));
		//		//String splitScop = "d1cph.1";
		//		String splitScop = "d1bxza1";
		//		ScopDomain dom = scop.getDomainByScopID(splitScop);
		//		System.out.println(dom);

		try {
			Structure s = cache.getStructure(pdpId);
			Atom[] ca = StructureTools.getAtomCAArray(s);
			System.out.println(ca.length);
			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(s);


			//			Structure s2 = cache.getStructure(splitScop);
			//			StructureAlignmentJmol jmol2 = new StructureAlignmentJmol();
			//			jmol2.setStructure(s2);

		} catch (Exception e){
			e.printStackTrace();
		}
	}
}
