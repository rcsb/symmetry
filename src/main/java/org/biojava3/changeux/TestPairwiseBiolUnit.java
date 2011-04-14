package org.biojava3.changeux;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.rcsb.fatcat.server.PdbChainKey;

public class TestPairwiseBiolUnit {

	public static void main(String[] args){

		PdbChainKey repre = PdbChainKey.fromName("3KFK.A");
		PdbChainKey member =PdbChainKey.fromName("3LOS.A");

		try {
			Structure s1 = AlignQuaternaryStructures.cache.getBiologicalUnit(repre.getPdbId());
			Structure s2 = AlignQuaternaryStructures.cache.getBiologicalUnit(member.getPdbId());

			s1.getPDBHeader().setTitle("Biological Unit of " + repre.getPdbId());
			s2.getPDBHeader().setTitle("Biological Unit of " + member.getPdbId());
			s1.setPDBCode(repre.getPdbId());
			s2.setPDBCode(member.getPdbId());

//
//			StructureAlignmentJmol jmol1 = new  StructureAlignmentJmol(null,null,null);
//			jmol1.setStructure(s1);
//
//			StructureAlignmentJmol jmol2 = new  StructureAlignmentJmol(null,null,null);
//			jmol2.setStructure(s2);
			
			StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(SmithWaterman3Daligner.algorithmName);
			 AlignQuaternaryStructures.align(repre, member, ce, null,true);
			 System.out.println("done");
			
		} catch (Exception e){
			e.printStackTrace();
		}

	}
}
