package org.biojava3.structure.sequenceui;

import java.util.SortedSet;

import javax.swing.JFrame;

import org.biojava3.structure.dbscan.GetRepresentatives;
import org.rcsb.fatcat.server.PdbChainKey;
import org.rcsb.sequence.biojavadao.BioJavaPubMedFactory;
import org.rcsb.sequence.biojavadao.BioJavaResidueInfoFactory;
import org.rcsb.sequence.biojavadao.BioJavaSequenceCollectionFactory;
import org.rcsb.sequence.core.PubMedProvider;
import org.rcsb.sequence.core.ResidueProvider;
import org.rcsb.sequence.core.SequenceCollectionFactory;
import org.rcsb.sequence.core.SequenceCollectionProvider;
import org.rcsb.sequence.model.ResidueInfoFactory;
import org.rcsb.sequence.util.ImageWrapper;


public class TestAllChains {

	public static void main(String[] args){

		initBioJavaView();

		SortedSet<PdbChainKey> reps = GetRepresentatives.getRepresentatives();

		int counter = 0;
		for (PdbChainKey rep : reps){
			counter++;
			try {
				System.out.println("#"  + counter + " " + rep);
				JFrame frame = showChain(rep);

				frame.dispose();
			} catch (Exception e){
				e.printStackTrace();
				System.err.println("Could not display " + rep);
				System.exit(0);
			}
		}


	}

	private static JFrame showChain(PdbChainKey rep) { 

		return ImageWrapper.showSeq(rep.getPdbId(), rep.getChainId());


	}

	/** provide a default view using BioJava.. could be done using some proper configuration managment...
	 * 
	 */
	public static void initBioJavaView(){

		// first the Residue Provider
		ResidueInfoFactory refactory = new BioJavaResidueInfoFactory();
		ResidueProvider.setResidueInfoFactory(refactory);

		// next the SequenceCollection Provider
		SequenceCollectionFactory sfactory = new BioJavaSequenceCollectionFactory();
		SequenceCollectionProvider.setSequenceCollectionFactory(sfactory);

		BioJavaPubMedFactory pfactory = new BioJavaPubMedFactory();
		PubMedProvider.setPubMedFactory(pfactory);


	}
}
