package org.biojava3.structure.align.symm.census2.benchmark;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.LinkedHashMap;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.io.FastaAFPChainConverter;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.CasePreservingProteinSequenceCreator;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.template.FastaHeaderParserInterface;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;

public class SymDFasta {

	public static void main(String[] args) throws IOException, StructureException {
		if (args.length != 1) System.err.println("Usage: " + SymDFasta.class.getSimpleName() + " pdbid-best.fasta");
		File file = new File(args[0]);
		display(file);
	}

	public static void display(File fastaFile) throws IOException, StructureException {
		AFPChain afpChain = getAlignment(fastaFile);
		Structure structure1 = StructureTools.getStructure(afpChain.getName1());
		Structure structure2 = StructureTools.getStructure(afpChain.getName2());
		Atom[] ca1 = StructureTools.getAtomCAArray(structure1);
		Atom[] ca2 = StructureTools.getAtomCAArray(structure2);
		String xml = AFPChainXMLConverter.toXML(afpChain);
		System.out.println(xml);
		StructureAlignmentDisplay.display(afpChain, ca1, ca2);
	}

	public static AFPChain getAlignment(File fastaFile) throws IOException, StructureException {
		
		InputStream inStream = new FileInputStream(fastaFile);
		SequenceCreatorInterface<AminoAcidCompound> creator = new CasePreservingProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet());
		FastaHeaderParserInterface<ProteinSequence, AminoAcidCompound> headerParser = new GenericFastaHeaderParser<ProteinSequence, AminoAcidCompound>();
		FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence, AminoAcidCompound>(inStream, headerParser, creator);
		LinkedHashMap<String, ProteinSequence> sequences = fastaReader.process();
		inStream.close();
		
		Iterator<ProteinSequence> seqIter = sequences.values().iterator();
		ProteinSequence firstSeq = seqIter.next();
		ProteinSequence secondSeq = seqIter.next();
		Iterator<String> namesIter = sequences.keySet().iterator();
		
		// first should be something like: 1nziA-original, is=21 (best)
		String firstName = namesIter.next();
		String secondName = namesIter.next();

		String pdbName = firstName.split("-")[0];
		Structure structure = StructureTools.getStructure(pdbName);
		// if we have a SCOP Id we're fine, but we might also have, e.g. 1nziA
		// in this latter case, we want 1nzi.A instead (note the dot)
		if (structure == null) {
			pdbName = pdbName.substring(0, pdbName.length()-1) + '.' + pdbName.charAt(pdbName.length()-1);
			structure = StructureTools.getStructure(pdbName);
		}
		if (structure == null) throw new IllegalArgumentException("No structure for " + pdbName + " was found");
		
		// second should be something like: 1W0P-permuted  is=-393 (best)
		Integer cpSite = Integer.parseInt(secondName.split("\\s+")[1].substring(3));
		
		AFPChain afpChain = FastaAFPChainConverter.cpFastaToAfpChain(firstSeq, secondSeq, structure, cpSite);
		return afpChain;
	}
	
}
