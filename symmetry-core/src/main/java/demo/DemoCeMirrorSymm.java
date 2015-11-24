package demo;

import java.io.IOException;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.ce.CECPParameters;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CEMirrorSymm;
import org.biojava.nbio.structure.align.util.AtomCache;

public class DemoCeMirrorSymm {

	public static void main(String[] args) throws IOException, StructureException {
		String name;
		boolean mirrorCoords, mirrorSeq;
		name = "1qys";
		mirrorSeq = true;
		mirrorCoords = true;
		// name = "d1kkeb2"; mirrorSeq = true; mirrorCoords = false;
		// name = "d2okua1"; mirrorSeq = true; mirrorCoords = false;
		name = "2cb2.A";
		mirrorSeq = true;
		mirrorCoords = false;
		// name = "d1in0a2"; mirrorSeq = true; mirrorCoords = true; //mirror
		// topology, but poor rmsd

		AtomCache cache = new AtomCache();

		Atom[] ca1 = cache.getAtoms(name);
		Atom[] ca2 = StructureTools.cloneAtomArray(ca1);

		// PDBFileReader pdbreader = new PDBFileReader();
		//
		// try{
		// name = "1qys_mirror";
		// mirrorCoords = true;
		// mirrorSeq = false;
		//
		// Structure struc =
		// pdbreader.getStructure("/Users/blivens/dev/bourne/1qys_mirror.pdb");
		// ca1 = StructureTools.getAtomCAArray(struc);
		// Structure struc2 =
		// pdbreader.getStructure("/Users/blivens/dev/bourne/1qys_mirror.pdb");
		// ca2 = StructureTools.getAtomCAArray(struc2);
		// } catch (Exception e){
		// e.printStackTrace();
		// return;
		// }

		CEMirrorSymm ce = new CEMirrorSymm(mirrorCoords, mirrorSeq);
		CECPParameters params = new CECPParameters();
		AFPChain afp = ce.align(ca1, ca2, params);
		StructureAlignmentDisplay.display(afp, ca1, ca2);
	}
	
}
