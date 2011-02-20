package org.biojava3.structure.utils;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;

public class SymmetryTools {
	
	
	public static Atom[] mirrorCoordinates(Atom[] ca2O) {
		for(int i=0;i<ca2O.length;i++) {
			//ca2O[i].setX(-ca2O[i].getX());
			Group g = ca2O[i].getGroup();
			for ( Atom a : g.getAtoms()){
				a.setX(-a.getX());
			}
		}

		return ca2O;
	}
	
	
	public static Atom[] duplicateMirrorCA2(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length*2];

		int pos = ca2clone.length - 1;

		Chain c = new ChainImpl();
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos--;
		}


		// Duplicate ca2!
		for (Atom a : ca2){
			Group g = (Group)a.getGroup().clone();
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos--;
		}

		return ca2clone;


	}
}
