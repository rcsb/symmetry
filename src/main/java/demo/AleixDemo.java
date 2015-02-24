// 24 February 2015 -- Aleix Lafita
// Biojava 3 Tutorial on Github

package demo;

import java.util.List;

import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;

public class AleixDemo {
	
	public static void main(String[] args) {
		try {
			
			Structure structure = StructureIO.getStructure("4HHB");
			StructureAlignmentJmol jmolPanel = new StructureAlignmentJmol();
			
			jmolPanel.setStructure(structure);
			
			jmolPanel.evalString("select *; color chain;");
			jmolPanel.evalString("select *; spacefill off; wireframe off; cartoon on; ");
			jmolPanel.evalString("select ligands; cartoon off; wireframe 0.3; spacefill 0.5; color cpk;");
			
			//print the number of atoms of the protein
			System.out.println(StructureTools.getNrAtoms(structure));
		} catch (Exception e){
			e.printStackTrace();
		}

	}

}
