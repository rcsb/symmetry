package demo;

import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.align.symm.CeSymm;

/**
 * Shows the centroid of a symmetric structure.
 * @author dmyerstu
 */
public class CalcAlignmentCentroid {

	public static void main(String[] args) throws Exception {
		String name = "d1srdb_";
		display(name);
	}
	
	public static void display(String name) throws Exception {
		AtomCache cache = new AtomCache();
		Atom[] ca1 = cache.getAtoms(name);
		Atom[] ca2 = cache.getAtoms(name);
		CeSymm ceSymm = new CeSymm();
		AFPChain afpChain = ceSymm.align(ca1, ca2);
		Atom alignedCentroid = calcCentroidFromAfpChain(afpChain, ca1);
		RotationAxis axis = new RotationAxis(afpChain);
		StructureAlignmentJmol jmolPanel = StructureAlignmentDisplay.display(afpChain, ca1, ca2);
		Atom domainCentroid = Calc.getCentroid(ca1);
		String centroidCmd = "draw ID alignCentroid color green diameter 3.0 CIRCLE {" + alignedCentroid.getX() + " " + alignedCentroid.getY() + " " + alignedCentroid.getZ() + "}";
		String domainCmd = "draw ID structCentroid color red diameter 3.0 CIRCLE {" + domainCentroid.getX() + " " + domainCentroid.getY() + " " + domainCentroid.getZ() + "}";
		String axisCmd = axis.getJmolScript(ca1);
		jmolPanel.evalString(domainCmd);
		jmolPanel.evalString(centroidCmd);
		jmolPanel.evalString(axisCmd);
	}


	private static Atom calcCentroidFromAfpChain(AFPChain afpChain, Atom[] ca) throws StructureException {
		Map<Integer,Integer> map = AlignmentTools.alignmentAsMap(afpChain);
		Atom[] alignedAtoms = new Atom[map.size()];
		int j = 0;
		for (int x : map.keySet()) {
			alignedAtoms[j] = ca[x];
			j++;
		}
		return Calc.getCentroid(alignedAtoms);
	}

}
