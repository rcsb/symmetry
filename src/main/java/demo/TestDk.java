package demo;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JFrame;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.gui.ScaleableMatrixPanel;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.nbio.structure.utils.SymmetryTools;



/**
 * 
 * @deprecated
 */
@Deprecated
public class TestDk {


	// 1btm vs 1btm

	public static void main(String[] args){
		TestDk dk = new TestDk();

		int k = 6;
		int fragmentLength = 10;
		
		
		try {
			AtomCache cache = new AtomCache();
			
			
			// intra and intermolecular symm
			//String chainId1 = "1jnr.D";
			//String chainId2 = "1jnr.B";
			
			// intramolecular symm
			//String chainId1 = "1mer.A";
			//String chainId2 = "1mer.A";
			
			
			// not related
			String chainId1 = "1cdg.A";
			String chainId2 = "1tim.A";
			
			
			// symm
			//String chainId1 = "1bp7.A";
			//String chainId2 = "1bp7.B";
			
			
			//String chainId1 = "4hhb.A";
			//String chainId2 = "4hhb.A";
			
			Structure s = cache.getStructure(chainId1);
			
			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(s);
			
			Atom[] ca1 = cache.getAtoms(chainId1);
			Atom[] ca2 = cache.getAtoms(chainId2);
			

			dk.align(ca1,ca2, k, fragmentLength);
			
		} catch (Exception e){
			e.printStackTrace();
		}
	}

	public void align(Atom[] ca1, Atom[] ca2, int k, int fragmentLength){
		
		//System.out.println("rows "  + rows + " " + cols +
		//      " ca1 l " + ca1.length + " ca2 l " + ca2.length);

		Matrix m2 = SymmetryTools.getDkMatrix(ca1, ca2, k, fragmentLength);

		ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
		JFrame frame = new JFrame();
		frame.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				JFrame f = (JFrame) e.getSource();
				f.setVisible(false);
				f.dispose();
				System.exit(0);
			}
		});
		smp.setMatrix(m2);
		//smp.getMatrixPanel().setScale(0.8f);
		
		frame.getContentPane().add(smp);

		frame.pack();
		frame.setVisible(true);





	}

	
}
