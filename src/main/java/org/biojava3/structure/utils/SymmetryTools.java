package org.biojava3.structure.utils;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JFrame;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ChainImpl;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.ce.CECalculator;
import org.biojava.bio.structure.align.helper.AlignTools;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.gui.ScaleableMatrixPanel;
import org.biojava.bio.structure.jama.Matrix;

public class SymmetryTools {

	// there won;t be an instance of this
	private SymmetryTools(){}


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


	public static Atom[] duplicateCA2(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length*2];

		int pos = 0;

		Chain c = new ChainImpl();
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}


		// Duplicate ca2!
		for (Atom a : ca2){
			Group g = (Group)a.getGroup().clone();
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}

		return ca2clone;


	}


	public static Atom[] duplicateMirrorCA2(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length];

		int pos = ca2clone.length - 1;

		Chain c = new ChainImpl();
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom
			c.addGroup(g);
			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos--;
		}


//		// Duplicate ca2!
//		for (Atom a : ca2){
//			Group g = (Group)a.getGroup().clone();
//			c.addGroup(g);
//			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);
//
//			pos--;
//		}

		return ca2clone;


	}

	public static void showMatrix(Matrix m, String string) {
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

		smp.setMatrix((Matrix)m.clone());
		//smp.getMatrixPanel().setScale(0.8f);

		frame.setTitle(string);
		frame.getContentPane().add(smp);

		frame.pack();
		frame.setVisible(true);

	}

	public Matrix  getDkMatrix(Atom[] ca1, Atom[] ca2,int fragmentLength,
			double[] dist1, double[] dist2, int rows, int cols) {
		Matrix diffDistMax =  Matrix.identity(ca1.length, ca2.length);

		for ( int i = 0 ; i< rows; i++){
			double score1 = 0;
			for ( int x=0 ; x < fragmentLength ; x++){
				score1 += dist1[i+x];
			}
			for ( int j = 0 ; j < cols ; j++){
				double score2 = 0;
				for ( int y=0 ; y < fragmentLength ; y++){
					score2 += dist2[j+y];
				}

				// if the intramolecular distances are very similar
				// the two scores should be similar, i.e. the difference is close to 0
				diffDistMax.set(i,j, Math.abs(score1-score2));
			}
		}


		// symmetry hack, disable main diagonale

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					diffDistMax.set(i,j, 99);
				}
				int diff2 = Math.abs(i-(j-ca2.length/2));
				if ( diff2 < 15 ){
					diffDistMax.set(i,j, 99);
				}
			}
		}
		return diffDistMax;

	}

	public static Matrix blankOutPreviousAlignment(AFPChain afpChain, Atom[] ca2,
			int rows, int cols, CECalculator calculator, Matrix max) {

		double[][] dist1 = calculator.getDist1();
		double[][] dist2 = calculator.getDist2();

		
		int[][][] optAln = afpChain.getOptAln();
		int blockNum = afpChain.getBlockNum();

		int[] optLen = afpChain.getOptLen();
		for(int bk = 0; bk < blockNum; bk ++)       {

			//Matrix m= afpChain.getBlockRotationMatrix()[bk];
			//Atom shift = afpChain.getBlockShiftVector()[bk];
			for ( int i=0;i< optLen[bk];i++){
				int pos1 = optAln[bk][0][i];
				int pos2 = optAln[bk][1][i];
				// blank out area around these positions...

				int dist = 10 ;
				int start1 = Math.max(pos1-dist,0);
				int start2 = Math.max(pos2-dist,0);
				int end1 = Math.min(pos1+dist, rows-1);
				int end2 = Math.min(pos2+dist, cols-1);
				//System.out.println(pos1 + "  " + pos2 + " " + start1 + " " + end1 + " " + start2 + " " + end2);
				for ( int i1 = start1; i1< end1 ; i1++){
					for ( int k=0; k < 8 ; k ++){
						if ( i-k >= 0)
							dist1[i-k][i-k] = 999;
						if ( i+k < rows)
							dist1[i+k][i+k] = 999;
						
					}

					for ( int j2 = start2 ; j2 < end2 ; j2++){
						//System.out.println(i1 + " " + j2);
						max.set(i1,j2,99);
						for ( int k=0; k < 8 ; k ++){							
							if ( j2-k >=0)
								dist2[j2-k][j2-k] = 0;
							if ( j2+k < cols)
								dist2[j2+k][j2+k] = 0;
						}
					}
				}

			}

		}
		return max;

	}

	public static Matrix blankOutCEOrig(Atom[] ca2, int rows, int cols,
			CECalculator calculator, Matrix origM) {
		origM = new Matrix( calculator.getMatMatrix());
		// symmetry hack, disable main diagonale

		//double[][] dist1 = calculator.getDist1();
		//double[][] dist2 = calculator.getDist2();

		for ( int i = 0 ; i< rows; i++){
			for ( int j = 0 ; j < cols ; j++){
				int diff = Math.abs(i-j);

				if ( diff < 15 ){
					origM.set(i,j, 99);
//					for ( int k=0; k < 5 ; k ++){
//						if ( i-k >= 0)
//							dist1[i][i-k] = 99;
//						if ( i+k < rows)
//							dist1[i][i+k] = 99;
//						if ( j-k >=0)
//							dist2[j][j-k] = 499;
//						if ( j+k < cols)
//							dist2[j][j+k] = 499;
//					}
				}
				int diff2 = Math.abs(i-(j-ca2.length/2));
				if ( diff2 < 15 ){
					origM.set(i,j, 99);
					
//					for ( int k=0; k < 5 ; k ++){
//
//						if ( i-k >= 0)
//							dist1[i][i-k] = 99;
//						if ( i+k < rows)
//							dist1[i][i+k] = 99;
//						if ( j-k >=0)
//							dist2[j][j-k] = 99;
//						if ( j+k < cols)
//							dist2[j][j+k] = 99;
//					}
				}
			}
		}
		return origM;
	}

	public static Atom[] cloneAtoms(Atom[] ca2) throws StructureException{
		// we don't want to rotate input atoms, do we?
		Atom[] ca2clone = new Atom[ca2.length];

		int pos = 0;
		for (Atom a : ca2){
			Group g = (Group) a.getGroup().clone(); // works because each group has only a CA atom

			ca2clone[pos] = g.getAtom(StructureTools.caAtomName);

			pos++;
		}

		return ca2clone;
	}

	public static Matrix getDkMatrix(Atom[] ca1, Atom[] ca2, int k, int fragmentLength) {
		double[] dist1 = AlignTools.getDiagonalAtK(ca1, k);

		double[] dist2 = AlignTools.getDiagonalAtK(ca2, k);

		int rows = ca1.length - fragmentLength - k + 1;
		int cols = ca2.length - fragmentLength - k + 1;

		// Matrix that tracks similarity of a fragment of length fragmentLength
		// starting a position i,j.

		Matrix m2 = new Matrix(rows,cols); 

		for ( int i = 0 ; i< rows; i++){
			double score1 = 0;
			for ( int x=0 ; x < fragmentLength ; x++){
				score1 += dist1[i+x];
			}
			for ( int j = 0 ; j < cols ; j++){
				double score2 = 0;
				for ( int y=0 ; y < fragmentLength ; y++){
					score2 += dist2[j+y];
				}	

				// if the intramolecular distances are very similar
				// the two scores should be similar, i.e. the difference is close to 0
				m2.set(i,j, Math.abs(score1-score2));
			}
		}
		return m2;
	}


}
