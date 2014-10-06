package demo;
/*
 *                  BioJava development code
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
 * Created on Feb 22, 2011
 *
 */




import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.jama.EigenvalueDecomposition;
import org.biojava.bio.structure.jama.Matrix;


/**
 *
 * @author Peter Rose
 * @author Andreas Prlic
 */
public class MomentsOfInertia {


	public static void main (String[] args){
		AtomCache cache = new AtomCache();


		int rotations = 3; 

		try {
			
			Atom[] ca = cache.getAtoms("1FE0.A");
			//MomentsOfInertia moi = new MomentsOfInertia(ca);
			//float[] momi = moi.getPrincipalMomentsOfInertia();
//			for ( int i = 0 ; i < momi.length; i++){
//				System.out.println(momi[i]);
//			}
			

			for ( int r = 0 ; r < rotations ; r++){
				
				Atom[] ca2 = StructureTools.cloneCAArray(ca);
				
				ca2 = Calc.centerAtoms(ca2);
				
				double angle = (360 /(double)rotations) * r;
				System.out.println("rotation: " + r + " angle: " + angle);
				
				Matrix m = Calc.matrixFromEuler( 0, angle, 0);

				// X axis
				for ( Atom a : ca2){
					Group g = a.getGroup();
					Calc.rotate(g,m);
				}
				
	
				MomentsOfInertia moi2 = new MomentsOfInertia(ca2);
				float[] momi2 = moi2.getPrincipalMomentsOfInertia();
				for ( int i = 0 ; i < momi2.length; i++){
					System.out.println(i+" "+ momi2[i]);
				}
				System.out.println("principle axes "+ moi2.getPrincipalAxes());
			}

		} catch (Exception e){
			e.printStackTrace();
		}
	}




	private boolean modified = true;

	private float[] principalMomentsOfInertia = new float[3];
	private Atom[] principalAxes = new Atom[3];

	Atom[] ca;


	/** Creates a new instance of MomentsOfInertia */
	public MomentsOfInertia(Atom[] ca) {
		this.ca = ca;
	}






	public float[] getPrincipalMomentsOfInertia() throws StructureException{
		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		return principalMomentsOfInertia;
	}

	public Atom[] getPrincipalAxes() throws StructureException {
		if (modified) {
			diagonalizeTensor();
			modified = false;
		}
		return principalAxes;
	}



	private double[][] getInteriaTensor() throws StructureException {
		//        http://en.wikipedia.org/wiki/Moment_of_inertia
		Atom p = new AtomImpl();

		double[][] tensor = new double[3][3];
		
		// calculate the interia tensor at center of mass
		Atom com = Calc.getCentroid( ca );
		System.out.println("Center of mass: " + com);
		for (int i = 0, n = ca.length; i < n; i++) {
			Atom a = ca[i];
			float mass = a.getElement().getAtomicMass();
			p = Calc.subtract(ca[i], com);


			tensor[0][0] += mass * (p.getY() * p.getY() + p.getZ() * p.getZ());
			tensor[1][1] += mass * (p.getX() * p.getX() + p.getZ() * p.getZ());
			tensor[2][2] += mass * (p.getX() * p.getX() + p.getY() * p.getY());

			tensor[0][1] -= mass * p.getX() * p.getY();
			tensor[0][2] -= mass * p.getX() * p.getZ();
			tensor[1][2] -= mass * p.getY() * p.getZ();
		}

		tensor[1][0] = tensor[0][1];
		tensor[2][0] = tensor[0][2];
		tensor[2][1] = tensor[1][2];

		return tensor;
	}

	private float[] diagonalizeTensor() throws StructureException{
		Matrix m = new Matrix(getInteriaTensor());

		EigenvalueDecomposition eig = m.eig();
		double[] eigenValues = eig.getRealEigenvalues();
		double[][] eigenVectors = eig.getV().getArray();

		for (int i = 0; i < 3; i++) {
			principalMomentsOfInertia[i] = (float)eigenValues[i];

			double x = eigenVectors[i][0];
			double y = eigenVectors[i][1];            
			double z = eigenVectors[i][2];

			principalAxes[i] = new AtomImpl();



			principalAxes[i] = new AtomImpl();
			principalAxes[i].setX(x);
			principalAxes[i].setY(y);
			principalAxes[i].setZ(z);
		}

		return principalMomentsOfInertia;
	}
}
