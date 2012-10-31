/**
 * 
 */
package org.biojava3.structure.align.symm.quaternary;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Quat4d;
import javax.vecmath.Vector3d;

/**
 * @author Peter
 *
 */
public class JmolScriptGenerator {
	private static String PRINCIPAL_AXIS_COLOR = "darkorange";
	private static String MINOR_AXIS_COLOR = "darkturquoise"; // royalblue
	
	private AxisTransformation axisTransformation = null;
	private Subunits subunits = null;
	private RotationGroup rotationGroup = null;
	
	public JmolScriptGenerator(AxisTransformation axisTransformation, Subunits subunits, RotationGroup rotationGroup) {
		this.axisTransformation = axisTransformation;
		this.subunits = subunits;
		this.rotationGroup = rotationGroup;
	}
	
	public String jmolDrawInertiaAxes() {
		StringBuilder s = new StringBuilder();
		Point3d centroid = axisTransformation.calcGeometricCenter();
		Vector3d[] axes = axisTransformation.getPrincipalAxesOfInertia();

		for (int i = 0; i < axes.length; i++) {
			s.append("draw l");
			s.append(200+i);
			s.append(" ");
			s.append("line");
			s.append(" ");
			Point3d v1 = new Point3d(axes[i]);
			if (i == 0) {
				v1.scale(1.2*axisTransformation.getYRadius());
			} else if (i == 1) {
				v1.scale(1.2*axisTransformation.getXRadius());
			} else if (i == 2) {
				v1.scale(1.2*axisTransformation.getZRadius());
			}
			Point3d v2 = new Point3d(v1);
			v2.negate();
			v1.add(centroid);
			v2.add(centroid);
			s.append(getJmolPoint(v1));
			s.append(getJmolPoint(v2));
			s.append(" width 0.5 ");
			s.append(" color white");
			s.append(";");
		}
        return s.toString();
	}

	/**
	 * Returns a Jmol script to rotate a structure into
	 * a canonical orientation. For symmetric structures,
	 * the highest order axis is aligned with the z-axis, 
	 * and for asymmetric structures, the principal axis of
	 * inertia is aligned with the y-axis.
	 * @return Jmol script
	 */
	public String getJmolTransformation() {
		Matrix4d t = axisTransformation.getTransformation();
		Quat4d q = new Quat4d();
		t.get(q);
		Point3d centroid = subunits.getCentroid();
		if (rotationGroup.getPointGroup().equals("C1")) {
			centroid = axisTransformation.calcGeometricCenter();
		}
		int zoom = 80;
		if (rotationGroup.getPointGroup().equals("T")) {
			zoom = 60;
		} else if (rotationGroup.getPointGroup().equals("O") || rotationGroup.getPointGroup().equals("I")) {
			zoom = 80;
		}
		return "rotate quaternion {" + jMolFloat(q.x) + " " + jMolFloat(q.y) + " " + jMolFloat(q.z) + " " + jMolFloat(q.w) + "}; zoom " + zoom + "; center {"+ centroid.x + " " + centroid.y + " " + centroid.z + "};";
	//	return "rotate quaternion {" + jMolFloat(q.x) + " " + jMolFloat(q.y) + " " + jMolFloat(q.z) + " " + jMolFloat(q.w) + "}; zoom " + zoom + ";";

	}

	public String getJmolSymmetryAxes() {
		StringBuilder s = new StringBuilder();

		int n = rotationGroup.getOrder();
		if (n == 0) {
			return s.toString();
		}

		float diameter = 0.5f;
		double radius = 0;
		String color = "";

		for (int i = 0; i < n; i++) {
			Rotation rotation = rotationGroup.getRotation(i);
			int direction =  rotation.getDirection();

			// don't draw redundant n-fold rotations around principal axis
			if (i > 0 && direction == 0) {
				continue;
			}

			if (rotationGroup.getPointGroup().startsWith("C") || rotationGroup.getPointGroup().startsWith("D")) {
				if (direction == 0) {
					radius = 1.2 * axisTransformation.getZRadius(); // principal axis uses z-dimension
				} else {
					radius = 1.1 * axisTransformation.getXYRadius();
				} 
			} else if (rotationGroup.getPointGroup().equals("T")) {
				radius = 0.9 * Math.max(axisTransformation.getZRadius(), axisTransformation.getXYRadius());
			} else { 
				radius = Math.max(axisTransformation.getZRadius(), axisTransformation.getXYRadius()); // for O, I point group
			}

			if (direction == 0) { // for principal axis
	//			radius = 1.2 * radius;
				color = PRINCIPAL_AXIS_COLOR;
				diameter = 0.5f;
			} else { // for all other axes
	//			radius = 1.1 * radius;
				color = MINOR_AXIS_COLOR;
				diameter = 0.25f;
			} 

			Point3d center = axisTransformation.calcGeometricCenter();
			AxisAngle4d axisAngle = rotation.getAxisAngle();
			Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);

			s.append(JmolScriptGenerator.getSymmetryAxis(i*100, direction, rotationGroup.getPointGroup(), rotation.getFold(), axisTransformation.getPrincipalRotationAxis(), axisTransformation.getRotationReferenceAxis(), radius, diameter, color, center, axis));

		}
		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure in Jmol
	 * into canonical orientations and animates the rotation
	 * around unique symmetry axes.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	public String getJmolAnimation(int delay) {
		String animation = "";
		String pointGroup = rotationGroup.getPointGroup();
		if (pointGroup.equals("C1")) {
			animation = getJmolAnimationC1(delay);
		} else if (pointGroup.equals("C2")) {
			animation = getJmolAnimationC2(delay);
		} else if (pointGroup.startsWith("C")) {
			animation = getJmolAnimationCyclic(delay);
		} else if (pointGroup.startsWith("D")) {
			animation = getJmolAnimationDihedral(delay);
		} else if (pointGroup.equals("T")) {
			animation = getJmolAnimationTetrahedral(delay);
		} else if (pointGroup.equals("O")) {
			animation = getJmolAnimationOctahedral(delay);
		} else if (pointGroup.equals("I")) {
			animation = getJmolAnimationIcosahedral(delay);
		}
		return animation;
	}

	/**
	 * Returns a Jmol script that rotates a structure without symmetry (C1)
	 * animates the rotation around the Cn axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	private String getJmolAnimationC1(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());	

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show inertia axes
		s.append(jmolDrawInertiaAxes());

		RectangularPrism p = new RectangularPrism(axisTransformation.getZRadius()*2, axisTransformation.getXRadius()*2, axisTransformation.getYRadius()*2);
		s.append(JmolScriptGenerator.getJmolWirePolyhedron(100, p, axisTransformation.calcPrismTransformation()));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show Side view 1
		s.append("set echo top center;");
		s.append("echo Top view;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Back view
		s.append("set echo top center;");
		s.append("echo Side view;");
		s.append("move  0 90 0 0 0 0 0 0 2;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view

		s.append("set echo top center;");
		s.append("echo;");
		s.append("move 0 -90 0 0 0 0 0 0 2;");
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("move -90 0 0 0 0 0 0 0 2;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}
	/**
	 * Returns a Jmol script that rotates a structure with
	 * cyclic (Cn) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the Cn axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	private String getJmolAnimationC2(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());	

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo ");
		s.append("2-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());

		// TODO the center for Cn should be the geometric center, not the centroid (part of revereseTransformationMatrix)

		RectangularPrism p = new RectangularPrism(axisTransformation.getZRadius()*2, axisTransformation.getXRadius()*2, axisTransformation.getYRadius()*2);
		s.append(JmolScriptGenerator.getJmolWirePolyhedron(100, p, axisTransformation.calcPrismTransformation()));

		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around 2-fold axis
		s.append("move 0 180 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Side view 1
		s.append("set echo top center;");
		s.append("echo Top view ;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Back view
		s.append("set echo top center;");
		s.append("echo Side view;");
		s.append("move 0 90 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo;");
		s.append("move 0 -90 0 0 0 0 0 0 2;");
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("move -90 0 0 0 0 0 0 0 2;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * cyclic (Cn) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the Cn axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	private String getJmolAnimationCyclic(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());	

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		int fold = getPrincipalRotationAxisFold();

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo ");
		s.append(fold);
		s.append("-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
//		s.append(jmolDrawInertiaAxes());
		s.append(getJmolSymmetryAxes());

		// TODO the center for Cn should be the geometric center, not the centroid (part of revereseTransformationMatrix)
	
		Prism p = new Prism(fold);
		p.setHeight(axisTransformation.getZRadius()*2);
		p.setInscribedRadius(axisTransformation.getXYRadius());
		s.append(JmolScriptGenerator.getJmolWirePolyhedron(100, p, axisTransformation.calcPrismTransformation()));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		fold  = getPrincipalRotationAxisFold();
		float angle = 360.0f/fold;
		s.append("move 0 0 ");
		s.append(angle);
		s.append(" 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Side view 1
		s.append("set echo top center;");
		s.append("echo Side view 1;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Back view
		s.append("set echo top center;");
		s.append("echo Back view;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Side view 2
		s.append("set echo top center;");
		s.append("echo Side view 2;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("move 90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * dihedral (Dn) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the Cn and perpendicular
	 * C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	public String getJmolAnimationDihedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

		// show front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append(rotationGroup.getPointGroup());
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		int fold  = getPrincipalRotationAxisFold();
		s.append("set echo top center;");
		s.append("echo ");
		s.append(fold);
		s.append("-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());

		Prism p = new Prism(fold);
		p.setHeight(axisTransformation.getZRadius()*2);
		p.setInscribedRadius(axisTransformation.getXYRadius());
		s.append(JmolScriptGenerator.getJmolWirePolyhedron(100, p, axisTransformation.getReverseTransformation()));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		float angle = 360.0f/fold;
		s.append("move 0 0 ");
		s.append(angle);
		s.append(" 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show side view
		s.append("set echo top center;");
		s.append("echo Side view;");
		s.append("color echo white;");
		// align minor axis with z-axis
		s.append("move 90 0 0 0 0 0 0 0 4;");

		// rotate around minor axis
		fold  = getMinorRotationAxisFold();
		angle = 360.0f/fold;
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo ");
		s.append(fold);
		s.append("-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 ");
		s.append(angle);
		s.append(" 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// rotate back
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move 0 0 ");
		s.append(360.0-angle);
		s.append(" 0 0 0 0 0 0;");

		// show Front view
		s.append("move -90 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * tetrahedral (T) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the C3 and C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	public String getJmolAnimationTetrahedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

		// show front view
		s.append("set echo top center;");
		s.append("echo Front view: C3-axis;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append("T (2,3)");
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
		s.append(getJmolSymmetryAxes());

		double radius = Math.max(axisTransformation.getZRadius(), axisTransformation.getXYRadius());
		Tetrahedron t = new Tetrahedron();
		t.setMidRadius(radius);
		s.append(JmolScriptGenerator.getJmolWirePolyhedron(100, t, axisTransformation.getReverseTransformation()));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		s.append("move 0 0 120 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		s.append("set echo top center;");
		s.append("echo Edge view: C2-axis;");
		s.append("color echo white;");
		double tetrahedralAngle = Math.toDegrees(Math.acos(-1.0/3.0));
		double angle = 180 - 0.5 * tetrahedralAngle;
		s.append("move ");
		s.append(jMolFloat(angle));
		s.append(" 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// rotate around 2-fold axis
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 180 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// rotate back to front view
		s.append("echo Front view: C3-axis;");
		s.append("color echo white;");
		s.append("move 0 0 -180 0 0 0 0 0 0;");
		s.append("move ");
		s.append(jMolFloat(-angle));
		s.append(" 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * octahedral (O) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the C4, C3, and C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	public String getJmolAnimationOctahedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

		// show front view
		s.append("set echo top center;");
		s.append("echo Front view: C4-axis;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append("O (4,3,2)");
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo 4-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
//		s.append(getJmolSymmetryAxes());
		double radius = Math.max(axisTransformation.getZRadius(), axisTransformation.getXYRadius());
		Octahedron o = new Octahedron();
		o.setMidRadius(radius);
		s.append(JmolScriptGenerator.getJmolWirePolyhedron(100, o, axisTransformation.getReverseTransformation()));
		s.append(drawAxis(300, o.getC4Axis(1.2)));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		s.append("move 0 0 90 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show side view (45 deg.)
		s.append("set echo top center;");
		s.append("echo Diagonal view: C3-axis;");
		s.append("color echo white;");
		s.append("move 0 45 0 0 0 0 0 0 0;");
		s.append("move 35.3 0 0 0 0 0 0 0 4;");
		s.append(drawAxis(301, o.getC3Axis(1.2)));
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 120 0 0 0 0 0 4;");

		// side view
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo Side view: C2-axis;");
		s.append("color echo white;");
		s.append("move -35.3 0 0 0 0 0 0 0 4;");
		s.append(drawAxis(302, o.getC2Axis(1.2)));
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// 2-fold rotation
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 180 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move 90 0 0 0 0 0 0 0 0;");
		s.append("move 0 0 45 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	/**
	 * Returns a Jmol script that rotates a structure with
	 * icosahedral (I) symmetry in Jmol into canonical orientations and 
	 * animates the rotation around the C5, C3, and C2 axis.
	 * @param delay period of time in seconds that the animation stops to show canonical views
	 * @return Jmol script
	 */
	public String getJmolAnimationIcosahedral(int delay) {
		StringBuilder s = new StringBuilder();
		s.append(getJmolTransformation());

		// show front view
		s.append("set echo top center;");
		s.append("echo Front view: C5-axis;");
		s.append("color echo white;");
		s.append("font echo 24 sanserif;");
		s.append("set echo bottom center;");
		s.append("echo Point group ");
		s.append("I (5,3,2)");
		s.append(";");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// show symmetry axes
		s.append("set echo top center;");
		s.append("echo 5-fold rotation;");
		s.append("color echo ");
		s.append(PRINCIPAL_AXIS_COLOR);
		s.append(";");
//		s.append(getJmolSymmetryAxes());
		double radius = Math.max(axisTransformation.getZRadius(), axisTransformation.getXYRadius());
		Icosahedron p = new Icosahedron();
		p.setMidRadius(radius);
		s.append(JmolScriptGenerator.getJmolWirePolyhedron(100, p, axisTransformation.getReverseTransformation()));
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		// rotate around principal axis
		s.append("move 0 0 72 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// align minor axis with z-axis

		// show side view (45 deg.)
		s.append("set echo top center;");
		s.append("echo Diagonal view: C3-axis;");
		s.append("color echo white;");
		s.append("move 100.81231696357165 0 0 0 0 0 0  0 4;");
		s.append("delay ");
		s.append(delay); 
		s.append(";");

		s.append("set echo top center;");
		s.append("echo 3-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 120 0 0 0 0 0 4;");

		// side view
		s.append("delay ");
		s.append(delay); 
		s.append(";");
		s.append("set echo top center;");
		s.append("echo Side view: C2-axis;");
		s.append("color echo white;");
		s.append("move 21 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// 2-fold rotation
		s.append("set echo top center;");
		s.append("echo 2-fold rotation;");
		s.append("color echo ");
		s.append(MINOR_AXIS_COLOR);
		s.append(";");
		s.append("move 0 0 180 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");

		// show Front view
		s.append("set echo top center;");
		s.append("echo Front view;");
		s.append("color echo white;");
		s.append("move -121.812 0 0 0 0 0 0 0 4;");
		s.append("delay ");
		s.append(delay);
		s.append(";");
		s.append("set echo top center;");
		s.append("echo;");

		return s.toString();
	}

	public static String getJmolWirePolyhedron(int index, Polyhedron p, Matrix4d transformation) {
		String color = "darkorange";
		StringBuilder s = new StringBuilder();

		Point3d[] vertices = p.getVertices();
		for (int i = 0; i < vertices.length; i++) {
			transformation.transform(vertices[i]);
		}

		for (int[] lineLoop: p.getLineLoops()) {
			s.append("draw l");
			s.append(index++);
			s.append(" line ");
			for (int i: lineLoop) {
				s.append(getJmolPoint(vertices[i]));
			}
			s.append(" width 0.75");
			s.append(" color ");
			s.append(color);
			s.append(";");
		}

		return s.toString();
	}
	
	private static String getJmolPoint(Point3d point) {
		StringBuilder s = new StringBuilder();
		s.append("{");
		s.append(jMolFloat(point.x));
		s.append(" ");
		s.append(jMolFloat(point.y));
		s.append(" ");
		s.append(jMolFloat(point.z));
		s.append("}");
		return s.toString();
	}

	public static String getSymmetryAxis(int index, int direction, String pointGroup, int n, Vector3d principalAxis, Vector3d referenceAxis, double radius, float diameter, String color, Point3d center, Vector3d axis) {
		if ((pointGroup.startsWith("C") ||pointGroup.startsWith("D")) && direction == 1) {
			Prism p = new Prism(n);
			p.setInscribedRadius(radius);
			radius = p.getCirumscribedRadius();
	    }
		if (pointGroup.equals("T")) {
			Tetrahedron t = new Tetrahedron();
			t.setMidRadius(radius);
			radius = t.getCirumscribedRadius();
		}
		if (pointGroup.equals("O")) {
			Octahedron o = new Octahedron();
			o.setMidRadius(radius);
			radius = o.getCirumscribedRadius();
		}
		
		Point3d p1 = new Point3d(axis);
		p1.scaleAdd(-radius, center);

		Point3d p2 = new Point3d(axis);
		p2.scaleAdd(radius, center);

		StringBuilder s = new StringBuilder();
		s.append("draw");
		s.append(" c");
		s.append(index);
		s.append(" cylinder ");
		s.append(getJmolPoint(p1));
		s.append(getJmolPoint(p2));
		s.append(" diameter ");
		s.append(diameter);
		//		s.append(" color translucent ");
		s.append(" color ");
		s.append(color);
		s.append(";");

		p1 = new Point3d(axis);
		p1.scaleAdd(-radius*0.95, center);

		p2 = new Point3d(axis);
		p2.scaleAdd(radius*0.95, center);

		if (n == 2) {
			s.append(getLineJmol(index, p1, axis, principalAxis, referenceAxis, color));
			s.append(getLineJmol(index + n + 1, p2, axis, principalAxis, referenceAxis, color));
		} else if (n > 2) {
			double polygonRadius = 3;
			s.append(getPolygonJmol(index, p1, principalAxis, referenceAxis, axis, n, color, polygonRadius));
			s.append(getPolygonJmol(index + n + 1, p2, principalAxis, referenceAxis, axis, n, color, polygonRadius));
		}

		return s.toString();
	}
	
	private String drawAxis(int index, Point3d axis) {
		StringBuilder s = new StringBuilder();
		s.append("draw l");
		s.append(index);
		s.append(" ");
		s.append("line");
		s.append(" ");
		Point3d v1 = new Point3d(axis);
		Point3d v2 = new Point3d(axis);
		v2.negate();
		Matrix4d m = axisTransformation.getReverseTransformation();
		m.transform(v1);
		m.transform(v2);
		s.append(getJmolPoint(v1));
		s.append(getJmolPoint(v2));
//		s.append(" width 0.25 ");
		s.append(" color purple");
		s.append(";");

		return s.toString();
	}
	
	private static String getLineJmol(int index, Point3d center, Vector3d axis, Vector3d principalAxis, Vector3d referenceAxis, String color) {
		StringBuilder s = new StringBuilder();
		s.append("draw l");
		s.append(index);
		s.append(" ");
		s.append("line");
		s.append(" ");

		Vector3d[] vertexes = getPolygonVertices(principalAxis, axis, referenceAxis, center, 2, 3);
		// create vertex list
		for (Vector3d v: vertexes) {
			s.append("{");
			s.append(jMolFloat(v.x));
			s.append(" ");
			s.append(jMolFloat(v.y));
			s.append(" ");
			s.append(jMolFloat(v.z));
			s.append("}");
		}

		s.append(" width 1.0 ");
		s.append(" color ");
		s.append(color);
		s.append(";");

		return s.toString();
	}

	private static String getPolygonJmol(int index, Point3d center, Vector3d principalAxis, Vector3d referenceAxis, Vector3d axis, int n, String color, double radius) {
		StringBuilder s = new StringBuilder();
		s.append("draw p");
		s.append(index);
		s.append(" ");
		s.append("polygon");
		s.append(" ");
		s.append(n+1); 
		s.append(" ");

		s.append("{");
		s.append(jMolFloat(center.x));
		s.append(" ");
		s.append(jMolFloat(center.y));
		s.append(" ");
		s.append(jMolFloat(center.z));
		s.append("}");

		Vector3d[] vertexes = getPolygonVertices(principalAxis, axis, referenceAxis, center, n, radius);
		//	System.out.println("AxisTransformation: polygon: " + Arrays.toString(vertexes));
		// create vertex list
		for (Vector3d v: vertexes) {
			s.append("{");
			s.append(jMolFloat(v.x));
			s.append(" ");
			s.append(jMolFloat(v.y));
			s.append(" ");
			s.append(jMolFloat(v.z));
			s.append("}");
		}

		// create face list
		s.append(" ");
		s.append(n);
		s.append(" ");

		for (int i = 1; i <= n; i++) {
			s.append("[");
			s.append(0);
			s.append(" ");
			s.append(i);
			s.append(" ");
			if (i < n) {
				s.append(i+1);
			} else {
				s.append(1);
			}
			s.append(" ");
			s.append(7);
			s.append("]");
		}

		s.append(" mesh");
		//	s.append(" color translucent ");
		s.append(" color ");
		s.append(color);
		s.append(";");

		return s.toString();
	}
	private static Vector3d[] getPolygonVertices(Vector3d principalAxis, Vector3d axis, Vector3d referenceAxis, Point3d center, int n, double radius) {
		Vector3d perp = new Vector3d(axis);
		// if axis coincides with principal axis, use the reference axis to orient polygon
		// TODO need to check alignment with y-axis for all cases of symmetry
		if (Math.abs(axis.dot(principalAxis)) > 0.9) {
			perp.set(referenceAxis);
		} else {
			perp.cross(perp, principalAxis);
		}
		//	System.out.println("AxisTransformation: perp. axis: " + referenceAxis);
		perp.scale(radius);		

		AxisAngle4d axisAngle = new AxisAngle4d(axis, 0);
		Vector3d[] vectors = new Vector3d[n];
		Matrix4d m = new Matrix4d();

		for (int i = 0; i < n; i++) {
			axisAngle.angle = i * 2 * Math.PI/n;
			vectors[i] = new Vector3d(perp);		
			m.set(axisAngle);
			m.transform(vectors[i]);
			vectors[i].add(center);
		}
		return vectors;
	}
	
	/**
	 * Returns the fold of the highest order rotation axis
	 * @return highest order rotation
	 */
	private int getPrincipalRotationAxisFold() {
		if (rotationGroup.getOrder() == 0) {
			return 1;
		}
		Rotation rotation = rotationGroup.getRotation(0); // the rotation around the principal axis is the first rotation
		return rotation.getFold();
	}
	
	/**
	 * Returns the fold of a minor order rotation axis
	 * @return minor order rotation
	 */
	private int getMinorRotationAxisFold() {
		// find axis that is not the rotation principal axis (direction = 1)
		for (int i = 0; i < rotationGroup.getOrder(); i++) {
			if (rotationGroup.getRotation(i).getDirection() == 1) {
				return rotationGroup.getRotation(i).getFold();
			}
		}
		return 1;
	}
	
	/**
	 * Returns a lower precision floating point number for Jmol
	 * @param f
	 * @return
	 */
	private static float jMolFloat(double f) {
		if (Math.abs(f) < 1.0E-7) {
			return 0.0f;
		}
		return (float)f;
	}
}
