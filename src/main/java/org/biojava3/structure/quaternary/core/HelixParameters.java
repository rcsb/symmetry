/**
 * 
 */
package org.biojava3.structure.quaternary.core;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

/**
 * @author Peter
 *
 */
public class HelixParameters {
	private List<List<Integer>> layerLines = null;
	private Matrix4d transformation = null;
	private double rmsd = 0;
	private int interactions = 0;
	private Matrix4d alignmentMatrix = null;

    public HelixParameters (Matrix4d transformation, double rmsd, int interactions, List<List<Integer>> layerLines) {
    	this.transformation = transformation;
    	this.rmsd = rmsd;
    	this.interactions = interactions;
    	this.layerLines = layerLines;
    }

	public int getStart() {
		return this.layerLines.size();
	}

	public List<List<Integer>> getLayerLines() {
		return this.layerLines;
	}

	public Matrix4d getTransformation() {
		return this.transformation;
	}
    
	public Vector3d getAxisVector() {
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(this.transformation);
		Vector3d axisVector = new Vector3d(axis.x, axis.y, axis.z);
		axisVector.normalize();
		return axisVector;
	}
	
	public double getAngle() {
		AxisAngle4d axis = new AxisAngle4d();
		axis.set(this.transformation);
		return axis.angle;
	}
	
	public double getRise() {
		Vector3d translation = new Vector3d();
		this.transformation.get(translation);
		return translation.length();
	}
	
	public double getRmsd() {
		return this.rmsd;
	}
	
	public int getInteractions() {
		return this.interactions;
	}
	
	public Matrix4d getAlignmentMatrix() {
		return this.alignmentMatrix;
	}
	
	public void setAlignmentMatrix(Matrix4d alignmentMatrix) {
		this.alignmentMatrix = alignmentMatrix;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getStart());
		sb.append("-helix: angle=");
		sb.append(Math.toDegrees(getAngle()));
		sb.append(" rise=");
		sb.append(getRise());
		sb.append(" rmsd=");
		sb.append(getRmsd());
		sb.append(" interactions=");
		sb.append(getInteractions());
		sb.append(" layer lines=");
		sb.append(getLayerLines());
		return sb.toString();
	}
}
