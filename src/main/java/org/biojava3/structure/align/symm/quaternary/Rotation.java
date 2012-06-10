
package org.biojava3.structure.align.symm.quaternary;

import java.util.List;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;

/**
 *
 * @author Peter
 */
public class Rotation {
    private double subunitRmsd = Double.MAX_VALUE;
    private double traceRmsd = Double.MAX_VALUE;
    private List<Integer> permutation;
    private Matrix4d transformation;
    private AxisAngle4d axisAngle;
    private int direction;
    private int fold;

    /**
     * @return the subunitRmsd
     */
    public double getSubunitRmsd() {
        return subunitRmsd;
    }

    /**
     * @param subunitRmsd the subunitRmsd to set
     */
    public void setSubunitRmsd(double subunitRmsd) {
        this.subunitRmsd = subunitRmsd;
    }

    /**
     * @return the traceRmsd
     */
    public double getTraceRmsd() {
        return traceRmsd;
    }

    /**
     * @param traceRmsd the traceRmsd to set
     */
    public void setTraceRmsd(double traceRmsd) {
        this.traceRmsd = traceRmsd;
    }

    /**
     * @return the permutation
     */
    public List<Integer> getPermutation() {
        return permutation;
    }

    /**
     * @param permutation the permutation to set
     */
    public void setPermutation(List<Integer> permutation) {
        this.permutation = permutation;
    }

    /**
     * @return the transformation
     */
    public Matrix4d getTransformation() {
        return transformation;
    }

    /**
     * @param transformation the transformation to set
     */
    public void setTransformation(Matrix4d transformation) {
        this.transformation = transformation;
    }

    /**
     * @return the fold
     */
    public int getFold() {
        return fold;
    }

    /**
     * @param fold the fold to set
     */
    public void setFold(int fold) {
        this.fold = fold;
    }

    /**
     * @return the direction
     */
    public int getDirection() {
        return direction;
    }

    /**
     * @param direction the direction to set
     */
    public void setDirection(int axis) {
        this.direction = axis;
    }

    /**
     * @return the axisAngle
     */
    public AxisAngle4d getAxisAngle() {
        return axisAngle;
    }

    /**
     * @param axisAngle the axisAngle to set
     */
    public void setAxisAngle(AxisAngle4d axisAngle) {
        this.axisAngle = axisAngle;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("fold: " + fold);
        sb.append(" orientation: " + direction);
        sb.append(" RMSD: " + subunitRmsd);
        sb.append(" permutation: " + permutation);
        return sb.toString();
    }
}
