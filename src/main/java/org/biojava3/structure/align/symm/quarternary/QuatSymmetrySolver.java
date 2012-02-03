
package org.biojava3.structure.align.symm.quarternary;

/**
 *
 * @author Peter
 */
public interface QuatSymmetrySolver {

    RotationGroup getSymmetryOperations();

    void setSubunitRmsdThreshold(double subunitRmsdThreshold);
    
    void setRmsdThreshold(double rmsdThreshold);

    void setGtsThreshold(double gtsThreshold);
    
    void setPseudoSymmetryAllowed(boolean pseudoSymmetryAllowed);
}
