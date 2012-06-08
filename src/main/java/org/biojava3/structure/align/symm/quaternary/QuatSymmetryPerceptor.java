/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.structure.align.symm.quaternary;

/**
 *
 * @author Peter
 */
public class QuatSymmetryPerceptor {
    private double subunitRmsdThreshold = 5.0f;
    private double rmsdThreshold = 5.0f;
    private double gtsThreshold = 0.0f;
    private boolean pseudoSymmetryAllowed = false;
    private int maxOrder = 60;
    private Subunits subunits = null;
    private RotationGroup symmetryOperations = new RotationGroup();
    private String method = "";

    public QuatSymmetryPerceptor(Subunits subunits) {
        this.subunits = subunits;
    }

    public void setSubunitRmsdThreshold(double subunitRmsdThreshold) {
        this.subunitRmsdThreshold = subunitRmsdThreshold;
    }
    
    public void setRmsdThreshold(double rmsdThreshold) {
        this.subunitRmsdThreshold = rmsdThreshold;
    }

    public void setGtsTreshold(double gtsThreshold) {
        this.gtsThreshold = gtsThreshold;
    }
    
    public void setPseudoSymmetryAllowed(boolean pseudoSymmetryAllowed) {
    	this.pseudoSymmetryAllowed = pseudoSymmetryAllowed;
    }
    
    public void setMaxOrder(int maxOrder) {
    	this.maxOrder = maxOrder;
    }
    
    public String getMethod() {
    	return method;
    }

    public RotationGroup getSymmetryOperations() {
        if (symmetryOperations.getOrder() == 0) {
            findSymmetryOperations();
            symmetryOperations.complete();
        }
        return symmetryOperations;
    }

    private void findSymmetryOperations() {
    	if (subunits.getSubunitCount() <= 1) {
    		symmetryOperations =  new RotationGroup();
    		return;
    	}
    	
        QuatSymmetrySolver solver = null;

        RotationGroup systematic = null;
        // For small systems use the systematic solver
        // which is more robust and faster, however, due
        // to its n! complexity becomes too slow for more
        // than 8 subunits.
        if (subunits.getSubunitCount() == 2) {
//        	System.out.println("C2 Rotation solver");
        	method = "C2rotation";
        	solver = new C2RotationSolver(subunits);
        	 solver.setRmsdThreshold(subunitRmsdThreshold);
             solver.setGtsThreshold(gtsThreshold);
             solver.setPseudoSymmetryAllowed(pseudoSymmetryAllowed);
             symmetryOperations = solver.getSymmetryOperations();
             addCalphaRmsd();
  //      } else if (subunits.getSubunitCount() < 9) {
        } else if (subunits.getSubunitCount() < 9) {
        	System.out.println("Systematic solver");
        	method = "systematic";
            solver = new SystematicSolver(subunits);
            solver.setRmsdThreshold(subunitRmsdThreshold);
            solver.setGtsThreshold(gtsThreshold);
            solver.setPseudoSymmetryAllowed(pseudoSymmetryAllowed);
            systematic = solver.getSymmetryOperations();
        } 
        if (subunits.getSubunitCount() > 2) {
        	System.out.println("Rotation solver");
        	method = "rotation";
            solver = new RotationSolver(subunits);
            ((RotationSolver) solver).setMaxOrder(maxOrder);
            solver.setRmsdThreshold(subunitRmsdThreshold);
            solver.setGtsThreshold(gtsThreshold);
            solver.setPseudoSymmetryAllowed(pseudoSymmetryAllowed);
            symmetryOperations = solver.getSymmetryOperations();
            addCalphaRmsd();
        }
        
        // uses the method that find the highest symmetry or use the method that doesn't over-predict symmetry (i.e,. 1YVK)
        if (systematic != null) {
  //      	if (systematic.getOrder() > symmetryOperations.getOrder()) {
        	if ((systematic.getOrder() > symmetryOperations.getOrder() && systematic.getOrder() <= subunits.getSubunitCount())
        			|| symmetryOperations.getOrder() > subunits.getSubunitCount()) {
        	   symmetryOperations = systematic;
        	   method = "systematic";
        	   addCalphaRmsd();
        	}
        }    
    }
    
    private void addCalphaRmsd() {
    	QuatSuperpositionScorer scorer = new QuatSuperpositionScorer(subunits);
    	for (int i = 0; i < symmetryOperations.getOrder(); i++) {
    		Rotation r = symmetryOperations.getRotation(i);
    		if (r.getTraceRmsd() == Double.MAX_VALUE) {
    			double rmsd = scorer.calcCalphaRMSD(r.getTransformation(), r.getPermutation());
    			if (rmsd < 0.0) {
    				System.out.println("Removing negative RMSD");
    				symmetryOperations.removeRotation(i);
    				i--;
    			}
    			r.setTraceRmsd(rmsd);	
    		}
    	}
    }
}
