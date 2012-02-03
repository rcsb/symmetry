/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava3.structure.align.symm.quarternary;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.Vector3d;

/**
 *
 * @author Peter
 */
public class RotationGroup {
     private List<Rotation> rotations = new ArrayList<Rotation>();
     private int principalAxisIndex = 0;
     private int higherOrderRotationAxis = 0;
     private int twoFoldsPerpendicular = 0;
     private int highestOrder = 0;
     private String pointGroup = "C1";
     private boolean complete = true;
     private boolean modified = true;

     public int getOrder() {
         return rotations.size();
     }

     public Rotation getRotation(int index) {
         return rotations.get(index);
     }

    public void addRotation(Rotation rotation) {
        rotations.add(rotation);
        modified = true;
    }
    
    public void removeRotation(int index) {
    	rotations.remove(index);
    	modified = true;
    }

    public String getPointGroup() {
        if (modified) {
            if (rotations.size() == 0) {
                return "C1";
            }
            findHighestOrderAxis();
            setEAxis();
            calcAxesDirections();
            findHigherOrderAxes();
            findTwoFoldsPerpendicular();
            calcPointGroup();
            modified = false;
        }
        return pointGroup;
    }
    
    public double getAverageSubunitRmsd() {
    	if (rotations.size() < 2) {
    		return 0.0;
    	}
    	double rmsd = 0;
    	// note, this loop starts at 1, because we don't take into account 
    	// RMSD of first operation (E)
    	for (int i = 1; i < rotations.size(); i++) {
    		rmsd += rotations.get(i).getSubunitRmsd();
    	}
    	return rmsd/(rotations.size()-1);
    }
    
    public double getAverageTraceRmsd() {
    	if (rotations.size() < 2) {
    		return 0.0;
    	}
    	double rmsd = 0;
    	// note, this loop starts at 1, because we don't take into account 
    	// RMSD of first operation (E)
    	for (int i = 1; i < rotations.size(); i++) {
    		double r = rotations.get(i).getTraceRmsd();
    		// if any invalid rmsd is found, stop
    		if (r < 0.0) {
    			return -1.0;
    		}
    		rmsd += rotations.get(i).getTraceRmsd();
    	}
    	return rmsd/(rotations.size()-1);
    }
    
    public double getAverageTraceGtsMin() {
    	if (rotations.size() < 2) {
    		return 100.0;
    	}
    	double gts = 0;
    	// note, this loop starts at 1, because we don't take into account 
    	// gts of first operation (E)
    	for (int i = 1; i < rotations.size(); i++) {
    		gts += rotations.get(i).getTraceGtsMin();
    	}
    	return gts/(rotations.size()-1);
    }
    
    public boolean isComplete() {
    	return complete;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Rotations: " + rotations.size() + "\n");
        for (Rotation s: rotations) {
            sb.append(s.toString() + "\n");
        }
        return sb.toString();
    }

    private void findHighestOrderAxis() {
        highestOrder = 1;
        principalAxisIndex = 0;
        double gts = 0;
        for (int i = 0; i < rotations.size(); i++) {
            Rotation s = rotations.get(i);
            if (s.getFold() > highestOrder) {
                highestOrder = s.getFold();
                principalAxisIndex = i;
                gts = s.getTraceGtsMin();
            } else if (s.getFold() >= highestOrder && s.getTraceGtsMin() > gts) {
                highestOrder = s.getFold();
                principalAxisIndex = i;
                gts = s.getTraceGtsMin();
            }
        }
//        System.out.println("highestOrderAxis: " + highestOrder);
//        System.out.println("principalAxisIndex: " + principalAxisIndex);
    }

    /**
     * Add E operation to the highest order rotation axis. By definition
     * E belongs to the highest order axis (see Papula p ??).
     */
    private void setEAxis() {
        Rotation e = rotations.get(0);
        Rotation h = rotations.get(principalAxisIndex);
        e.setAxisAngle(h.getAxisAngle());
    }

    private void findHigherOrderAxes() {
        higherOrderRotationAxis = 0;
        for (Rotation s: rotations) {
            if (s.getFold() > 2 && s.getDirection() == 1) {
               higherOrderRotationAxis++;
            }
        }
//        System.out.println("higherOrderRotationAxis: " + higherOrderRotationAxis);
    }

    private void calcAxesDirections() {
        if (highestOrder == 1) {
              for (Rotation s: rotations) {
                  s.setDirection(0);
              }
              return;
        }
        
        AxisAngle4d pa = rotations.get(principalAxisIndex).getAxisAngle();
        Vector3d pv = new Vector3d(pa.getX(), pa.getY(), pa.getZ());

        for (Rotation s: rotations) {
           AxisAngle4d axis = s.getAxisAngle();
           Vector3d av = new Vector3d(axis.getX(), axis.getY(), axis.getZ());
           if (Math.abs(pv.dot(av)) > 0.9f) {
               // co-linear with principal axis
               s.setDirection(0);
//               System.out.println("Axis co-linear: " + axis + " - " + pa);
           } else {
               // not co-linear or perpendicular to principal axis
               s.setDirection(1);
  //             System.out.println("Axis perpendicular: " + axis + " - " + pa);
           }
        }
        rotations.get(0).setDirection(0); // set the E axis to the principal axis (by definition)
    }

    private void findTwoFoldsPerpendicular() {
        twoFoldsPerpendicular = 0;
        // s.getFold() == 0 -> include E
        for (Rotation s: rotations) {
            if (s.getFold() == 2 && s.getDirection() == 1) {
                twoFoldsPerpendicular++;
            }
        }
//        System.out.println("twoFoldsPerpendicular: " + twoFoldsPerpendicular);
    }

    private void calcPointGroup() {
        if (higherOrderRotationAxis > 1) {
            // cubic groups
            if (highestOrder == 5) {
                // rotational icosahedral symmetry or chiral icosahedral symmetry
                pointGroup = "I";
                complete = rotations.size() == 60;
                return;
            } else if (highestOrder == 4) {
                // rotational octahedral symmetry or chiral octahedral symmetry
                pointGroup = "O";
                complete = rotations.size() == 24;
                return;
            } else if (highestOrder == 3) {
                // rotational tetrahedral symmetry or chiral tetrahedral symmetry
                pointGroup = "T";
                complete = rotations.size() == 12;
                return;
            }
        } else {
            // Cn and Dn groups
            // if E is not counted, subtract 1
            if (Math.abs(twoFoldsPerpendicular - highestOrder) <= 1 && highestOrder > 1) {
                pointGroup = "D" + highestOrder;
                complete = rotations.size() == 2 * highestOrder;
                return;
            } else {
                pointGroup = "C" + highestOrder;
                complete = rotations.size() == highestOrder;
                return;
            }
        }

        pointGroup = "C1";
        return;
    }

    public void sortByFoldDecending() {
        Collections.sort(rotations, new Comparator<Rotation>() {
			public int compare(Rotation o1, Rotation o2) {
				// check this ???
				return Math.round(Math.signum(o1.getFold() - o2.getFold()));
			}
        });
    }
}
