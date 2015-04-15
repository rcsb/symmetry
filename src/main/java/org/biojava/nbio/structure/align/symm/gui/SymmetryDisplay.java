package org.biojava.nbio.structure.align.symm.gui;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.jcolorbrewer.ColorBrewer;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

/**
 * Code adapted from StructureAlignmentDisplay class of biojava.
 * 
 * @author lafita
 * 
 */
public class SymmetryDisplay extends StructureAlignmentDisplay{

   
   /** Display the alignment
    * 
    * @param afpChain
    * @param ca1
    * @param ca2
    * @return a StructureAlignmentJmol instance
    * @throws StructureException
    */
   public static SymmetryJmol display(AFPChain afpChain, Atom[] ca1, Atom[] ca2, Color[] subunitColors) throws StructureException {
      
      if ( ca1.length < 1 || ca2.length < 1){
         throw new StructureException("length of atoms arrays is too short! " + ca1.length + "," + ca2.length);
      }
      
      Group[] twistedGroups = prepareGroupsForDisplay(afpChain, ca1, ca2);
            
      List<Group> hetatms  = new ArrayList<Group>();
      List<Group> nucs1    = new ArrayList<Group>();
      Group g1 = ca1[0].getGroup();
      Chain c1 = null;
      if ( g1 != null) {
         c1 = g1.getChain();
         if ( c1 != null){
            hetatms = c1.getAtomGroups(GroupType.HETATM);;
            nucs1  = c1.getAtomGroups(GroupType.NUCLEOTIDE);
         }
      }
      List<Group> hetatms2 = new ArrayList<Group>();
      List<Group> nucs2    = new ArrayList<Group>();
      Group g2 = ca2[0].getGroup();
      Chain c2 = null;
      if ( g2 != null){
         c2 = g2.getChain();
         if ( c2 != null){
            hetatms2 = c2.getAtomGroups(GroupType.HETATM);
            nucs2 = c2.getAtomGroups(GroupType.NUCLEOTIDE);
         }
      }
         
      return DisplaySymmAFP.display(afpChain, twistedGroups, ca1, ca2,hetatms, nucs1, hetatms2, nucs2, subunitColors);
      
   }
   
   public static SymmetryJmol display(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureException {
	   
	    //Set the color to the DEFAULT CeSymm Colors
	    Color[] subunitColors = null;
	   	CESymmParameters.SubunitColors COLOR = CESymmParameters.SubunitColors.DEFAULT;
	   	
		switch(COLOR){
		case COLOR_SET: 
			subunitColors = ColorBrewer.Set1.getColorPalette(afpChain.getBlockNum());
			break;
		case SPECTRAL:
			subunitColors = ColorBrewer.Spectral.getColorPalette(afpChain.getBlockNum());
			break;
		case PAIRED:
			subunitColors = ColorBrewer.Paired.getColorPalette(afpChain.getBlockNum());
		case GRADUAL:
			break;
		}
	   
	   return display(afpChain,ca1,ca2,subunitColors);
   }
}
