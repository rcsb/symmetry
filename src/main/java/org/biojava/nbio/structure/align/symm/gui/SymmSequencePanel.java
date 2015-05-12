package org.biojava.nbio.structure.align.symm.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.biojava.bio.structure.StructureTools;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.JPrintPanel;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.aligpanel.AFPChainCoordManager;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;


/** 
 * A JPanel that can display the subunit colored sequence of a symmetry analzsis in a nice way and interact with Jmol.
 * 
 * @author Aleix Lafita
 *
 */
public class SymmSequencePanel  extends JPrintPanel implements AlignmentPositionListener, WindowListener{

   private static final long serialVersionUID = -6892229111166263764L;

   private AFPChain afpChain;
   private AFPChainCoordManager coordManager;
   private Font seqFont;
   private SymmetryJmol jmol;
   private SymmSequencePanelMouseMotionListener mouseMoLi;

   private BitSet selection;

   private boolean selectionLocked;
   private Atom[] ca1;
   private Color[] subunitColors;

   public SymmSequencePanel(){
      super();
      this.setBackground(Color.white);
      coordManager = new AFPChainCoordManager();
      seqFont = new Font("SansSerif",Font.PLAIN,12);

      mouseMoLi = new SymmSequencePanelMouseMotionListener(this);
      this.addMouseMotionListener(mouseMoLi);
      this.addMouseListener(mouseMoLi);
      mouseMoLi.addAligPosListener(this);

      selection = new BitSet();
   }

   public AFPChainCoordManager getCoordManager() {
      return coordManager;
   }

   public void addAlignmentPositionListener(AlignmentPositionListener li){
      mouseMoLi.addAligPosListener(li);
   }

   public void destroy(){

      afpChain = null;;
      mouseMoLi.destroy();	
      jmol = null;
      ca1 = null;
      selection = null;
   }

   public AFPChain getAFPChain(){
      return afpChain;
   }

   public void setAFPChain(AFPChain afpChain) {

      this.afpChain = (AFPChain) afpChain.clone();
      //Change the alignment sequence to display (delete gaps basically, we only want the sequence)
      int length = afpChain.getOptAln()[afpChain.getBlockNum()-1][0][afpChain.getOptLen()[0]-1] - afpChain.getOptAln()[0][0][0];
      char[] seq = new char[length];
      for (int i=0; i<length; i++)
    	  seq[i] = StructureTools.get1LetterCode(ca1[afpChain.getOptAln()[0][0][0]+i].getGroup().getPDBName());
      this.afpChain.setAlnLength(length);
      this.afpChain.setAlnseq1(seq);
      this.afpChain.setAlnseq2(seq);
      this.afpChain.setAlnsymb(seq);
      
      coordManager.setAFPChain(this.afpChain);
      if (afpChain != null) {
         selection = new BitSet(afpChain.getAlnLength());
      }
   }


@Override
public void paintComponent(Graphics g){

      super.paintComponent(g);

      Graphics2D g2D = (Graphics2D) g;
      g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
            RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

      char[] seq = afpChain.getAlnseq1();
      
      String summary = afpChain.getName1();
      g2D.drawString(summary, 20, coordManager.getSummaryPos());
      
      // draw a darker background
      Rectangle sig = new Rectangle(10,10,10,10);				
      g2D.fill(sig);
      
      int blockNum = afpChain.getBlockNum();

      List<Integer> boundaries = new ArrayList<Integer>();
      List<Integer> alignedPos = new ArrayList<Integer>();
      for (int bk=0; bk<blockNum; bk++){
     	 for (int res=0; res<afpChain.getOptLen()[bk]; res++){
     		 alignedPos.add(afpChain.getOptAln()[bk][0][res]);
     		 if (res==afpChain.getOptLen()[bk]-1) boundaries.add(afpChain.getOptAln()[bk][0][res]);
     	 }
      }
      
      int colorPos = 0;
      for ( int i = 0 ; i < afpChain.getOptLength() ;i++){

    	 int res = afpChain.getOptAln()[0][0][0]+i;
    	 if (res > boundaries.get(colorPos)) colorPos++;
    	 
         char c = seq[i];
         g2D.setFont(seqFont);

         Point p1 = coordManager.getPanelPos(0,i);
         int xpos1 = p1.x;
         int ypos1 = p1.y;
         
         Color bg = Color.white;

         if (!alignedPos.contains(res)) bg = Color.white;
         else bg = subunitColors[colorPos];
            
        // draw a darker background
        g2D.setPaint(bg);
        Rectangle rec = new Rectangle(p1.x-1,p1.y-11, 12, 12);
        g2D.fill(rec);

         if ( isSelected(i)){
            // draw selection
            Color bgn = Color.YELLOW;
            g2D.setPaint(bgn);
            // draw a darker backgroun
            Rectangle recn = new Rectangle(p1.x-1,p1.y-11, 12, 12);
            g2D.fill(recn);
         }

         // draw the AA sequence
         g2D.setColor(Color.black);
         g2D.drawString(c+"",xpos1,ypos1);
      }

      int nrLines = (afpChain.getAlnLength() -1) / AFPChainCoordManager.DEFAULT_LINE_LENGTH;


      for (int i = 0; i <= nrLines; i++){

         try {
            // draw legend at i
            Point p1 = coordManager.getLegendPosition(i,0);


            int aligPos = i * AFPChainCoordManager.DEFAULT_LINE_LENGTH;
            Atom a1 = DisplayAFP.getAtomForAligPos(afpChain, 0,aligPos, ca1,false);
            String label1 = JmolTools.getPdbInfo(a1,false);				
            g2D.drawString(label1, p1.x,p1.y);

            Point p3 = coordManager.getEndLegendPosition(i,0);

            aligPos = i * AFPChainCoordManager.DEFAULT_LINE_LENGTH + AFPChainCoordManager.DEFAULT_LINE_LENGTH -1 ;
            if (aligPos > afpChain.getAlnLength())
               aligPos = afpChain.getAlnLength() - 1;
            Atom a3 = DisplayAFP.getAtomForAligPos(afpChain, 0,aligPos, ca1,true);

            String label3 = JmolTools.getPdbInfo(a3,false);

            g2D.drawString(label3, p3.x,p3.y);


         } catch (StructureException e){
            e.printStackTrace();
         }
      }
   }

 


   protected boolean isSelected(int alignmentPosition) {

      return selection.get(alignmentPosition);

   }


   @Override
public void mouseOverPosition(AlignedPosition p) {
      //System.out.println("AligPanel: mouse over position " + p.getPos1() );

      if ( ! selectionLocked)
         selection.clear();
      selection.set(p.getPos1());

      updateJmolDisplay();

      this.repaint();

   }

   private void updateJmolDisplay() {

      if ( jmol == null)
         return;

      int size = ca1.length;

      StringBuffer cmd = new StringBuffer("select ");

      int nrSelected = 0;
      try {

         for (int i = 0 ; i< size ; i++){
            if ( selection.get(i)){

               Atom a1 = DisplayAFP.getAtomForAligPos(afpChain, 0,i, ca1, false);

               String select1 = "";

               if ( a1 != null ) 
                  select1 = JmolTools.getPdbInfo(a1);

               // nothing to display
               if ( select1.equals(""))
                  continue;

               if ( nrSelected > 0)
                  cmd.append(", ");

               cmd.append(select1);
               nrSelected++;
            }
         }


      } catch (StructureException e){
         e.printStackTrace();
      }
      if ( nrSelected == 0)
         cmd.append(" none;");
      else
         cmd.append("; set display selected;");

      jmol.evalString(cmd.toString());


   }


   @Override
public void positionSelected(AlignedPosition p) {
      mouseOverPosition(p);

   }

   @Override
public void rangeSelected(AlignedPosition start, AlignedPosition end) {
      //System.out.println("AligPanel: range selected " + start.getPos1() + " - " + end.getPos1() + " selectionLockedL " + selectionLocked);
      if ( ! selectionLocked )
         selection.clear();
      selection.set(start.getPos1(), end.getPos1()+1);
      updateJmolDisplay();
      this.repaint();

   }

   @Override
public void selectionLocked() {
      selectionLocked = true;

   }

   @Override
public void selectionUnlocked() {
      selectionLocked = false;
      selection.clear();
      this.repaint();

   }


   @Override
public void toggleSelection(AlignedPosition p) {
      selection.flip(p.getPos1());
      //System.out.println("AligPanel: toggle selection " + p.getPos1() + " " + selection.get(p.getPos1()));
      updateJmolDisplay();
      this.repaint();

   }



   public void setStructureAlignmentJmol(SymmetryJmol jmol) {
      this.jmol = jmol;

   }


   @Override
public void windowActivated(WindowEvent e) {

      // TODO Auto-generated method stub

   }


   @Override
public void windowClosed(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowClosing(WindowEvent e) {
      destroy();

   }


   @Override
public void windowDeactivated(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowDeiconified(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowIconified(WindowEvent e) {
      // TODO Auto-generated method stub

   }


   @Override
public void windowOpened(WindowEvent e) {
      // TODO Auto-generated method stub

   }

   @Override
public void actionPerformed(ActionEvent e) {
      String cmd = e.getActionCommand();
      // print is handled by superclass
      if ( cmd.equals(MenuCreator.PRINT)) {
         super.actionPerformed(e);
      } else if ( cmd.equals(MenuCreator.SELECT_EQR)){
         selectEQR();
      }else {
         System.err.println("Unknown command:" + cmd);
      }

   }

   private void selectEQR() {

      selection.clear();

      List<Integer> pos1 = DisplayAFP.getEQRAlignmentPos(afpChain);

      for (int pos : pos1){
         selection.flip(pos);
      }
      mouseMoLi.triggerSelectionLocked(true);
      updateJmolDisplay();
      this.repaint();

   }

   public Atom[] getCa1() {
      return ca1;
   }


   public void setCa1(Atom[] ca1) {
      this.ca1 = ca1;
   }

	public Color[] getSubunitColors() {
		return subunitColors;
	}
	
	public void setSubunitColors(Color[] subunitColors) {
		this.subunitColors = subunitColors;
	}


}


