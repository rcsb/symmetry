package org.biojava3.ramachandran;
import java.awt.RenderingHints;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.StructureName;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.secstruc.SecStruc;
import org.biojava.bio.structure.secstruc.SecStrucGroup;
import org.biojava.bio.structure.secstruc.SecStrucState;
import org.biojava3.structure.StructureIO;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.FastScatterPlot;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

/**
 * A demo of the fast scatter plot.
 *
 */
public class FastScatterPlotDemo extends ApplicationFrame {

	/** A constant for the number of items in the sample dataset. */
	//private static final int COUNT = 500000;

	/** The data. */
	private float[][] data = new float[2][0];

	/**
	 * Creates a new fast scatter plot demo.
	 *
	 * @param title  the frame title.
	 */
	public FastScatterPlotDemo(final String title) {

		super(title);
		populateData();
		final NumberAxis domainAxis = new NumberAxis("X");
		domainAxis.setAutoRangeIncludesZero(false);
		final NumberAxis rangeAxis = new NumberAxis("Y");
		rangeAxis.setAutoRangeIncludesZero(false);
		final FastScatterPlot plot = new FastScatterPlot(this.data, domainAxis, rangeAxis);
		final JFreeChart chart = new JFreeChart("Fast Scatter Plot", plot);
		//        chart.setLegend(null);

		// force aliasing of the rendered content..
		chart.getRenderingHints().put
		(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

		final ChartPanel panel = new ChartPanel(chart, true);
		panel.setPreferredSize(new java.awt.Dimension(500, 270));
		//      panel.setHorizontalZoom(true);
		//    panel.setVerticalZoom(true);
		panel.setMinimumDrawHeight(10);
		panel.setMaximumDrawHeight(2000);
		panel.setMinimumDrawWidth(20);
		panel.setMaximumDrawWidth(2000);

		setContentPane(panel);

	}

	// ****************************************************************************
	// * JFREECHART DEVELOPER GUIDE                                               *
	// * The JFreeChart Developer Guide, written by David Gilbert, is available   *
	// * to purchase from Object Refinery Limited:                                *
	// *                                                                          *
	// * http://www.object-refinery.com/jfreechart/guide.html                     *
	// *                                                                          *
	// * Sales are used to provide funding for the JFreeChart project - please    * 
	// * support us so that we can continue developing free software.             *
	// ****************************************************************************

	/**
	 * Populates the data array with random values.
	 */
	private void populateData() {

		//        for (int i = 0; i < this.data[0].length; i++) {
		//            final float x = (float) i + 100000;
		//            this.data[0][i] = x;
		//            this.data[1][i] = 100000 + (float) Math.random() * COUNT;
		//        }

		AtomCache cache = new AtomCache();
		
		SortedSet<String> repres = JFatCatClient.getRepresentatives("http://emmy.rcsb.org/jfatcatserver/align/", 40);

		List<float[]> values = new ArrayList<float[]>();
		
		List<String> knownStructures = new ArrayList<String>();
		int count = -1;
		for (String repre : repres){
			count++;
//			if ( count > 800)
//				break;
			System.out.println("#" +count + "/" + repres.size() + " " + repre + " " + knownStructures.size());

			try {
				StructureName key = new StructureName(repre);
				if ( knownStructures.contains(key.getPdbId()))
					continue;
				knownStructures.add(key.getPdbId());
				Structure s = StructureIO.getStructure(key.getPdbId());

				SecStruc ss = new SecStruc();
				ss.assign(s);
				//System.out.println(ss);
				SecStrucGroup[] groups = ss.getGroups(); 


				
				for ( SecStrucGroup g : groups){
				
					SecStrucState state = (SecStrucState) g.getProperty("secstruc");

					float[] val = new float[2];
					val[0] = (float)state.getPhi();
					val[1] = (float) state.getPsi();
					values.add(val);
				}
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		data = new float[2][values.size()];
		int i = 0;
		for ( float[] v : values){ 
			data[0][i]=v[0];
			data[1][i]=v[1];
			i++;
		}




	}



	/**
	 * Starting point for the demonstration application.
	 *
	 * @param args  ignored.
	 */
	public static void main(final String[] args) {

		final FastScatterPlotDemo demo = new FastScatterPlotDemo("Fast Scatter Plot Demo");
		demo.pack();
		RefineryUtilities.centerFrameOnScreen(demo);
		demo.setVisible(true);

	}

}
