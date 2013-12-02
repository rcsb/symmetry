package org.biojava3.ramachandran;

import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import javax.swing.JPanel;

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
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.xy.XYShapeRenderer;
import org.jfree.chart.title.PaintScaleLegend;
import org.jfree.data.Range;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XYZDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RefineryUtilities;

/**
 * A simple demo for the {@link XYShapeRenderer} class.
 */
public class XYShapeRendererDemo1 extends ApplicationFrame {

	/**
	 * 
	 */
	private static final long serialVersionUID = -8048818543969127101L;

	/**
	 * A demonstration application showing a bubble chart.
	 *
	 * @param title  the frame title.
	 */
	public XYShapeRendererDemo1(String title) {
		super(title);
		JPanel chartPanel = createDemoPanel();
		chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
		setContentPane(chartPanel);
	}

	/**
	 * Creates a chart.
	 *
	 * @param dataset  the dataset.
	 *
	 * @return The chart.
	 */
	private static JFreeChart createChart(XYZDataset dataset) {

		// set up the Axis

		NumberAxis xAxis = new NumberAxis("Phi");
		xAxis.setAutoRange(false);
		xAxis.setRange(new Range(-180, 180));


		//xAxis.setAutoRangeIncludesZero(false);

		NumberAxis yAxis = new NumberAxis("Psi");
		yAxis.setAutoRange(false);
		yAxis.setRange(new Range(-180, 180));
		//yAxis.setAutoRangeIncludesZero(false);

		NumberAxis zAxis = new NumberAxis("Score");
		zAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());


		//XYItemRenderer renderer = new XYDotRenderer();
		XYShapeRenderer renderer = new XYShapeRenderer();

		System.out.println(renderer.getBaseShape().getBounds());
		Ellipse2D el = new Ellipse2D.Float();
		//System.out.println(s.);
		el.setFrame(-2,-2,4,4);
		renderer.setBaseShape(el);
		System.out.println(renderer.getBaseShape().getBounds());

		LookupPaintScale ps = new LookupPaintScale(0, 100, new Color(200, 200, 255));
		for ( int i=0;i<=80;i++){
			
			int r = (int)Math.round(255/100.0*i);
			int g = (int)Math.round(255/100.0*i);
			int b = 255;
			//System.out.println(i+ " " + r + " " + g + " " + b);
			
			ps.add(100-i, new Color(r,g,b));
					

		}

		//GrayPaintScale ps = new GrayPaintScale(0,200);


		renderer.setPaintScale(ps);





		XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);

		//plot.setDomainPannable(true);
		//plot.setRangePannable(true);




		JFreeChart chart = new JFreeChart("Ramachandran plot", plot);
		chart.removeLegend();


		PaintScaleLegend psl = new PaintScaleLegend(ps, zAxis);
		psl.setPosition(RectangleEdge.RIGHT);
		psl.setMargin(4, 4, 40, 4);
		psl.setAxisLocation(AxisLocation.BOTTOM_OR_RIGHT);
		chart.addSubtitle(psl);

		//ChartUtilities.applyCurrentTheme(chart);



		return chart;
	}

	/**
	 * Creates a sample dataset.
	 *
	 * @return A sample dataset.
	 */
	public static XYZDataset createDataset() {

		AtomCache cache = new AtomCache();
		

		SortedSet<String> repres = JFatCatClient.getRepresentatives("http://emmy.rcsb.org/jfatcatserver/align/", 40);

		List<String> knownStructures = new ArrayList<String>();
		int count = -1;

		TreeSet<ScoreMap> maps = new TreeSet<ScoreMap>();

		double maxCount = 0;

		for (String repre : repres){

			count++;

			//if ( count > 2000)
				//break;

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

					ScoreMap map = new ScoreMap();

					map.setPhi( (int) Math.round(state.getPhi()*100/100.0d));
					map.setPsi( (int) Math.round(state.getPsi()*100/100.0d));
					map.incrementCount();
					if ( map.getPhi() ==0 && map.getPsi()==0)
						continue;
					
					if ( maps.contains(map)){
						
						ScoreMap m = maps.floor(map);
						m.incrementCount();
						//System.out.println("re-using old map " + m);
						if ( m.getCount() > maxCount)
							maxCount = m.getCount();
					} else {
						maps.add(map);
					}
				}
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		
		//System.out.println(maps);
		
		DefaultXYZDataset dataset = new DefaultXYZDataset();
		double[]x = new double[maps.size()];
		double[]y = new double[maps.size()];
		double[]z = new double[maps.size()];

		System.out.println("maxCount:" + maxCount );
		int i =-1;
		for (ScoreMap m : maps){
			i++;
			x[i] = m.getPhi();
			y[i] = m.getPsi();
			z[i] = m.getCount()/(double)maxCount *100.0;

			//if ( m.getCount() > maxCount)
			//	maxCount = m.getCount();

		}

		System.out.println("maxCount:" + maxCount );


		double[][] series = new double[][] { x, y, z };
		dataset.addSeries("Series 1", series);

		return dataset;
	}

	/**
	 * Creates a panel for the demo (used by SuperDemo.java).
	 *
	 * @return A panel.
	 */
	public static JPanel createDemoPanel() {
		JFreeChart chart = createChart(createDataset());
		ChartPanel chartPanel = new ChartPanel(chart);
		// chartPanel.setMouseWheelEnabled(true);
		return chartPanel;
	}

	/**
	 * Starting point for the demonstration application.
	 *
	 * @param args  ignored.
	 */
	public static void main(String[] args) {
		XYShapeRendererDemo1 demo = new XYShapeRendererDemo1(
				"Ramachandran Plot");
		demo.pack();
		RefineryUtilities.centerFrameOnScreen(demo);
		demo.setVisible(true);
		demo.requestFocus();
	}

}
