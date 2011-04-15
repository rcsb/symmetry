package results;

import java.awt.Dimension;
import java.awt.RenderingHints;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.biojava.bio.structure.align.webstart.JNLPProxy;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.entity.ChartEntity;
import org.jfree.chart.labels.CustomXYToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.urls.CustomXYURLGenerator;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.rcsb.fatcat.server.DbUtils;
import org.rcsb.fatcat.server.PdbChainKey;
import org.rcsb.fatcat.server.dao.DBAlignment;
import org.rcsb.fatcat.server.dao.SplitDatabase;

public class PlotRMSDvsProbability {
	private static final String BASE_URL="http://beta.rcsb.org/pdb/workbench/showPrecalcAlignment.do?action=pw_fatcat&pdb1=%s&chain1=%s&pdb2=%s&chain2=%s";

	CustomXYToolTipGenerator toolTipGenerator;
	CustomXYURLGenerator urlGenerator;
	ChartPanel chartPanel;
	
	
	public static void main(String[] args){
		String name = "12AS.B";
		
		PdbChainKey repre = PdbChainKey.fromName(name);
		
		PlotRMSDvsProbability me = new PlotRMSDvsProbability();
		me.showData(repre);
	}
	
	public void showData(PdbChainKey repre){
		
		XYDataset series = getData4Repre(repre);
		JFreeChart chart = createChart(series);
		chartPanel = new ChartPanel(chart);
		chartPanel.setFillZoomRectangle(true);
		
		chartPanel.setPreferredSize(new Dimension(500, 270));
		
		chartPanel.addChartMouseListener(new ChartMouseListener() {
			

			public void chartMouseMoved(ChartMouseEvent event) {
				// TODO Auto-generated method stub

			}

			public void chartMouseClicked(ChartMouseEvent event) {

				JFreeChart chart = event.getChart();
				int x = event.getTrigger().getX();
				int y = event.getTrigger().getY();
				//event.get
				ChartEntity entity = event.getEntity();
				if ( entity == null)
					return;
				System.out.println(entity.getURLText());
				try {
					
					JNLPProxy.showDocument(new URL(entity.getURLText()));
					ChartRenderingInfo info = chartPanel.getChartRenderingInfo();
					chart.handleClick(x, y, info);
				} catch (Exception e){
					e.printStackTrace();
				}

			}
		});
		
		
		JFrame f = new JFrame("Results for " + repre.toName());
		
		f.setContentPane(chartPanel);
		f.pack();
		f.setVisible(true);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	
	}
	
	private  XYDataset getData4Repre(PdbChainKey repre) {
		
		
		DefaultXYDataset datas = new DefaultXYDataset();
		toolTipGenerator = new CustomXYToolTipGenerator();

		urlGenerator = new CustomXYURLGenerator();

		List<String> toolTips = new ArrayList<String>();
		List<String> urls = new ArrayList<String>();

//        for (int i = 0; i < data[0].length; i++) {
//            final double x = (double) i + 100000;
//            data[0][i] = x;
//            data[1][i] = 100000 + (double) Math.random() * COUNT;
//        }

	    SplitDatabase db = new SplitDatabase();
	    List<DBAlignment> aligs = db.getAllAlignments(repre);
	
	    double[][] data = new double[2][aligs.size()];
	    int x = -1;
	    for ( DBAlignment alig: aligs){
	    	x++;
	    	data[0][x] = alig.getRmsdOpt();
	    	data[1][x] = alig.getProbability();
	    	
	    	String toolTip = alig.getName1() + " : " + alig.getName2() + " | " + alig.toString();

			toolTips.add(toolTip);

			PdbChainKey n1 = PdbChainKey.fromName(alig.getName1());
			PdbChainKey n2 = PdbChainKey.fromName(alig.getName2());
			String url = String.format(BASE_URL,n1.getPdbId(),n1.getChainId(),n2.getPdbId(),n2.getChainId());
			urls.add(url);
	    }
	    
		datas.addSeries("repre", data);
		toolTipGenerator.addToolTipSeries(toolTips);
		urlGenerator.addURLSeries(urls);
		return datas;
	}

	private  JFreeChart createChart(XYDataset dataset) {

		//JFreeChart chart = ChartFactory.createBarChart("Symmetry in SCOP superfamilies",
			//	"SCOP Class", "% symmetry", dataset, PlotOrientation.VERTICAL,true,true,true);
		
		JFreeChart chart = ChartFactory.createScatterPlot("RMSD vs Probability", "RMSD", "Probability", dataset, 
				PlotOrientation.HORIZONTAL,true,true,true);

		XYPlot plot = (XYPlot) chart.getPlot();

		// ******************************************************************
		//  More than 150 demo applications are included with the JFreeChart
		//  Developer Guide...for more information, see:
		//
		//  >   http://www.object-refinery.com/jfreechart/guide.html
		//
		// ******************************************************************

		// set the range axis to display integers only...
		NumberAxis rangeAxis = (NumberAxis) plot.getDomainAxis();
		//rangeAxis.setStandardTickUnits(NumberAxis.);
		rangeAxis.setLowerBound(0);
		rangeAxis.setUpperBound(25);
		
		
		// disable bar outlines...
		XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) plot.getRenderer();
		
		chart.getRenderingHints().put
        (RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		if ( toolTipGenerator != null)
			renderer.setToolTipGenerator(toolTipGenerator);
		if ( urlGenerator != null)
			renderer.setURLGenerator(urlGenerator);
		
		// set up gradient paints for series...
//		GradientPaint gp0 = new GradientPaint(0.0f, 0.0f, Color.blue,
//				0.0f, 0.0f, new Color(0, 0, 64));
//		GradientPaint gp1 = new GradientPaint(0.0f, 0.0f, Color.green,
//				0.0f, 0.0f, new Color(0, 64, 0));
//		GradientPaint gp2 = new GradientPaint(0.0f, 0.0f, Color.red,
//				0.0f, 0.0f, new Color(64, 0, 0));
		//renderer.setSeriesPaint(0, gp0);
		//renderer.setSeriesPaint(1, gp1);
		//renderer.setSeriesPaint(2, gp2);

		
		
		
		
		
		// rotate legend
//		CategoryAxis domainAxis = plot.getDomainAxis();
//		domainAxis.setCategoryLabelPositions(
//				CategoryLabelPositions.createUpRotationLabelPositions(
//						Math.PI / 6.0));

		return chart;

	}
}
