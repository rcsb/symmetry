package results;

import java.util.Locale;


import org.jfree.chart.ChartFactory; 
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.general.DefaultPieDataset;
import org.jfree.data.statistics.SimpleHistogramBin;
import org.jfree.data.statistics.SimpleHistogramDataset;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GradientPaint;
import java.awt.Paint;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarPainter;
import org.jfree.chart.renderer.category.BarRenderer;

import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;

/**
 * A simple demonstration application showing how to create a bar chart.
 */
public class SymmetryInSCOPSuperfamilies    {

	/**
	 * Creates a new demo instance.
	 *
	 * @param title  the frame title.
	 */
	public SymmetryInSCOPSuperfamilies() {
		
		CategoryDataset dataset = createDataset();
		JFreeChart chart = createChart(dataset);
		ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setFillZoomRectangle(true);
		
		chartPanel.setPreferredSize(new Dimension(500, 270));
		
		JFrame f = new JFrame("SCOP Classes");
		
		f.setContentPane(chartPanel);
		f.pack();
		f.setVisible(true);
	}

	/**
	 * Returns a sample dataset.
	 *
	 * @return The dataset.
	 */
	private static CategoryDataset createDataset() {

		// row keys...
		String series1 = "score >3.5, 0.35";
		String series2 = "score >3.5, 0.35, RMSD < 5";
		

		// column keys...
		String category0 = "Overall";
		String category1 = "A - Alpha";
		String category2 = "B - Beta";
		String category3 = "C - A / B";
		String category4 = "D - A+B";
		String category5 = "E - Multi domain";
		String category6 = "F - Membrane";
		String category7 = "G - Small";
		
		// create the dataset...
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();

		dataset.addValue(0.42, series1, category0);
		dataset.addValue(0.54, series1, category1);
		dataset.addValue(0.50, series1, category2);
		dataset.addValue(0.45, series1, category3);
		dataset.addValue(0.41, series1, category4);
		dataset.addValue(0.29, series1, category5);
		dataset.addValue(0.6,  series1, category6);
		dataset.addValue(0.08, series1, category7);
		
		dataset.addValue(0.33, series2, category0);
		dataset.addValue(0.39, series2, category1);
		dataset.addValue(0.44, series2, category2);
		dataset.addValue(0.40, series2, category3);
		dataset.addValue(0.26, series2, category4);
		dataset.addValue(0.14, series2, category5);
		dataset.addValue(0.43, series2, category6);
		dataset.addValue(0.06, series2, category7);
		
		return dataset;

	}

	/**
	 * Creates a sample chart.
	 *
	 * @param dataset  the dataset.
	 *
	 * @return The chart.
	 */
	private static JFreeChart createChart(CategoryDataset dataset) {

		JFreeChart chart = ChartFactory.createBarChart("Symmetry in SCOP superfamilies",
				"SCOP Class", "% symmetry", dataset, PlotOrientation.VERTICAL,true,true,true);                   

		CategoryPlot plot = (CategoryPlot) chart.getPlot();

		// ******************************************************************
		//  More than 150 demo applications are included with the JFreeChart
		//  Developer Guide...for more information, see:
		//
		//  >   http://www.object-refinery.com/jfreechart/guide.html
		//
		// ******************************************************************

		// set the range axis to display integers only...
		NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
		rangeAxis.setLowerBound(0);
		rangeAxis.setUpperBound(1);
		//rangeAxis.setStandardTickUnits(NumberAxis.);

		// disable bar outlines...
		BarRenderer renderer = (BarRenderer) plot.getRenderer();
		renderer.setDrawBarOutline(true);

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

		
		renderer.setShadowVisible(false);
		
		
		
		// rotate legend
//		CategoryAxis domainAxis = plot.getDomainAxis();
//		domainAxis.setCategoryLabelPositions(
//				CategoryLabelPositions.createUpRotationLabelPositions(
//						Math.PI / 6.0));

		return chart;

	}

	/**
	 * Starting point for the demonstration application.
	 *
	 * @param args  ignored.
	 */
	public static void main(String[] args) {
		SymmetryInSCOPSuperfamilies demo = new SymmetryInSCOPSuperfamilies();
				
		
	}



}
