package org.biojava3.structure.align.symm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.statistics.HistogramDataset;
import org.jfree.data.statistics.HistogramType;

@Deprecated
public class SymmRefiner {


	/**
	 * Takes a self alignment and applies it <tt>n</tt> times. Returns a histogram
	 * of the distance between each residue and its n-image.
	 * 
	 * @param afpChain A self alignment, ie one where both proteins are the same
	 * @param n The number of times to apply the alignment
	 * @return
	 * @throws StructureException
	 */
	public static Map<Integer,Integer> applyAlignment(AFPChain afpChain, int n) throws StructureException {
		if(n < 1) {
			throw new IllegalArgumentException("n must be at least 1");
		}
		
		// quick check that afpChain represents a self-alignment
		assert(afpChain.getCa1Length() == afpChain.getCa2Length());

		Map<Integer,Integer> deltas = new HashMap<Integer,Integer>();
		
		// convert alignment to a map between residue indices
		Map<Integer,Integer> alignment = CeSymm.alignmentAsMap(afpChain);
		alignment.put(null, null);
		
		//iterate over the alignment
		for(Integer preimage : alignment.keySet()) {
			Integer image = preimage;
			
			//Apply the mapping n times
			for(int i=0;i<n;i++) {
				image = alignment.get(image);
			}
			
			// calculate the number of residues between image and preimage
			Integer delta = image==null? null : image - preimage;
			
			if(preimage != null) {
				deltas.put(preimage, delta);
			}
		}
		
		return deltas;
	}

	
	public static void main(String[] args) {
		try {
			String name;
			
			name = "1itb.A"; // b-trefoil, C3
			//name = "1tim.A"; // tim-barrel, C8
			name = "d1p9ha_"; // not rotational symmetry
			name = "3HKE.A"; // very questionable alignment
			name = "d1jlya1";
			
			AtomCache cache = new AtomCache();
			Atom[] ca1 = cache.getAtoms(name);
			Atom[] ca2 = cache.getAtoms(name);
			
			StructureAlignmentFactory.addAlgorithm(new CeSymm());
			CeSymm ce = (CeSymm) StructureAlignmentFactory.getAlgorithm(CeSymm.algorithmName);
			
			AFPChain afpChain = ce.align(ca1, ca2);
			/*
			displayHist(afpChain,1,8);
			for(int n=1;n<=8;n++) {
				displayHist(afpChain,n,n);
			}*/
			int symm = CeSymm.getSymmetryOrder(afpChain);
			System.out.println("Symmetry="+symm);
			
			StructureAlignmentDisplay.display(afpChain, ca1, ca2);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/*
	private static void printHist(Map<Integer, Integer> hist) {
		List<Integer> preimages = new ArrayList<Integer>(hist.keySet());
		
		Collections.sort(preimages, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				//sort nulls to the end
				if(o1 == null) {
					if(o2 == null) { return 0; }
					else { return -1; }
				} else {
					if(o2 == null) { return 1; }
					else { return o1.compareTo(o2); }
				}
			}
		} );
		
		for(Integer pre : preimages) {
			if(pre != null) {
				System.out.format("%d\t%d\n", pre,hist.get(pre));
			}
		}
	}
	*/
	
	private static void displayHist(AFPChain afpChain,int minSymmetry, int maxSymmetry) throws StructureException {
		List<List<Integer>> values = new ArrayList<List<Integer>>(maxSymmetry);
		int minValue = Integer.MAX_VALUE;
		int maxValue = Integer.MIN_VALUE;
		
		for(int n=minSymmetry;n<=maxSymmetry;n++) {
			// Apply afpChain n times to each residue
			Map<Integer,Integer> deltas = applyAlignment(afpChain, n);
			
			List<Integer> nValues = new ArrayList<Integer>(deltas.size());
			for(Integer v : deltas.values()) {
				if(v != null) {
					nValues.add(v);

					if(v<minValue) {
						minValue = v;
					}
					if(v>maxValue) {
						maxValue = v;
					}
				}
			}
			values.add(nValues);
		}
		
		
		HistogramDataset dataset = new HistogramDataset();
		//SimpleHistogramDataset dataset = new SimpleHistogramDataset("histogram");
		dataset.setType(HistogramType.FREQUENCY);
		
		int bins = maxValue-minValue+1;
		System.out.println(bins+" bins");
		for(int n=minSymmetry;n<=maxSymmetry;n++) {
			List<Integer> nValues = values.get(n-minSymmetry);
			double[] valuesArr = new double[nValues.size()];
			for(int i=0;i<valuesArr.length;i++)
				valuesArr[i] = nValues.get(i);
			dataset.addSeries(""+n, valuesArr, bins,(double)minValue,(double)maxValue+1);
		}
			
		for(int item = 0;item < dataset.getItemCount(0);item++) {
			System.out.format("%f\t%f\t%f\n",dataset.getStartXValue(0, item),dataset.getEndXValue(0, item),dataset.getStartYValue(0, item));
		}
		
		
		
		String plotTitle = "Residue differences after rotating n times"; 
		String xaxis = "Number of residues different";
		String yaxis = "Frequency"; 
		PlotOrientation orientation = PlotOrientation.VERTICAL; 
		boolean legend = true; 
		boolean toolTips = false;
		boolean urls = false; 
		JFreeChart chart = ChartFactory.createHistogram( plotTitle, xaxis, yaxis, 
				dataset, orientation, legend, toolTips, urls);
		ChartFrame frame = new ChartFrame(plotTitle, chart);
		frame.pack();
		frame.setVisible(true);
		/*
		String filename = "/tmp/chart.png";
		int width = 500;
		int height = 300; 
		try {
			ChartUtilities.saveChartAsPNG(new File(filename), chart, width, height);
		} catch (IOException e) {}
		*/
	}

}
