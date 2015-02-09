package demo;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.client.PdbPair;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.benchmark.Benchmark;
import org.biojava.nbio.benchmark.FisherBenchmark;

public class RunFisherBenchmark {

	public void main(String[] args){

		Map<String, Integer> best = new HashMap<String, Integer>();
		AtomCache cache = new AtomCache();

		cache.setFetchCurrent(true);
		cache.setFetchFileEvenIfObsolete(true);		
		cache.setPath("/Volumes/Macintosh HD2/PDB/");

		Benchmark benchmark = new FisherBenchmark();

		List<PdbPair> pairs = benchmark.getPairs();
		int counter =0;
		System.out.print("#");
		for (StructureAlignment alig : StructureAlignmentFactory.getAllAlgorithms()){
			if ( alig.getAlgorithmName().contains("Optimal"))
				continue;
			System.out.print(alig.getAlgorithmName() +"\t");
		}
		System.out.println("");
		for ( PdbPair pair : pairs){
			int currentBest = 0;
			String bestA = "";

			try {
				counter++;
				Atom[] ca1 = cache.getAtoms(pair.getName1());
				Atom[] ca2 = cache.getAtoms(pair.getName2());

				System.out.print("#" + counter + " " + pair.getName1() + " " + pair.getName2() + " ");
				for (StructureAlignment alig : StructureAlignmentFactory.getAllAlgorithms()){
					if ( alig.getAlgorithmName().contains("Optimal"))
						continue;
					AFPChain afpChain = alig.align(ca1, ca2);
					System.out.print(afpChain.getNrEQR() + "(" + String.format("%.1f",afpChain.getTotalRmsdOpt()) + ")"+ "\t" );
					if ( afpChain.getNrEQR() >= currentBest) {
						currentBest = afpChain.getNrEQR();
						bestA = alig.getAlgorithmName();
					}
				}
				if ( !best.containsKey(bestA)){
					best.put(bestA, 0);
				}
				int val = best.get(bestA);
				val++;
				best.put(bestA,val);

				System.out.println("");

				//System.out.println("#" + counter+" " + AfpChainWriter.toDBSearchResult(afpChain));


			} catch (Exception e){
				e.printStackTrace();			
			}

			System.out.println(best);
		}
	}
}
