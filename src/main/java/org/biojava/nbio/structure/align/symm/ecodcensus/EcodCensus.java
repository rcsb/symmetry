package org.biojava.nbio.structure.align.symm.ecodcensus;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomPositionMap;
import org.biojava.nbio.structure.ResidueRangeAndLength;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.symm.CESymmParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.OrderDetectorMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CeSymm;
import org.biojava.nbio.structure.align.symm.order.OrderDetectionFailedException;
import org.biojava.nbio.structure.align.symm.order.OrderDetector;
import org.biojava.nbio.structure.align.symm.order.SequenceFunctionOrderDetector;
import org.biojava.nbio.structure.align.symm.protodomain.Protodomain;
import org.biojava.nbio.structure.align.symm.protodomain.ProtodomainCreationException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.ecod.EcodDatabase;
import org.biojava.nbio.structure.ecod.EcodDomain;
import org.biojava.nbio.structure.ecod.EcodFactory;
import org.biojava.nbio.structure.io.LocalPDBDirectory.ObsoleteBehavior;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class EcodCensus {
	private static final Logger logger = LoggerFactory.getLogger(EcodCensus.class);

	public static void main(String[] args) {
		String outFilename;
		outFilename = "/tmp/EcodCensus.tsv";
//		outFilename = "-";
		PrintStream out = System.out;
		if( !outFilename.equalsIgnoreCase("-") ) {
			try {
				out = new PrintStream(new FileOutputStream(outFilename),true);
			} catch (FileNotFoundException e) {
				logger.error("Can't write to "+outFilename);
			}
		}
		
		String ecodVersion = "develop83";
		EcodDatabase ecod = EcodFactory.getEcodDatabase(ecodVersion);
		AtomCache cache = new AtomCache();
		cache.setObsoleteBehavior(ObsoleteBehavior.FETCH_OBSOLETE);
		
		CeSymm cesymm = new CeSymm();
		CESymmParameters param = (CESymmParameters) cesymm.getParameters();
		param.setRefineMethod(RefineMethod.REFINE);
		
		OrderDetector detector;
		//detector = new RotationOrderDetector(8, RotationOrderMethod.SINGLE_CUSP_FIXED_SSE);
		detector = new SequenceFunctionOrderDetector();
		param.setOrderDetectorMethod(OrderDetectorMethod.SEQUENCE_FUNCTION);
		
		EcodRepresentatives ecodreps = new EcodRepresentatives(ecod);
		List<EcodDomain> reps;
		try {
			reps = ecodreps.getDomains(3);
		} catch (IOException e) {
			logger.error("Error fetching ECOD domains",e);
			System.exit(1); return;
		}
		
		for(EcodDomain d: reps) {
			String rangeStr = String.format("%s.%s",d.getPdbId(),d.getRange());
			
			Atom[] ca1, ca2;
			try {
				Structure struct = cache.getStructure(rangeStr);
				ca1 = StructureTools.getRepresentativeAtomArray(struct);
				ca2 = StructureTools.getRepresentativeAtomArray(struct.clone());
			} catch (IOException e) {
				logger.error("Error getting structure for "+d.getDomainId(),e);
				continue;
			} catch (StructureException e) {
				logger.error("Error getting structure for "+d.getDomainId(),e);
				continue;
			}
			
			AFPChain afpChain;
			try {
				
				afpChain = cesymm.align(ca1, ca2);
//				afpChain.setName1(d.getDomainId());
//				afpChain.setName2(d.getDomainId());
				afpChain.setName1(rangeStr);
				afpChain.setName2(rangeStr);
			} catch (StructureException e) {
				logger.error("Error running CE-Symm on "+d.getDomainId(),e);
				continue;
			}
		
			int order;
			try {
				order = detector.calculateOrder(afpChain, ca1);
			} catch (OrderDetectionFailedException e) {
				logger.error("Error getting order for "+d.getDomainId(),e);
				order = 1;
			}
			
//			Protodomain protodomain;
//			try {
//				protodomain = Protodomain.fromSymmetryAlignment(afpChain, ca1, order, cache);
//			} catch (ProtodomainCreationException e) {
//				logger.error("Error getting protodomain for "+d.getDomainId(),e);
//				continue;
//			}

			String protodomainName = d.getDomainId();
			
//			Iterator<ResidueRangeAndLength> it = protodomain.getRanges().iterator();
			
			AtomPositionMap map = new AtomPositionMap(ca1);
			List<ResidueRangeAndLength> unsplicedRanges = ResidueRangeAndLength.parseMultiple(d.getRange(),map);
			List<ResidueRangeAndLength> splicedRanges = Protodomain.spliceApproxConsecutive(map, unsplicedRanges, 4);
			
			Iterator<ResidueRangeAndLength> it = splicedRanges.iterator();
			
			StringBuilder protodomainRangeStr = new StringBuilder();
			if(it.hasNext()) {
				protodomainRangeStr.append(d.getPdbId());
				protodomainRangeStr.append("_");
				ResidueRangeAndLength range = it.next();
				protodomainRangeStr.append(range.toString());
			}
			while(it.hasNext()) {
				ResidueRangeAndLength range = it.next();
				protodomainRangeStr.append(range.toString());
			}
			
			out.format("%s\t%s%n",protodomainName,protodomainRangeStr);
		}
		if( out != System.out) {
			out.close();
		}
	}
	
}
