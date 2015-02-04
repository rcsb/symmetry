package org.biojava.nbio.structure.align.symm;

import java.util.List;

import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeParameters;

/**
 * Provides parameters to {@link CeCPMain}
 * 
 * @author Spencer Bliven
 *
 */
public class CESymmParameters extends CeParameters {

	private boolean refineResult;
	private OrderDetectorMethod orderDetectorMethod;
	private int maxNrAlternatives; // Not exposed in UI

	
	public static enum OrderDetectorMethod {
		SEQUENCE_FUNCTION;
		public static OrderDetectorMethod DEFAULT = SEQUENCE_FUNCTION;
	}
	
	public CESymmParameters() {
		super();
		refineResult = false;
		orderDetectorMethod = OrderDetectorMethod.DEFAULT;
		maxNrAlternatives = 1;
	}

	@Override
	public String toString() {
		return "CECPParameters [scoringStrategy=" + scoringStrategy 
		+ ", maxGapSize=" + maxGapSize 
		+ ", rmsdThr=" + rmsdThr 
		+ ", rmsdThrJoin="+ rmsdThrJoin 
		+ ", winSize=" + winSize 
		+ ", showAFPRanges=" + showAFPRanges 
		+ ", maxOptRMSD=" + maxOptRMSD
		+ ", seqWeight=" + seqWeight
		+ "]";
	}


	@Override
	public void reset(){
		super.reset();
		refineResult = false;
		orderDetectorMethod = OrderDetectorMethod.DEFAULT;
		maxNrAlternatives = 1;
	}


	@Override
	public List<String> getUserConfigHelp() {
		List<String> params = super.getUserConfigHelp();
		StringBuilder orderTypes = new StringBuilder("Order detection method: ");
		OrderDetectorMethod[] vals = OrderDetectorMethod.values();
		if(vals.length == 1) {
			orderTypes.append(vals[0].name());
		} else if(vals.length > 1 ) {
			for(int i=0;i<vals.length-1;i++) {
				orderTypes.append(vals[i].name());
				orderTypes.append(", ");
			}
			orderTypes.append("or ");
			orderTypes.append(vals[vals.length].name());
		}
		params.add(orderTypes.toString());
		params.add("Refine the result to a multiple alignment");
		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = super.getUserConfigParameters();
		params.add("OrderDetectorMethod");
		params.add("RefineResult");
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames(){
		List<String> params = super.getUserConfigParameterNames();
		
		params.add("Order detection method");
		params.add("Refine Result");
		return params;
	}

	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = super.getUserConfigTypes();
		params.add(OrderDetectorMethod.class);
		params.add(Boolean.class);
		return params;
	}

	public boolean isRefineResult() {
		return refineResult;
	}

	public void setRefineResult(Boolean refineResult) {
		this.refineResult = refineResult;
	}

	public OrderDetectorMethod getOrderDetectorMethod() {
		return orderDetectorMethod;
	}

	public void setOrderDetectorMethod(OrderDetectorMethod orderDetectorMethod) {
		this.orderDetectorMethod = orderDetectorMethod;
	}
	
	public void setMaxNrAlternatives(int max) {
		maxNrAlternatives = max;
	}

	public int getMaxNrAlternatives() {
		return maxNrAlternatives;
	}
}
