package org.biojava.nbio.structure.align.symm;

import java.util.List;

import org.biojava.nbio.structure.align.ce.CeCPMain;
import org.biojava.nbio.structure.align.ce.CeParameters;

/**
 * Provides parameters to {@link CeCPMain}
 * 
 * @author Spencer Bliven
 *
 */
public class CESymmParameters extends CeParameters {

	private int maxSymmOrder; //Renamed, old variable maxNrAlternatives (now means max nr. of iterations/order of symmetry)
	private OrderDetectorMethod orderDetectorMethod;
	private RefineMethod refineMethod;
	
	public static enum OrderDetectorMethod {
		SEQUENCE_FUNCTION;
		public static OrderDetectorMethod DEFAULT = SEQUENCE_FUNCTION;
	}
	
	public static enum RefineMethod {
		NOT_REFINED,
		SINGLE,
		MULTIPLE,
		SINGLE_OPTIMIZE,
		MULTIPLE_OPTIMIZE;
		public static RefineMethod DEFAULT = SINGLE;
	}
	
	public CESymmParameters() {
		super();
		maxSymmOrder = 8;
		refineMethod = RefineMethod.DEFAULT;
		orderDetectorMethod = OrderDetectorMethod.DEFAULT;
	}

	@Override
	public String toString() {
		return "CESymmParameters [maxSymmOrder=" + maxSymmOrder
				+ ", orderDetectorMethod=" + orderDetectorMethod
				+ ", refineMethod=" + refineMethod + ", winSize=" + winSize
				+ ", rmsdThr=" + rmsdThr + ", rmsdThrJoin=" + rmsdThrJoin
				+ ", maxOptRMSD=" + maxOptRMSD + ", scoringStrategy="
				+ scoringStrategy + ", maxGapSize=" + maxGapSize
				+ ", showAFPRanges=" + showAFPRanges
				+ ", sideChainScoringType=" + sideChainScoringType
				+ ", gapOpen=" + gapOpen + ", gapExtension=" + gapExtension
				+ ", distanceIncrement=" + distanceIncrement + ", oRmsdThr="
				+ oRmsdThr + ", maxNrIterationsForOptimization="
				+ maxNrIterationsForOptimization + ", substitutionMatrix="
				+ substitutionMatrix + ", seqWeight=" + seqWeight + "]";
	}


	@Override
	public void reset(){
		super.reset();
		maxSymmOrder = 8;
		orderDetectorMethod = OrderDetectorMethod.DEFAULT;
		refineMethod = RefineMethod.DEFAULT;
	}


	@Override
	public List<String> getUserConfigHelp() {
		List<String> params = super.getUserConfigHelp();
		
		//maxSymmOrder help explanation
		params.add("Sets the maximum order of symmetry of the protein.");
		
		StringBuilder orderTypes = new StringBuilder("Order Detection Method: ");
		OrderDetectorMethod[] vals = OrderDetectorMethod.values();
		if(vals.length == 1) {
			orderTypes.append(vals[0].name());
		} else if(vals.length > 1 ) {
			for(int i=0;i<vals.length-1;i++) {
				orderTypes.append(vals[i].name());
				orderTypes.append(", ");
			}
			orderTypes.append("or ");
			orderTypes.append(vals[vals.length-1].name());
		}
		params.add(orderTypes.toString());
		
		StringBuilder refineTypes = new StringBuilder("Refinement Method: ");
		RefineMethod[] values = RefineMethod.values();
		if(values.length == 1) {
			refineTypes.append(values[0].name());
		} else if(values.length > 1 ) {
			for(int i=0;i<values.length-1;i++) {
				refineTypes.append(values[i].name());
				refineTypes.append(", ");
			}
			refineTypes.append("or ");
			refineTypes.append(values[values.length-1].name());
		}
		params.add(refineTypes.toString());
		
		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = super.getUserConfigParameters();
		params.add("MaxSymmOrder");
		params.add("OrderDetectorMethod");
		params.add("RefineMethod");
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames(){
		List<String> params = super.getUserConfigParameterNames();
		params.add("Maximum Order of Symmetry");
		params.add("Order Detection Method");
		params.add("Refinement Method");
		return params;
	}

	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = super.getUserConfigTypes();
		params.add(Integer.class);
		params.add(OrderDetectorMethod.class);
		params.add(RefineMethod.class);
		return params;
	}

	public RefineMethod getRefineMethod() {
		return refineMethod;
	}

	public void setRefineMethod(RefineMethod refineMethod) {
		this.refineMethod = refineMethod;
	}
	
	@Deprecated
	public void setRefineResult(boolean doRefine) {
		if (!doRefine){
			refineMethod = RefineMethod.NOT_REFINED;
		}
		else{
			refineMethod = RefineMethod.DEFAULT;
		}
	}

	public OrderDetectorMethod getOrderDetectorMethod() {
		return orderDetectorMethod;
	}

	public void setOrderDetectorMethod(OrderDetectorMethod orderDetectorMethod) {
		this.orderDetectorMethod = orderDetectorMethod;
	}
	
	public void setMaxSymmOrder(Integer max) {
		maxSymmOrder = max;
	}

	public int getMaxSymmOrder() {
		return maxSymmOrder;
	}
}
