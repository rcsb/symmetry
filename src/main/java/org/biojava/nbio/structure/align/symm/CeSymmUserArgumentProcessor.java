package org.biojava.nbio.structure.align.symm;

import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.CeUserArgumentProcessor;
import org.biojava.nbio.structure.align.ce.StartupParameters;
import org.biojava.nbio.structure.align.symm.CESymmParameters.OrderDetectorMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.align.symm.CESymmParameters.SubunitColors;

public class CeSymmUserArgumentProcessor extends CeUserArgumentProcessor{
	
	protected static class CeSymmStartupParams extends CeUserArgumentProcessor.CeStartupParams {
		
		protected int maxSymmOrder;
		protected OrderDetectorMethod orderDetectorMethod;
		protected RefineMethod refineMethod;
		private SubunitColors subunitColors;
		

		public CeSymmStartupParams() {
			super();
			maxSymmOrder = 8;
			orderDetectorMethod = OrderDetectorMethod.DEFAULT;
			refineMethod = RefineMethod.DEFAULT;
			subunitColors = SubunitColors.DEFAULT;
		}

		public OrderDetectorMethod getOrderDetectorMethod() {
			return orderDetectorMethod;
		}
		
		public void setOrderDetectorMethod(OrderDetectorMethod orderDetectorMethod) {
			this.orderDetectorMethod = orderDetectorMethod;
		}
		
		public RefineMethod getRefineMethod() {
			return refineMethod;
		}
		
		public void setRefineMethod(RefineMethod refineMethod) {
			this.refineMethod = refineMethod;
		}
		
		public SubunitColors getSubunitColors() {
			return subunitColors;
		}

		public void setSubunitColors(SubunitColors colors) {
			this.subunitColors = colors;
		}

		public int getMaxSymmOrder() {
			return maxSymmOrder;
		}

		public void setMaxSymmOrder(Integer max) {
			this.maxSymmOrder = max;
		}

		@Override
		public String toString() {
			return "CeSymmStartupParams [orderDetectorMethod="
					+ orderDetectorMethod + ", refineMethod=" + refineMethod
					+ ", maxSymmOrder=" + maxSymmOrder + ", subunitColors=" + subunitColors
					+ ", getWinSize()=" + getWinSize()
					+ ", getScoringStrategy()=" + getScoringStrategy()
					+ ", getGapOpen()=" + getGapOpen() + ", getGapExtension()="
					+ getGapExtension() + ", getMaxGapSize()="
					+ getMaxGapSize() + ", isShowAFPRanges()="
					+ isShowAFPRanges() + ", getMaxOptRMSD()="
					+ getMaxOptRMSD() + ", toString()=" + super.toString()
					+ ", getSearchFile()=" + getSearchFile()
					+ ", getAlignPairs()=" + getAlignPairs()
					+ ", getSaveOutputDir()=" + getSaveOutputDir()
					+ ", isShowMenu()=" + isShowMenu() + ", isPrintCE()="
					+ isPrintCE() + ", getPdb1()=" + getPdb1() + ", getPdb2()="
					+ getPdb2()
					+ ", isPrintXML()=" + isPrintXML() + ", isPrintFatCat()="
					+ isPrintFatCat() + ", getPdbFilePath()="
					+ getPdbFilePath() + ", getCacheFilePath()="
					+ getCacheFilePath() + ", isShow3d()=" + isShow3d()
					+ ", getOutFile()=" + getOutFile() + ", isAutoFetch()="
					+ isAutoFetch() + ", getShowDBresult()="
					+ getShowDBresult() + ", getNrCPU()=" + getNrCPU()
					+ ", getFile1()=" + getFile1() + ", getFile2()="
					+ getFile2() + ", isOutputPDB()=" + isOutputPDB()
					+ ", isDomainSplit()=" + isDomainSplit() + ", getClass()="
					+ getClass() + ", hashCode()=" + hashCode() + "]";
		}

	}
	
	@Override
	protected StartupParameters getStartupParametersInstance() {
		return new CeSymmStartupParams();
	}
	
	@Override
	public StructureAlignment getAlgorithm() {
		return new CeSymm();
	}
	

	@Override
	public Object getParameters() {
		
		StructureAlignment alignment = getAlgorithm();
		
		CESymmParameters aligParams = (CESymmParameters) alignment.getParameters();
		CeSymmStartupParams startParams = (CeSymmStartupParams) params;
		
		if ( aligParams == null)
			aligParams = new CESymmParameters();
		
		// Copy relevant parameters from the startup parameters
		aligParams.setMaxGapSize(startParams.getMaxGapSize());
		aligParams.setWinSize(startParams.getWinSize());
		aligParams.setScoringStrategy(startParams.getScoringStrategy());
		aligParams.setMaxOptRMSD(startParams.getMaxOptRMSD());
		aligParams.setGapOpen(startParams.getGapOpen());
		aligParams.setGapExtension(startParams.getGapExtension());
		aligParams.setShowAFPRanges(startParams.isShowAFPRanges());
		aligParams.setMaxSymmOrder(startParams.getMaxSymmOrder());
		aligParams.setOrderDetectorMethod(startParams.getOrderDetectorMethod());
		aligParams.setRefineMethod(startParams.getRefineMethod());
		aligParams.setSubunitColors(startParams.getSubunitColors());
		
		return aligParams;
	}
	

	@Override
	public String getDbSearchLegend(){
		//String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
		//return legend;
		
		return "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\tDescription\t " ;
		
	}
}
