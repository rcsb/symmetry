package org.biojava3.structure.align.symm.census;

import java.io.Serializable;

import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AfpChainWriter;



public class ProtoDomainDBSearchResult implements Serializable{

	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private static String newline  = AfpChainWriter.newline;
	
	String name1;
	String name2;
	Float alignScore;
	String algorithmName;
	Double probability;
	Float rmsd;
	Integer ca1Length;
	Integer ca2Length;
	Integer coverage1;
	Integer coverage2;
	Float identity;
	Float similarity;
	Integer eqr;
	Boolean isCP;
	String classificationID2;
	String originalScopID2;
	
	public String getName1() {
		return name1;
	}
	public void setName1(String name1) {
		this.name1 = name1;
	}
	public String getName2() {
		return name2;
	}
	public void setName2(String name2) {
		this.name2 = name2;
	}
	public Float getAlignScore() {
		return alignScore;
	}
	public void setAlignScore(Float alignScore) {
		this.alignScore = alignScore;
	}
	public String getAlgorithmName() {
		return algorithmName;
	}
	public void setAlgorithmName(String algorithmName) {
		this.algorithmName = algorithmName;
	}
	public Double getProbability() {
		return probability;
	}
	public void setProbability(Double probability) {
		this.probability = probability;
	}
	public Float getRmsd() {
		return rmsd;
	}
	public void setRmsd(Float rmsd) {
		this.rmsd = rmsd;
	}
	public Integer getCa1Length() {
		return ca1Length;
	}
	public void setCa1Length(Integer ca1Length) {
		this.ca1Length = ca1Length;
	}
	public Integer getCa2Length() {
		return ca2Length;
	}
	public void setCa2Length(Integer ca2Length) {
		this.ca2Length = ca2Length;
	}
	public Integer getCoverage1() {
		return coverage1;
	}
	public void setCoverage1(Integer coverage1) {
		this.coverage1 = coverage1;
	}
	public Integer getCoverage2() {
		return coverage2;
	}
	public void setCoverage2(Integer coverage2) {
		this.coverage2 = coverage2;
	}
	public Float getIdentity() {
		return identity;
	}
	public void setIdentity(Float identity) {
		this.identity = identity;
	}
	public Float getSimilarity() {
		return similarity;
	}
	public void setSimilarity(Float similarity) {
		this.similarity = similarity;
	}
	public Integer getEqr() {
		return eqr;
	}
	public void setEqr(Integer eqr) {
		this.eqr = eqr;
	}


	public Boolean getIsCP() {
		return isCP;
	}
	public void setIsCP(Boolean isCP) {
		this.isCP = isCP;
	}
	public String toString(){

		StringBuffer str = new StringBuffer();

		str.append(getName1());
		str.append("\t");
		str.append(getName2());
		str.append("\t");
		str.append(String.format("%.2f",getAlignScore()));
		str.append("\t");     
		if ( getAlgorithmName().equalsIgnoreCase(CeMain.algorithmName)){
			str.append(String.format("%.2f",getProbability()));
		} else {
			str.append(String.format("%.2e",getProbability()));
		}
		str.append("\t");
		str.append(String.format("%.2f",getRmsd()));
		str.append("\t");
		str.append(getEqr());
		str.append("\t");
		str.append(getCa1Length());
		str.append("\t");
		str.append(getCa2Length());      
		str.append("\t");
		str.append(getCoverage1());
		str.append("\t");
		str.append(getCoverage2());
		str.append("\t");
		str.append(String.format("%.2f",getIdentity()));
		str.append("\t");
		str.append(String.format("%.2f",getSimilarity()));
		str.append("\t");
		str.append(isCP);
		str.append("\t");
		str.append(newline);

		return str.toString();
	}
	
	
	
	public String getClassificationID2() {
		return classificationID2;
	}
	public void setClassificationID2(String classificationID2) {
		this.classificationID2 = classificationID2;
	}
	public String getOriginalScopID2() {
		return originalScopID2;
	}
	public void setOriginalScopID2(String originalScopID2) {
		this.originalScopID2 = originalScopID2;
	}
	public static String getHTMLHeader(){
		StringBuffer str = new StringBuffer();
		str.append("<tr><th>");
		str.append("Name1");
		str.append("</th><th>");
		str.append("Name2");
		str.append("</th><th>");
		str.append("AlignScore");
		str.append("</th><th>");     
		
		str.append("Z-Score");
		
		str.append("</th><th>");
		str.append("RMSD");
		str.append("</th><th>");
		str.append("EQR");
		str.append("</th><th>");
		str.append("Length1");
		str.append("</th><th>");
		str.append("Length2");      
		str.append("</th><th>");
		str.append("Coverage1");
		str.append("</th><th>");
		str.append("Coverage2");
		str.append("</th><th>");
		str.append("Identity");
		str.append("</th><th>");
		str.append("Similarity");
		str.append("</th><th>");
		str.append("CP?");
		str.append("</th><th>");
		str.append("scopID2");
		str.append("</th><th>");
		str.append("classificationID2");
		
		str.append("</th>");
		str.append("</tr>");
		str.append(newline);

		return str.toString();
	}
	public String toHTML(){
		StringBuffer str = new StringBuffer();
		str.append("<tr><td>");
		str.append(getName1());
		str.append("</td><td>");
		str.append(getName2());
		str.append("</td><td>");
		str.append(String.format("%.2f",getAlignScore()));
		str.append("</td><td>");     
		
		str.append(String.format("%.2f",getProbability()));
		
		str.append("</td><td>");
		str.append(String.format("%.2f",getRmsd()));
		str.append("</td><td>");
		str.append(getEqr());
		str.append("</td><td>");
		str.append(getCa1Length());
		str.append("</td><td>");
		str.append(getCa2Length());      
		str.append("</td><td>");
		str.append(getCoverage1());
		str.append("</td><td>");
		str.append(getCoverage2());
		str.append("</td><td>");
		str.append(String.format("%.2f",getIdentity()));
		str.append("</td><td>");
		str.append(String.format("%.2f",getSimilarity()));
		str.append("</td><td>");
		str.append(isCP);
		str.append("</td><td>");
		str.append(originalScopID2);
		str.append("</td><td>");
		str.append(classificationID2);
		str.append("</td>");
		str.append("</tr>");
		str.append(newline);

		return str.toString();
	}

}
