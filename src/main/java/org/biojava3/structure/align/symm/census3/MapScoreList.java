package org.biojava3.structure.align.symm.census3;

import java.util.HashMap;
import java.util.Map;

/**
 * Simple map-based implementation of an AdditionalScoreList
 */
public class MapScoreList extends AdditionalScoreList { 
	private static final long serialVersionUID = -3342892699452329060L;
	
	private Map<String,Number> scores;
	public MapScoreList() {
		this.scores = new HashMap<String, Number>();
	}
	public MapScoreList(Map<String,Number> scores) {
		this.scores = scores;
	}
	public Map<String, Number> getScores() {
		return scores;
	}
	public void setScores(Map<String, Number> scores) {
		this.scores = scores;
	}
	@Override
	public String[] getScoreNames() {
		String[] names = new String[scores.size()];
		int i= 0;
		for(String name: scores.keySet()) {
			names[i++] = name;
		}
		return names;
	}
	@Override
	public Number getScore(String scoreName) {
		return scores.get(scoreName);
	}
}