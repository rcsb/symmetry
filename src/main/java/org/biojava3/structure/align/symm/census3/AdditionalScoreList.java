package org.biojava3.structure.align.symm.census3;

import java.io.Serializable;

/**
 * @author dmyersturnbull
 * TODO
 */
public abstract class AdditionalScoreList implements Serializable {

	private static final long serialVersionUID = 155894445075321750L;
	
	public abstract Number getScore(String scoreName);
	public abstract String[] getScoreNames();
	
}
