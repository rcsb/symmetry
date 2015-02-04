package org.biojava.nbio.structure.align.symm.census3;

import java.io.Serializable;

import javax.xml.bind.annotation.XmlTransient;


/**
 * @author dmyersturnbull
 * TODO
 */
@XmlTransient
public abstract class AdditionalScoreList implements Serializable {

	private static final long serialVersionUID = 155894445075321750L;
	
	public abstract Number getScore(String scoreName);
	public abstract String[] getScoreNames();
	
}
