package org.biojava.nbio.structure.align.symm.census3.analysis.mespeus;


/**
 * Something that returns true or false for a {@link MespeusEntry}.
 * @author dmyersturnbull
 */
public interface MespeusEntryMatcher {
	boolean matches(MespeusEntry entry);
}
