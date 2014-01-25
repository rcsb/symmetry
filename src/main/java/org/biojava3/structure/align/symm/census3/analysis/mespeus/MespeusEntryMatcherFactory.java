package org.biojava3.structure.align.symm.census3.analysis.mespeus;

/**
 * A factory for {@link MespeusEntryMatcher MespeusEntryMatchers}.
 * @author dmyersturnbull
 */
public class MespeusEntryMatcherFactory {

	public static MespeusEntryMatcher everything() {
		return new MespeusEntryMatcher() {
			@Override
			public boolean matches(MespeusEntry entry) {
				return true;
			}
		};
	}
	public static MespeusEntryMatcher coordinationNumberEquals(final int coordinationNumber) {
		return new MespeusEntryMatcher() {
			@Override
			public boolean matches(MespeusEntry entry) {
				return entry.getCoordinationNumber() == coordinationNumber;
			}
		};
	}
	public static MespeusEntryMatcher coordinationEquals(final int coordinationNumber) {
		return new MespeusEntryMatcher() {
			@Override
			public boolean matches(MespeusEntry entry) {
				return entry.getCoordinationNumber() == coordinationNumber;
			}
		};
	}
	public static MespeusEntryMatcher coordinationNumberAtLeast(final int coordinationNumber) {
		return new MespeusEntryMatcher() {
			@Override
			public boolean matches(MespeusEntry entry) {
				return entry.getCoordinationNumber() >= coordinationNumber;
			}
		};
	}
	public static MespeusEntryMatcher coordinationNumberAtMost(final int coordinationNumber) {
		return new MespeusEntryMatcher() {
			@Override
			public boolean matches(MespeusEntry entry) {
				return entry.getCoordinationNumber() <= coordinationNumber;
			}
		};
	}
	public static MespeusEntryMatcher metalName(final String[] metalNames) {
		return new MespeusEntryMatcher() {
			@Override
			public boolean matches(MespeusEntry entry) {
				String s = entry.getMetalName().split(" ")[0];
				for (String name : metalNames) {
					if (name.equalsIgnoreCase(s)) return true;
				}
				return false;
			}
		};
	}
}
