package org.biojava3.structure.align.symm.benchmark.comparison;

import org.biojava3.structure.align.symm.benchmark.KnownInfo;

/**
 * A factory for {@link GroupGuesser GroupGuessers}.
 * @author dmyerstu
 */
public class GroupGuesserFactory {

	public static GroupGuesser and(final GroupGuesser A, final GroupGuesser B) {
		return new GroupGuesser() {

			@Override
			public boolean hasEquivalentGroup(String a, String b) {
				return A.hasEquivalentGroup(a, b) && B.hasEquivalentGroup(a, b);
			}

			@Override
			public boolean hasEquivalentOrder(int a, int b) {
				return A.hasEquivalentOrder(a, b) && B.hasEquivalentOrder(a, b);
			}

		};
	}

	public static GroupGuesser exact() {

		return new GroupGuesser() {

			@Override
			public boolean hasEquivalentGroup(String A, String B) {
				int a = KnownInfo.getOrderFromGroup(A);
				int b = KnownInfo.getOrderFromGroup(B);
				return hasEquivalentOrder(a, b);
			}

			@Override
			public boolean hasEquivalentOrder(int a, int b) {
				return a == b;
			}

		};
	}

	public static GroupGuesser or(final GroupGuesser A, final GroupGuesser B) {
		return new GroupGuesser() {

			@Override
			public boolean hasEquivalentGroup(String a, String b) {
				return A.hasEquivalentGroup(a, b) || B.hasEquivalentGroup(a, b);
			}

			@Override
			public boolean hasEquivalentOrder(int a, int b) {
				return A.hasEquivalentOrder(a, b) || B.hasEquivalentOrder(a, b);
			}

		};
	}

	public static GroupGuesser withDivisorsOk() {

		return new GroupGuesser() {

			@Override
			public boolean hasEquivalentGroup(String A, String B) {
				int a = KnownInfo.getOrderFromGroup(A);
				int b = KnownInfo.getOrderFromGroup(B);
				return hasEquivalentOrder(a, b);
			}

			@Override
			public boolean hasEquivalentOrder(int a, int b) {
				if (a == 1 || b == 1) {
					return a == 1 && b == 1;
				}
				return b % a == 0;
			}

		};
	}

	public static GroupGuesser withMultiplesAndDivisorsOk() {
		return or(withMultiplesOk(), withDivisorsOk());
	}

	public static GroupGuesser withMultiplesOk() {

		return new GroupGuesser() {

			@Override
			public boolean hasEquivalentGroup(String A, String B) {
				int a = KnownInfo.getOrderFromGroup(A);
				int b = KnownInfo.getOrderFromGroup(B);
				return hasEquivalentOrder(a, b);
			}

			@Override
			public boolean hasEquivalentOrder(int a, int b) {
				if (a == 1 || b == 1) {
					return a == 1 && b == 1;
				}
				return a % b == 0;
			}

		};
	}

}
