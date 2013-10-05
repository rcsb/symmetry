package org.biojava3.structure.align.symm.benchmark.comparison.order;

import org.biojava3.structure.align.symm.benchmark.KnownInfo;

/**
 * A factory for {@link GroupComparator GroupGuessers}.
 * @author dmyerstu
 */
public class GroupComparisonFactory {

	public static GroupComparator and(final GroupComparator A, final GroupComparator B) {
		return new GroupComparator() {

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

	public static GroupComparator exact() {

		return new GroupComparator() {

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

	public static GroupComparator or(final GroupComparator A, final GroupComparator B) {
		return new GroupComparator() {

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

	public static GroupComparator withDivisorsOk() {

		return new GroupComparator() {

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

	public static GroupComparator withMultiplesAndDivisorsOk() {
		return or(withMultiplesOk(), withDivisorsOk());
	}

	public static GroupComparator withMultiplesOk() {

		return new GroupComparator() {

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
