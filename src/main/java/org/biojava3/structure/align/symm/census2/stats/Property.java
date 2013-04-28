package org.biojava3.structure.align.symm.census2.stats;

import org.biojava3.structure.align.symm.census2.Result;

public abstract class Property {

	public abstract double getProperty(Result result) throws PropertyUndefinedException;
	
	public abstract String getName();
	
	@Override
	public String toString() {
		return getName();
	}

	public static Property identity() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAlignment().getIdentity();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "identity";
			}
		};
	}
	
	public static Property similarity() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAlignment().getSimilarity();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "similarity";
			}
		};
	}

	public static Property order() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getOrder();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "order";
			}
		};
	}

	public static Property hasOrder() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				if (result.getOrder() == null) return 0;
				return 1;
			}
			@Override
			public String getName() {
				return "order";
			}
		};
	}
	
	public static Property tmScore() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAlignment().getTmScore();
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "TM-score";
			}
		};
	}

	public static Property epsilon() {
		return new Property() {
			@Override
			public double getProperty(Result result) throws PropertyUndefinedException {
				try {
					return result.getAxis().evaluateEpsilon(result.getOrder());
				} catch (NullPointerException e) {
					throw new PropertyUndefinedException(e);
				}
			}
			@Override
			public String getName() {
				return "epsilon";
			}
		};
	}
	
}
