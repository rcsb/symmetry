package org.biojava3.structure.align.symm.census3.analysis.ligands;

import java.util.LinkedHashMap;
import java.util.Map;

import org.biojava.bio.structure.Element;
import org.biojava.bio.structure.ElementType;

/**
 * Something to decide "yes" or "no" given a {@link CensusLigand}.
 * @author dmyersturnbull
 */
public abstract class LigandMatcher {

	public static LigandMatcher and(final LigandMatcher... as) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				for (LigandMatcher a : as) {
					if (!a.matches(ligand)) return false;
				}
				return true;
			}
		};
	}

	public static LigandMatcher everything() {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				return true;
			}
		};
	}

	public static LigandMatcher hasAtomicMass(final double min) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (e.getAtomicMass() >= min) return true;
				}
				return false;
			}
		};
	}

	public static LigandMatcher hasElectronegativity(final float min) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (e.getPaulingElectronegativity() >= min) return true;
				}
				return false;
			}
		};
	}

	public static LigandMatcher hasElement(final Element element) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				// System.out.println(ligand.getFormula());
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (e.isMetal()) {
					}
					if (e.equals(element)) return true;
				}
				return false;
			}
		};
	}

	public static LigandMatcher hasElementType(final ElementType elementType) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (e.getElementType().equals(elementType)) return true;
				}
				return false;
			}
		};
	}

	public static LigandMatcher hasMetal() {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				return ligand.isMetallic();
			}
		};
	}

	public static LigandMatcher hasTotalOxidationMagnitude(final int min) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (e.equals(Element.R)) continue;
					int sum = 0;
					for (int o : e.getAllOxidationStates()) {
						sum += o;
					}
					if (Math.abs(sum) >= min) return true;
				}
				return false;
			}
		};
	}

	public static LigandMatcher hasOxidationMagnitude(final int min) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (e.equals(Element.R)) continue;
					for (int o : e.getAllOxidationStates()) {
						if (Math.abs(o) >= min) return true;
					}
				}
				return false;
			}
		};
	}

	public static LigandMatcher hasPeriod(final int min) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (e.getPeriod() >= min) return true;
				}
				return false;
			}
		};
	}

	public static LigandMatcher isOnlyMetal() {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				for (Element e : elements.keySet()) {
					if (!e.isMetal()) return false;
				}
				return true;
			}
		};
	}

	public static LigandMatcher isOrganicMetal() {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				boolean containsMetal = false;
				boolean containsCarbon = false;
				for (Element e : elements.keySet()) {
					if (e.isMetal()) containsMetal = true;
					if (e.equals(Element.C)) containsCarbon = true;
				}
				return containsMetal && containsCarbon;
			}
		};
	}

	public static LigandMatcher isTrueSalt() {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				boolean containsMetal = false;
				boolean containsNonmetal = false;
				for (Element e : elements.keySet()) {
					if (e.getElementType().equals(ElementType.ALKALI_METAL)
							|| e.getElementType().equals(ElementType.ALKALINE_EARTH_METAL)) containsMetal = true;
					if (e.isNonMetal()) containsNonmetal = true;
				}
				return containsMetal && containsNonmetal;
			}
		};
	}

	public static LigandMatcher nElements(final int n) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				Map<Element, Integer> elements = parse(ligand);
				return elements.size() == n;
			}
		};
	}

	public static LigandMatcher nothing() {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				return false;
			}
		};
	}

	public static LigandMatcher or(final LigandMatcher... as) {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				for (LigandMatcher a : as) {
					if (a.matches(ligand)) return true;
				}
				return false;
			}
		};
	}

	public static Map<Element, Integer> parse(CensusLigand ligand) {
		Map<Element, Integer> compounds = new LinkedHashMap<Element, Integer>();
		StringBuilder currentCompound = null;
		int currentNumber = 1;
		for (int i = 0; i < ligand.getFormula().length(); i++) {
			char c = ligand.getFormula().charAt(i);
			if (Character.isUpperCase(c)) {
				if (currentCompound != null) {
					compounds.put(Element.valueOfIgnoreCase(currentCompound.toString()), currentNumber);
				}
				currentCompound = new StringBuilder();
				currentNumber = 1;
				currentCompound.append(c);
			} else if (Character.isDigit(c)) {
				currentNumber = Character.getNumericValue(c);
			} else {
				currentCompound.append(c);
			}
		}
		if (currentCompound != null) {
			compounds.put(Element.valueOfIgnoreCase(currentCompound.toString()), currentNumber);
		}
		return compounds;
	}

	public boolean containsAMatch(LigandsOfStructure ligands) {
		for (CensusLigand ligand : ligands.getLigands()) {
			if (matches(ligand)) return true;
		}
		return false;
	}

	public abstract boolean matches(CensusLigand ligand);

	public LigandMatcher not() {
		return new LigandMatcher() {
			@Override
			public boolean matches(CensusLigand ligand) {
				return !LigandMatcher.this.matches(ligand);
			}
		};
	}

}
