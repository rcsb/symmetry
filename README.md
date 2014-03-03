symmetry-benchmark
==================

This project was designed to benchmark results of CE-Symm, accompanying [symmetry](https://github.com/rcsb/symmetry). Please see that project for more general information.

Users may be interested in this project to:
 - Benchmark alternative data sets using CE-Symm.
 - Reproduce the benchmarks in the CE-Symm paper.
 - Benchmark protein symmetry-detection algorithms other than CE-Symm.

symmetry-benchmark has many dependencies and is best built using [Maven](http://maven.apache.org/).

### Benchmarking an alternative data set

The following sections describe how to benchmark the performance of CE-Symm against a benchmark set. (Another algorithm could also be used, but this is not discussed here.) The user should first obtain a checkout of the sourcecode and build it using Maven.

#### Obtaining a benchmark Sample file

The main method of ```SampleBuilder``` takes a CE-Symm results XML file and a tab-delimited file describing the correct, or expected, results. It produces an output XML file containing both of these data.

Usage: ```SampleBuilder results_file.xml benchmark_file.xml [correct_results.tsv]```

The results file (```results_file.xml```) is an XML serialization of the class CensusResultList in [symmetry](https://github.com/rcsb/symmetry). CE-Symm's command-line produces a file of this format when the ```-xml``` option is specified.

An example of the tab-delimited file (```correct_results.tsv```) can be found at https://github.com/rcsb/symmetry-benchmark/blob/master/src/main/resources/domain_symm_benchmark.tsv If this argument is not specified, that file  will be used. Only the numerical order is important for most uses.

The output file (```benchmark_file.xml```) is an XML serialization of the class Sample in symmetry-benchmark.

#### Overall accuracy and ROC curves

The class ```AccuracyFinder``` can be used to determine the accuracy of CE-Symm in the binary decision of symmetry/asymmetric.

Receiver Operating Characteristic curves can be generated using ROCCurves.
Usage: ```ROCCurves benchmark_file.xml rocs.png```

Here, ```benchmark_file.xml``` is the file obtained from the previous step. A poorly rendered PNG of the ROCs will be rendered to ```rocs.png```.

The full results will also be printed to standard output.
For each metric (including TM-score and TM-score and order, as used in the paper), two lines are printed:
The first line contains the X-coordinates (false positive rate), and the second line contains the Y-coordinates (true positive rate).

#### Order-detection accuracy

The correctness of the order-detection can be measured using ```OrderAccuracy```.
Usage: ```OrderAccuracy benchmark_file.xml```
Alternatively, this Java snippet is appropriate where multiples of the correct order are deemed acceptable:
```
OrderAccuracy finder = new OrderAccuracy(benchmarkFile, SignificanceFactory.symmetric(), GroupGuesserFactory.withMultiplesOk());
System.out.println(finder);
```

### Package documentation
The core classes and subpackages are contained in ```org.biojava3.structure.align.symm.benchmark```, which contains ```Sample``` and ```SampleBuilder```.
 - ```org.biojava3.structure.align.symm.benchmark.comparison``` compares predictions against known results.
 - ```org.biojava3.structure.align.symm.benchmark.comparison.order``` compares predictions of symmetry order against known orders.
 - ```org.biojava3.structure.align.symm.benchmark.external``` was used to benchmark the CE-Symm competitor SymD.
 - ```org.biojava3.structure.align.symm.benchmark.folds``` was used to reproduce the secondary benchmark.
 - ```org.biojava3.structure.align.symm.benchmark.set``` was used to generate and verify the benchmark set itself, irrespective of predictions.


### Build Status
[![Build Status](https://travis-ci.org/rcsb/symmetry-benchmark.png)](https://travis-ci.org/rcsb/symmetry-benchmark)
