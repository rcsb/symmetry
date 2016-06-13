package writers;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Parent class for all output formats. Writers need to implement the
 * writeHeader method and define the a method to write results to the file.
 * <p>
 * All the output writers need to be thread safe, because parallelization is
 * expected in the calculation.
 * 
 * @author Aleix Lafita
 * 
 */
public abstract class OutputWriter {
	
	protected static final Logger logger = LoggerFactory
			.getLogger(OutputWriter.class);
	
	protected PrintWriter writer;

	/**
	 * Constructor with a 'filename'. Opens the file and initializes a
	 * PrintWriter.
	 * 
	 * @param filename
	 * @throws IOException
	 */
	public OutputWriter(String filename) throws IOException {
		this.writer = openOutputFile(filename);
	}

	/**
	 * Writes the first line with headers for each column of results.
	 * Implementations of this method have to be synchronized.
	 * 
	 * @throws IOException
	 */
	abstract public void writeHeader() throws IOException;

	/**
	 * Flush and close the writer.
	 */
	public synchronized void close() {
		if (writer != null) {
			writer.flush();
			writer.close();
		}
	}

	/**
	 * Opens 'filename' for writing.
	 * 
	 * @param filename
	 *            Name of output file, or '-' for standard out
	 * @throws IOException
	 */
	private static PrintWriter openOutputFile(String filename)
			throws IOException {
		if (filename.equals("-")) {
			return new PrintWriter(System.out, true);
		}
		return new PrintWriter(new BufferedWriter(new FileWriter(filename)));
	}

}