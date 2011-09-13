package org.biojava3.changeux;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.beans.*;

@Deprecated
public class ProgressBar extends JPanel
                             implements ActionListener, 
                                        PropertyChangeListener {

    /**
	 * 
	 */
	private static final long serialVersionUID = 5760932164981676730L;
	private JProgressBar progressBar;
    
    private JTextArea taskOutput;
    private Task task;

    static String newline = System.getProperty("line.separator");
    class Task extends SwingWorker<Void, Void> {
        /*
         * Main task. Executed in background thread.
         */
        @Override
        public Void doInBackground() {
            
            int progress = 0;
            //Initialize progress property.
            setProgress(0);
            while (progress < 100) {
                //Sleep for up to one second.
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException ignore) {}
                //Make random progress.
                
                setProgress(Math.min(progress, 100));
            }
            return null;
        }

        /*
         * Executed in event dispatching thread
         */
        @Override
        public void done() {
        	//Toolkit.getDefaultToolkit().beep();
        	setProgress(100);
            setCursor(null); //turn off the wait cursor
            //taskOutput.append("Done!" + newline );
        }
    }
    
    public void addStatus(String txt){
    	taskOutput.append(txt + "" +newline);
    }
    
    public void done() {
    	
        Toolkit.getDefaultToolkit().beep();
        task.done();
        progressBar.setValue(100);
        progressBar.setIndeterminate(false);
        
        setCursor(null); //turn off the wait cursor
        
    }

    public ProgressBar() {
        super(new BorderLayout());

        //Create the demo's UI.
        //startButton = new JButton("Start");
        //startButton.setActionCommand("start");
        //startButton.addActionListener(this);

        progressBar = new JProgressBar(0, 50);
        progressBar.setValue(0);
        //progressBar.setStringPainted(false);
        progressBar.setIndeterminate(true);
        taskOutput = new JTextArea(5, 50);
        taskOutput.setMargin(new Insets(5,5,5,5));
        taskOutput.setEditable(false);

        JPanel panel = new JPanel();

        panel.add(progressBar);

        add(panel, BorderLayout.PAGE_START);
        add(new JScrollPane(taskOutput), BorderLayout.CENTER);
        setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));

    }

    /**
     * Invoked when the user presses the start button.
     */
    public void actionPerformed(ActionEvent evt) {
        
        setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        //Instances of javax.swing.SwingWorker are not reusuable, so
        //we create new instances as needed.
        task = new Task();
        task.addPropertyChangeListener(this);
        task.execute();
    }

    /**
     * Invoked when task's progress property changes.
     */
    public void propertyChange(PropertyChangeEvent evt) {
        if ("progress" == evt.getPropertyName()) {
            int progress = (Integer) evt.getNewValue();
            progressBar.setValue(progress);
            //taskOutput.append(String.format(
            //        "Completed %d%% of task.\n", task.getProgress()));
        } 
    }


    /**
     * Create the GUI and show it. As with all GUI code, this must run
     * on the event-dispatching thread.
     */
    public static ProgressBar createAndShowGUI() {
        //Create and set up the window.
        JFrame frame = new JFrame("Detecting Symmetry");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Create and set up the content pane.
        ProgressBar newContentPane = new ProgressBar();
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
        return newContentPane;
    }

    public static void main(String[] args) {
        //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGUI();
            }
        });
    }
}
