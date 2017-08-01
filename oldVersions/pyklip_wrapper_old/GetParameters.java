//Elijah Spiro
//Version 1.0 - 9/17/16

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class GetParameters{
    
    static JTextField filepath2;
    static JTextField annuli2;
    static JTextField IWA2;
    static JTextField movement2;
    static JTextField output2;
    static JTextField klmodes2;
    static JTextField subsections2;

    public static void main(String[] args) throws FileNotFoundException{
	
	final JFrame frame = new JFrame("KLIP Launcher");
	frame.setSize(400,400);
	frame.setLocationRelativeTo(null);
	frame.setVisible(true);
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	frame.getContentPane().setBackground(Color.black);

	addFields(frame);
	
	frame.repaint();

    }

    public static void addFields(JFrame frame) throws FileNotFoundException{

	JLabel annuli1 = new JLabel("Annuli");
        annuli1.setForeground(Color.white);
        annuli1.setSize(100,40);
        annuli1.setLocation(20,20);

        annuli2 = new JTextField("9");
        annuli2.setBackground(Color.white);
        annuli2.setSize(50,30);
        annuli2.setLocation(16,50);

	JLabel IWA1 = new JLabel("IWA");
	IWA1.setForeground(Color.white);
        IWA1.setSize(100,40);
	IWA1.setLocation(85,20);

	IWA2 = new JTextField("10");
	IWA2.setBackground(Color.white);
	IWA2.setSize(50,30);
        IWA2.setLocation(76,50);

	JLabel movement1 = new JLabel("Movement");
        movement1.setForeground(Color.white);
        movement1.setSize(100,40);
        movement1.setLocation(132,20);

        movement2 = new JTextField("2.5");
        movement2.setBackground(Color.white);
        movement2.setSize(50,30);
        movement2.setLocation(135,50);
	
	JLabel klmodes1 = new JLabel("KL Modes");
        klmodes1.setForeground(Color.white);
        klmodes1.setSize(100,40);
        klmodes1.setLocation(250,20);

        klmodes2 = new JTextField("1,2,3,4,5,10,20,50,100");
        klmodes2.setBackground(Color.white);
        klmodes2.setSize(180,30);
        klmodes2.setLocation(200,50);

	JLabel subsections1 = new JLabel("Subsections");
        subsections1.setForeground(Color.white);
        subsections1.setSize(100,40);
        subsections1.setLocation(60,90);

        subsections2 = new JTextField("1");
        subsections2.setBackground(Color.white);
        subsections2.setSize(50,30);
        subsections2.setLocation(67,120);

	JLabel output1 = new JLabel("Output Filename");
        output1.setForeground(Color.white);
        output1.setSize(200,40);
        output1.setLocation(200,90);

        output2 = new JTextField("KlippedImage");
        output2.setBackground(Color.white);
        output2.setSize(180,30);
        output2.setLocation(170,120);

	JLabel filepath1 = new JLabel("Path to Desired Directory");
        filepath1.setForeground(Color.white);
        filepath1.setSize(200,40);
        filepath1.setLocation(150,170);

        filepath2 = new JTextField("example/star_name/date/sliced/");
        filepath2.setBackground(Color.white);
        filepath2.setSize(325,30);
        filepath2.setLocation(70,200);

	final JButton searcher = new JButton("!");
	searcher.setSize(20,20);
	searcher.setLocation(30, 204);
	searcher.setOpaque(true);
	searcher.setForeground(Color.blue);
	searcher.addActionListener(new ActionListener() 
	{
	    public void actionPerformed(ActionEvent e){
		final JFileChooser fc = new JFileChooser();
		fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		fc.showOpenDialog(searcher);
		String path = (fc.getSelectedFile().getAbsolutePath());
		filepath2.setText(path);
	    }
	});
	

	JButton launcher = new JButton("Run KLIP");
	launcher.setSize(100,40);
	launcher.setLocation(150,320);
	launcher.setOpaque(true);
	launcher.setForeground(Color.blue);
	launcher.addActionListener(new ActionListener()
	{
	    public void actionPerformed(ActionEvent e) {
		try{
		PrintWriter writer = new PrintWriter("parameters.txt", "UTF-8");
		writer.println(filepath2.getText());
		writer.println(annuli2.getText());
		writer.println(IWA2.getText());
		writer.println(movement2.getText());
		writer.println(output2.getText());
		writer.println(klmodes2.getText());
		writer.println(subsections2.getText());
		writer.close();
		System.exit(0);
		} catch (Exception e1){}
	    }
	});

	JLabel divider = new JLabel("____________________________________________________________");
	divider.setSize(500,20);
	divider.setLocation(0,280);
	divider.setVisible(true);
	divider.setForeground(Color.blue);
	
	frame.getContentPane().add(annuli1);
        frame.getContentPane().add(annuli2);
	frame.getContentPane().add(IWA1);
	frame.getContentPane().add(IWA2);
        frame.getContentPane().add(movement1);
        frame.getContentPane().add(movement2);
	frame.getContentPane().add(klmodes1);
	frame.getContentPane().add(klmodes2);
	frame.getContentPane().add(subsections1);
	frame.getContentPane().add(subsections2);
	frame.getContentPane().add(output1);
	frame.getContentPane().add(output2);
	frame.getContentPane().add(launcher);
	frame.getContentPane().add(divider);
	frame.getContentPane().add(filepath1);
	frame.getContentPane().add(filepath2);
	frame.getContentPane().add(searcher);
       
    }    


}

