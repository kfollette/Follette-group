//Elijah Spiro
//Version 1.0 - 9/17/16
//Version 2.0 - 3/23/17

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.IOException;
import javax.swing.ImageIcon;


public class GetParameters{
    
    static JTextField filepath2;
    static JTextField annuli2;
    static JTextField IWA2;
    static JTextField movement2;
    static JTextField output2;
    static JTextField klmodes2;
    static JTextField subsections2;
    static JCheckBox snr;
    static JTextField a1;
    static JTextField a2;
    static JTextField a3;
    static JTextField m1;
    static JTextField m2;
    static JTextField m3;
    static JTextField s1;
    static JTextField s2;
    static JTextField s3;
    static JTextField wid;
    static JTextField pa;
    static JTextField rad;
    
    public static void main(String[] args) throws FileNotFoundException{
	
	relaunchMain();

    }

    public static void relaunchMain() throws FileNotFoundException{
    
	final JFrame frame = new JFrame("Master KLIP");
        frame.setSize(400,150);
        frame.setLocationRelativeTo(null);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().setBackground(Color.black);
        frame.setLayout(null);

        addFields1(frame);
        frame.setVisible(true);
        frame.repaint();

	
    }

    public static void singleReduction() throws FileNotFoundException{
	final JFrame frame = new JFrame("Single KLIP Reduction");
        frame.setSize(400,400);
        frame.setLocationRelativeTo(null);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().setBackground(Color.black);
        frame.setLayout(null);

        addFields2(frame);

        frame.repaint();
        frame.setVisible(true);
    }

    public static void automateGUI() throws FileNotFoundException{
        final JFrame frame = new JFrame("Automated Parameter Search");
        frame.setSize(400,400);
        frame.setLocationRelativeTo(null);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().setBackground(Color.black);
        frame.setLayout(null);

        addFields3(frame);

        frame.repaint();
        frame.setVisible(true);
    }

     
    
    public static void addFields1(final JFrame frame) throws FileNotFoundException{
    
	JButton single = new JButton("Single Reduction");
        single.setSize(150,40);
        single.setLocation(25,40);
        single.setOpaque(true);
        single.setForeground(Color.blue);
        single.addActionListener(new ActionListener()
	    {
		public void actionPerformed(ActionEvent e) {
		    try{
			frame.dispose();
			singleReduction();
		    } catch (Exception e1){}
		}
	    });

	JButton automated = new JButton("Automate Parameters");
	automated.setSize(150,40);
        automated.setLocation(225,40);
        automated.setOpaque(true);
        automated.setForeground(Color.blue);
        automated.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent e) {
                    try{
			frame.dispose();
			automateGUI();
                    } catch (Exception e1){}
                }
            });
	
	frame.getContentPane().add(single);
	frame.getContentPane().add(automated);
	
    }

    public static void addFields2(final JFrame frame) throws FileNotFoundException{

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

        output2 = new JTextField("Klipped_Image");
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

	snr = new JCheckBox("SNR Analysis");
	snr.setLocation(67,235);
	snr.setSize(120,50);
	snr.setVisible(true);
	snr.setForeground(Color.white);
        
   	final JCheckBox save = new JCheckBox("Save 4D KLIP data");
	save.setLocation(187,235);
	save.setSize(170,50);
	save.setVisible(true);
	save.setForeground(Color.white);    

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
	
	JButton back = new JButton("<--");
        back.setSize(50,30);
        back.setLocation(20, 325);
        back.setOpaque(true);
        back.setForeground(Color.red);
        back.addActionListener(new ActionListener()
	    {
		public void actionPerformed(ActionEvent e){
		    try{
			frame.dispose();
			relaunchMain();
		    }
		    catch (Exception e1){
			e1.printStackTrace();
		    }
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
		    PrintWriter writer = new PrintWriter("single_reduction_parameters.txt", "UTF-8");
		    writer.println(filepath2.getText());
		    writer.println(annuli2.getText());
		    writer.println(IWA2.getText());
		    writer.println(movement2.getText());
		    writer.println(output2.getText());
		    writer.println(klmodes2.getText());
		    writer.println(subsections2.getText());
		    writer.println(snr.isSelected());
            writer.println(save.isSelected());
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
	frame.getContentPane().add(snr);
    frame.getContentPane().add(save);
	frame.getContentPane().add(back);
	
    }    

    public static void addFields3(final JFrame frame) throws FileNotFoundException{
	
	final JLabel annuli = new JLabel("Annuli");
        annuli.setForeground(Color.white);
        annuli.setSize(140,40);
        annuli.setLocation(40, 2+15);
	annuli.setFont(new Font("Serif", Font.BOLD, 20));
	

	final JLabel movement = new JLabel("Movement");
        movement.setForeground(Color.white);
        movement.setSize(140,40);
        movement.setLocation(40, 30+5+15);
        movement.setFont(new Font("Serif", Font.BOLD, 20));
        

	final JLabel subsections = new JLabel("Subsections");
        subsections.setForeground(Color.white);
        subsections.setSize(140,40);
        subsections.setLocation(40, 58+10+15);
        subsections.setFont(new Font("Serif", Font.BOLD, 20));


	JLabel divider = new JLabel("____________________________________________________________");
        divider.setSize(500,20);
        divider.setLocation(0,280+10+20);
        divider.setVisible(true);
        divider.setForeground(Color.blue);
	
	JButton back = new JButton("<--");
        back.setSize(50,30);
        back.setLocation(20, 325+12);
        back.setOpaque(true);
        back.setForeground(Color.red);
        back.addActionListener(new ActionListener()
            {
                public void actionPerformed(ActionEvent e){
                    try{
                        frame.dispose();
                        relaunchMain();
                    }
                    catch (Exception e1){
                        e1.printStackTrace();
                    }
                }
            });

	JLabel start = new JLabel("Start");
        start.setForeground(Color.white);
        start.setSize(100,50);
        start.setLocation(178,-10);
        
    JLabel stop = new JLabel("Stop");
        stop.setForeground(Color.white);
        stop.setSize(100,50);
        stop.setLocation(178+50,-10);
        
    JLabel inc = new JLabel("Increment");
        inc.setForeground(Color.white);
        inc.setSize(100,50);
        inc.setLocation(178+100,-10);
        
        
	a1 = new JTextField("3");
        a1.setBackground(Color.green);
        a1.setSize(50,30);
        a1.setLocation(200-30-10,7+17);

	a2 = new JTextField("7");
        a2.setBackground(Color.red);
        a2.setSize(50,30);
        a2.setLocation(250-30-10,7+17);

	a3 = new JTextField("1");
        a3.setBackground(Color.white);
        a3.setSize(100,30);
        a3.setLocation(300-30-10,7+17);

	m1 = new JTextField("5");
        m1.setBackground(Color.green);
        m1.setSize(50,30);
        m1.setLocation(200-30-10,35+5+17);
	

        m2 = new JTextField("7");
        m2.setBackground(Color.red);
        m2.setSize(50,30);
        m2.setLocation(250-30-10,35+5+17);
	

        m3 = new JTextField("1");
        m3.setBackground(Color.white);
        m3.setSize(100,30);
        m3.setLocation(300-30-10,35+5+17);
	

	s1 = new JTextField("1");
        s1.setBackground(Color.green);
        s1.setSize(50,30);
        s1.setLocation(200-30-10,63+10+17);
	

        s2 = new JTextField("1");
        s2.setBackground(Color.red);
        s2.setSize(50,30);
        s2.setLocation(250-30-10,63+10+17);
	

        s3 = new JTextField("1");
        s3.setBackground(Color.white);
        s3.setSize(100,30);
        s3.setLocation(300-30-10,63+10+17);
	

	JLabel output1 = new JLabel("Output Filename");
        output1.setForeground(Color.white);
        output1.setSize(225,40);
        output1.setLocation(230,173);

        output2 = new JTextField("star_date_set");
        output2.setBackground(Color.white);
        output2.setSize(180,30);
        output2.setLocation(206,203);

        JLabel filepath1 = new JLabel("Path to Desired Directory");
        filepath1.setForeground(Color.white);
        filepath1.setSize(200,40);
        filepath1.setLocation(118,117);

        filepath2 = new JTextField("example/star_name/date/sliced/");
        filepath2.setBackground(Color.white);
        filepath2.setSize(325,30);
        filepath2.setLocation(58,147);

        snr = new JCheckBox("SNR Analysis");
        snr.setLocation(135,235);
        snr.setSize(120,50);
        snr.setVisible(true);
        snr.setForeground(Color.white);

        final JButton searcher = new JButton("!");
        searcher.setSize(20,20);
        searcher.setLocation(25, 152);
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

    wid = new JTextField("10,15");
        wid.setBackground(Color.white);
        wid.setSize(100,30);
        wid.setLocation(300-30,281+2);
        
    pa = new JTextField("120");
        pa.setBackground(Color.white);
        pa.setSize(100,30);
        pa.setLocation(147,281+2);
        
    rad = new JTextField("13");
        rad.setBackground(Color.white);
        rad.setSize(100,30);
        rad.setLocation(25,281+2);
        
        
    JLabel widL = new JLabel("Width");
        widL.setForeground(Color.white);
        widL.setSize(100,40);
        widL.setLocation(300,254+2);

    JLabel paL = new JLabel("Position Angle");
        paL.setForeground(Color.white);
        paL.setSize(100,40);
        paL.setLocation(150,254+2);
        
    JLabel radL = new JLabel("Radius");
        radL.setForeground(Color.white);
        radL.setSize(100,40);
        radL.setLocation(55,254+2);


	JLabel IWA1 = new JLabel("IWA");
        IWA1.setForeground(Color.white);
        IWA1.setSize(100,40);
        IWA1.setLocation(27,173);

        IWA2 = new JTextField("10");
        IWA2.setBackground(Color.white);
        IWA2.setSize(50,30);
        IWA2.setLocation(15,203);

	JLabel klmodes1 = new JLabel("KL Modes");
        klmodes1.setForeground(Color.white);
        klmodes1.setSize(100,40);
        klmodes1.setLocation(100,173);

        klmodes2 = new JTextField("1,5,10,50,100");
        klmodes2.setBackground(Color.white);
        klmodes2.setSize(130,30);
        klmodes2.setLocation(71,203);

	JLabel planets = new JLabel("Planet Locations:");
        planets.setForeground(Color.white);
        planets.setSize(190,30);
        planets.setLocation(117,238);
	planets.setFont(new Font("Serif", Font.BOLD, 20));
        
  
   final JCheckBox save2 = new JCheckBox("Save 4D KLIP data");
	save2.setLocation(115,328);
	save2.setSize(170,50);
	save2.setVisible(true);
	save2.setForeground(Color.white);    
   
    

        
	JButton launcher = new JButton("Run KLIP");
        launcher.setSize(80,35);
        launcher.setLocation(298,335);
        launcher.setOpaque(true);
        launcher.setForeground(Color.blue);
        launcher.addActionListener(new ActionListener()
	    {
		public void actionPerformed(ActionEvent e) {
		    try{
			PrintWriter writer = new PrintWriter("automation_parameters.txt", "UTF-8");
			writer.println(filepath2.getText());
			writer.println(IWA2.getText());
			writer.println(klmodes2.getText());
			writer.println(output2.getText());
			writer.println(a1.getText());
			writer.println(a2.getText());
			writer.println(a3.getText());
			writer.println(m1.getText());
			writer.println(m2.getText());
			writer.println(m3.getText());
			writer.println(s1.getText());
			writer.println(s2.getText());
			writer.println(s3.getText());
            writer.println(rad.getText());
            writer.println(pa.getText());
            writer.println(wid.getText());
            writer.println(save2.isSelected());
			writer.close();
			System.exit(0);
		    } catch (Exception e1){}
		}
	    });

       
	frame.getContentPane().add(launcher);
	frame.getContentPane().add(planets);
	frame.getContentPane().add(klmodes1);
	frame.getContentPane().add(klmodes2);
	frame.getContentPane().add(annuli);
	frame.getContentPane().add(divider);
	frame.getContentPane().add(back);
	frame.getContentPane().add(a1);
	frame.getContentPane().add(a2);
	frame.getContentPane().add(a3);
	frame.getContentPane().add(movement);
        frame.getContentPane().add(m1);
        frame.getContentPane().add(m2);
        frame.getContentPane().add(m3);
	frame.getContentPane().add(output1);
        frame.getContentPane().add(output2);
        frame.getContentPane().add(filepath1);
        frame.getContentPane().add(filepath2);
	frame.getContentPane().add(searcher);
	frame.getContentPane().add(s1);
        frame.getContentPane().add(s2);
        frame.getContentPane().add(s3);
        frame.getContentPane().add(wid);
        frame.getContentPane().add(widL);
        frame.getContentPane().add(paL);
        frame.getContentPane().add(radL);
        frame.getContentPane().add(pa);
        frame.getContentPane().add(rad);
         frame.getContentPane().add(start);
        frame.getContentPane().add(stop);
        frame.getContentPane().add(inc);
	frame.getContentPane().add(subsections);
	frame.getContentPane().add(IWA1);
        frame.getContentPane().add(IWA2);
        frame.getContentPane().add(save2);


    }

}

