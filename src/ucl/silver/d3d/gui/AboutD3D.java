package ucl.silver.d3d.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import ucl.silver.d3d.core.*;

/**
 * <p>
 * Title: D3D</p>
 * <p>
 * Description: 3D Reaction-Diffusion Simulator</p>
 * <p>
 * Copyright: Copyright (c) 2022</p>
 * <p>
 * Company: The Silver Lab at University College London</p>
 *
 * @author Jason Rothman
 * @version 2.1
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * 
 * Contact j.rothman@ucl.ac.uk
 * https://silverlab.org/
 * https://github.com/SilverLabUCL/D3D
 */
public class AboutD3D extends JDialog implements ActionListener {

  JPanel panel1 = new JPanel();
  JPanel panel2 = new JPanel();
  JPanel insetsPanel1 = new JPanel();
  JPanel insetsPanel2 = new JPanel();
  JPanel insetsPanel3 = new JPanel();
  JButton button1 = new JButton();
  JLabel imageLabel = new JLabel();
  JLabel label1 = new JLabel();
  JLabel label2 = new JLabel();
  JLabel label3 = new JLabel();
  JLabel label4 = new JLabel();
  ImageIcon image1 = new ImageIcon();
  BorderLayout borderLayout1 = new BorderLayout();
  BorderLayout borderLayout2 = new BorderLayout();
  FlowLayout flowLayout1 = new FlowLayout();
  GridLayout gridLayout1 = new GridLayout();
  String product = "D3D, a 3D Reaction-Diffusion simulator";
  String version = "Version 2.0";
  String copyright = "Copyright (c) 2022";
  String comments = "Created by The Silver Lab at University College London";
  String nl = System.lineSeparator();
  String gnu1 = "This program is free software: you can redistribute it and/or modify it " + nl + "under the terms of the GNU General Public License as published by the Free Software Foundation, " + nl + "either version 3 of the License, or (at your option) any later version.";
  String gnu2 = "This program is distributed in the hope that it will be useful, " + nl + "but WITHOUT ANY WARRANTY; without even the implied warranty " + nl + "of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  " + nl + "See the GNU General Public License for more details.";
  String gnu3 = "You should have received a copy of the GNU General Public License along with this program.  " + nl + "If not, see <https://www.gnu.org/licenses/>.";
  
  public AboutD3D(Frame parent) {
    super(parent);
    MainFrame mf = (MainFrame) parent;
    enableEvents(AWTEvent.WINDOW_EVENT_MASK);
    try {
      jbInit(Master.D3D_VERSION);
    }
    catch(Exception e) {
    }
  }
  //Component initialization
  private void jbInit(String D3Dv) throws Exception  {

      //image1 = new ImageIcon(ucl.silver.d3d.gui.MainFrame.class.getResource("about.png"));
    //imageLabel.setIcon(image1);

    version = "Version " + D3Dv;

    setTitle("About");
    panel1.setLayout(borderLayout1);
    panel2.setLayout(borderLayout2);
    insetsPanel1.setLayout(flowLayout1);
    insetsPanel2.setLayout(flowLayout1);
    insetsPanel2.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
    gridLayout1.setRows(4);
    gridLayout1.setColumns(1);
    label1.setText(product);
    label2.setText(version);
    label3.setText(copyright);
    label4.setText(comments);
    insetsPanel3.setLayout(gridLayout1);
    insetsPanel3.setBorder(BorderFactory.createEmptyBorder(10, 60, 10, 10));
    button1.setText("Ok");
    button1.addActionListener(this);
    insetsPanel2.add(imageLabel, null);
    panel2.add(insetsPanel2, BorderLayout.WEST);
    this.getContentPane().add(panel1, null);
    insetsPanel3.add(label1, null);
    insetsPanel3.add(label2, null);
    insetsPanel3.add(label3, null);
    insetsPanel3.add(label4, null);
    panel2.add(insetsPanel3, BorderLayout.CENTER);
    insetsPanel1.add(button1, null);
    panel1.add(insetsPanel1, BorderLayout.SOUTH);
    panel1.add(panel2, BorderLayout.NORTH);
    setResizable(true);
    Master.log(nl);
    Master.log(product);
    Master.log(version);
    Master.log(copyright);
    Master.log(comments);
    Master.log(nl);
    Master.log(gnu1);
    Master.log(nl);
    Master.log(gnu2);
    Master.log(nl);
    Master.log(gnu3);
    Master.log(nl);
  }
  //Overridden so we can exit when window is closed
  @Override
  protected void processWindowEvent(WindowEvent e) {
    if (e.getID() == WindowEvent.WINDOW_CLOSING) {
      cancel();
    }
    super.processWindowEvent(e);
  }
  //Close the dialog
  void cancel() {
    dispose();
  }
  //Close the dialog on a button event
  public void actionPerformed(ActionEvent e) {
    if (e.getSource() == button1) {
      cancel();
    }
  }
}
