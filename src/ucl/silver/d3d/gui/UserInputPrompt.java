package ucl.silver.d3d.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

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
 */
public class UserInputPrompt
    extends JPanel implements ActionListener {

  UserInputAction action;

  JFrame frame;

  public String callerID;

  public UserInputPrompt(String id, String title, UserInputField fArray[], UserInputAction ua) {

    super(new BorderLayout(10, 10));

    action = ua;
    callerID = id;

    GridLayout g = new GridLayout(0, 1);
    g.setVgap(10);

    JPanel labelPane = new JPanel(g);
    JPanel fieldPane = new JPanel(g);
    JPanel buttonPane = new JPanel();

    for (int i = 0; i < fArray.length; i += 1) {
      labelPane.add(fArray[i].iLabel);
      fieldPane.add(fArray[i].iField);
    }

    setBorder(BorderFactory.createEmptyBorder(20, 40, 20, 40));

    this.add(labelPane, BorderLayout.CENTER);
    this.add(fieldPane, BorderLayout.LINE_END);

    JButton ok = new JButton("OK");
    ok.addActionListener(this);
    buttonPane.add(ok);

    JButton cancel = new JButton("Cancel");
    cancel.addActionListener(this);
    buttonPane.add(cancel);

    add(buttonPane, BorderLayout.SOUTH);

    createAndShowGUI(title);

  }

  public void createAndShowGUI(String title) {

    // make sure we have nice window decorations.
    JFrame.setDefaultLookAndFeelDecorated(true);

    // create and set up the window.
    frame = new JFrame(title);

    // create and set up the content pane.
    JComponent newContentPane = this;
    newContentPane.setOpaque(true); //content panes must be opaque
    frame.setContentPane(newContentPane);
    frame.pack();

    // center window and display
    Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    Dimension frameSize = frame.getSize();
    frame.setLocation( (screenSize.width - frameSize.width) / 2,
                      (screenSize.height - frameSize.height) / 2);
    frame.setVisible(true);

  }

  /**
   * called when the user clicks button
   */
  public void actionPerformed(ActionEvent e) {
    action.actionButton(e);
  }

  public void exit() {
    frame.dispose();
  }

}
