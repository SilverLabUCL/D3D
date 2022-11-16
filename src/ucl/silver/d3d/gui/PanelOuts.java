package ucl.silver.d3d.gui;

import java.awt.*;
import javax.swing.JPanel;
import ucl.silver.d3d.core.*;
import ucl.silver.d3d.utils.*;

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

public class PanelOuts extends JPanel {

  Detector[] detect;

  private double out[][] = null;

  private int width = 0, height = 0; // panel width, height

  BorderLayout borderLayout1 = new BorderLayout();
  FlowLayout FlowLayout1 = new FlowLayout();

  JPanel jPanelLeft = new JPanel();

  public PanelOuts() {
    try {
      detect = Master.project.detectors;
      jbInit();
    }
    catch(Exception ex) {
      //ex.printStackTrace();
    }
  }

  final void jbInit() throws Exception {
    this.setLayout(borderLayout1);
    jPanelLeft.setMinimumSize(new Dimension(120, 250));
    jPanelLeft.setPreferredSize(new Dimension(120, 250));
    jPanelLeft.setLayout(FlowLayout1);
  }

  @Override
  public void paintComponent(Graphics g) {

    int x0, y0, x1, y1;
    double dmin = 0, dmax = 0, dh;
    int i, j, k, iw, jinc;

    int dataPoint = 0;

    if (detect == null) {
      return;
    }

    width = this.getWidth();
    height = this.getHeight();
    height *= 0.9;

    g.setColor(Color.white);
    g.fillRect(0, 0, width, height); // clear panel

    for (Detector d : detect) {

      out = d.save.data;

      if (out == null) {
        continue;
      }

      dmin = Math.min(dmin, Utility.getMin(out));
      dmax = Math.max(dmax, Utility.getMax(out));

    }

    dh = Math.abs(dmax - dmin);
    i = 0;

    for (Detector d : detect) {

      out = d.save.data;

      if (out == null) {
        continue;
      }

      iw = 0;
      jinc = 0;

      while (iw <= 1) {
        jinc++;
        iw = width * jinc / out.length;
      }

      x1 = 0;
      y1 = (int) ((height * 1.0) * (out[dataPoint][0] - dmin) / dh);
      y1 = -1 * y1 + height;

      k = 1;

      for (j = 1; j < out.length; j+=jinc) {

        x0 = x1;
        y0 = y1;

        x1 = k * iw;
        y1 = (int) ((height * 1.0) * (out[dataPoint][j] - dmin) / dh);
        y1 = -1 * y1 + height;

        g.setColor(getColor(i));
        g.drawLine(x0, y0, x1, y1);

        k++;

      }
      
      i++;

    }

  }

  public Color getColor(int i) {
    switch(i) {
      case 0:
        return Color.red;
      case 1:
        return Color.blue;
      case 2:
        return Color.green;
      case 3:
        return Color.black;
      case 4:
        return Color.pink;
      case 5:
        return Color.magenta;
      case 6:
        return Color.orange;
      case 7:
        return Color.cyan;
      case 8:
        return Color.lightGray;
      case 9:
        return Color.yellow;
    }
    return Color.red;
  }

}
