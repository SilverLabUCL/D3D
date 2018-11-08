package ucl.silver.d3d.gui;

import java.awt.*;
import javax.swing.*;

/**
 * <p>
 * Title: D3D</p>
 * <p>
 * Description: 3D Reaction-Diffusion Simulator</p>
 * <p>
 * Copyright: Copyright (c) 2018</p>
 * <p>
 * Company: The Silver Lab at University College London</p>
 *
 * @author Jason Rothman
 * @version 1.0
 */
public class PanelLog extends JPanel {

    JPanel jPanel = new JPanel();
    BorderLayout borderLayout = new BorderLayout();
    FlowLayout flowLayout = new FlowLayout();

    protected JTextArea textArea;
    protected JScrollPane scrollPane;

    private final static String newline = "\n";

    public PanelLog() {
        try {
            jbInit();
        } catch (Exception ex) {
            //ex.printStackTrace();
        }
    }

    final void jbInit() throws Exception {

        setLayout(borderLayout);

        jPanel.setMinimumSize(new Dimension(120, 250));
        jPanel.setPreferredSize(new Dimension(120, 250));
        jPanel.setLayout(flowLayout);

        textArea = new JTextArea();
        textArea.setEditable(false);

        scrollPane = new JScrollPane(textArea);

        add(scrollPane);
        
    }

    public void append(String str, boolean endWithNewLine) {
        if (endWithNewLine) {
            textArea.append(str + newline);
        } else {
            textArea.append(str);
        }
    }

    public void setText(String str, boolean endWithNewLine) {
        if (endWithNewLine) {
            textArea.setText(str + newline);
        } else {
            textArea.setText(str);
        }
    }

}
