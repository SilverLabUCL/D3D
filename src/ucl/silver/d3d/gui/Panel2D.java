package ucl.silver.d3d.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
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
 */
public class Panel2D
        extends JPanel {

    public Grid grid2D = null;

    JPanel jPanelTop = new JPanel();
    JPanel jPanelLeft = new JPanel();
    JPanel jPanelLeft1 = new JPanel();
    JPanel jPanelLeft2 = new JPanel();
    JPanel jPanelBottom = new JPanel();
    JPanel jPanelGrid = new JPanel();
    JPanel jPanelSpacer1 = new JPanel();
    JPanel jPanelSpacer2 = new JPanel();

    JLabel labelDisplay = new JLabel();
    JComboBox<String> comboDisplayMode = new JComboBox<String>();
    JRadioButton radioEditMode = new JRadioButton();
    JButton buttonColor = new JButton();

    JLabel labelPlane = new JLabel();
    JSpinner spinnerPlane = new JSpinner();
    JSlider sliderPlane = new JSlider();
    JRadioButton radioAxisMode = new JRadioButton();
    JButton buttonSwitch = new JButton();
    JButton buttonCopy = new JButton();
    JButton buttonPaste = new JButton();
    JButton buttonClear = new JButton();

    JLabel labelZoom = new JLabel();
    JSlider sliderZoom = new JSlider();
    JRadioButton cross = new JRadioButton();
    JRadioButton voxelGrid = new JRadioButton();
    JRadioButton symmetry4 = new JRadioButton();

    JTextField xCoord = new JTextField();
    JTextField yCoord = new JTextField();
    JTextField zCoord = new JTextField();
    JTextField vCoord = new JTextField();
    JRadioButton coordDim = new JRadioButton();
    JRadioButton coordClick = new JRadioButton();

    JScrollPane jScrollPaneGrid = new JScrollPane();

    private String comboEditString = "All";

    private boolean controlsOff = true;

    public Panel2D() {
        try {
            jbInit();
        } catch (Exception ex) {
            //ex.printStackTrace();
        }
    }

    final void jbInit() throws Exception {

        int ydim = 25;

        jPanelLeft.setMinimumSize(new Dimension(170, 250));
        jPanelLeft.setPreferredSize(new Dimension(170, 250));

        jPanelLeft1.setBorder(BorderFactory.createEtchedBorder());
        jPanelLeft1.setMinimumSize(new Dimension(130, 255));
        jPanelLeft1.setPreferredSize(new Dimension(130, 255));
        jPanelLeft1.setLayout(new FlowLayout());

        jPanelLeft2.setLayout(new FlowLayout());
        jPanelLeft2.setBorder(BorderFactory.createEtchedBorder());
        jPanelLeft2.setMinimumSize(new Dimension(130, 260));
        jPanelLeft2.setPreferredSize(new Dimension(130, 260));

        jPanelBottom.setLayout(new FlowLayout(FlowLayout.CENTER));

        jPanelSpacer1.setMinimumSize(new Dimension(100, 5));
        jPanelSpacer1.setPreferredSize(new Dimension(100, 5));

        jPanelSpacer2.setMinimumSize(new Dimension(140, 5));
        jPanelSpacer2.setPreferredSize(new Dimension(140, 5));
        
        jPanelGrid.setAlignmentY((float) 0.5);
        //jPanelGrid.setMinimumSize(new Dimension(100, 100));
        //jPanelGrid.setPreferredSize(new Dimension(100, 100));

        labelDisplay.setText("DISPLAY");
        labelDisplay.setFont(new Font("Arial", Font.BOLD, 12));
        labelDisplay.setPreferredSize(new Dimension(100, ydim));
        labelDisplay.setMaximumSize(new Dimension(100, ydim));
        labelDisplay.setHorizontalAlignment(JLabel.CENTER);

        comboDisplayMode.setToolTipText("select what to edit");
        comboDisplayMode.addItemListener(new Panel2D_comboDisplayMode_itemAdapter(this));
        comboDisplayMode.setPreferredSize(new Dimension(100, ydim));

        radioEditMode.setToolTipText("edit mode selector");
        radioEditMode.setText("edit off");
        radioEditMode.addActionListener(new Panel2D_radioEditMode_actionAdapter(this));
        radioEditMode.setMaximumSize(new Dimension(100, ydim));
        radioEditMode.setPreferredSize(new Dimension(100, ydim));
        radioEditMode.setEnabled(false);

        buttonColor.setToolTipText("set display color");
        buttonColor.setPreferredSize(new Dimension(100, ydim));
        buttonColor.setText("Color");
        buttonColor.addMouseListener(new Panel2D_buttonColor_mouseAdapter(this));

        labelPlane.setText("PLANE");
        labelPlane.setFont(new Font("Arial", Font.BOLD, 12));
        labelPlane.setPreferredSize(new Dimension(100, ydim));
        labelPlane.setMaximumSize(new Dimension(100, ydim));
        labelPlane.setHorizontalAlignment(JLabel.CENTER);

        spinnerPlane.setToolTipText("plane number");
        spinnerPlane.addChangeListener(new Panel2D_spinnerPlane_changeAdapter(this));
        spinnerPlane.setPreferredSize(new Dimension(60, ydim));

        sliderPlane.setToolTipText("plane number");
        sliderPlane.addChangeListener(new Panel2D_sliderPlane_changeAdapter(this));
        sliderPlane.setPreferredSize(new Dimension(100, ydim));

        radioAxisMode.setToolTipText("xyz plane orientation selector");
        radioAxisMode.setHorizontalAlignment(SwingConstants.CENTER);
        radioAxisMode.setText("xy plane");
        radioAxisMode.addActionListener(new Panel2D_radioAxisMode_actionAdapter(this));
        radioAxisMode.setMaximumSize(new Dimension(100, ydim));
        radioAxisMode.setPreferredSize(new Dimension(100, ydim));

        buttonSwitch.setToolTipText("rotate/switch axes");
        buttonSwitch.setPreferredSize(new Dimension(100, ydim));
        buttonSwitch.setText("Rotate");
        buttonSwitch.addMouseListener(new Panel2D_buttonSwitch_mouseAdapter(this));

        buttonCopy.setToolTipText("copy selected plane");
        buttonCopy.setPreferredSize(new Dimension(100, ydim));
        buttonCopy.setText("Copy");
        buttonCopy.addMouseListener(new Panel2D_buttonCopy_mouseAdapter(this));

        buttonPaste.setToolTipText("paste copied plane");
        buttonPaste.setPreferredSize(new Dimension(100, ydim));
        buttonPaste.setText("Paste");
        buttonPaste.addMouseListener(new Panel2D_buttonPaste_mouseAdapter(this));

        buttonClear.setToolTipText("clear selected plane");
        buttonClear.setPreferredSize(new Dimension(100, ydim));
        buttonClear.setText("Clear");
        buttonClear.addMouseListener(new Panel2D_buttonClear_mouseAdapter(this));

        labelZoom.setText("Zoom");
        labelZoom.setPreferredSize(new Dimension(30, ydim));

        sliderZoom.setToolTipText("zoom");
        sliderZoom.setName("Zoom");
        sliderZoom.addChangeListener(new Panel2D_sliderZoom_changeAdapter(this));
        sliderZoom.setPreferredSize(new Dimension(70, ydim));

        cross.setText("0 cross");
        cross.addActionListener(new Panel2D_cross_actionAdapter(this));
        cross.setMaximumSize(new Dimension(100, ydim));
        cross.setPreferredSize(new Dimension(100, ydim));

        voxelGrid.setText("grid");
        voxelGrid.addActionListener(new Panel2D_voxelGrid_actionAdapter(this));
        voxelGrid.setMaximumSize(new Dimension(100, ydim));
        voxelGrid.setPreferredSize(new Dimension(100, ydim));

        symmetry4.setText("4 symmetry");
        symmetry4.addActionListener(new Panel2D_symmetry4_actionAdapter(this));
        symmetry4.setMaximumSize(new Dimension(100, ydim));
        symmetry4.setPreferredSize(new Dimension(100, ydim));

        xCoord.setText("");
        xCoord.setPreferredSize(new Dimension(100, ydim));

        yCoord.setText("");
        yCoord.setPreferredSize(new Dimension(100, ydim));

        zCoord.setText("");
        zCoord.setPreferredSize(new Dimension(100, ydim));

        vCoord.setText("");
        vCoord.setPreferredSize(new Dimension(100, ydim));

        coordDim.addActionListener(new Panel2D_coordDim_actionAdapter(this));
        coordDim.setMaximumSize(new Dimension(80, ydim));
        coordDim.setPreferredSize(new Dimension(80, ydim));
        
        coordClick.addActionListener(new Panel2D_coordClick_actionAdapter(this));
        coordClick.setMaximumSize(new Dimension(120, ydim));
        coordClick.setPreferredSize(new Dimension(120, ydim));

        setLayout(new BorderLayout());
        addComponentListener(new Panel2D_this_componentAdapter(this));

        jPanelLeft1.add(labelDisplay, null);
        jPanelLeft1.add(comboDisplayMode, null);
        jPanelLeft1.add(buttonColor, null);
        jPanelLeft1.add(radioEditMode, null);
        jPanelLeft1.add(cross, null);
        jPanelLeft1.add(voxelGrid, null);
        jPanelLeft1.add(symmetry4, null);
        jPanelLeft1.add(labelZoom, null);
        jPanelLeft1.add(sliderZoom, null);

        jPanelLeft2.add(labelPlane, null);
        jPanelLeft2.add(spinnerPlane, null);
        jPanelLeft2.add(sliderPlane, null);
        jPanelLeft2.add(radioAxisMode, null);
        jPanelLeft2.add(buttonSwitch, null);
        jPanelLeft2.add(buttonClear, null);
        jPanelLeft2.add(buttonCopy, null);
        jPanelLeft2.add(buttonPaste, null);

        jPanelLeft.add(jPanelLeft1, null);
        jPanelLeft.add(jPanelSpacer1, null);
        jPanelLeft.add(jPanelLeft2, null);

        jPanelBottom.add(jPanelSpacer2, null);
        jPanelBottom.add(xCoord, null);
        jPanelBottom.add(yCoord, null);
        jPanelBottom.add(zCoord, null);
        jPanelBottom.add(vCoord, null);
        jPanelBottom.add(coordDim, null);
        jPanelBottom.add(coordClick, null);

        jPanelLeft.setOpaque(true);
        jPanelLeft.setLayout(new FlowLayout());

        grid2D = new Grid(1000, 800, this);

        jPanelGrid.add(jScrollPaneGrid, null);

        add(jPanelGrid, BorderLayout.CENTER);
        add(jPanelTop, BorderLayout.NORTH);
        add(jPanelLeft, BorderLayout.WEST);
        add(jPanelBottom, BorderLayout.SOUTH);

        jScrollPaneGrid.getViewport().add(grid2D);

    }

    public void updateControls() {

        grid2D.adjustGridSize();

        cross.setSelected(grid2D.cross);
        voxelGrid.setSelected(grid2D.voxelGrid);
        symmetry4.setSelected(grid2D.sym4);
        
        displayModeUpdate();
        objectColorUpdate();
        editModeUpdate();
        zoomSliderUpdate();

        gridPlaneUpdate();
        axesModeUpdate();
        updateEditButtons();
        
        coordinatesDimUpdate();
        coordinatesClickUpdate();

        controlsOff = false;

    }

    void displayModeUpdate() {

        int numDiffusants = Master.project.numDiffusants();
        int numSources = Master.project.numSources();
        int numDetectors = Master.project.numDetectors();

        int items = comboDisplayMode.getItemCount();

        controlsOff = true; // turn off because the following changes cause a state change

        if (items == 0) {
            comboDisplayMode.addItem("All");
            comboDisplayMode.addItem("---");
            comboDisplayMode.addItem("Geometry");
        } else if (items > 0) {
            for (int i = items - 1; i > 2; i--) {
                comboDisplayMode.removeItemAt(i);
            }
        }

        if (numDiffusants > 0) {

            comboDisplayMode.addItem("---");

            if (numDiffusants > 1) {
                comboDisplayMode.addItem("Diffusant.all");
            }

            for (int i = 0; i < numDiffusants; i++) {
                comboDisplayMode.addItem("Diffusant." + i);
            }

        }

        if (numSources > 0) {

            comboDisplayMode.addItem("---");

            if (numSources > 1) {
                comboDisplayMode.addItem("Source.all");
            }

            for (int i = 0; i < numSources; i++) {
                comboDisplayMode.addItem("Source." + i);
            }

        }

        if (numDetectors > 0) {

            comboDisplayMode.addItem("---");

            if (numDetectors > 1) {
                comboDisplayMode.addItem("Detector.all");
            }

            for (int i = 0; i < numDetectors; i++) {
                comboDisplayMode.addItem("Detector." + i);
            }

        }

        comboDisplayMode.setSelectedItem(comboEditString);

        controlsOff = false;

    }

    public void displayModeSet(String eSelect) {

        if (controlsOff) {
            return;
        }

        int emode = displayModeGet(eSelect);
        int anum = displayModeNum(eSelect);

        if ((emode >= 0) && (emode <= 4)) {
            comboEditString = eSelect;
            grid2D.setViewMode(emode);
            grid2D.setViewArrayNum(anum);
        } else {
            comboEditString = "All";
            grid2D.setViewMode(4);
            grid2D.setViewArrayNum(-1);
        }

        objectColorUpdate();
        editModeUpdate();
        updateEditButtons();

        grid2D.repaint();

    }

    public int displayModeGet(String s) {

        s = displayModePrefix(s);

        if (s.compareToIgnoreCase("Geometry") == 0) {
            return 0;
        } else if (s.compareToIgnoreCase("Diffusant") == 0) {
            return 1;
        } else if (s.compareToIgnoreCase("Source") == 0) {
            return 2;
        } else if (s.compareToIgnoreCase("Detector") == 0) {
            return 3;
        } else if (s.compareToIgnoreCase("All") == 0) {
            return 4;
        }

        return -1;

    }

    public String displayModePrefix(String item) {
        int i = item.lastIndexOf(".");
        if (i > 0) {
            return item.substring(0, i);
        } else {
            return item;
        }
    }

    public int displayModeNum(String item) {
        int i = item.lastIndexOf(".");
        int j;
        String n;

        if ((i > 0) && (i + 1 < item.length())) {
            n = item.substring(i + 1, item.length());
            try {
                j = Integer.parseInt(n);
                return j;
            } catch (NumberFormatException e) {
                return -1;
            }
        } else {
            return -1;
        }

    }

    public void objectColorUpdate() {

        String colorStr = grid2D.getCurrentColor();

        if (ColorD3D.isColorScale(colorStr)) {
            buttonColor.setText("Color Scale");
        } else {
            buttonColor.setText("Color");
        }

        if (colorStr == null) {
            buttonColor.setBackground(Color.white);
            buttonColor.setEnabled(false);
        } else {
            buttonColor.setBackground(ColorD3D.string2color(colorStr));
            buttonColor.setEnabled(true);
        }

    }

    public void objectColorSet() {

        String colorStr = grid2D.getCurrentColor();

        if ((colorStr == null) || (colorStr.length() == 0)) {
            return; // no color to change
        }

        colorStr = ColorD3D.promptColorOrScale(colorStr);

        if ((colorStr == null) || (colorStr.length() == 0)) {
            return; // cancel
        }

        grid2D.setCurrentColor(colorStr);

        objectColorUpdate();

    }

    public void editModeUpdate(){

        radioEditMode.setSelected(true);

        if (grid2D.displayMode == 0) {
            radioEditMode.setEnabled(true);
        } else {
            radioEditMode.setEnabled(false);
        }

        switch(grid2D.editMode) {
            case 0:
                radioEditMode.setText("edit off");
                break;
            case 1:
                radioEditMode.setText("non-space");
                break;
            case 2:
                radioEditMode.setText("space");
                break;
            default:
                radioEditMode.setText("marquee");
        }

    }

    public void editModeSet(int m) {

        if (controlsOff) {
            return;
        }

        switch (m) {
            case 0:
            case 1:
            case 2:
            case 3:
                break;
            default:
                m = 0;
        }

        grid2D.setEditMode(m);
        objectColorUpdate();
        editModeUpdate();

    }

    public void zoomSliderUpdate() {

        controlsOff = true;

        sliderZoom.setMinimum(grid2D.pixelsPerVoxelMin);
        sliderZoom.setMaximum(grid2D.pixelsPerVoxelMax);
        sliderZoom.setValue(grid2D.pixelsPerVoxel);

        controlsOff = false;

    }

    public void resetVoxelWidth() {
        grid2D.initVoxelWidth(-1, -1);
        zoomSliderUpdate();
        grid2D.repaint();
    }

    public void gridPlaneUpdate() {

        controlsOff = true;

        spinnerPlane.setValue(grid2D.kVoxel);

        sliderPlane.setMinimum(0);
        sliderPlane.setMaximum(grid2D.kVoxels() - 1);
        sliderPlane.setValue(grid2D.kVoxel);

        controlsOff = false;

    }

    public void gridPlaneSet(int i) {

        if (controlsOff) {
            return;
        }

        int max = grid2D.kVoxels() - 1;

        if (i < 0) {
            i = 0;
        }

        if (i > max) {
            i = max;
        }

        spinnerPlane.setValue(i);
        sliderPlane.setValue(i);
        grid2D.setGridPlane(i);

        

    }

    public void axesModeUpdate() {

        String xyz = "";

        switch(grid2D.axesMode) {
            case 0:
                if (grid2D.ijSwitch) {
                    xyz = "yx";
                } else {
                    xyz = "xy";
                }
                break;
            case 1:
                if (grid2D.ijSwitch) {
                    xyz = "zy";
                } else {
                    xyz = "yz";
                }
                break;
            case 2:
                if (grid2D.ijSwitch) {
                    xyz = "xz";
                } else {
                    xyz = "zx";
                }
                break;
        }

        radioAxisMode.setText(xyz + " plane");
        radioAxisMode.setSelected(true);

    }

    public void axesModeSet(int m) {

        if (controlsOff) {
            return;
        }

        switch (m) {
            case 0:
            case 1:
            case 2:
                break;
            default:
                m = 0;
        }

        grid2D.setAxesMode(m);
        axesModeUpdate();
        gridPlaneUpdate();

    }

    public void updateEditButtons() {

        if (grid2D.displayMode == 0) {
            buttonCopy.setEnabled(true);
            buttonPaste.setEnabled(true);
            buttonClear.setEnabled(true);
        } else {
            buttonCopy.setEnabled(false);
            buttonPaste.setEnabled(false);
            buttonClear.setEnabled(false);
        }

    }

    public void coordinatesSet(String xc, String yc, String zc, String vc) {

        if (controlsOff) {
            return;
        }

        xCoord.setText(xc);
        yCoord.setText(yc);
        zCoord.setText(zc);
        vCoord.setText(vc);

    }

    public void coordinatesDimUpdate() {

        if (grid2D.coordVoxels) {
            coordDim.setText("voxels");
        } else {
            coordDim.setText(Master.project.spaceUnits);
        }

        coordDim.setSelected(true);

    }

    public void coordinatesDimToggle() {

        if (controlsOff) {
            return;
        }
        
        grid2D.coordVoxels = !grid2D.coordVoxels;

        coordinatesDimUpdate();

    }

    public void coordinatesClickUpdate() {

        if (grid2D.coordClick) {
            coordClick.setText("mouse click");
        } else {
            coordClick.setText("mouse move");
        }

        coordClick.setSelected(true);

    }

    public void coordinatesClickToggle() {

        if (controlsOff) {
            return;
        }
        
        grid2D.coordClick = !grid2D.coordClick;

        coordinatesClickUpdate();

    }

    void spinnerPlane_stateChanged(ChangeEvent e) {
        Integer v = (Integer) spinnerPlane.getValue();
        gridPlaneSet(v);
    }

    void sliderPlane_stateChanged(ChangeEvent e) {
        Integer v = (Integer) sliderPlane.getValue();
        gridPlaneSet(v);
    }

    void radioAxisMode_actionPerformed(ActionEvent e) {
        axesModeSet(grid2D.axesMode+1);
    }

    void buttonSwitch_mouseClicked(MouseEvent e) {
        grid2D.switchAxesToggle();
        axesModeUpdate();
    }

    void comboDisplayMode_itemStateChanged(ItemEvent e) {
        JComboBox cb = (JComboBox) e.getSource();
        String item = (String) cb.getSelectedItem();
        displayModeSet(item);
    }

    void buttonColor_mouseClicked(MouseEvent e) {
        objectColorSet();
    }

    void buttonCopy_mouseClicked(MouseEvent e) {
        if (grid2D.displayMode == 0) {
            grid2D.copyPlane();
        }
    }

    void buttonPaste_mouseClicked(MouseEvent e) {
        if (grid2D.displayMode == 0) {
            grid2D.pastePlane();
        }
    }

    void buttonClear_mouseClicked(MouseEvent e) {
        if (grid2D.displayMode == 0) {
            grid2D.clearPlane();
        }
    }

    void radioEditMode_actionPerformed(ActionEvent e) {
        editModeSet(grid2D.editMode+1);
    }

    void coordDim_actionPerformed(ActionEvent e) {
        coordinatesDimToggle();
    }

    void coordClick_actionPerformed(ActionEvent e) {
        coordinatesClickToggle();
    }

    void this_componentResized(ComponentEvent e) {
        grid2D.revalidate();
        jScrollPaneGrid.setPreferredSize(new Dimension(jPanelGrid.getWidth(), jPanelGrid.getHeight()));
        jScrollPaneGrid.revalidate();
    }

    void sliderZoom_stateChanged(ChangeEvent e) {
        Integer v = (Integer) sliderZoom.getValue();
        grid2D.setVoxelWidth(v);
    }

    void crossLinesToggle() {
        grid2D.setCross(!grid2D.cross);
        cross.setSelected(grid2D.cross);
    }

    void voxelGridToggle() {
        grid2D.setVoxelGrid(!grid2D.voxelGrid);
        voxelGrid.setSelected(grid2D.voxelGrid);
    }

    void sym4toggle() {
        grid2D.setSym4(!grid2D.sym4);
        symmetry4.setSelected(grid2D.sym4);
    }
}

class Panel2D_comboDisplayMode_itemAdapter
        implements java.awt.event.ItemListener {

    Panel2D adaptee;

    Panel2D_comboDisplayMode_itemAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void itemStateChanged(ItemEvent e) {
        if (e.getStateChange() == 2) {
            adaptee.comboDisplayMode_itemStateChanged(e);
        }
    }
}

class Panel2D_radioEditMode_actionAdapter
        implements java.awt.event.ActionListener {

    Panel2D adaptee;

    Panel2D_radioEditMode_actionAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.radioEditMode_actionPerformed(e);
    }
}

class Panel2D_buttonColor_mouseAdapter
        extends java.awt.event.MouseAdapter {

    Panel2D adaptee;

    Panel2D_buttonColor_mouseAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.buttonColor_mouseClicked(e);
    }
}

class Panel2D_spinnerPlane_changeAdapter
        implements javax.swing.event.ChangeListener {

    Panel2D adaptee;

    Panel2D_spinnerPlane_changeAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        adaptee.spinnerPlane_stateChanged(e);
    }
}

class Panel2D_sliderPlane_changeAdapter
        implements javax.swing.event.ChangeListener {

    Panel2D adaptee;

    Panel2D_sliderPlane_changeAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        adaptee.sliderPlane_stateChanged(e);
    }
}

class Panel2D_radioAxisMode_actionAdapter
        implements java.awt.event.ActionListener {

    Panel2D adaptee;

    Panel2D_radioAxisMode_actionAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.radioAxisMode_actionPerformed(e);
    }
}

class Panel2D_buttonSwitch_mouseAdapter
        extends java.awt.event.MouseAdapter {

    Panel2D adaptee;

    Panel2D_buttonSwitch_mouseAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.buttonSwitch_mouseClicked(e);
    }
}

class Panel2D_buttonCopy_mouseAdapter
        extends java.awt.event.MouseAdapter {

    Panel2D adaptee;

    Panel2D_buttonCopy_mouseAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.buttonCopy_mouseClicked(e);
    }
}

class Panel2D_buttonPaste_mouseAdapter
        extends java.awt.event.MouseAdapter {

    Panel2D adaptee;

    Panel2D_buttonPaste_mouseAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.buttonPaste_mouseClicked(e);
    }
}

class Panel2D_buttonClear_mouseAdapter
        extends java.awt.event.MouseAdapter {

    Panel2D adaptee;

    Panel2D_buttonClear_mouseAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.buttonClear_mouseClicked(e);
    }
}

class Panel2D_sliderZoom_changeAdapter
        implements javax.swing.event.ChangeListener {

    Panel2D adaptee;

    Panel2D_sliderZoom_changeAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        adaptee.sliderZoom_stateChanged(e);
    }
}

class Panel2D_cross_actionAdapter
        implements java.awt.event.ActionListener {

    Panel2D adaptee;

    Panel2D_cross_actionAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.crossLinesToggle();
    }
}

class Panel2D_voxelGrid_actionAdapter
        implements java.awt.event.ActionListener {

    Panel2D adaptee;

    Panel2D_voxelGrid_actionAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.voxelGridToggle();
    }
}

class Panel2D_symmetry4_actionAdapter
        implements java.awt.event.ActionListener {

    Panel2D adaptee;

    Panel2D_symmetry4_actionAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.sym4toggle();
    }
}

class Panel2D_coordDim_actionAdapter
        implements java.awt.event.ActionListener {

    Panel2D adaptee;

    Panel2D_coordDim_actionAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.coordDim_actionPerformed(e);
    }
}

class Panel2D_coordClick_actionAdapter
        implements java.awt.event.ActionListener {

    Panel2D adaptee;

    Panel2D_coordClick_actionAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.coordClick_actionPerformed(e);
    }
}

class Panel2D_this_componentAdapter
        extends java.awt.event.ComponentAdapter {

    Panel2D adaptee;

    Panel2D_this_componentAdapter(Panel2D adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void componentResized(ComponentEvent e) {
        adaptee.this_componentResized(e);
    }
}
