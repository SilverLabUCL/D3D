package ucl.silver.d3d.gui;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import java.awt.*;
import java.awt.event.*;
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
public class PanelParams
        extends JPanel {

    public PanelParamsTable paramTable = new PanelParamsTable();
    JComboBox<String> paramSelect = new JComboBox<String>();
    JSpinner paramArray = new JSpinner();
    JButton newArrayItem = new JButton();
    JButton killArrayItem = new JButton();
    JPanel jPanelLeft = new JPanel();
    JPanel jPanelLeftCenter = new JPanel();
    JPanel jPanelTable = new JPanel();
    JPanel jPanelSpacer1 = new JPanel();
    JScrollPane scrollPane;
    
    

    private int paramSelectNum = 1; // (0) InitProject (1) Project (2) Geometry (3) Diffusant (4) Source (5) Detector (6) Batch (7) Error

    private int paramArrayNum = 0; // display array number
    private int diffusantNum = 0;
    private int sourceNum = 0;
    private int detectorNum = 0;
    private int batchNum = 0;
    private int errorNum = 0;

    public PanelParams() {

        paramSelect.addItemListener(new PanelParams_ParamSelect_itemAdapter(this));
        paramSelect.setPreferredSize(new Dimension(85, 21));
        paramSelect.setToolTipText("Select what to edit");
        paramSelect.addItem("InitProject");
        paramSelect.addItem("Project");
        paramSelect.addItem("Geometry");
        paramSelect.addItem("Diffusant");
        paramSelect.addItem("Source");
        paramSelect.addItem("Detector");
        paramSelect.addItem("Batch");
        paramSelect.addItem("Error");

        paramArray.setPreferredSize(new Dimension(85, 25));
        paramArray.addChangeListener(new PanelParams_objectArrayNum_changeAdapter(this));

        newArrayItem.setPreferredSize(new Dimension(60, 25));
        newArrayItem.setToolTipText("Create a new object.");
        newArrayItem.setText("New");
        newArrayItem.addMouseListener(new PanelParams_buttonNew_mouseAdapter(this));

        killArrayItem.setPreferredSize(new Dimension(60, 25));
        killArrayItem.setToolTipText("Kill selected object.");
        killArrayItem.setText("Kill");
        killArrayItem.addMouseListener(new PanelParams_buttonKill_mouseAdapter(this));

        jPanelLeft.setMinimumSize(new Dimension(120, 250));
        jPanelLeft.setPreferredSize(new Dimension(120, 250));

        jPanelLeftCenter.setMinimumSize(new Dimension(100, 250));
        jPanelLeftCenter.setPreferredSize(new Dimension(100, 250));
        jPanelLeftCenter.setLayout(new FlowLayout());

        //jPanelTable.setLayout(new BorderLayout());
        //jPanelTable.setMinimumSize(new Dimension(paramTable.width(), 2000));
        //jPanelTable.setPreferredSize(new Dimension(paramTable.width(), 2000));
        //jPanelTable.add(paramTable, BorderLayout.CENTER);
        //jPanelTable.add(paramTable.getTableHeader(), BorderLayout.NORTH);

        //scrollPane = new JScrollPane(jPanelTable);
        scrollPane = new JScrollPane(paramTable);

        jPanelSpacer1.setMinimumSize(new Dimension(100, 200));
        jPanelSpacer1.setPreferredSize(new Dimension(100, 200));

        setLayout(new BorderLayout());

        jPanelLeftCenter.add(paramSelect);
        jPanelLeftCenter.add(paramArray);
        jPanelLeftCenter.add(newArrayItem);
        jPanelLeftCenter.add(killArrayItem);

        jPanelLeft.add(jPanelSpacer1);
        jPanelLeft.add(jPanelLeftCenter);

        add(scrollPane, BorderLayout.CENTER);
        add(jPanelLeft, BorderLayout.WEST);

        setParamSelect(-1, -1);

    }

    public int paramArrayNum() {
        Integer arrayNum = (Integer) paramArray.getValue();
        return arrayNum;
    }

    public int paramVectorLength() {

        switch (paramSelectNum) {
            case 0:
            case 1:
            case 2:
                return 0;
            case 3:
                if (Master.project.diffusants == null) {
                    return 0;
                } else {
                    return Master.project.diffusants.length;
                }
            case 4:
                if (Master.project.sources == null) {
                    return 0;
                } else {
                    return Master.project.sources.length;
                }
            case 5:
                if (Master.project.detectors == null) {
                    return 0;
                } else {
                    return Master.project.detectors.length;
                }
            case 6:
                if (Master.project.batches == null) {
                    return 0;
                } else {
                    return Master.project.batches.length;
                }
            case 7:
                if (Master.project.errors == null) {
                    return 0;
                } else {
                    return Master.project.errors.length;
                }

        }

        return 0;

    }

    public void updateControls() {
        
        if (paramSelectNum != paramSelect.getSelectedIndex()) {
            paramSelect.setSelectedIndex(paramSelectNum);
        }
        
        if (paramArrayNum != paramArrayNum()) {
            paramArray.setValue(new Integer(paramArrayNum));
        }

        paramTable.setFields();
        objectArrayNum_OnOff();
        newKillButtons_OnOff();

    }

    private void objectArrayNum_OnOff() {

        int arrayLength = paramVectorLength();

        switch (paramSelectNum) {
            case 0:
            case 1:
            case 2:
                paramArray.setVisible(false);
                break;
            default:
                paramArray.setVisible(true);
        }

        if (arrayLength == 0) {
            paramArray.setEnabled(false);
        } else {
            paramArray.setEnabled(true);
        }

    }

    private void newKillButtons_OnOff() {

        switch (paramSelectNum) {
            case 0:
            case 1:
            case 2:
                newArrayItem.setVisible(false);
                killArrayItem.setVisible(false);
                break;
            default:
                newArrayItem.setVisible(true);
                killArrayItem.setVisible(true);
        }
    }

    public void setParamSelect(int pSelect, int arrayNum) {

        if (pSelect < 0) {
            pSelect = paramSelectNum;
        } else {
            paramSelectNum = pSelect;
        }

        paramTable.paramVectorSelect = null;

        paramArrayNum = 0;

        switch (pSelect) {
            case 0: // InitProject
                paramTable.paramVectorSelect = Master.project.initProject;
                break;
            default: // Project
                paramTable.paramVectorSelect = Master.project;
                paramSelectNum = 1;
                break;
            case 2: // Geometry
                paramTable.paramVectorSelect = Master.project.geometry;
                break;
            case 3: // Diffusant
                diffusantNum = setParamArraySelect(Master.project.diffusants, arrayNum, diffusantNum);
                break;
            case 4: // Source
                sourceNum = setParamArraySelect(Master.project.sources, arrayNum, sourceNum);
                break;
            case 5: // Detector
                detectorNum = setParamArraySelect(Master.project.detectors, arrayNum, detectorNum);
                break;
            case 6: // Batch
                batchNum = setParamArraySelect(Master.project.batches, arrayNum, batchNum);
                break;
            case 7: // Errors
                errorNum = setParamArraySelect(Master.project.errors, arrayNum, errorNum);
                break;
        }

        updateControls();

    }

    private int setParamArraySelect(ParamVector[] pVector, int arrayNum, int defaultNum) {

        if (pVector == null) {
            return 0;
        }

        if (arrayNum == -1) {
            arrayNum = defaultNum;
        }

        if (arrayNum < 0) {
            arrayNum = 0;
        }

        if (arrayNum >= pVector.length) {
            arrayNum = pVector.length - 1;
        }

        paramArrayNum = arrayNum;
        paramTable.paramVectorSelect = pVector[arrayNum];

        return arrayNum;

    }

    void paramVectorNew() {

        int i = paramSelect.getSelectedIndex();

        switch (i) {
            case 3: // Diffusant
                Master.addDiffusantPrompt();
                break;
            case 4: // Source
                break;
            case 5: // Detector
                break;
            case 6: // Batch
                Master.addBatch(0, "", 0);
                break;
            case 7: // Error
                Master.project.addError("");
                break;
        }

        setParamSelect(-1, -1);

    }

    void paramVectorKill(int pSelect, int arrayNum) {

        String title = "D3D Kill";
        String aStr;

        if (pSelect < 0) {
            pSelect = paramSelectNum;
        }

        if (arrayNum < 0) {
            arrayNum = paramArrayNum;
        }

        aStr = Integer.toString(arrayNum);

        switch (pSelect) {
            case 3: // Diffusant
                if ((Master.project.diffusants == null) || (Master.project.diffusants.length == 0)) {
                    return;
                }
                if (Master.promptYes(title, "Kill Diffusant #" + aStr + "?")) {
                    Master.project.killDiffusant(arrayNum);
                }
                break;
            case 4: // Source
                if ((Master.project.sources == null) || (Master.project.sources.length == 0)) {
                    return;
                }
                if (Master.promptYes(title, "Kill Source #" + aStr + "?")) {
                    Master.project.killSource(arrayNum);
                }
                break;
            case 5: // Detector
                if ((Master.project.detectors == null) || (Master.project.detectors.length == 0)) {
                    return;
                }
                if (Master.promptYes(title, "Kill Detector #" + aStr + "?")) {
                    Master.project.killDetector(arrayNum);
                }
                break;
            case 6: // Batch
                if ((Master.project.batches == null) || (Master.project.batches.length == 0)) {
                    return;
                }
                if (Master.promptYes(title, "Kill Barch #" + aStr + "?")) {
                    Master.project.killBatch(arrayNum);
                }
                break;
            case 7: // Error
                if ((Master.project.errors == null) || (Master.project.errors.length == 0)) {
                    return;
                }
                if (Master.promptYes(title, "Kill Error #" + aStr + "?")) {
                    Master.project.killError(arrayNum);
                }
                break;
        }

        setParamSelect(-1, -1);

    }

    void paramSelect_itemStateChanged(ItemEvent e) {
        int i = paramSelect.getSelectedIndex();
        setParamSelect(i, -1);
    }

    public void paramArrayNum_stateChanged(ChangeEvent e) {
        setParamSelect(-1, paramArrayNum());
    }

    void buttonNew_mouseClicked(MouseEvent e) {
        paramVectorNew();
    }

    void buttonKill_mouseClicked(MouseEvent e) {
       paramVectorKill(-1, -1);
    }

}

class PanelParams_ParamSelect_itemAdapter
        implements java.awt.event.ItemListener {

    PanelParams adaptee;

    PanelParams_ParamSelect_itemAdapter(PanelParams adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void itemStateChanged(ItemEvent e) {
        adaptee.paramSelect_itemStateChanged(e);
    }
}

class PanelParams_objectArrayNum_changeAdapter
        implements javax.swing.event.ChangeListener {

    PanelParams adaptee;

    PanelParams_objectArrayNum_changeAdapter(PanelParams adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void stateChanged(ChangeEvent e) {
        adaptee.paramArrayNum_stateChanged(e);
    }
}

class PanelParams_buttonNew_mouseAdapter
        extends java.awt.event.MouseAdapter {

    PanelParams adaptee;

    PanelParams_buttonNew_mouseAdapter(PanelParams adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.buttonNew_mouseClicked(e);
    }
}

class PanelParams_buttonKill_mouseAdapter
        extends java.awt.event.MouseAdapter {

    PanelParams adaptee;

    PanelParams_buttonKill_mouseAdapter(PanelParams adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.buttonKill_mouseClicked(e);
    }
}

