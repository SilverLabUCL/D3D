package ucl.silver.d3d.gui;

import ucl.silver.d3d.core.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
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
public class MainFrame
        extends JFrame {

    JPanel contentPane;

    JFileChooser fc = new JFileChooser();

    JMenuBar jManuMain = new JMenuBar();

    JMenu jMenuD3D = new JMenu();
    JMenuItem jMenuD3DAbout = new JMenuItem();
    JMenuItem jMenuD3DExit = new JMenuItem();

    JMenu jMenuProject = new JMenu();
    JMenuItem jMenuProjectNew = new JMenuItem();
    JMenuItem jMenuProjectOpen = new JMenuItem();
    JMenuItem jMenuProjectSave = new JMenuItem();
    JMenuItem jMenuProjectInit = new JMenuItem();
    JMenuItem jMenuProjectDivider = new JMenuItem();

    JMenu jMenuGeometry = new JMenu();
    JMenuItem jMenuGeometryClear = new JMenuItem();
    JMenuItem jMenuGeometryEllipsoid = new JMenuItem();
    JMenuItem jMenuGeometryCylinder = new JMenuItem();
    JMenuItem jMenuGeometryBouton = new JMenuItem();
    JMenuItem jMenuGeometryGlomerulus = new JMenuItem();
    JMenuItem jMenuGeometryOpen = new JMenuItem();
    JMenuItem jMenuGeometrySave = new JMenuItem();
    JMenuItem jMenuGeometryExport = new JMenuItem();
    JMenuItem jMenuGeometryImport = new JMenuItem();

    JMenu jMenuDiffusant = new JMenu();
    JMenu jMenuAddDiffusant = new JMenu();
    JMenuItem jMenuAddSimple = new JMenuItem();
    JMenuItem jMenuaddBuffer = new JMenuItem();
    JMenuItem jMenuAddPhotolysis = new JMenuItem();
    JMenuItem jMenuKillDiffusant = new JMenuItem();

    JMenu jMenuSource = new JMenu();
    JMenu jMenuAddSource = new JMenu();
    JMenuItem jMenuAddSourceQuanta = new JMenuItem();
    JMenuItem jMenuAddSourceGauss = new JMenuItem();
    JMenuItem jMenuKillSource = new JMenuItem();

    JMenu jMenuDetector = new JMenu();
    JMenu jMenuAddDetector = new JMenu();
    JMenuItem jMenuAddDetectorAvg = new JMenuItem();
    JMenuItem jMenuAddDetectorWeighted = new JMenuItem();
    JMenuItem jMenuAddDetectorSnap = new JMenuItem();
    JMenuItem jMenuKillDetector = new JMenuItem();

    JMenu jMenuPSF = new JMenu("PSF");
    JMenuItem jMenuPSFopen = new JMenuItem();
    JMenuItem jMenuPSFsave = new JMenuItem();
    JMenuItem jMenuPSFclose = new JMenuItem();
    JMenuItem jMenuPSFexport = new JMenuItem();

    JMenu jMenuHelp = new JMenu();
    JMenuItem jMenuHelpTopics = new JMenuItem();

    JTabbedPane tabbedPane = new JTabbedPane();
    JPanel tabParam = new JPanel();
    JPanel tab2Dview = new JPanel();
    JPanel tabLog = new JPanel();
    JPanel tabOuts = new JPanel();
    //JPanel tab3Dview = new JPanel();

    public PanelParams panelParams;
    public Panel2D panel2D = null;
    //public Panel3D panel3D;
    public PanelLog panelLog;
    public PanelOuts panelOuts;
    public JPanel panelCompute = new JPanel();

    JButton jButtonComputeJava = new JButton();
    JButton jButtonComputeNative = new JButton();
    JButton jButtonPause = new JButton();
    JButton jButtonPreview = new JButton();
    JButton jButtonCancel = new JButton();
    JButton jButtonInit = new JButton();

    FlowLayout flowLayoutCompute = new FlowLayout();
    BorderLayout borderLayoutParam = new BorderLayout();
    BorderLayout borderLayout2Dview = new BorderLayout();
    BorderLayout borderLayoutLog = new BorderLayout();
    BorderLayout borderLayoutOuts = new BorderLayout();
    BorderLayout borderLayout3Dview = new BorderLayout();
    BorderLayout borderLayoutContent = new BorderLayout();
    FlowLayout flowLayout1 = new FlowLayout();
    
    private UserInput user = new UserInput();

    // construct the frame
    public MainFrame() {

        enableEvents(AWTEvent.WINDOW_EVENT_MASK);

        try {
            jbInit();
            validate();
        } catch (Exception e) {
            //e.printStackTrace();
        }

    }

    // Component initialization
    public final void jbInit() throws Exception {

        contentPane = (JPanel) this.getContentPane();
        contentPane.setLayout(borderLayoutContent);
        contentPane.setBorder(BorderFactory.createLoweredBevelBorder());
        setSize(new Dimension(1000, 800));
        //this.setTitle("D3D");

        jMenuD3D.setText("D3D");

        jMenuD3DAbout.setText("About");
        jMenuD3DAbout.addActionListener(new MainFrame_jMenuD3DAbout_ActionAdapter(this));
        jMenuD3DExit.setText("Quit");
        jMenuD3DExit.addActionListener(new MainFrame_jMenuD3DExit_ActionAdapter(this));

        jMenuProject.setText("Project");

        jMenuProjectNew.setText("New");
        jMenuProjectNew.addActionListener(new MainFrame_jMenuProject_actionAdapter(this));
        jMenuProjectOpen.setText("Open");
        jMenuProjectOpen.addActionListener(new MainFrame_jMenuProject_actionAdapter(this));
        jMenuProjectSave.setText("Save");
        jMenuProjectSave.addActionListener(new MainFrame_jMenuProject_actionAdapter(this));
        jMenuProjectInit.setText("Init");
        jMenuProjectInit.addActionListener(new MainFrame_jMenuProject_actionAdapter(this));
        jMenuProjectDivider.setText("---");

        jMenuGeometry.setText("Geometry");

        jMenuGeometryClear.setText("Clear");
        jMenuGeometryClear.addActionListener(new MainFrame_jMenuGeometryClear_actionAdapter(this));

        jMenuGeometryEllipsoid.setText("Ellipsoid");
        jMenuGeometryEllipsoid.addActionListener(new MainFrame_jMenuGeometryEllipsoid_actionAdapter(this));

        jMenuGeometryCylinder.setText("Cylinder");
        jMenuGeometryCylinder.addActionListener(new MainFrame_jMenuGeometryCylinder_actionAdapter(this));

        jMenuGeometryBouton.setText("Bouton");
        jMenuGeometryBouton.addActionListener(new MainFrame_jMenuGeometryBouton_actionAdapter(this));

        jMenuGeometryGlomerulus.setText("Glomerulus");
        jMenuGeometryGlomerulus.addActionListener(new MainFrame_jMenuGeometryGlomerulus_actionAdapter(this));

        jMenuGeometryOpen.setText("Open");
        jMenuGeometryOpen.addActionListener(new MainFrame_jMenuGeometryOpen_actionAdapter(this));

        jMenuGeometrySave.setText("Save");
        jMenuGeometrySave.addActionListener(new MainFrame_jMenuGeometrySave_actionAdapter(this));

        jMenuGeometryImport.setText("Import");
        jMenuGeometryImport.addActionListener(new MainFrame_jMenuGeometryImport_actionAdapter(this));

        jMenuGeometryExport.setText("Export");
        jMenuGeometryExport.addActionListener(new MainFrame_jMenuGeometryExport_actionAdapter(this));

        jMenuDiffusant.setText("Diffusant");
        jMenuAddDiffusant.setText("Add");
        jMenuAddSimple.setText("Simple");
        jMenuAddSimple.addActionListener(new MainFrame_jMenuAddSimple_actionAdapter(this));

        jMenuAddPhotolysis.setText("Photolysis");
        jMenuAddPhotolysis.addActionListener(new MainFrame_jMenuAddPhotolysis_actionAdapter(this));

        jMenuaddBuffer.setText("Buffer");
        jMenuaddBuffer.addActionListener(new MainFrame_jMenuaddBuffer_actionAdapter(this));
        jMenuKillDiffusant.setText("Kill");
        jMenuKillDiffusant.addActionListener(new MainFrame_jMenuKillDiffusant_actionAdapter(this));

        jMenuSource.setText("Source");
        jMenuAddSource.setText("Add");
        jMenuAddSourceQuanta.setText("Quanta");
        jMenuAddSourceQuanta.addActionListener(new MainFrame_jMenuAddSourceQuanta_actionAdapter(this));
        jMenuAddSourceGauss.setText("Gaussian");
        jMenuAddSourceGauss.addActionListener(new MainFrame_jMenuAddSourceGauss_actionAdapter(this));
        jMenuKillSource.setText("Kill");
        jMenuKillSource.addActionListener(new MainFrame_jMenuKillSource_actionAdapter(this));

        jMenuDetector.setText("Detector");
        jMenuAddDetector.setText("Add");
        jMenuAddDetectorAvg.setText("Average");
        jMenuAddDetectorAvg.addActionListener(new MainFrame_jMenuAddDetectorAvg_actionAdapter(this));
        jMenuAddDetectorWeighted.setText("Weighted");
        jMenuAddDetectorWeighted.addActionListener(new MainFrame_jMenuAddDetectorWeighted_actionAdapter(this));
        jMenuAddDetectorSnap.setText("Snapshot");
        jMenuAddDetectorSnap.addActionListener(new MainFrame_jMenuAddDetectorSnap_actionAdapter(this));
        jMenuKillDetector.setText("Kill");
        jMenuKillDetector.addActionListener(new MainFrame_jMenuKillDetector_actionAdapter(this));

        jMenuPSFopen.setText("Open");
        jMenuPSFopen.addActionListener(new MainFrame_jMenuPSFopen_actionAdapter(this));

        jMenuPSFsave.setText("Save");
        jMenuPSFsave.addActionListener(new MainFrame_jMenuPSFsave_actionAdapter(this));

        jMenuPSFclose.setText("Close");
        jMenuPSFclose.addActionListener(new MainFrame_jMenuPSFclose_actionAdapter(this));

        jMenuPSFexport.setText("Export");
        jMenuPSFexport.addActionListener(new MainFrame_jMenuPSFexport_actionAdapter(this));

        jMenuHelp.setText("Help");
        jMenuHelpTopics.setText("Help Topics");

        jButtonComputeJava.setPreferredSize(new Dimension(100, 25));
        jButtonComputeJava.setActionCommand("Compute");
        jButtonComputeJava.setText("Run");
        jButtonComputeJava.addActionListener(new MainFrame_jButtonComputeJava_actionAdapter(this));

        jButtonComputeNative.setPreferredSize(new Dimension(100, 25));
        jButtonComputeNative.setText("NativeC");
        jButtonComputeNative.addActionListener(new MainFrame_jButtonComputeNative_actionAdapter(this));

        panelParams = new PanelParams();
        panel2D = new Panel2D();
        //panel3D = new Panel3D(proj);
        panelLog = new PanelLog();
        panelOuts = new PanelOuts();

        //panelParams.setMinimumSize(new Dimension(800, 500));
        //panelParams.setPreferredSize(new Dimension(800, 500));
        //panelParams.setSize(new Dimension(800, 500));
        //panelParams.setLayout(flowLayout1);
        //panelParams.validate();

        tabParam.setLayout(borderLayoutParam);

        jButtonPause.setPreferredSize(new Dimension(100, 25));
        jButtonPause.setActionCommand("Pause");
        jButtonPause.setText("Pause");
        jButtonPause.addActionListener(new MainFrame_jButtonPause_actionAdapter(this));

        jButtonPreview.setPreferredSize(new Dimension(100, 25));
        jButtonPreview.setActionCommand("Preview");
        jButtonPreview.setText("Preview");
        jButtonPreview.addActionListener(new MainFrame_jButtonPreview_actionAdapter(this));

        jButtonCancel.setAlignmentY((float) 0.5);
        jButtonCancel.setPreferredSize(new Dimension(100, 25));
        jButtonCancel.setActionCommand("Cancel");
        jButtonCancel.setText("Cancel");
        jButtonCancel.addActionListener(new MainFrame_jButtonCancel_actionAdapter(this));

        jButtonInit.setAlignmentY((float) 0.5);
        jButtonInit.setPreferredSize(new Dimension(100, 25));
        jButtonInit.setActionCommand("Init");
        jButtonInit.setText("Init");
        jButtonInit.addActionListener(new MainFrame_jButtonInit_actionAdapter(this));

        tabParam.add(panelParams, BorderLayout.CENTER);

        tab2Dview.setLayout(borderLayout2Dview);
        tab2Dview.add(panel2D, BorderLayout.CENTER);

        //tab3Dview.setBorder(BorderFactory.createLoweredBevelBorder());
        //tab3Dview.setLayout(borderLayout3Dview);
        //tab3Dview.add(panel3D, BorderLayout.CENTER);

        tabLog.setLayout(borderLayoutLog);
        tabLog.add(panelLog, BorderLayout.CENTER);

        tabOuts.setLayout(borderLayoutOuts);
        tabOuts.add(panelOuts, BorderLayout.CENTER);

        tabbedPane.add(tabParam, "Parameters");
        tabbedPane.add(tab2Dview, "2D view");
        //tabbedPane.add(tab3Dview, "3D view");
        tabbedPane.add(tabLog, "Log");
        //tabbedPane.add(tabOuts, "Outputs");
        tabbedPane.addMouseListener(new MainFrame_tabbedPane_mouseAdapter(this));

        panelCompute.setLayout(flowLayoutCompute);
        panelCompute.add(jButtonInit, null);
        panelCompute.add(jButtonComputeJava, null);
        //panelCompute.add(jButtonComputeNative, null);
        panelCompute.add(jButtonPause, null);
        panelCompute.add(jButtonPreview, null);
        panelCompute.add(jButtonCancel, null);
        contentPane.add(panelCompute, BorderLayout.SOUTH);

        jMenuD3D.add(jMenuD3DAbout);
        jMenuD3D.add(jMenuD3DExit);

        //jMenuProject.add(jMenuProjectNew);
        //jMenuProject.add(jMenuProjectOpen);
        //jMenuProject.add(jMenuProjectSave);
        //jMenuProject.add(jMenuProjectInit);

        jMenuGeometry.add(jMenuGeometryClear);
        jMenuGeometry.add(jMenuGeometryEllipsoid);
        jMenuGeometry.add(jMenuGeometryCylinder);
        jMenuGeometry.add(jMenuGeometryBouton);
        jMenuGeometry.add(jMenuGeometryGlomerulus);
        jMenuGeometry.add(jMenuGeometryOpen);
        jMenuGeometry.add(jMenuGeometrySave);
        jMenuGeometry.add(jMenuGeometryImport);
        jMenuGeometry.add(jMenuGeometryExport);

        jMenuPSF.add(jMenuPSFopen);
        jMenuPSF.add(jMenuPSFsave);
        jMenuPSF.add(jMenuPSFexport);

        jMenuAddDiffusant.add(jMenuAddSimple);
        jMenuAddDiffusant.add(jMenuaddBuffer);
        jMenuAddDiffusant.add(jMenuAddPhotolysis);

        jMenuDiffusant.add(jMenuAddDiffusant);
        jMenuDiffusant.add(jMenuKillDiffusant);

        jMenuAddSource.add(jMenuAddSourceQuanta);
        jMenuAddSource.add(jMenuAddSourceGauss);

        jMenuSource.add(jMenuAddSource);
        jMenuSource.add(jMenuKillSource);

        jMenuAddDetector.add(jMenuAddDetectorAvg);
        jMenuAddDetector.add(jMenuAddDetectorWeighted);
        jMenuAddDetector.add(jMenuAddDetectorSnap);

        jMenuDetector.add(jMenuAddDetector);
        jMenuDetector.add(jMenuKillDetector);

        jMenuHelp.add(jMenuHelpTopics);

        jManuMain.add(jMenuD3D);
        jManuMain.add(jMenuProject);
        jManuMain.add(jMenuGeometry);
        jManuMain.add(jMenuDiffusant);
        jManuMain.add(jMenuSource);
        jManuMain.add(jMenuDetector);
        jManuMain.add(jMenuPSF);
        jManuMain.add(jMenuHelp);

        setJMenuBar(jManuMain);
        contentPane.add(tabbedPane, BorderLayout.CENTER);

        initTabs();

    }

    public void updateMenuProject() {

        String projectList[] = Master.projectList();

        int numProjects = projectList.length;

        jMenuProject.removeAll();

        jMenuProject.add(jMenuProjectNew);
        jMenuProject.add(jMenuProjectOpen);
        jMenuProject.add(jMenuProjectSave);
        jMenuProject.add(jMenuProjectInit);

        if (numProjects > 0) {

            JMenuItem jMenuProjectList[] = new JMenuItem[numProjects];

            jMenuProject.add(jMenuProjectDivider);
            
            for (int i = 0; i < numProjects; i++) {

                jMenuProjectList[i] = new JMenuItem();
                jMenuProjectList[i].setText(projectList[i]);
                jMenuProjectList[i].addActionListener(new MainFrame_jMenuProject_actionAdapter(this));

                if (projectList[i].equalsIgnoreCase(Master.project.name)){
                    jMenuProjectList[i].setForeground(Color.red);
                }

                jMenuProject.add(jMenuProjectList[i]);

            }

        }

    }

    public void initTabs() {

        if (panelParams != null) {
            panelParams.updateControls();
        }

        updateMenuProject();

        repaint();

    }

    public void tabbedPaneClicked(MouseEvent e) {
        
        switch(tabbedPane.getSelectedIndex()) {
            case 0:
                break;
            case 1:
                panel2D.updateControls();
                break;
            case 2:
                break;
        }

    }

    // File | Exit action performed
    public void jMenuD3DExit_actionPerformed(ActionEvent e) {
        System.exit(0);
    }

    // Help | About action performed
    public void jMenuD3DAbout_actionPerformed(ActionEvent e) {
        AboutD3D jDialog = new AboutD3D(this);
        Dimension dlgSize = jDialog.getPreferredSize();
        Dimension frmSize = getSize();
        Point loc = getLocation();
        jDialog.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x,
                (frmSize.height - dlgSize.height) / 2 + loc.y);
        jDialog.setModal(true);
        jDialog.pack();
        jDialog.setVisible(true);
        tabbedPane.setSelectedIndex(2); // change to log tab
    }

    // Overridden so we can exit when window is closed
    @Override
    protected void processWindowEvent(WindowEvent e) {
        super.processWindowEvent(e);
        if (e.getID() == WindowEvent.WINDOW_CLOSING) {
            jMenuD3DExit_actionPerformed(null);
        }
    }

    void jButtonComputeJava_actionPerformed(ActionEvent e) {
        Master.project.simulationStart(false);
    }

    void jButtonComputeNative_actionPerformed(ActionEvent e) {
        // no longer working
    }

    void jButtonPause_actionPerformed(ActionEvent e) {
        Master.project.simulationPause();
    }

    void jButtonPreview_actionPerformed(ActionEvent e) {
        Master.project.simulationStart(true);
    }

    void jButtonCancel_actionPerformed(ActionEvent e) {
        Master.project.simulationCancel();
    }

    void jButtonInit_actionPerformed(ActionEvent e) {

        if (Master.project.monteCarlo != null) {
            Master.project.monteCarlo.initAll();
        }

        Master.project.init();

    }

    void jMenuProject_actionPerformed(ActionEvent e) {

        String pname = "Project0";
        String fn;
        String command = e.getActionCommand();

        if (command.equalsIgnoreCase("New")) {

            pname = (String) JOptionPane.showInputDialog(
                    "Please enter new project name:", pname);

            if (pname != null) {
                Master.createProject("FD", pname, null);
                initTabs();
            }

        } else if (command.equalsIgnoreCase("Open")) {

            fn = openFileName();

            if (fn.length() > 0) {
                Master.projectOpen(fn);
                initTabs();
            }

        } else if (command.equalsIgnoreCase("Save")) {

            fn = saveFileName(Master.project.name + ".d3d");

            if (fn.length() > 0) {
                Master.projectSave(fn);
            }

        } else if (command.equalsIgnoreCase("Init")) {
            Master.promptInitProject();
        } else {
            Master.projectSelect(command);
            initTabs();
        }

    }

    void jMenuGeometryClear_actionPerformed(ActionEvent e) {
        Master.project.geometry.clear();
    }

    void jMenuGeometryEllipsoid_actionPerformed(ActionEvent e) {
        Master.project.geometry.ellipsoid();
        Master.project.geometry.name = "ellipsoid";
    }

    void jMenuGeometryCylinder_actionPerformed(ActionEvent e) {
        Master.project.geometry.cylinder();
        Master.project.geometry.name = "cylinder";
    }

    void jMenuGeometryBouton_actionPerformed(ActionEvent e) {
        user.createBoutonPrompt();
    }

    void jMenuGeometryGlomerulus_actionPerformed(ActionEvent e) {
        user.createGlomerulusPrompt();
    }

    void jMenuGeometrySave_actionPerformed(ActionEvent e) {

        String fileName = saveFileName(Master.project.name + "_Geometry.d3d");

        if (fileName.length() > 0) {
            Master.saveGeometry(fileName);
        }

    }

    void jMenuGeometryOpen_actionPerformed(ActionEvent e) {

        String fileName = openFileName();

        if (fileName.length() > 0) {
            Master.openGeometry(fileName);
            initTabs();
        }

    }

    void jMenuGeometryExport_actionPerformed(ActionEvent e) {

        String fileName = saveFileName(Master.project.name + "_Geometry.txt");

        if (fileName.length() > 0) {
            Master.project.geometry.exportAB(fileName);
        }

    }

    void jMenuGeometryImport_actionPerformed(ActionEvent e) {

        String fileName = openFileName();

        boolean fileContainsNonSpaceBorder = false;

        if (fileName.length() > 0) {
            Master.project.geometry.importAB(fileName, fileContainsNonSpaceBorder);
            initTabs();
        }

    }

    void jMenuAddSimple_actionPerformed(ActionEvent e) {
        user.addDiffusantPrompt();
    }

    void jMenuaddBuffer_actionPerformed(ActionEvent e) {
        user.addDiffusantBufferPrompt();
    }

    void jMenuKillDiffusant_actionPerformed(ActionEvent e) {
        user.killDiffusantPrompt(0);
    }

    void jMenuAddSourceQuanta_actionPerformed(ActionEvent e) {
        user.addSourceImpulsePrompt();
    }

    void jMenuAddSourceGauss_actionPerformed(ActionEvent e) {
        user.addSourceGaussPrompt();
    }

    void jMenuKillSource_actionPerformed(ActionEvent e) {
        user.killSourcePrompt(0);
    }

    void jMenuAddDetectorAvg_actionPerformed(ActionEvent e) {
        user.addDetectorPrompt();
    }

    void jMenuAddDetectorWeighted_actionPerformed(ActionEvent e) {
        user.addDetectorWeightedPrompt();
    }

    void jMenuAddDetectorSnap_actionPerformed(ActionEvent e) {
        user.addDetectorSnapPrompt();
    }

    void jMenuKillDetector_actionPerformed(ActionEvent e) {
        user.killDetectorPrompt(0);
    }

    void jMenuAddPhotolysis_actionPerformed(ActionEvent e) {
        user.addDiffusantPhotolysisPrompt();
    }

    void jMenuPSFopen_actionPerformed(ActionEvent e) {
        String fn = openFileName();
        if (fn.length() > 0) {
            //Master.openPSF(fn);
            //panel2D.repaint();
        }
    }

    void jMenuPSFsave_actionPerformed(ActionEvent e) {
        String fn = saveFileName(Master.project.name + "_PSF.d3d");
        if (fn.length() > 0) {
            //user.savePSF(fn);
        }
    }

    void jMenuPSFclose_actionPerformed(ActionEvent e) {
    }

    void jMenuPSFexport_actionPerformed(ActionEvent e) {
        String fn = saveFileName(Master.project.name + "_PSF.d3d");
        if (fn.length() > 0) {
            user.writePSF(fn);
        }
    }

    // get file name
    String openFileName() {

        fc.setSelectedFile(new File(""));
        int approve = fc.showOpenDialog(this);
        String fileName = "";

        if (approve == JFileChooser.APPROVE_OPTION) {
            fileName = fc.getCurrentDirectory() + "\\" + fc.getSelectedFile().getName();
        }

        return fileName;

    }

    // get file name
    String saveFileName(String fn) {
        
        fc.setSelectedFile(new File(fn));
        int approve = fc.showSaveDialog(this);
        String fileName = "";

        if (approve == JFileChooser.APPROVE_OPTION) {
            fileName = fc.getCurrentDirectory() + "\\" + fc.getSelectedFile().getName();
        }

        return fileName;

    }
}

class MainFrame_jMenuD3DExit_ActionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuD3DExit_ActionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuD3DExit_actionPerformed(e);
    }
}

class MainFrame_jMenuD3DAbout_ActionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuD3DAbout_ActionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuD3DAbout_actionPerformed(e);
    }
}

class MainFrame_jButtonComputeJava_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonComputeJava_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonComputeJava_actionPerformed(e);
    }
}

class MainFrame_jButtonComputeNative_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonComputeNative_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonComputeNative_actionPerformed(e);
    }
}

class MainFrame_jMenuProject_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuProject_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuProject_actionPerformed(e);
    }
}

class MainFrame_jMenuGeometryClear_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryClear_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryClear_actionPerformed(e);
    }
}

class MainFrame_jMenuGeometryEllipsoid_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryEllipsoid_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryEllipsoid_actionPerformed(e);
    }
}

class MainFrame_jMenuGeometryCylinder_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryCylinder_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryCylinder_actionPerformed(e);
    }
}

class MainFrame_jMenuGeometryBouton_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryBouton_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryBouton_actionPerformed(e);
    }
}

class MainFrame_jMenuGeometryGlomerulus_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryGlomerulus_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryGlomerulus_actionPerformed(e);
    }
}

class MainFrame_jMenuGeometrySave_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometrySave_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometrySave_actionPerformed(e);
    }

}

class MainFrame_jMenuGeometryOpen_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryOpen_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryOpen_actionPerformed(e);
    }
}

class MainFrame_jMenuGeometryExport_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryExport_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryExport_actionPerformed(e);
    }

}

class MainFrame_jMenuGeometryImport_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuGeometryImport_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuGeometryImport_actionPerformed(e);
    }
}

class MainFrame_jMenuAddSimple_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuAddSimple_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuAddSimple_actionPerformed(e);
    }
}

class MainFrame_jMenuaddBuffer_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuaddBuffer_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuaddBuffer_actionPerformed(e);
    }
}

class MainFrame_jMenuKillDiffusant_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuKillDiffusant_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuKillDiffusant_actionPerformed(e);
    }
}

class MainFrame_jMenuAddSourceQuanta_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuAddSourceQuanta_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuAddSourceQuanta_actionPerformed(e);
    }
}

class MainFrame_jMenuAddSourceGauss_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuAddSourceGauss_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuAddSourceGauss_actionPerformed(e);
    }
}

class MainFrame_jMenuKillSource_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuKillSource_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuKillSource_actionPerformed(e);
    }
}

class MainFrame_jMenuAddDetectorAvg_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuAddDetectorAvg_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuAddDetectorAvg_actionPerformed(e);
    }
}

class MainFrame_jMenuAddDetectorWeighted_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuAddDetectorWeighted_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuAddDetectorWeighted_actionPerformed(e);
    }
}

class MainFrame_jMenuAddDetectorSnap_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuAddDetectorSnap_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuAddDetectorSnap_actionPerformed(e);
    }
}

class MainFrame_jMenuKillDetector_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuKillDetector_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuKillDetector_actionPerformed(e);
    }
}

class MainFrame_jMenuAddPhotolysis_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuAddPhotolysis_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuAddPhotolysis_actionPerformed(e);
    }
}

class MainFrame_jMenuPSFopen_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuPSFopen_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuPSFopen_actionPerformed(e);
    }
}

class MainFrame_jMenuPSFsave_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuPSFsave_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuPSFsave_actionPerformed(e);
    }
}

class MainFrame_jMenuPSFclose_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuPSFclose_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuPSFclose_actionPerformed(e);
    }
}

class MainFrame_jMenuPSFexport_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuPSFexport_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuPSFexport_actionPerformed(e);
    }
}

class MainFrame_tabbedPane_mouseAdapter
        extends java.awt.event.MouseAdapter {

    MainFrame adaptee;

    MainFrame_tabbedPane_mouseAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.tabbedPaneClicked(e);
    }

}

class MainFrame_jButtonPause_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonPause_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonPause_actionPerformed(e);
    }
}

class MainFrame_jButtonPreview_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonPreview_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonPreview_actionPerformed(e);
    }
}

class MainFrame_jButtonCancel_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonCancel_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonCancel_actionPerformed(e);
    }
}

class MainFrame_jButtonInit_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonInit_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonInit_actionPerformed(e);
    }
}
