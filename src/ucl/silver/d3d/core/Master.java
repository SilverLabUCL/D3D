package ucl.silver.d3d.core;

import ucl.silver.d3d.gui.*;
import ucl.silver.d3d.init.*;
import ucl.silver.d3d.utils.*;
import java.awt.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import javax.swing.*;
import java.io.*;

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
 * @version 2.0
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
public final class Master {

    public static final String D3D_VERSION = "2.0";

    public static Date startDate = null;

    public static Project project = null;
    public static MainFrame mainframe = null;

    private static Project[] projects = null;

    public static String[] initProjectList = {"InitProject", "InitFD_Demo", "InitMC_Demo"};
    public static String[] diffusantList = {"Diffusant", "DiffusantCalretinin", "DiffusantCalretininCa", "DiffusantParvalbuminCa", "DiffusantPhoto", "DiffusantReactant", "DiffusantVesicles", "DiffusantVesiclesAZ"};
    public static String[] detectorList = {"Detector", "DetectorAvg", "DetectorDye", "DetectorPSF", "DetectorSnapshot"};
    public static String[] sourceList = {"Source", "SourceGauss", "SourceGamma", "SourceUptake"};
    
    public static String[] psfList = {"PSF", "PSFgauss", "PSFtorok", "PSFwilson"};

    public static String[] colorScaleList = {"BTC", "BTY", "Gray", "LinGray", "Heat", "Optimal", "LinOptimal", "Magenta", "Rainbow"};

    public static ColorScaleBTC colorScaleBTC = null;
    public static ColorScaleBTY colorScaleBTY = null;
    public static ColorScaleGray colorScaleGray = null;
    public static ColorScaleLinGray colorScaleLinGray = null;
    public static ColorScaleHeat colorScaleHeat = null;
    public static ColorScaleOptimal colorScaleOptimal = null;
    public static ColorScaleLinOptimal colorScaleLinOptimal = null;
    public static ColorScaleMagenta colorScaleMagenta = null;
    public static ColorScaleRainbow colorScaleRainbow = null;

    public static BufferedWriter bw1 = null;
    private static BufferedWriter bw2 = null;

    private static String[] logList = null;

    public static boolean writeToFiles = false; // allows writing to external file

    public static boolean foundMainStartUpArguments = false;
    
    //private static String initClassAndFunction = "InitProject.initCube";
    
    //private static String initClassAndFunction = "InitFD_Demo.initCalciumChannel";
    //private static String initClassAndFunction = "InitFD_Demo.initGlutamateSource";
    //private static String initClassAndFunction = "InitFD_Demo.initGlutamateSourceCleft";
    //private static String initClassAndFunction = "InitFD_Demo.initCrankCylinder";
    private static String initClassAndFunction = "InitFD_Demo.initCrankCylinderMovie";
    //private static String initClassAndFunction = "InitFD_Demo.initCrankPlane";
    //private static String initClassAndFunction = "InitFD_Demo.initGlomerulus";
    //private static String initClassAndFunction = "InitFD_Demo.initFrapAxelrod";
    //private static String initClassAndFunction = "InitFD_Demo.init_MFT_Scott_Rusakov";
    
    //private static String initClassAndFunction = "InitMC_Demo.initCichocki";
    //private static String initClassAndFunction = "InitMC_Demo.initFrapMovie";
    //private static String initClassAndFunction = "InitMC_Demo.initFrapMFT";
    //private static String initClassAndFunction = "InitMC_Demo.initActiveZone";

    //private static String initClassAndFunction = "InitFD_Ca_David.initCube";
    //private static String initClassAndFunction = "InitFD_Ca_David.initCubeHighRes";
    //private static String initClassAndFunction = "InitFD_Ca_David.initBouton";

    //private static String initClassAndFunction = "InitFD_Ca_Federico.initBouton";
    
    //private static String initClassAndFunction = "InitFD_Ca_Sandrine.initGiant";
    //private static String initClassAndFunction = "InitFD_Ca_Sandrine.initBouton";
    
    //private static String initClassAndFunction = "InitMC_Frap_Jason.initMonteCarloFrap";
    //private static String initClassAndFunction = "InitMC_Frap_Jason.initFiniteDifferenceFrap";

    //private static String initClassAndFunction = "InitMC_AZ_Jason.initMonteCarloAZ";
    //private static String initClassAndFunction = "InitMC_AZEM_Jason.initMonteCarloAZ";

    //private static String initClassAndFunction = "InitFD_GlutamateUncaging.initLargeSpot";
    //private static String initClassAndFunction = "InitFD_GlutamateUncaging.initGlomerulus";

    private Master(){
        // cannot be instantiated
    }
    
    private static boolean initProjectClass(String initClass) {

        if (initClass.equalsIgnoreCase("InitProject")) {
            project.initProject = new InitProject(project);
        } else if (initClass.equalsIgnoreCase("InitFD_Demo")) {
            project.initProject = new InitFD_Demo(project);
        } else if (initClass.equalsIgnoreCase("InitMC_Demo")) {
            project.initProject = new InitMC_Demo(project);
        } else {
            log("Master.initProject error: failed to find init Class " + initClass);
            return true; // error
        }

        return false;

    }

    public static boolean initStartupArgs(String[][] startupArgs, boolean onlyProject, boolean onlyInitProject) {

        String longName, strValue;

        if ((startupArgs == null) || (startupArgs.length == 0)) {
            return false;
        }

        if (startupArgs[0].length != 2) {
            return true; // array has wrong dimensions
        }

        for (int i = 0; i < startupArgs.length; i++) {
            
            longName = startupArgs[i][0];
            strValue = startupArgs[i][1];

            if ((longName == null) || (strValue == null)) {
                continue;
            }

            if (project.setObjectParam(longName, strValue, onlyProject, onlyInitProject)) {
                log("executed the following StartUp command: " + longName + " = " + strValue);
            } else {
                System.err.println("Fatal Error: the following StartUp command could not be executed: " + longName + " = " + strValue);
                return true;
            }

        }

        return false;

    }

    public static boolean initMainFrame(boolean packFrame) {

        mainframe = new MainFrame();

        if (packFrame) {
            mainframe.pack();
        } else {
            mainframe.validate();
        }

        // center the window
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension frameSize = mainframe.getSize();

        if (frameSize.height > screenSize.height) {
            frameSize.height = screenSize.height;
        }

        if (frameSize.width > screenSize.width) {
            frameSize.width = screenSize.width;
        }

        mainframe.setLocation((screenSize.width - frameSize.width) / 2,
                (screenSize.height - frameSize.height) / 2);

        mainframe.initTabs();

        updateMainFrameTitle();

        return true;

    }

    public static void updatePanel2D() {

        if ((mainframe == null) || (mainframe.panel2D == null)) {
            return;
        }

        mainframe.panel2D.updateControls();

    }

    public static Grid grid() {

        if ((mainframe == null) || (mainframe.panel2D == null)) {
            return null;
        }

        return mainframe.panel2D.grid2D;

    }

    public static void exit(String errorStr) {
        log("encountered fatal error: " + errorStr);
        System.exit(0);
    }

    public static void log(String message) {

        if (message == null) {
            return;
        }

        System.out.println(message);

        if (project != null) {
            writeToLogFile(message, false);
        }

        if (mainframe != null) {
            mainframe.panelLog.append(message, true);
        }

    }

    public void fatalError(String message, boolean exit) {

        if (message != null) {
            System.err.println(message + " (" + toString() + ")");
        }

        if (project != null) {
            writeToLogFile(message, false);
        }

        if (mainframe != null) {
            mainframe.panelLog.append(message, true);
        }

        if (exit) {

            if ((project != null) && (mainframe != null)) {
                mainframe.dispose();
            }

            System.exit(0);

        }

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Log Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean writeToLogFile(String s, boolean newFile) {

        String logFile;
        BufferedWriter bw;
        
        if (!writeToFiles) {
            return writeToLogList(s);
        }

        logFile = project.fullDirectory() + project.logFileName();

        try {

            if (newFile){
                bw = new BufferedWriter(new FileWriter(logFile)); // new file
            } else {
                bw = new BufferedWriter(new FileWriter(logFile, true)); // append
            }

            if (logList != null) {

                for (String log : logList) {

                    if ((log == null) || log.equalsIgnoreCase("")) {
                        break;
                    }

                    bw.write(log); // write saved logs first
                    bw.newLine();

                }

                logList = null;

            }

            bw.write(s);
            bw.newLine();
            bw.close();

        } catch (IOException e) {
            Master.exit("fatal error: unable to write to log file: " + logFile);
        }

        return true;

    }

    public static boolean writeToLogList(String s) {

        if (logList == null) {
            logList = new String[1000];
        }

        for (int i = 0; i < logList.length; i++) {
            if ((logList[i] == null) || logList[i].equalsIgnoreCase("")) {
                logList[i] = s;
                break;
            }
        }

        return true;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Open/Save Object Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean writeObject(String fileName, Object o) {

        FileOutputStream f_out;
        ObjectOutputStream obj_out;

        try {
            f_out = new FileOutputStream(fileName);
            obj_out = new ObjectOutputStream(f_out);
            obj_out.writeObject(o);
            log("saved " + o.getClass() + " to file " + fileName);
            return true;
        } catch (NotSerializableException e) {
            log("ReadWrite Error: NotSerializableException " + o.getClass());
            return false;
        } catch (IOException e) {
            log("ReadWrite Error: IOException " + o.getClass());
            return false;
        }

    }

    public static Object readObject(String fileName) {

        FileInputStream f_in;
        ObjectInputStream obj_in;
        Object o = null;

        try {
            f_in = new FileInputStream(fileName);
            obj_in = new ObjectInputStream(f_in);
            o = obj_in.readObject();
            log("opened file " + fileName + " to class: " + o.getClass());
        } catch (ClassNotFoundException e) {
        } //catch (OptionalDataException e) {}
        catch (IOException e) {
            log("ReadWrite Error: unable to read object from file: " + fileName);
        }

        return o;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Project Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static int addProject(Project newProject) {

        int i = 0;

        if (projects != null) {
            i = projects.length;
        }

        Project[] newArray = new Project[i+1];

        if (i > 0) {
            System.arraycopy(projects, 0, newArray, 0, i);
        }

        newArray[i] = newProject;

        projects = newArray; // replace old array with new one

        Master.log("created Project #" + i + " : " + newProject.name);

        return i;

    }

    public static boolean killProject(int i) {

        int k = 0;
        String name = "";

        if ((i < 0) || (i >= projects.length)) {
            return false;
        }

        if (projects.length == 1) {
            projects = null;
            name = project.name;
            project = null;
            Master.log("killed Project #0: " + name);
            return true;
        }

        Project[] newArray = new Project[projects.length - 1]; // new array

        for (int j = 0; j < projects.length; j++) {

            if (j == i) {
                name = projects[j].name;
                continue;
            }

            newArray[k] = projects[j];
            k++;

        }

        projects = newArray; // replace old array with new one

        if (name.length() > 0) {
            Master.log("killed Project #" + i + ": " + name);
        }

        return true;

    }

    public static String[] projectList() {

        int numProjects = projects.length;

        String projectList[] = new String[numProjects];

        for (int i = 0; i < numProjects; i++) {
            if (projects[i] != null) {
                projectList[i] = projects[i].name;
            }
        }

        return projectList;

    }

    public static boolean projectSelect(String projectName) {

        if (projects == null) {
            return false;
        }

        if (project.name.equalsIgnoreCase(projectName)) {
            return false; // already selected
        }

        for (int i = 0; i < projects.length; i++) {
            if (projects[i] == null) {
                continue;
            }
            if (projects[i].name.equalsIgnoreCase(projectName)) {
                project = projects[i];
                mainframe.panelParams.setParamSelect(-1, 0);
                Master.log("changed to Project #" + i + " : " + projectName);
                return true;
            }
        }

        return false;

    }

    public static boolean projectSelect(int projectNum) {

        if ((projects == null) || (projectNum < 0) || (projectNum >= projects.length)) {
            return false;
        }

        if (projects[projectNum] == null) {
            return false;
        }

        project = projects[projectNum];

        mainframe.panelParams.setParamSelect(-1, 0);

        return true;

    }
    
    public static void createProject(String projectType, String projectName, String[][] startupArgsProject) {

        project = new Project(projectType);
        
        if ((projectName != null) && (projectName.length() > 0)) {
            project.name = projectName;
        }

        if (startupArgsProject != null) {
            initStartupArgs(startupArgsProject, true, false);
        }

        addProject(project);

        if (mainframe != null) {
            mainframe.panelParams.setParamSelect(-1, 0);
            mainframe.initTabs();
            //updateMainFrameTitle();
        }

    }

    public static void projectOpen(String fn) {

        Object o = readObject(fn);

        if (o instanceof Project) {

            project = (Project) o;
            project.file = fn;
            addProject(project);

            resetAllParamVectors();

            if (mainframe != null) {
                mainframe.panelParams.setParamSelect(-1, 0);
                mainframe.initTabs();
            }

            updateMainFrameTitle();

            log("opened project: " + fn);

        } else {
            log("open project error: " + fn + " is not of appropriate format.");
        }

    }

    public static boolean projectSave(String fn) {
        return writeObject(fn, project);
    }

    public static void updateMainFrameTitle() {
        if (mainframe != null) {
            //mainframe.setTitle("D3D " + project.name);
            mainframe.setTitle(project.name);
        }
    }

    public static String initProjectClassStr(String initClassAndFunction) {

        if (initClassAndFunction == null) {
            return null;
        }

        String[] splitStr = initClassAndFunction.split("\\.");

        if (splitStr.length == 2) {
            return splitStr[0];
        }

        return null;

    }

    public static String initProjectFunctionStr(String initClassAndFunction) {

        if (initClassAndFunction == null) {
            return null;
        }

        String[] splitStr = initClassAndFunction.split("\\.");

        if (splitStr.length == 2) {
            return splitStr[1];
        }

        return null;

    }

    public static boolean promptInitProject() {

        boolean error;
        String initClass, initFunction;

        String[] possibilities = initProjectList;

        ImageIcon icon = null;

        String s = (String) JOptionPane.showInputDialog(
                mainframe,
                "Choose Class:",
                "Initialize Project",
                JOptionPane.PLAIN_MESSAGE,
                icon,
                possibilities,
                initProjectClassStr(initClassAndFunction));

        if ((s == null) || (s.length() == 0)) {
            return true; // error
        }

        initClass = s;

        if (initProjectClass(initClass)) {
            return true; // error
        }

        possibilities = project.initProject.initFuncList;

        if (possibilities.length < 1) {
            return true; // error
        }

        s = (String) JOptionPane.showInputDialog(
                mainframe,
                "Choose Function:",
                "Initialize Project",
                JOptionPane.PLAIN_MESSAGE,
                icon,
                project.initProject.initFuncList,
                possibilities[0]);

        if ((s == null) || (s.length() == 0)) {
            return true; // error
        }

        initFunction = s;

        project.nullArrays();
        project.geometry.resizeAndFill(1, 1, 1, Geometry.SPACEVALUE);

        log("executing " + initClass + "." + initFunction + "...");

        error = project.initProject.initFunction(initFunction);

        if (error) {

            log("failed to execute " + initClass + "." + initFunction);

        } else {

            initClassAndFunction = initClass + "." + initFunction;
            project.initClassAndFunction = initClass + "." + initFunction;

            resetAllParamVectors();
            project.init();

            if (mainframe != null) {
                mainframe.panel2D.resetVoxelWidth();
                mainframe.panelParams.setParamSelect(-1, -1);
            }

            log("finished " + initClassAndFunction);

        }

        return error;

    }

    public static boolean initProject(String[][] startupArgs) {
        return initProject(initClassAndFunction, startupArgs);
    }

    public static boolean initProject(String classAndFunction, String[][] startupArgsInitProject) {

        boolean error;

        String initClass = initProjectClassStr(classAndFunction);
        String initFunction = initProjectFunctionStr(classAndFunction);

        if (initProjectClass(initClass)) {
            return true; // error
        }

        project.nullArrays();

        if (startupArgsInitProject != null) {
            initStartupArgs(startupArgsInitProject, false, true);
        }

        log("executing " + classAndFunction + "...");

        error = project.initProject.initFunction(initFunction);

        if (error) {

            log("failed to execute " + classAndFunction);
            
        } else {

            initClassAndFunction = classAndFunction;
            project.initClassAndFunction = classAndFunction;

            log("finished " + classAndFunction);

        }

        if (mainframe != null) {
            mainframe.panel2D.resetVoxelWidth();
        }

        return error;

    }

    public static String currentDate() {

        Date date = new Date(System.currentTimeMillis());

        SimpleDateFormat dateFormat = new SimpleDateFormat("EEE, d MMM yyyy HH:mm:ss Z");

        return dateFormat.format(date);

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Export Text file functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean writeParamVectorOpen(String directory, String fileName) {

        String fullName;

        if ((fileName == null) || (fileName.length() == 0)) {
            log("bad output parameter file name");
            return false;
        }

        fullName = directory + fileName + ".dat";

        //log("attempting to create output parameter file: " + fullName);

        try {

            bw1 = new BufferedWriter(new FileWriter(fullName));

            bw1.write("#D3D Parameter File v" + D3D_VERSION);
            bw1.newLine();
            bw1.write("#" + currentDate());
            bw1.newLine();

            log("created output parameter file: " + fullName);

            return true;

        } catch (IOException e) {
            log("unable to write file: " + fullName);
            return false;
        }

    }

    public static boolean writeParamVectorClose() {

        if (bw1 == null) {
            return false;
        }

        try {
            bw1.close();
            bw1 = null;
            return true;
        } catch (IOException e) {
            log("unable to close file: " + bw1.toString());
            return false;
        }

    }

    public static boolean writeParamVectorArray(Object[] pvArray) {

        ParamVector pv;

        if (bw1 == null) {
            return false;
        }

        if (pvArray == null) {
            return false;
        }

        for (Object o : pvArray) {
            pv = (ParamVector) o;
            writeParamVector(pv);
        }

        return true;

    }

    public static boolean writeParamVector(ParamVector pv) {

        String oname;
        String otype;

        if (bw1 == null) {
            return false;
        }

        if (pv == null) {
            return false;
        }

        ParamObject o;
        ParamObject[] v = pv.getVector();

        try {

            for (ParamObject po : v) {

                o = po;
                oname = o.getName();
                otype = o.getType();

                if (otype.compareToIgnoreCase("Class") == 0) {
                    bw1.newLine();
                    bw1.write(otype + "=" + o.getValue());
                } else if (oname.length() > 0) {
                    bw1.write(oname + "=" + o.getValue());
                } else {
                    continue;
                }

                bw1.newLine();

            }

            bw1.newLine();
            bw1.write("***************************************");
            bw1.newLine();

        } catch (IOException e) {
            log("unable to write to file: " + bw1.toString());
            return false;
        }

        return true;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     XML functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean writeParamVectorXMLopen(String directory, String fileName) {

        String fullName1 = directory + fileName + ".xml";
        String fullName2 = directory + fileName + ".xsl";

        if ((fileName == null) || (fileName.length() == 0)) {
            return false;
        }

        try {

            bw1 = new BufferedWriter(new FileWriter(fullName1));

            bw1.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            bw1.newLine();
            bw1.write("<?xml-stylesheet type=\"text/xsl\" href=\"" + fileName + ".xsl\"?>");
            bw1.newLine();
            bw1.write("<D3D version=\"" + D3D_VERSION + "\">");

            log("created output parameter file: " + fullName1);

        } catch (IOException e) {
            log("unable to write file: " + fullName1);
            return false;
        }

        try {

            bw2 = new BufferedWriter(new FileWriter(fullName2));

            bw2.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            bw2.newLine();
            bw2.write("<xsl:stylesheet version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">");
            bw2.newLine();
            bw2.write("<xsl:template match=\"/\">");
            bw2.newLine();
            bw2.write("<html>");
            bw2.newLine();
            bw2.write("<body>");
            bw2.newLine();
            bw2.write("<center>");
            bw2.newLine();
            bw2.write("<h2>D3D Simulation Parameters</h2>");

            log("created output parameter file: " + fullName2);

        } catch (IOException e) {
            log("unable to write file: " + fullName2);
            return false;
        }

        return true;

    }

    public static boolean writeParamVectorXMLclose() {

        if ((bw1 == null) || (bw2 == null)){
            bw1 = null;
            bw2 = null;
            return false;
        }

        try {
            bw1.newLine();
            bw1.write("</D3D>");
            bw1.close();
            bw1 = null;
        } catch (IOException e) {
            log("unable to close file: " + bw1.toString());
            return false;
        }

        try {
            bw2.newLine();
            bw2.write("</center>");
            bw2.newLine();
            bw2.write("</body>");
            bw2.newLine();
            bw2.write("</html>");
            bw2.newLine();
            bw2.write("</xsl:template>");
            bw2.newLine();
            bw2.write("</xsl:stylesheet>");
            bw2.close();
            bw2 = null;
        } catch (IOException e) {
            log("unable to close file: " + bw2.toString());
            return false;
        }

        return true;

    }

    public static boolean writeParamVectorArrayXML(Object[] pvarray) {

        ParamVector pv;

        if ((bw1 == null) || (bw2 == null)){
            return false;
        }

        if (pvarray == null) {
            return false;
        }

        for (Object o : pvarray) {
            pv = (ParamVector) o;
            writeParamVectorXML(pv);
        }

        return true;

    }

    public static boolean writeParamVectorXML(ParamVector pv) {

        String oname;
        String otype;
        String classIndent = "  ";
        String[] classes = new String[10];
        String classTRcolor = "<tr bgcolor=\"#99ccff\">";

        boolean tableOpen = false;

        int classCounter = 0;

        classes[0]="D3D";

        if ((bw1 == null) || (bw2 == null)){
            return false;
        }

        if (pv == null) {
            return false;
        }

        ParamObject o;
        ParamObject[] v = pv.getVector();

        try {

            bw2.newLine();
            bw2.newLine();
            bw2.write("<table border=\"1\">");

            for (ParamObject po : v) {

                o = po;
                oname = o.getName();
                otype = o.getType();

                if (otype.compareToIgnoreCase("Class") == 0) {

                    if (tableOpen) {
                        bw2.newLine();
                        bw2.write("</table>");
                    }

                    tableOpen = true;

                    bw2.newLine();
                    bw2.newLine();
                    bw2.write("<table border=\"0\">");
                    bw2.newLine();
                    bw2.write(classTRcolor);
                    bw2.newLine();
                    bw2.write("<th width=\"300\"></th>");
                    bw2.newLine();
                    bw2.write("<th width=\"400\">" + o.getValue() + "</th>");
                    bw2.newLine();
                    bw2.write("<th width=\"100\"></th>");
                    bw2.newLine();
                    bw2.write("</tr>");

                    classTRcolor = "<tr bgcolor=\"#cccc99\">";

                    bw1.newLine();

                    for (int j = 0; j < classCounter+1; j++) {
                        bw1.write(classIndent);
                    }

                    bw1.write("<" + o.getValue() + ">"); // open new class
                    classCounter += 1;

                    classes[classCounter]=o.getValue();

                } else if (otype.compareToIgnoreCase("EndClass") == 0) {

                    bw1.newLine();
                    classCounter -= 1;

                    for (int j = 0; j < classCounter+1; j++) {
                        bw1.write(classIndent);
                    }

                    bw1.write("</" + o.getValue() + ">"); // close class

                    if (tableOpen){
                        bw2.newLine();
                        bw2.write("</table>");
                        tableOpen = false;
                    }

                } else if (oname.length() > 0) {

                    bw2.newLine();
                    bw2.write("<tr>");
                    bw2.newLine();
                    bw2.write("<td>" + oname + "</td>");
                    bw2.newLine();
                    bw2.write("<td>");
                    bw2.write("<xsl:value-of select=\"");

                    for (int j = 0; j < classCounter+1; j++) {
                        bw2.write(classes[j]+"/");
                    }

                    bw2.write(oname+"/@value");
                    bw2.write("\"/>");
                    bw2.write("</td>");
                    bw2.newLine();
                    bw2.write("<td>");
                    bw2.write("<xsl:value-of select=\"");

                    for (int j = 0; j < classCounter+1; j++) {
                        bw2.write(classes[j]+"/");
                    }

                    bw2.write(oname+"/@units");
                    bw2.write("\"/>");
                    bw2.write("</td>");
                    bw2.newLine();
                    bw2.write("</tr>");

                    bw1.newLine();

                    for (int j = 0; j < classCounter+1; j++) {
                        bw1.write(classIndent);
                    }

                    //bw1.write("<" + oname + " units=\"" + o.getUnits() + "\">" + o.getValue() + "</" + oname + ">");

                    bw1.write("<" + oname + " value=\"" + o.getValue() + "\" units=\"" + o.getUnits() + "\" />");

                }

            }

            bw2.newLine();
            bw2.newLine();
            bw2.write("</table>");
            bw2.newLine();
            bw2.write("<p></p>");

        } catch (IOException e) {
            return false;
        }

        return true;

    }

    // export 3D array to a text file
    public static void export3DArray(String file, CoordinatesVoxels c, double[][][] p) {

        int kcount = 50;

        if (p == null) {
            return;
        }

        DataOutputStream dout;
        StopWatch time = new StopWatch();

        if (file.length() == 0) {
            return;
        }

        try {

            dout = new DataOutputStream(new FileOutputStream(file));

            log("writing to file: " + file);

            time.start();

            for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {

                if ((k > 0) && Math.IEEEremainder(k, kcount) == 0) {
                    time.stop();
                    log("writing k: " + k + ", t: " + time.toString());
                    time.start();
                }

                for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                    for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                        if ((i >= 0) && (i < p.length) && (j >= 0) && (j < p[0].length) &&
                                (k >= 0) && (k < p[0][0].length)) {
                            dout.writeBytes("\n");
                            dout.writeBytes(Integer.toString(k));
                            dout.writeBytes(",");
                            dout.writeBytes(Integer.toString(j));
                            dout.writeBytes(",");
                            dout.writeBytes(Integer.toString(i));
                            dout.writeBytes(",");
                            dout.writeBytes(Double.toString(p[i][j][k]));
                        }
                    }
                }
            }

            time.stop();

            dout.close();

            log("finished exporting PSF");

        } catch (IOException e) {
            log("unable to write file: " + file);
        }

    }

    // export 3D array to a text file
    public static void export3DArrayLaci(String file, CoordinatesVoxels c, double[][][] p) {
        int kcount = 50;

        if (p == null) {
            return;
        }

        DataOutputStream dout;
        StopWatch time = new StopWatch();

        if (file.length() == 0) {
            return;
        }

        try {

            dout = new DataOutputStream(new FileOutputStream(file));

            log("writing to file: " + file);

            time.start();

            for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {

                if ((k > 0) && Math.IEEEremainder(k, kcount) == 0) {
                    time.stop();
                    log("writing k: " + k + ", t: " + time.toString());
                    time.start();
                }

                for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                    for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                        if ((i >= 0) && (i < p.length) && (j >= 0) && (j < p[0].length) &&
                                (k >= 0) && (k < p[0][0].length)) {

                            dout.writeBytes(Double.toString(p[i][j][k]));
                            dout.writeBytes("\n"); // tab

                        }
                    }
                }

                dout.writeBytes("\n"); // new line

            }

            time.stop();

            dout.close();

            log("finished exporting PSF");

        } catch (IOException e) {
            log("unable to write file: " + file);
        }

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Geometry Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static void saveGeometry(String fn) {

        boolean save = writeObject(fn, project.geometry);

        if (save) {
            log("saved geometry to: " + fn);
        } else {
            log("failed to save geometry to: " + fn);
        }

    }

    public static Geometry openGeometry(String fn) {

        Object o = readObject(fn);
        Geometry g = null;

        if (o instanceof Geometry) {

            g = (Geometry) o;
            project.geometry = g;
            project.geomFile = fn;
            g.file = fn;

            resetAllParamVectors();

            if (mainframe != null) {
                mainframe.initTabs();
            }

            log("opened geometry: " + fn);

        } else {
            log("open geometry error: " + fn + " is not of appropriate format.");
        }

        return g;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Diffusant Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static String diffusantName() {
        return "Diffusant" + Integer.toString(project.numDiffusants());
    }

    public static int addDiffusant() {
        double C = 0;
        double D = 1;
        return addDiffusant(C, D);
    }

    public static int addDiffusant(double C, double D) {
        return addDiffusant(diffusantName(), C, D);
    }

    // create simple diffusant
    public static int addDiffusant(String name, double C, double D) {
        int dnum = project.addDiffusant(new Diffusant(project, name, C, D, null));
        updatePanel2D();
        return dnum;
    }

    public static int addDiffusantReactant(double cTotal, double D, int DiffusantNum, double Kon,
            double Koff) {
        return addDiffusantReactant(diffusantName(), cTotal, D, DiffusantNum, Kon, Koff);
    }

    // create a buffer diffusant that reacts with another diffusant
    public static int addDiffusantReactant(String name, double cTotal, double D, int DiffusantNum,
            double Kon, double Koff) {
        int dd = project.addDiffusant(new DiffusantReactant(project, name, cTotal, D, null, DiffusantNum, Kon, Koff));
        updatePanel2D();
        return dd;
    }

    public static int addDiffusantDye(double cTotal, double D, int DiffusantNum, double Kon,
            double Koff, double R) {
        return addDiffusantDye(diffusantName(), cTotal, D, DiffusantNum, Kon, Koff, R);
    }

    // create a dye diffusant that reacts with another diffusant
    public static int addDiffusantDye(String name, double cTotal, double D, int DiffusantNum,
            double Kon, double Koff, double R) {
        int dd = project.addDiffusant(new DiffusantDye(project, name, cTotal, D, null, DiffusantNum, Kon, Koff, R));
        updatePanel2D();
        return dd;
    }

    public static Diffusant addDiffusantPrompt() {

        Object[] possibilities = diffusantList;

        ImageIcon icon = null;

        String s = (String) JOptionPane.showInputDialog(
                mainframe,
                "Choose Class:",
                "Add Diffusant",
                JOptionPane.PLAIN_MESSAGE,
                icon,
                possibilities,
                initProjectClassStr(initClassAndFunction));

        if ((s == null) || (s.length() == 0)) {
            return null;
        }

        return null;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Source Functions
    //
    ///////////////////////////////////////////////////////////////////
    
    public static String sourceName() {
        return "Source" + Integer.toString(project.numSources());
    }

    public static Source addSource() {

        int diffusantNum = 0;
        double onset = 1;
        double duration = 1;
        double ctotal = 1; // mM
        
        PulseTimer pt = new PulseTimer(project, onset, duration, ctotal);

        return addSource(diffusantNum, project.geometry, pt);

    }
    
    public static Source addSource(int DiffusantNum, CoordinatesVoxels c, double onset, double duration, double ctotal) {
        PulseTimer pt = new PulseTimer(project, onset, duration, ctotal);
        return addSource(sourceName(), DiffusantNum, c, pt);
    }
    
    public static Source addSourceNtotal(int DiffusantNum, CoordinatesVoxels c, double onset, double duration, double ntotal) {
        double ctotal = Utility.N2C(ntotal, project.dx, c.spaceVoxels);
        PulseTimer pt = new PulseTimer(project, onset, duration, ctotal);
        return addSource(sourceName(), DiffusantNum, c, pt);
    }
    
    public static Source addSource(int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        return addSource(sourceName(), DiffusantNum, c, pt);
    }
    
    public static Source addSource(String name, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        Source s = new Source(project, name, DiffusantNum, c, pt);
        project.addSource(s);
        updatePanel2D();
        return s;
    }
    
    public static Source addSourceClamp(int DiffusantNum, CoordinatesVoxels c, double clampValue) {
        return addSourceClamp(sourceName(), DiffusantNum, c, clampValue);
    }
    
    public static Source addSourceClamp(String name, int DiffusantNum, CoordinatesVoxels c, double clampValue) {
        Source s = new Source(project, name, DiffusantNum, c, clampValue);
        project.addSource(s);
        updatePanel2D();
        return s;
    }

    public static SourceGauss addSourceGauss() {

        int diffusantNum = 0;
        double tpeak = 1;
        double stdv = 1;
        double ctotal = 1;

        return addSourceGauss(diffusantNum, project.geometry, tpeak, stdv, ctotal);

    }
    
    public static SourceGauss addSourceGauss(int DiffusantNum, CoordinatesVoxels c,
            double tpeak, double stdv, double cTotal) {
        PulseTimer pt = new PulseTimer(project, tpeak, stdv, cTotal);
        return addSourceGauss(sourceName(), DiffusantNum, c, pt);
    }
    
    public static SourceGauss addSourceGauss(int DiffusantNum, CoordinatesVoxels c,
            PulseTimer pt) {
        return addSourceGauss(sourceName(), DiffusantNum, c, pt);
    }

    public static SourceGauss addSourceGauss(String name, int DiffusantNum, CoordinatesVoxels c,
            PulseTimer pt) {
        SourceGauss s = new SourceGauss(project, name, DiffusantNum, c, pt);
        project.addSource(s);
        updatePanel2D();
        return s;
    }
    
    public static SourceGamma addSourceGamma() {

        int diffusantNum = 0;
        double onset = 1;
        double tau = 1;
        double ctotal = 1;

        return addSourceGamma(diffusantNum, project.geometry, onset, tau, ctotal);

    }
    
    public static SourceGamma addSourceGamma(int DiffusantNum, CoordinatesVoxels c,
            double onset, double tau, double cTotal) {
        PulseTimer pt = new PulseTimer(project, onset, tau, cTotal);
        return addSourceGamma(sourceName(), DiffusantNum, c, pt);
    }
    
    public static SourceGamma addSourceGamma(int DiffusantNum, CoordinatesVoxels c,
            PulseTimer pt) {
        return addSourceGamma(sourceName(), DiffusantNum, c, pt);
    }

    public static SourceGamma addSourceGamma(String name, int DiffusantNum, CoordinatesVoxels c,
            PulseTimer pt) {
        SourceGamma s = new SourceGamma(project, name, DiffusantNum, c, pt);
        project.addSource(s);
        updatePanel2D();
        return s;
    }
    
    public static SourceUptake addSourceUptake() {
        int diffusantNum = 0;
        double krate = 1;
        return addSourceUptake(diffusantNum, project.geometry, krate);
    }
    
    public static SourceUptake addSourceUptake(int DiffusantNum, CoordinatesVoxels c,
            double krate) {
        return addSourceUptake(sourceName(), DiffusantNum, c, krate);
    }
    
    public static SourceUptake addSourceUptake(String name, int DiffusantNum, CoordinatesVoxels c,
            double krate) {
        SourceUptake s = new SourceUptake(project, name, DiffusantNum, c, krate);
        project.addSource(s);
        updatePanel2D();
        return s;
    }
    
    public static Source addSourceImpulse(int DiffusantNum, CoordinatesVoxels c, double t, double ctotal) {
        PulseTimer pt = new PulseTimer(project, t, 0, ctotal);
        return addSource(sourceName(), DiffusantNum, c, pt);
    }

    public static void addSourceImpulseRandom(int DiffusantNum, int num, int Ntotal) {

        int i, j, k, im, jm, km;
        double t;
        int numVoxels = 1;
        double ctotal = Utility.N2C(Ntotal, project.dx, numVoxels);

        CoordinatesVoxels coordinates = new CoordinatesVoxels(project);

        int xVoxels = project.geometry.xVoxels;
        int yVoxels = project.geometry.yVoxels;
        int zVoxels = project.geometry.zVoxels;

        for (int ii = 0; ii < num; ii++) {
            im = xVoxels;
            jm = yVoxels;
            km = zVoxels;
            i = (int) (im * Math.random());
            j = (int) (jm * Math.random());
            k = (int) (km * Math.random());
            t = project.simTime * Math.random();
            if ((i > 0) && (j > 0) && (k > 0)) {
                if ((i < im - 1) && (j < jm - 1) && (k < km - 1)) {
                    coordinates.setVoxelPoint(i, j, k);
                    addSourceImpulse(DiffusantNum, coordinates, t, ctotal);
                }
            }
        }

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Detector Functions (Average, Weighted Average, Snapshot)
    //
    ///////////////////////////////////////////////////////////////////

    public static String detectorName() {
        return "Detector" + Integer.toString(project.numDetectors());
    }

    public static Detector addDetector() {
        int diffusantNum = 0;
        return addDetector(diffusantNum);
    }

    public static Detector addDetector(int DiffusantNum) {
        return addDetector(detectorName(), DiffusantNum);
    }

    public static Detector addDetector(String name, int DiffusantNum) {
        Detector dd = new Detector(project, name, DiffusantNum, null);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static Detector addDetector(int DiffusantNum, CoordinatesVoxels c) {
        return addDetector(detectorName(), DiffusantNum, c);
    }

    public static Detector addDetector(String name, int DiffusantNum, CoordinatesVoxels c) {
        Detector dd = new Detector(project, name, DiffusantNum, c);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorAvg addDetectorAvg() {
        int diffusantNum = 0;
        return addDetectorAvg(diffusantNum);
    }

    public static DetectorAvg addDetectorAvg(int DiffusantNum) {
        return addDetectorAvg(detectorName(), DiffusantNum);
    }

    // create average detector over entire shape
    public static DetectorAvg addDetectorAvg(String name, int DiffusantNum) {
        DetectorAvg dd = new DetectorAvg(project, name, DiffusantNum, null);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorAvg addDetectorAvg(int DiffusantNum, CoordinatesVoxels c) {
        return addDetectorAvg(detectorName(), DiffusantNum, c);
    }

    // create average detector within coordinates
    public static DetectorAvg addDetectorAvg(String name, int DiffusantNum, CoordinatesVoxels c) {
        DetectorAvg dd = new DetectorAvg(project, name, DiffusantNum, c);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorPSF addDetectorPSF(int DiffusantNum, PSF psf) {
        return addDetectorPSF(detectorName(), DiffusantNum, psf);
    }

    // create a weighted detector with PSF (over entire shape)
    public static DetectorPSF addDetectorPSF(String name, int DiffusantNum, PSF psf) {
        DetectorPSF dd = new DetectorPSF(project, name, DiffusantNum, null, psf);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorPSF addDetectorPSF(int DiffusantNum, CoordinatesVoxels c, PSF psf) {
        return addDetectorPSF(detectorName(), DiffusantNum, c, psf);
    }

    // create a weighted detector with PSF
    public static DetectorPSF addDetectorPSF(String name, int DiffusantNum, CoordinatesVoxels c, PSF psf) {
        DetectorPSF dd = new DetectorPSF(project, name, DiffusantNum, c, psf);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorPSF addDetectorPSF_Gauss() {
        int diffusantNum = 0;
        double xSTDV = 1;
        double ySTDV = 1;
        double zSTDV = 1;
        return addDetectorPSF_Gauss(diffusantNum, xSTDV, ySTDV, zSTDV);
    }

    public static DetectorPSF addDetectorPSF_Gauss(int DiffusantNum, double xSTDV, double ySTDV, double zSTDV) {
        return addDetectorPSF_Gauss(detectorName(), DiffusantNum, xSTDV, ySTDV, zSTDV);
    }

    public static DetectorPSF addDetectorPSF_Gauss(String name, int DiffusantNum, double xSTDV, double ySTDV, double zSTDV) {

        PSF psf = new PSFgauss(project, null, xSTDV, ySTDV, zSTDV);

        return addDetectorPSF(name, DiffusantNum, psf);

    }

    public static DetectorPSF addDetectorPSF_Gauss(int DiffusantNum, double xSTDV, double ySTDV, double zSTDV,
            double i0, double j0, double k0) {
        return addDetectorPSF_Gauss(detectorName(), DiffusantNum, xSTDV, ySTDV, zSTDV, i0, j0, k0);
    }

    public static DetectorPSF addDetectorPSF_Gauss(String name, int DiffusantNum, double xSTDV, double ySTDV, double zSTDV,
            double i0, double j0, double k0) {

        PSF psf = new PSFgauss(project, null, xSTDV, ySTDV, zSTDV);

        psf.setVoxelCenter(i0, j0, k0);

        return addDetectorPSF(name, DiffusantNum, psf);

    }

    public static DetectorPSF addDetectorPSF_Gauss(int DiffusantNum, CoordinatesVoxels c,
            double xSTDV, double ySTDV, double zSTDV) {
        return addDetectorPSF_Gauss(detectorName(), DiffusantNum, c, xSTDV, ySTDV, zSTDV);
    }

    public static DetectorPSF addDetectorPSF_Gauss(String name, int DiffusantNum, CoordinatesVoxels c,
            double xSTDV, double ySTDV, double zSTDV) {

        PSF psf = new PSFgauss(project, c, xSTDV, ySTDV, zSTDV);

        return addDetectorPSF(name, DiffusantNum, c, psf);

    }

    public static DetectorDye addDetectorDye(int DiffusantNum, double xSTDV, double ySTDV, double zSTDV,
            double i0, double j0, double k0) {
        return Master.addDetectorDye(detectorName(), DiffusantNum, xSTDV, ySTDV, zSTDV, i0, j0, k0);
    }

    public static DetectorDye addDetectorDye(String name, int DiffusantNum, double xSTDV, double ySTDV, double zSTDV,
            double i0, double j0, double k0) {

        PSF psf = new PSFgauss(project, null, xSTDV, ySTDV, zSTDV);

        psf.setVoxelCenter(i0, j0, k0);

        return Master.addDetectorDye(name, DiffusantNum, psf);

    }

    public static DetectorDye addDetectorDye(int DiffusantNum, PSF psf) {
        return Master.addDetectorDye(detectorName(), DiffusantNum, psf);
    }

    public static DetectorDye addDetectorDye(String name, int DiffusantNum, PSF psf) {
        DetectorDye dd = new DetectorDye(project, name, DiffusantNum, null, psf);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorDye addDetectorDye(int DiffusantNum, CoordinatesVoxels c, PSF psf) {
        return Master.addDetectorDye(detectorName(), DiffusantNum, c, psf);
    }

    public static DetectorDye addDetectorDye(String name, int DiffusantNum, CoordinatesVoxels c, PSF psf) {
        DetectorDye dd = new DetectorDye(project, name, DiffusantNum, c, psf);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorSnapshot addDetectorSnapshot() {
        int diffusantNum = 0;
        double t = 0;
        return addDetectorSnapshot(diffusantNum, t);
    }
    
    public static DetectorSnapshot addDetectorSnapshot(int DiffusantNum, double t) {
        return addDetectorSnapshot(detectorName(), DiffusantNum, t);
    }

    // create a snapshot detector at time t over entire shape
    public static DetectorSnapshot addDetectorSnapshot(String name, int DiffusantNum, double t) {
        DetectorSnapshot dd = new DetectorSnapshot(project, name, DiffusantNum, null, t);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    public static DetectorSnapshot addDetectorSnapshot(int DiffusantNum, String xyzSelect, int plane, double t) {
        return addDetectorSnapshot(detectorName(), DiffusantNum, xyzSelect, plane, t);
    }

    public static DetectorSnapshot addDetectorSnapshot(String name, int DiffusantNum, String xyzSelect, int plane, double t) {

        CoordinatesVoxels coordinates = new CoordinatesVoxels(project);

        if ((xyzSelect.compareToIgnoreCase("xy") == 0) ||
                (xyzSelect.compareToIgnoreCase("yx") == 0)) { // snapshot the whole xy plane
            coordinates.zVoxel1 = plane;
            coordinates.zVoxel2 = plane;
        } else if ((xyzSelect.compareToIgnoreCase("yz") == 0) ||
                (xyzSelect.compareToIgnoreCase("zy") == 0)) { // snapshot the whole yz plane
            coordinates.xVoxel1 = plane;
            coordinates.xVoxel2 = plane;
        } else if ((xyzSelect.compareToIgnoreCase("zx") == 0) ||
                (xyzSelect.compareToIgnoreCase("xz") == 0)) { // snapshot the whole zx plane
            coordinates.yVoxel1 = plane;
            coordinates.yVoxel2 = plane;
        }

        coordinates.update();

        return addDetectorSnapshot(name, DiffusantNum, coordinates, t);

    }

    public static DetectorSnapshot addDetectorSnapshot(int DiffusantNum, CoordinatesVoxels c, double t) {
        return addDetectorSnapshot(detectorName(), DiffusantNum, c, t);
    }

    // create a snapshot detector at time t
    public static DetectorSnapshot addDetectorSnapshot(String name, int DiffusantNum, CoordinatesVoxels c, double t) {
        DetectorSnapshot dd = new DetectorSnapshot(project, name, DiffusantNum, c, t);
        project.addDetector(dd);
        updatePanel2D();
        return dd;
    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Photolysis Functions (i.e. MNI uncaging or FRAP, PSFs)
    //
    ///////////////////////////////////////////////////////////////////
    
    public static boolean photoExists() {
        Diffusant[] diffusant = project.diffusants;
        if (diffusant != null) {
            for (Diffusant d : diffusant) {
                if (d instanceof DiffusantPhoto) {
                    return true;
                }
            }
        }
        return false;
    }

    public static DiffusantPhoto getPhoto() {
        Diffusant[] diffusant = project.diffusants;
        if (diffusant != null) {
            for (Diffusant d : diffusant) {
                if (d instanceof DiffusantPhoto) {
                    return (DiffusantPhoto) d;
                }
            }
        }
        return null;
    }

    public static boolean openPSFtorok(PSFtorok psf, String fn) {

        Object o = readObject(fn);

        if (o instanceof PSF) {

            //psf = (PSFtorok) o;
            project.psfFile = fn;
            resetAllParamVectors();

            if (mainframe != null) {
                mainframe.initTabs();
            }

            log("opened PSF: " + fn);

            return true;

        } else if (psf.readRandomAccessHeader(fn)) {

            psf.array = null; // delete old PSF if it exists
            psf.randomAccess = true;
            psf.fileName = fn;
            project.psfFile = fn;
            resetAllParamVectors();

            if (mainframe != null) {
                mainframe.initTabs();
            }

            return true;

        } else {
            log("open PSF error: " + fn + " is not of appropriate format.");
        }
        return false;
    }

    public static boolean savePSFtorok(PSFtorok psf, String fName) {
        return writeObject(fName, psf);
    }

    // create PSF files
    public static void writePSF(String fName, int imax, int jmax, int kmax) {
        if (photoExists()) {
            PSFtorok psf = (PSFtorok) getPhoto().psf;
            psf.writeRandomAccess(fName, imax, jmax, kmax);
        }
    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Batch Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static Batch addBatch(int num, String varName, double value) {
        Batch b = new Batch(num, varName, value, true);
        project.addBatch(b);
        return b;
    }

    public static Batch addBatch(int num, String varName, String strValue) {
        Batch b = new Batch(num, varName, strValue, true);
        project.addBatch(b);
        return b;
    }

    

    public static Batch addBatch(int num, String varName, double value, boolean execute) {
        Batch b = new Batch(num, varName, value, execute);
        project.addBatch(b);
        return b;
    }

    public static Batch addBatch(int num, String varName, String strValue, boolean execute) {
        Batch b = new Batch(num, varName, strValue, execute);
        project.addBatch(b);
        return b;
    }

    public static void addBatchList(String varName, String valueList) {
        addBatchList(0, varName, valueList);
    }

    public static void addBatchList(int startNum, String varName, String valueList) {

        double dbl;
        int items = Utility.itemsInList(valueList, ',');

        String str;
        Batch b;

        for (int i = 0; i < items; i += 1) {

            str = Utility.itemFromList(valueList, i);

            try {
                dbl = Double.parseDouble(str);
                b = new Batch(startNum + i, varName, dbl, true);
                project.addBatch(b);
            } catch (NumberFormatException e) {
            }

        }

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     Parameter Vector Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static void updateAllParamVectors() {

        if (project != null) {

            System.out.println("updateAllParamVectors");

            project.updateVectors();

            if (project.geometry != null) {
                project.geometry.updateVectors();
            }

            updateParamVectorArray(project.diffusants);
            updateParamVectorArray(project.sources);
            updateParamVectorArray(project.detectors);
            updateParamVectorArray(project.batches);
            updateParamVectorArray(project.errors);

        }

    }

    public static void updateParamVectorArray(Object[] pv) {

        ParamVector temp;

        if (pv != null) {
            for (Object o : pv) {
                if (o instanceof ParamVector) {
                    temp = (ParamVector) o;
                    temp.updateVectors();
                }
            }
        }

    }

    public static void resetAllParamVectors() {

        if (project != null) {

            project.clearVector();
            project.createVector(true);

            if (project.geometry != null) {
                project.geometry.clearVector();
                project.geometry.createVector(true);
            }

            resetParamVectorArray(project.diffusants);
            resetParamVectorArray(project.sources);
            resetParamVectorArray(project.detectors);
            resetParamVectorArray(project.batches);
            resetParamVectorArray(project.errors);
            
        }

    }

    public static void resetParamVectorArray(Object[] pv) {

        ParamVector temp;

        if (pv != null) {
            for (Object o : pv) {
                if (o instanceof ParamVector) {
                    temp = (ParamVector) o;
                    temp.clearVector();
                    temp.createVector(true);
                }
            }
        }

    }

    public static boolean promptYes(String title, String message) {

        int response = JOptionPane.showConfirmDialog(null, message, title,
                JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);

        switch(response) {
            case JOptionPane.NO_OPTION:
                return false;
            case JOptionPane.YES_OPTION:
                return true;
            case JOptionPane.CLOSED_OPTION:
                return false;
        }

        return false;

    }

    // get file name
    public static String getFileName(String pathname) {
        
        String pstr;
        
        if (pathname == null) {
            pathname = "";
        }

        pstr = pathname;
        
        if (pstr.isEmpty()) {
            pstr = System.getProperty("user.home");
        }

        JFileChooser fc = new JFileChooser();
        fc.setSelectedFile(new File(pstr));
        fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        int approve = fc.showOpenDialog(mainframe);

        if (approve == JFileChooser.APPROVE_OPTION) {
            pstr = fc.getCurrentDirectory() + "\\" + fc.getSelectedFile().getName();
            //pstr = pstr.substring(2);
            pstr = pstr.replace('\\', '/');
            return pstr;
        }
        
        return pathname;

    }

    public static double[] readTextFile(String file) {

        RandomAccessFile raFile;

        int counter = 0;
        long numBytes;
        String strValue;
        int chunks = 1000;
        boolean readFile = true;

        boolean foundNonZero = false;
        int counter2 = 10;

        double[] temp;

        try {

            raFile = new RandomAccessFile(new File(file), "r"); // open file

            Master.log("reading File: " + file);

            numBytes = raFile.length();

            temp = new double[chunks];

            do {

                strValue = raFile.readLine();

                if (strValue == null) {
                    break;
                }

                if (counter >= temp.length) {
                    temp = Utility.resizeArray(temp, temp.length + chunks);
                }

                if ((counter >= 0) && (counter < temp.length)) {

                    temp[counter] = Double.parseDouble(strValue);

                    //if (counter < 10) {
                    //    Master.log("" + counter + "," + temp[counter]+ "," + strValue);
                    //}

                    //if (Math.IEEEremainder(counter, 1000) == 0) {
                   //     Master.log("" + counter + "," + temp[counter] + "," + strValue);
                    //}

                    //if (!foundNonZero && (temp[counter] != 0)) {
                    //    Master.log("" + (counter -1) + "," + temp[counter-1]);
                    //    foundNonZero = true;
                    //}

                    //if (foundNonZero && (counter2 >= 0)) {
                    //    Master.log("" + counter + "," + temp[counter] + "," + strValue);
                    //    counter2--;
                    //}

                }

                counter++;

                if (counter >= numBytes) {
                    readFile = false;
                }

            } while (readFile);

        } catch (IOException e) {
            System.err.println("unable to read file: " + file);
            return null;
        }

        if (counter < temp.length) {
            temp = Utility.resizeArray(temp, counter);
        }

        Master.log("lines read: " + temp.length);

        try {
            raFile.close();
            Master.log("closed file: " + file);
        } catch (IOException e) {
        }

        return temp;

    }

    public static boolean saveToBinaryFile(String file, double[] array) {

        DataOutputStream binaryWriter;

        try {
            binaryWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, false))); // open and append
            Master.log("opened binary output file " + file);
        } catch (IOException e) {
            System.err.println("failed to open binary output file " + file);
            return false;
        }

        try {
            for (int i = 0; i < array.length; i++) {
                binaryWriter.writeDouble(array[i]);
            }
        } catch (IOException e) {
            System.err.println("unable to write to binary output file " + file);
            return false;
        }

        try {
            binaryWriter.close();
            Master.log("closed binary output file " + file);
            return true;
        } catch (IOException e) {
            System.err.println("unable to close binary output file " + file);
            return false;
        }

    }

    public static double[] readBinaryFile(String fileName) {

        DataInputStream binaryReader;

        int counter = 0;
        double dvalue;
        int chunks = 1000;

        double[] temp;

        boolean foundNonZero = false;
        int counter2 = 10;

        try {

            binaryReader = new DataInputStream(new BufferedInputStream(new FileInputStream(fileName)));

            Master.log("reading File: " + fileName + ", " + binaryReader.available() + " bytes");

            temp = new double[chunks];

            while (binaryReader.available() > 0) {

                dvalue = binaryReader.readDouble();

                if (counter >= temp.length) {
                    temp = Utility.resizeArray(temp, temp.length + chunks);
                }

                temp[counter] = dvalue;

                if (Math.IEEEremainder(counter, 5000000) == 0) {
                    Master.log("read " + counter + " doubles");
                }

                if (!foundNonZero && (counter > 1) && (temp[counter] != 0)) {
                    Master.log("" + (counter - 1) + "," + temp[counter - 1]);
                    foundNonZero = true;
                }

                if (foundNonZero && (counter2 >= 0)) {
                    Master.log("" + counter + "," + temp[counter]);
                    counter2--;
                }

                counter++;

            }

        } catch (IOException e) {
            System.err.println("unable to read file: " + fileName);
            return null;
        }

        if (counter < temp.length) {
            temp = Utility.resizeArray(temp, counter);
        }

        Master.log("lines read: " + temp.length);

        try {
            binaryReader.close();
            Master.log("closed file: " + fileName);
        } catch (IOException e) {
        }

        return temp;

    }

}


