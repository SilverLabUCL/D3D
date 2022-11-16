package ucl.silver.d3d.core;

import javax.swing.UIManager;
import ucl.silver.d3d.utils.*;
import java.util.Date;

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
public class StartD3D {

    public StartD3D(String[] args) {

        boolean runSimulation = false;
        boolean packFrame = false;
        long seed = -1;

        Master.startDate = new Date(System.currentTimeMillis());

        Master.log("D3D started " + Master.startDate.toString());

        String projectType = "FD";

        String initClassAndFunction = "";

        String[][] startupArgs = null;
        String[][] startupArgsProject;
        String[][] startupArgsInitProject;

        String[] splitStr;

        String[] testing = new String[7];
        testing[0] = "-MonteCarlo";
        testing[1] = "-Project.simTime";
        testing[2] = "1.2";
        testing[3] = "-Project.Directory";
        testing[4] = "/cluster/project3/silverlab1/d3d/jason/frap/";
        testing[5] = "-Project.EM_iseries";
        testing[6] = "3";
        //testing[4] = "-initProject.source_numchannels";
        //testing[5] = "16";
        //testing[6] = "-EGTA.cTotal";
        //testing[7] = "10";

        //testing[6] = "-Geometry.cubeWidth";
        //testing[7] = "2.4";
        //testing[5] = "-Particles.D";
        //testing[6] = "0.2";
        //testing[6] = "-Detector0.PSF.numericalAperture";
        //testing[7] = "0.7";

        //testing[0] = "-init";
        //testing[1] = "InitMC_Frap_Jason.initMonteCarloActiveZone";
        //testing[2] = "-Particles.D";
        //testing[3] = "8.250E-5";
        //testing[4] = "-Geometry.cubeWidth";
        //testing[5] = "2.0";
        //testing[6] = "-Particles.setImmobilePercent";
        //testing[7] = "0.1875";
        //testing[5] = "-MonteCarlo.minParticleStep";
        //testing[6] = "0.003";
        //testing[7] = "-MonteCarlo.MSDspatialTbgn";
        //testing[8] = "100";
        //testing[5] = "-MonteCarlo.azHemisphere";
        //testing[6] = "1";
        //testing[10] = "-Project.Directory";
        //testing[11] = "/cluster/project3/silverlab1/d3d/jason/frap/";
        //testing[12] = "-Project.Folder";
        //testing[13] = "c20_d080_im18p75_0";

        //args = testing; // THIS ALLOWS TESTING args

        if ((args != null) && (args.length > 0)) {

            Master.foundMainStartUpArguments = true;

            startupArgs = new String[args.length][2];

            for (int i = 0; i < args.length; i++) {

                if (args[i] == null) {
                    break;
                }

                if (args[i].equals("-?") || args[i].equalsIgnoreCase("-help")) {
                    printUsageAndQuit();
                } else if (args[i].equalsIgnoreCase("-noGUI")) {
                    System.out.println("Starting without GUI");
                    Master.showGUI = false;
                } else if (args[i].equalsIgnoreCase("-run")) {
                    System.out.println("Auto-run simulation");
                    runSimulation = true;
                } else if (args[i].equalsIgnoreCase("-MonteCarlo")) {
                    projectType = "MC";
                } else if (args[i].equalsIgnoreCase("-MonteCarloActiveZone")) {
                    projectType = "MCAZ";
                } else if (args[i].equalsIgnoreCase("-init")) {

                    if ((args[i + 1].length() == 0) || args[i + 1].startsWith("-")) {
                        printUsageAndQuit();
                    }

                    splitStr = args[i + 1].split("\\.");

                    if (splitStr.length != 2) { // should be split into 2 Strings
                        System.out.println("-init format error");
                        printUsageAndQuit();
                    }

                    if ((splitStr[0] == null) || (splitStr[0].length() == 0)) {
                        System.out.println("-init format error");
                        printUsageAndQuit();
                    }

                    if ((splitStr[1] == null) || (splitStr[1].length() == 0)) {
                        System.out.println("-init format error");
                        printUsageAndQuit();
                    }

                    initClassAndFunction = args[i + 1];

                } else if (args[i].startsWith("-")) {

                    if (i + 1 < args.length) {
                        startupArgs[i][0] = args[i].substring(1); // long name without first dash
                        startupArgs[i][1] = args[i + 1]; // value
                    } else {
                        System.out.println("-parameter format error: missing final parameter");
                        printUsageAndQuit();
                    }

                }
            }

        }

        startupArgsProject = getClassStartupArgs(startupArgs, "Project");
        startupArgsInitProject = getClassStartupArgs(startupArgs, "InitProject");

        Master.createProject(projectType, "", startupArgsProject);
        Master.initMersenneTwister();

        if (Master.showGUI) {
            Master.initMainFrame(packFrame);
        }

        if (initClassAndFunction.length() > 0) {
            if (Master.initProject(initClassAndFunction, startupArgsInitProject)) {
                System.exit(0);
            }
        } else if (Master.initProject(startupArgsInitProject)) {
            System.exit(0);
        }

        if (startupArgs != null) {
            if (Master.initStartupArgs(startupArgs, false, false)) {
                //System.out.println("failed initStartupArgs");
                printUsageAndQuit();
            }
        }
        
        Master.resetAllParamVectors();
        Master.project.init();

        if (Master.showGUI && (Master.mainframe != null)) {

            if (Master.mainframe.panelParams != null) {
                if (Master.project.errors != null) {
                    Master.mainframe.panelParams.setParamSelect(7, -1); // error
                } else {
                    Master.mainframe.panelParams.setParamSelect(1, -1);
                }
                Master.mainframe.panelParams.updateControls();
            }

            if (Master.mainframe.panel2D != null) {
                Master.mainframe.panel2D.resetVoxelWidth();
            }

            Master.mainframe.setVisible(true);

        }

        if (runSimulation) {
            //System.out.println("starting project simulation...");
            Master.project.simulationStart(false);
        }

    }

    private static String[][] getClassStartupArgs(String[][] startupArgs, String className) {

        int j = 0;
        String parent, parameterName;

        if ((startupArgs == null) || (startupArgs.length == 0)) {
            return null;
        }

        String[][] s = new String[startupArgs.length][2];

        for (int i = 0; i < startupArgs.length; i++) {

            parent = Utility.parentClass(startupArgs[i][0]);
            parameterName = Utility.parameterName(startupArgs[i][0]);

            if ((parent == null) || (parameterName == null)) {
                continue;
            }

            if (parent.equalsIgnoreCase(className)) {
                s[j][0] = startupArgs[i][0];
                s[j][1] = startupArgs[i][1];
                startupArgs[i][0] = null;
                startupArgs[i][1] = null;
                j++;
            }
            
        }

        return s;

    }

    private static void printUsageAndQuit() {
        System.out.println("\nUsage: java  ucl.silver.d3d.core.StartD3D [-options]");
        System.out.println("where options include: \n");
        System.out.println("-? -help        print this help message");
        System.out.println("-MonteCarlo     run Monte Carlo simulator (default finite difference)");
        System.out.println("-noGUI          run without GUI");
        System.out.println("-run            start simulation");
        System.out.println("-init Class.function   init class and function (-init InitProject.initCube");
        System.out.println("-name value     parameter name value pairs (-Project.simTime 2.5 -Geometry.xLength 0.8 -MonteCarlo.MSD.save2TextFile 1)");
        System.exit(0);
    }

    // main method, everything begins here
    public static void main(String[] args) {

        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            //e.printStackTrace();
        }

        StartD3D d3d = new StartD3D(args);

    }

}
