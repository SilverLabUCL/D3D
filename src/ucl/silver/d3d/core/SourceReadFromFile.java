/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ucl.silver.d3d.core;

import java.awt.Color;
import java.io.RandomAccessFile;
import ucl.silver.d3d.utils.Utility;

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
public class SourceReadFromFile extends Source {
    
    public boolean useProjectArray = false;
    private boolean useArray = false;
    
    public String importFileName = ""; // import source from external data file
    public String importFileType = "binary"; // "text" or "binary"
    
    public double Cvoxel_conversionFactor = 1;
    
    private final Color defaultColor = new Color(153, 0, 51);
    
    private final RandomAccessFile raFile = null;
    public double[] array = null;
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("fileName")) {
            return "DIR";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("Cvoxel_conversionFactor")) {
            return false;
        }
        //if (name.equalsIgnoreCase("Ivoxel_conversionFactor")) {
            //return false;
        //}
        return true;
    }
    
    public SourceReadFromFile(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt, String file) {

        super(p, NAME, DiffusantNum, c, pt);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;
        importFileName = file;
        pulseTimer = pt;
        save = new Save(p);

        if (c == null) {
            coordinates = null;
        } else {
            coordinates = new CoordinatesVoxels[1];
            coordinates[0] = new CoordinatesVoxels(p);
            coordinates[0].matchVoxels(c);
        }

        createVector(true); // sets ParamVector

    }
    
    @Override
    public void updateStats() {

        int numVoxels;
        double molesPerVoxelPerDT, charge, litersPerVoxel;
        double Qtotal, Qvoxel;

        if (totalVoxelsForSpaceVoxelsOnly) {
            numVoxels = spaceVoxels;
        } else {
            numVoxels = voxels;
        }

        charge = project.diffusants[diffusantNum].charge;
        litersPerVoxel = Utility.litersPerVoxel(project.dx);
        //C2I_conversionFactor = 1e12 * charge * Utility.FARADAY * litersPerVoxel / project.dt; // convert mM to pA

        // compute factor for converting pA/totalVoxels to mM/voxel
        if (charge > 0) {
            Qtotal = 1 / (charge * Utility.FARADAY * 1e3 * 1e12); // mols/msec
        } else {
            Qtotal = Double.NaN;
        }

        Qvoxel = Qtotal / numVoxels; // mols/msec
        molesPerVoxelPerDT = Qvoxel * project.dt;
        Cvoxel_conversionFactor = 1.0e3 * molesPerVoxelPerDT / litersPerVoxel; // mM

        // compute factor for converting mM/voxel to pA/voxel
        molesPerVoxelPerDT = 1 * litersPerVoxel * 1e-3;
        Qvoxel = molesPerVoxelPerDT / project.dt; // mols/msec

        //if (charge > 0) {
            //Ivoxel_conversionFactor = Qvoxel * charge * Utility.FARADAY * 1e3 * 1e12; // convert molar flux to pA
        //} else {
            //Ivoxel_conversionFactor = Double.NaN;
        //}

    }
    
    @Override
    public double release(RunFiniteDifference fd, Geometry geometry) {

        int index, numVoxels;
        double cvoxel, avgC = 0;
        boolean timerHigh = false;
        
        if (!useArray) {
            return 0;
        }

        if (spaceVoxels == 0) {
            return 0;
        }

        if ((fd.diffus == null) || (fd.it >= fd.itmax)) {
            return 0;
        }

        if ((pulseTimer == null) || (pulseTimer.timer == null)) {
            timerHigh = true; // on for all time
        } else if (pulseTimer.timer[fd.it] > 0) {
            timerHigh = true;
        }

        if (useProjectArray) {

            if (fd.it < project.sourceArray.length) {

                cvoxel = project.sourceArray[(int) fd.it];

                if (cvoxel < 0) {
                    cvoxel = 0;
                } else {
                    cvoxel *= Cvoxel_conversionFactor;
                }

                //ivoxel = cvoxel * Ivoxel_conversionFactor;
                
            } else {
                cvoxel = 0;
                //ivoxel = 0;
            }

        } else {

            if (fd.it < array.length) {

                cvoxel = array[(int) fd.it];

                if (cvoxel < 0) {
                    cvoxel = 0;
                } else {
                    cvoxel *= Cvoxel_conversionFactor;
                }

                //ivoxel = cvoxel * Ivoxel_conversionFactor;
            } else {
                cvoxel = 0;
                //ivoxel = 0;
            }

        }

        if (timerHigh) {

            if (fd.diffus[0].length == 1) { // single compartment

                fd.diffus[diffusantNum][0] += cvoxel;
                avgC = cvoxel;
                //avgI = ivoxel;

            } else {

                if (indexRFD == null) {
                    //Master.exit("Source error: coordinates index array has not be initialized.");
                    return 0;
                }

                numVoxels = indexRFD.length;

                for (int i = 0; i < numVoxels; i++) {

                    index = indexRFD[i];

                    fd.diffus[diffusantNum][index] += cvoxel;
                    avgC += cvoxel;
                    //avgI += ivoxel;

                }

                //if (numVoxels > 1) {
                    //avgC /= 1.0 * numVoxels;
                    //avgI /= 1.0 * numVoxels;
                //}

            }

        }
        
        saveValue(avgC);

        return avgC;

    }
    
    @Override
    public boolean check() {

        boolean ok = true;

        if (!super.check()) {
            return false;
        }

        if ((importFileName != null) && (importFileName.length() > 0)) {
            ok = readFile();
        }
        
        if (useProjectArray && project.sourceArray != null) {
            useArray = true;
        } else {
            useArray = (array != null);
        }
        
        return ok;
        
    }
    
    public boolean readFile() {

        if (importFileName.length() == 0) {
            return false;
        }

        if (useProjectArray) {

            if (project.sourceArray == null) {

                if (importFileType.equalsIgnoreCase("text")) {
                    project.sourceArray = Master.readTextFile(importFileName);
                } else if (importFileType.equalsIgnoreCase("binary")) {
                    project.sourceArray = Master.readBinaryFile(importFileName);
                }

                if (project.sourceArray == null) {
                    error("failed to read Source file: " + importFileName);
                }
                
            }

        } else {

            if (importFileType.equalsIgnoreCase("text")) {
                array = Master.readTextFile(importFileName);
            } else if (importFileType.equalsIgnoreCase("binary")) {
                array = Master.readBinaryFile(importFileName);
            }

            if (array == null) {
                error("failed to read Source file: " + importFileName);
            }

        }

        return true;

    }
    
}
