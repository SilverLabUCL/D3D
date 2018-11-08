package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;
import ucl.silver.d3d.utils.*;

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
public class InitMC_Demo extends InitProject {
    
    private final String directory = "/Jason/D3D/Simulations/";
    //private final String directory = "/Users/jason/Documents/D3D/Simulations/"
    private final String folder = "Testing";
    
    public double vesicle_radius = 0.025; // um
    public double Dshort = 0.060e-3; // um^2/ms // vesicle short-time diffusion constant
    
    public boolean FRAP_on = false;
    public double FRAP_kPhoto = 2.0; // photolysis k factor
    public double FRAP_onset = 100; // ms
    public double FRAP_bleach_length = 0.5; // ms
    public boolean FRAP_bleach_on = true; // bleach pulse
    public double FRAP_smallprobe_length = 2.0; // ms
    public double FRAP_smallprobe_amp_ratio = 0.0014; // ratio measured laser power after objective
    public boolean FRAP_smallprobe_on = true;

    public double excite_fwhm_xy = 0.343; // Gauss iPSF
    public double excite_fwhm_z = 1.285; // Gauss iPSF
    public double detect_fwhm_xy = 0.2551; // measured, bead deconvolved
    public double detect_fwhm_z = 0.9157; // measured, bead deconvolved

    public InitMC_GeometryWithMito_Jason igm;
    public PulseTimer timer = null;
    
    public InitMC_Demo(Project p) {
        super(p);
        igm = new InitMC_GeometryWithMito_Jason(p);
        String[] flist = {"initCichocki", "initFrapMovie", "initFrapMFT", "initActiveZone"};
        initFuncList = flist;
        project.finiteDifference = null; // turn off FD
        createVector(true);
    }

    @Override
    public boolean initFunction(String initSelect) {

        int i = initFunctionNum(initSelect);

        switch (i) {
            
            case 0:
                return initCichocki();
            case 1:
                return initFrapMovie();
            case 2:
                return initFrapMFT();
            case 3:
                return initActiveZone();
                
            default:
                error("initFunction", "initSelect", "failed to find init function");
                return true; // error
        }

    }
    
    private boolean initCichocki() {
        // Dynamic computer simulation of concentrated hard sphere suspensions: 
        // I. Simulation technique and mean square displacement data
        // B Cichocki K Hinsen (1990)
        // Physica A 166:473–491.
        // https://doi.org/10.1016/0378-4371(90)90068-4
        //
        // output is MSD(t)
        // compute D(t) = MSD/6t
        // compute D(t)/Dshort
        // steady-state value is Dlong/Dshort
        // run multiple reps (5-10) for each lambda and average
        //
        // Results for 0.3 vesicle density:
        // for lambda = 0.035, Dlong/Dshort = 0.
        // for lambda = 0.025, Dlong/Dshort = 0.
        // for lambda = 0.018, Dlong/Dshort = 0.
        // for lambda = 0, Dlong/Dshort = 0.481 (3-point extrap)
        // Cichocki Table I. Dlong/Dshort = 0.486

        RunMonteCarlo mc = project.newMonteCarlo();
        
        project.name = "Monte Carlo Cichocki";
        project.directory = directory;
        project.folder = "Cichocki";
        
        //double lambda = 0.035;
        double lambda = 0.025;
        //double lambda = 0.018;
        
        project.dx = 2 * vesicle_radius;
        project.set("printDT", 1);
        project.set("saveDT", 0.005);
        project.stability = 0.005;
        
        double cubeWidth = 1;
        
        geometry.resizeWithSpace(cubeWidth, cubeWidth, cubeWidth);

        double vesicleVolumeFraction = 0.30;
        //double vesicleVolumeFraction = 0.40; // Figure 1
        Dshort = 0.001; // um^2/ms

        DiffusantVesicles dv = new DiffusantVesicles(project, "Vesicles", 0, Dshort, null, vesicle_radius);

        dv.setVolumeFraction = vesicleVolumeFraction;
        dv.setImmobilePercent = 0;
        dv.saveXYZ = false;
        dv.saveMSD = true;
        
        int dnum = project.addDiffusant(dv);
        
        double vesicleStep = lambda * vesicle_radius;

        System.out.println("vesicle step = " + vesicleStep + " um");

        double t0 = vesicle_radius * vesicle_radius / ( 6 * Dshort);

        System.out.println("t0 = " + t0 + " ms");

        project.simTime = 16 * t0; // Figure 1
        project.simTime *= 10; // longer to reach steady-state
        
        mc.minVesicleStep = vesicleStep;
        mc.freeDiffusion = false;
        mc.PBC = true; // use periodic boundary conditions
        
        if (!Master.foundMainStartUpArguments) {
            mc.initAll(); // do not run when running batches
        }
        
        //Master.addBatchList(0, "geometry.cubeWidth", "1,1,1,1,1,1,1,1,1,1");
        
        return false;

    }
    
    private boolean initFrapMovie() {

        RunMonteCarloPhoto mc = project.newMonteCarloPhoto();
        
        project.name = "Monte Carlo FRAP";
        project.directory = directory;
        project.folder = folder;
        
        FRAP_onset = 2;
        project.simTime = FRAP_onset + 60;
        project.dx = 2 * vesicle_radius;
        project.set("printDT", 5);
        project.set("saveDT", 0.05);
        project.stability = 0.005;
  
        //double cubeWidth = 2.0;
        double cubeWidthxy = 1.1;
        double cubeWidthz = cubeWidthxy * 1.5; // longer in z-axis
        
        geometry.forceEvenVoxels = true;
        geometry.resizeWithSpace(cubeWidthxy, cubeWidthxy, cubeWidthz);
        
        double sdxy = Utility.gaussFWHM2STDV(excite_fwhm_xy);
        double sdz = Utility.gaussFWHM2STDV(excite_fwhm_z);

        PSFgauss iPSF = new PSFgauss(project, null, sdxy, sdxy, sdz);
        iPSF.xySymmetric = true;
        iPSF.set("xVoxelCenter", geometry.xVoxelCenter);
        iPSF.set("yVoxelCenter", geometry.yVoxelCenter);
        iPSF.set("zVoxelCenter", geometry.zVoxelCenter);
        iPSF.name = "Illumination PSF";
        
        int iPSF_index = geometry.addPSF(iPSF);

        timer = new PulseTimer(project, FRAP_onset, FRAP_bleach_length, 1.0);

        DiffusantVesiclesPhoto vs = new DiffusantVesiclesPhoto(project, "Vesicles", 0, Dshort, null, vesicle_radius, timer, iPSF);

        int dnum = project.addDiffusant(vs);
        
        vs.setVolumeFraction = 0.26; // Zoltan
        vs.setImmobilePercent = 0; //0.25
        vs.kPhoto = FRAP_kPhoto;
        vs.colorReady.setColor("[r=102,g=255,b=0]");
        vs.colorReady.setGradientColor("[r=255,g=0,b=0]");
        
        sdxy = Utility.gaussFWHM2STDV(detect_fwhm_xy);
        sdz = Utility.gaussFWHM2STDV(detect_fwhm_z);

        PSFgauss dPSF = new PSFgauss(project, null, sdxy, sdxy, sdz);
        dPSF.xySymmetric = true;
        dPSF.set("xVoxelCenter", geometry.xVoxelCenter);
        dPSF.set("yVoxelCenter", geometry.yVoxelCenter);
        dPSF.set("zVoxelCenter", geometry.zVoxelCenter);
        dPSF.name = "Detection PSF";
        
        DetectorPSF dw1 = Master.addDetectorPSF(dnum, dPSF);
        dw1.save.save2BinaryFile = true;
        dw1.save.ydim = "Fluorescence";
        dw1.save.autoDimensions = false;
        
        int dPSF_index = geometry.addPSF(dPSF);
        
        mc.minVesicleStep = 0.002;

        mc.frapOn = true;
        mc.saveFluorescence = false;
        mc.PBC = true;
        
        mc.setPSFSelect(iPSF_index, dPSF_index);
        
        mc.initAll();

        return false;

    }
    
    private boolean initFrapMFT() {
        // Physical determinants of vesicle mobility and supply at a central synapse
        // Rothman JS, Kocsis L, Herzog E, Nusser Z, Silver RA
        // Elife. 2016 Aug 19;5. pii: e15133. doi: 10.7554/eLife.15133.
        // Figures 3 and 4
        
        RunMonteCarloPhoto mc = project.newMonteCarloPhoto();
        
        project.name = "Monte Carlo MFT FRAP";
        project.directory = directory;
        project.folder = folder;

        double baselinetime = 500;

        FRAP_onset = baselinetime + 70; // first pulse is 70 ms before bleaching pulse
        project.simTime = 10000 + FRAP_onset;
        project.dx = 0.05;
        project.stability = 0.005;
        project.set("printDT", 200);
        project.set("saveDT", 0.05);
        
        double vesicle_width = 0.049; // um //  measured in MFT
        vesicle_radius = vesicle_width / 2.0;
        
        double cubeWidth = 3.0;

        igm.xdim = cubeWidth;
        igm.ydim = cubeWidth;
        igm.zdim = cubeWidth;
        
        igm.mitoNonSpaceVoxels = true; // mito are non-space voxels
        igm.mitoVolumeFraction = 0.28; // measured in MFT
        //igm.mito_ijkselect = -2; // orientation of mito chosen at random
        igm.mito_ijkselect = 0; // xy plane // Figure 3 orientation

        if (igm.initGeometryMC(-1)) {
            return true;
        }
        
        initPulseTimer();
        
        int iPSFselect = 0; // Gauss (fast to compute)
        //int iPSFselect = 1; // Torok (slow to compute) // this was used in eLife paper
        
        int iPSF_index;
        double sdxy, sdz;
        
        DiffusantVesiclesPhoto vs;
        
        switch (iPSFselect) {
            
            case 0:
                
                sdxy = Utility.gaussFWHM2STDV(excite_fwhm_xy);
                sdz = Utility.gaussFWHM2STDV(excite_fwhm_z);

                PSFgauss gPSF = new PSFgauss(project, null, sdxy, sdxy, sdz);

                gPSF.xySymmetric = true;
                gPSF.set("xVoxelCenter", geometry.xVoxelCenter);
                gPSF.set("yVoxelCenter", geometry.yVoxelCenter);
                gPSF.set("zVoxelCenter", geometry.zVoxelCenter);
                gPSF.name = "Illumination PSF";

                iPSF_index = geometry.addPSF(gPSF);
                
                vs = new DiffusantVesiclesPhoto(project, "Vesicles", 0, Dshort, null, vesicle_radius, timer, gPSF);

                break;

            case 1:

                PSFtorok tPSF = new PSFtorok(project, null);
                tPSF.setFitZoltan3();
                tPSF.xySymmetric = true;
                tPSF.set("xVoxelCenter", geometry.xVoxelCenter);
                tPSF.set("yVoxelCenter", geometry.yVoxelCenter);
                tPSF.set("zVoxelCenter", geometry.zVoxelCenter);
                tPSF.name = "Illumination PSF";
                tPSF.rscale = 1;

                iPSF_index = geometry.addPSF(tPSF);
                
                vs = new DiffusantVesiclesPhoto(project, "Vesicles", 0, Dshort, null, vesicle_radius, timer, tPSF);
                
                break;
                
            default:
                return true;
        }

        vs.colorReady.setColor("[r=147,g=210,b=21]");
        vs.colorImmobile.setColor("[r=150,g=150,b=150]");
        vs.color.setColor("[r=0,g=0,b=153]");
        vs.color.setGradientColor("[r=255,g=255,b=255]");
        vs.setVolumeFraction = 0.17; // measured in MFT
        vs.setImmobilePercent = 0.25; // estimated in MFT
        vs.kPhoto = FRAP_kPhoto;

        int dnum = project.addDiffusant(vs);
        
        sdxy = Utility.gaussFWHM2STDV(detect_fwhm_xy);
        sdz = Utility.gaussFWHM2STDV(detect_fwhm_z);

        PSFgauss dPSF = new PSFgauss(project, null, sdxy, sdxy, sdz);

        dPSF.xySymmetric = true;
        dPSF.set("xVoxelCenter", geometry.xVoxelCenter);
        dPSF.set("yVoxelCenter", geometry.yVoxelCenter);
        dPSF.set("zVoxelCenter", geometry.zVoxelCenter);
        dPSF.name = "Detection PSF";
        
        DetectorPSF d = Master.addDetectorPSF(dnum, dPSF);
        d.save.save2BinaryFile = true;
        d.save.save2TextFile = false;
        d.save.ydim = "Fluorescence";
        d.save.autoDimensions = false;
        
        int dPSF_index = geometry.addPSF(dPSF);

        //mc.minVesicleStep = 0.001;
        //mc.minVesicleStep = 0.002;
        mc.minVesicleStep = 0.005;

        mc.frapOn = true;
        mc.saveFluorescence = true;

        mc.setPSFSelect(iPSF_index, dPSF_index);

        mc.initAll();

        return false;

    }
    
    public void initPulseTimer() {

        if (!FRAP_on) {
            return;
        }

        if (FRAP_bleach_on) {
            timer = new PulseTimer(project, FRAP_onset, FRAP_bleach_length, 1.0);
        } else {
            timer = new PulseTimer(project, FRAP_onset, FRAP_bleach_length, 0);
        }

        if (FRAP_smallprobe_on) {
            timer.add(FRAP_onset - 70, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset - 40, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset - 10, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 20, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 50, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 80, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 110, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 150, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 200, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 260, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 310, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 390, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 490, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 810, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 1110, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 1710, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 2710, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 3910, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 4880, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 5910, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 6910, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 7910, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 8910, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
            timer.add(FRAP_onset + 9880, FRAP_smallprobe_length, FRAP_smallprobe_amp_ratio);
        }

    }
    
    public boolean initActiveZone() {

        RunMonteCarloAZ mc = project.newMonteCarloAZ();
        project.name = "Monte Carlo Active Zone";
        project.directory = directory;
        project.folder = folder;
        
        project.simTime = 1000;
        project.dx = 0.05;
        project.set("printDT", 100);
        project.set("saveDT", 0.1);

        //double cubeWidth = 2;
        double cubeWidth = 1.5;
        
        geometry.resizeWithSpace(cubeWidth, cubeWidth, cubeWidth);
        
        vesicle_radius = 0.049/2;

        DiffusantVesiclesAZ vs = new DiffusantVesiclesAZ(project, "Vesicles", 0, Dshort, null, vesicle_radius);

        project.addDiffusant(vs);

        vs.setVolumeFraction = 0.17;
        vs.setImmobilePercent = 0.25;
        vs.colorReady.setColor("[r=147,g=210,b=21]");
        vs.colorImmobile.setColor("[r=150,g=150,b=150]");
        vs.colorDocked.setColor("[r=255,g=102,b=0]");
        vs.colorReserve.setColor("[r=102,g=102,b=255]");

        //mc.minVesicleStep = 0.0020;
        mc.minVesicleStep = 0.0050;

        //mc.releaseProb = 1.0;
        //mc.releaseRate = 0.1; // kHz
        mc.releaseRate = Double.POSITIVE_INFINITY; // kHz
        //mc.releaseLimit = 50; // limit of number of release events
        //mc.releaseStartTime = 20;
        mc.releaseTimeOutLimit = 2000;

        //mc.replenishRate = Double.POSITIVE_INFINITY; // kHz
        mc.replenishFromAZdistance = 0.3; // um

        mc.dockingOn = true;

        mc.azWidth = 0.147; // cMFT // Zoltan EM

        mc.saveRelease = true;

        mc.initAll();

        return false;

    }
    
}
