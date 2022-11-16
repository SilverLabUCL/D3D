package ucl.silver.d3d.init;

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
public class InitFD_Demo extends InitProject {
    
    //private final String directory = "/Jason/D3D/Simulations/";
    //private final String directory = "/Users/jason/Documents/D3D/Simulations/";
    
    private final String folder = "Testing";
    
    public InitFD_Demo(Project p) {
        super(p);
        String[] flist = {"initCalciumChannel", "initGlutamateSource", "initGlutamateSourceCleft", "initCrankCylinder", "initCrankCylinderMovie", "initCrankPlane", "initGlomerulus_Nielsen2004", "initFrapAxelrod", "init_MFT_Scott_Rusakov"};
        initFuncList = flist;
        project.newFiniteDifference(); // turn on FD
        createVector(true);
    }

    @Override
    public boolean initFunction(String initSelect) {

        int i = initFunctionNum(initSelect);

        switch (i) {
            case 0:
                return initCalciumChannel();
            case 1:
                return initGlutamateSource();
            case 2:
                return initGlutamateSourceCleft();
            case 3:
                return initCrankCylinder();
            case 4:
                return initCrankCylinderMovie();
            case 5:
                return initCrankPlane();
            case 6:
                return initGlomerulus_Nielsen2004();
            case 7:
                return initFrapAxelrod();
            case 8:
                return init_MFT_Scott_Rusakov();
            default:
                error("initFunction", "initSelect", "failed to find init function");
                return true; // error
        }

    }
    
    private boolean initGlutamateSource() {
        //
        // To watch, go to 2D View, display Diffusant.0 and click Preview
        //
        double cubeWidth = 0.5;
        double xv1, yv1, zv1, xv2, yv2, zv2;

        project.name = "Cube with Glutamate Source";
        project.directory = directory;
        project.folder = folder;

        project.simTime = 2;
        project.dx = 0.01;
        project.set("printDT", 0.2);
        project.set("saveDT", 0.01);

        geometry.resizeWithSpace(cubeWidth, cubeWidth, project.dx);
        geometry.clear();

        xv1 = geometry.xVoxel1 + 1;
        yv1 = geometry.yVoxel1 + 1;
        zv1 = geometry.zVoxel2;
        xv2 = geometry.xVoxel1 + 1;
        yv2 = geometry.yVoxel2 - 1;
        zv2 = geometry.zVoxel2;

        coordinates.setVoxels(xv1, yv1, zv1, xv2, yv2, zv2);
        geometry.setSpace(coordinates, -1); // create a diffusion barrier

        double glut_C0 = 0;

        int idiff = Master.addDiffusant("glutamate", glut_C0, Utility.D_GLUTAMATE_37C); // add glutamate particle
        
        coordinates.setVoxelPoint(geometry.xVoxel1, geometry.yVoxelCenter, geometry.zVoxel1);
        coordinates.yVoxel2++;
        coordinates.update();
        
        int molecules = 4000; // # glutamate
        double tonset = 0.1; // ms
        double ctotal = Utility.N2C(molecules, project.dx, coordinates.voxels);
        double ctotal_final = Utility.N2C(molecules, project.dx, geometry.spaceVoxels); // end of simulation
        
        PulseTimer timer = new PulseTimer(project, tonset, 0, ctotal, "mM");
        
        timer.add(tonset + 1, 0, ctotal); // second pulse
        ctotal_final *= 2; // second pulse
        
        Source source = Master.addSource(idiff, coordinates, timer, "C");

        xv1 = geometry.xVoxel2;
        yv1 = geometry.yVoxel1;
        zv1 = geometry.zVoxel2;
        xv2 = geometry.xVoxel2;
        yv2 = geometry.yVoxel2;
        zv2 = geometry.zVoxel2;

        coordinates.setVoxels(xv1, yv1, zv1, xv2, yv2, zv2);

        Master.addDetectorAvg(idiff, coordinates);
        
        project.diffusants[idiff].displayMaxC = ctotal_final * 100;
        
        Master.log("start Ctotal=" + ctotal);
        Master.log("end Ctotal=" + ctotal_final);

        return false;

    }
    
    private boolean initGlutamateSourceCleft() {

        double xv1, yv1, zv1, xv2, yv2, zv2;

        int molecules = 4000; // # glutamate
        double release_onset = 0.5; // ms
        double release_duration = 0; // ms // impulse
        //double release_duration = 0.25; // ms
        
        Source source;
        PulseTimer timer;

        project.name = "Glutamate Source";
        project.directory = directory;
        project.folder = folder;

        project.simTime = 50;
        project.dx = 0.01;
        project.set("printDT", 5);
        project.set("saveDT", 0.01);
        
        double xwidth = 0.05;
        double ywidth = 1.5;

        geometry.resizeWithSpace(xwidth, ywidth, project.dx);
        geometry.clear();

        double glut_C0 = 0;

        int idiff = Master.addDiffusant("glutamate", glut_C0, Utility.D_GLUTAMATE_37C);

        xv1 = geometry.xVoxel1;
        yv1 = (int) geometry.yVoxelCenter;
        zv1 = geometry.zVoxel1;
        xv2 = geometry.xVoxel1;
        yv2 = yv1 + 1;
        zv2 = geometry.zVoxel1;
        
        coordinates.setVoxels(xv1, yv1, zv1, xv2, yv2, zv2);
        
        double ctotal = Utility.N2C(molecules, project.dx, coordinates.spaceVoxels);
        double ctotal_final = Utility.N2C(molecules, project.dx, geometry.spaceVoxels); // end of simulation
        
        timer = new PulseTimer(project, release_onset, release_duration, ctotal, "mM");
        timer.add(release_onset + 10, release_duration, ctotal); // second pulse
        source = Master.addSource(idiff, coordinates, timer, "C");
        //source.save.save2TextFile = true;
        
        ctotal_final *= 2; // second pulse
        
        Master.log("start Ctotal=" + ctotal);
        Master.log("end Ctotal=" + ctotal_final);
        
        double uptake_krate = 1; // 1/ms
        
        if (uptake_krate > 0) {
            source = Master.addSourceUptake(idiff, geometry, uptake_krate);
            //source.save.save2BinaryFile = true;
            //source.saveSelect = "pA";
        }

        xv1 = geometry.xVoxel2;
        yv1 = geometry.yVoxelCenter - 4;
        zv1 = geometry.zVoxel1;
        xv2 = geometry.xVoxel2;
        yv2 = geometry.yVoxelCenter + 5;
        zv2 = geometry.zVoxel2;

        coordinates.setVoxels(xv1, yv1, zv1, xv2, yv2, zv2);

        Master.addDetectorAvg(idiff, coordinates);
        
        if (false) {

            yv1 -= 40;
            yv2 -= 40;

            coordinates.setVoxels(xv1, yv1, zv1, xv2, yv2, zv2);

            Master.addDetectorAvg(idiff, coordinates); // spillover

        }
        
        //DetectorSnapshot snap = Master.addDetectorSnapshot(idiff, "xy", 0, release_onset + 0.2);
        
        project.diffusants[idiff].displayMaxC = ctotal_final * 100;

        return false;

    }

    private boolean initCalciumChannel() {

        double x0, y0, z0, source_width, cTotal = Double.NaN, cTotal_final;
        double cubeWidth = 1.0; // um

        project.name = "Cube with Calcium Source";
        project.directory = directory;
        project.folder = folder;

        project.simTime = 8; // ms
        project.dx = 0.02; // um
        project.set("printDT", 0.5);
        project.set("saveDT", 0.01);
        //project.set("saveRate", 99999);

        geometry.resizeWithSpace(cubeWidth, cubeWidth, cubeWidth);
        geometry.clear();

        double Ca_C0 = 0; // mM // initial calcium concentration

        int dnum = Master.addDiffusant("Ca", Ca_C0, Utility.D_CALCIUM_20C);

        x0 = geometry.xVoxelCenter;
        y0 = geometry.yVoxelCenter;
        z0 = geometry.zVoxel1;

        source_width = 4; // width in voxels

        coordinates.xyVoxelSquare(x0, y0, z0, source_width);
        //coordinates.setVoxelPoint(x0, y0, z0);

        double tpeak = 1; // ms
        double gaussSTDV = 0.2; // ms
        double gammaTau = 0.2; // ms
        double caCurrent = 2; // pA, peak
        
        //String sourceType = "gauss";
        String sourceType = "gamma";
        
        Source s = null;
        
        if (sourceType.equalsIgnoreCase("gauss")) {
            cTotal = Utility.gaussIpeak2Ctotal(caCurrent, gaussSTDV, caCurrent, coordinates.volume);
            s = Master.addSourceGauss(dnum, coordinates, tpeak, gaussSTDV, cTotal);
        } else if (sourceType.equalsIgnoreCase("gamma")) {
            cTotal = Utility.gammaIpeak2Ctotal(caCurrent, 2, 1.0 / gammaTau, caCurrent, coordinates.volume);
            s = Master.addSourceGamma(dnum, coordinates, tpeak, gaussSTDV, cTotal);
        }
        
        if (s != null) {

            s.save.save2BinaryFile = true;
            s.saveSelect = "pA";
            s.saveSelect = "mM/ms"; // Q

            cTotal_final = Utility.N2C(s.Ntotal, project.dx, geometry.spaceVoxels); // end of simulation

            Master.log("start Ctotal=" + cTotal);
            Master.log("end Ctotal=" + cTotal_final);

        }
        
        double uptake_krate = 1; // 1/ms
        
        if (uptake_krate > 0) {
            coordinates.matchVoxels(geometry);
            coordinates.zVoxel1 = coordinates.zVoxel2;
            s = Master.addSourceUptake(dnum, coordinates, uptake_krate);
            s.save.save2BinaryFile = true;
            s.saveSelect = "pA";
        }

        double detect_gauss_SDxy = 0.11;
        double detect_gauss_SDz = 0.36;

        x0 = geometry.xVoxelCenter;
        y0 = geometry.yVoxelCenter;
        z0 = geometry.zVoxelCenter;

        Master.addDetectorAvg(dnum);
        Master.addDetectorPSF_Gauss(dnum, detect_gauss_SDxy, detect_gauss_SDxy, detect_gauss_SDz, x0, y0, z0);
        
        geometry.resizeWithSpace(cubeWidth/2, cubeWidth/2, cubeWidth);
        // geometry is symmetric, so cut along symmetric lines

        return false;

    }

    private boolean initCrankCylinder() {
        // J Crank (1975) The Mathematics of Diffusion, Eq 2.6 or 2.7
        // Oxford UK: Clarendon Press
        // Simulation of a plane source in an infinite cylinder with unit cross section
        // Note, to match Eq 2.6 or 2.7, the output of detectors (mM) 
        // needs to be divided by the the volume of the source:
        // volume = dx * baseArea = dx * 1.0
        
        double D = 1; // um^2/s // diffusion coefficient
        double baseArea = 1.0; // um^2 // unit cross section
        double radius = Math.sqrt(baseArea / Math.PI);
        double diameter = 2 * radius;
        double length = 20; // long enough to mimic infinite length
        double dx = 0.05; // um
        
        Master.log("cylinder radius = " + radius + " um");
        
        //boolean halfLength = false; // eq. 2.6
        boolean halfLength = true; // eq. 2.7

        double source_cTotal = 1.0;

        int num_detectors = 10;
        double distance_between_detectors = 0.5;
        
        if (halfLength) {
            length *= 0.5;
        }

        int xVoxels = (int) (diameter / dx);
        int yVoxels = xVoxels;
        int zVoxels = (int) (length / dx);
        int k, dk;
        
        if (!halfLength) {
            zVoxels = Utility.makeOdd(zVoxels); // make odd, so source is centered
        }
        
        project.name = "Crank Cylinder Theoretical Solution";
        project.directory = directory;
        project.folder = folder;

        project.simTime = 2;
        project.dx = dx;
        project.set("printDT", 0.1);
        project.set("saveDT", 0.002);

        geometry.forceEvenVoxels = true;
        geometry.resizeWithSpace(xVoxels, yVoxels, zVoxels);
        geometry.cylinder();
        
        double zv0 = geometry.zVoxelCenter;

        if (halfLength) {
            zv0 = geometry.zVoxel2;
        }
        
        double xv1 = geometry.xVoxel1;
        double yv1 = geometry.yVoxel1;
        double xv2 = geometry.xVoxel2;
        double yv2 = geometry.yVoxel2;

        double C0 = 0; // mM
        
        int dnum = Master.addDiffusant("p", C0, D);

        coordinates.setVoxels(xv1, yv1, zv0, xv2, yv2, zv0);
        
        double source_spaceVoxels = coordinates.spaceVoxels;
        
        double releaseTime = 0.01;
        
        Master.addSourceImpulse(dnum, coordinates, releaseTime, source_cTotal, "mM");

        dk = (int) (distance_between_detectors / dx);

        k = (int) zv0;

        for (int kcnt = 0; kcnt < num_detectors; kcnt++) {
            coordinates.setVoxels(xv1, yv1, k, xv2, yv2, k);
            Master.addDetectorAvg(dnum, coordinates);
            k -= dk;
        }

        //Master.addDetectorAvg(dnum); // total concentration should be step function
        
        double cTotal_final = source_cTotal * source_spaceVoxels / geometry.spaceVoxels;
        
        Master.log("final C = " + cTotal_final);

        geometry.forceEvenVoxels = false;
        geometry.resizeWithSpace(diameter/2, diameter/2, -1);
        // because of symmetry, we can reduce simulation space
        
        return false;

    }

    private boolean initCrankCylinderMovie() {
        // To watch, go to 2D View, display Diffusant.0 and click Preview
        
        double source_ctotal = 1.0; // mM
        double D0 = 1; // um^2/ms
        double C0 = 0; // mM

        double simTime = 2.0; // ms

        int xVoxels = 141;
        int yVoxels = 30;
        int zVoxels = 1;

        project.name = "Crank Long Cylinder";
        project.directory = directory;
        project.folder = folder;

        project.simTime = simTime;
        project.dx = 0.02; // um
        project.set("printDT", 0.2);
        project.set("saveDT", 0.01);

        geometry.resizeWithSpace(xVoxels, yVoxels, zVoxels);

        int dnum = Master.addDiffusant("p", C0, D0);
        
        Diffusant d = project.diffusants[dnum];
        d.displayMinC = 0;
        d.displayMaxC = 0.1;
        
        double x0 = geometry.xVoxelCenter;
        double y1 = geometry.yVoxel1;
        double y2 = geometry.yVoxel2;

        coordinates.setVoxels(x0, y1, 0, x0, y2, 0);
        
        PulseTimer pt = new PulseTimer(project, 0.2, 0, source_ctotal, "mM");
        pt.add(0.7, 0, source_ctotal);
        pt.add(1.2, 0, source_ctotal);
        
        Master.addSource(dnum, coordinates, pt, "C");
        Master.addDetectorAvg(dnum, coordinates);

        return false;

    }
    
    private boolean initCrankPlane() {
        // J Crank (1975) The Mathematics of Diffusion, Eq 3.4
        // Oxford UK: Clarendon Press
        // Simulation of point source in infinite plane
        // Note, to match Eq 3.4, the output of detectors (mM) 
        // needs to be divided by the the xy-area of the source:
        // SA = dx * 2 * dx * 2
        // SA = 0.0016 um^2, for dx = 0.02
        
        double D = 1; // um^2/s // diffusion coefficient
        double width = 14; // wide enough to mimic infinite plane

        double source_cTotal = 1.0;

        int num_detectors = 10;
        double distance_between_detectors = 0.5;

        project.name = "Crank Plane with Point Source";
        project.directory = directory;
        project.folder = folder;

        project.simTime = 0.5;
        project.dx = 0.02;
        project.set("printDT", 0.1);
        project.set("saveDT", 0.002);
        
        int iWidth = (int) (width / project.dx);

        geometry.forceOddVoxels = true;
        geometry.resizeWithSpace(iWidth, iWidth, 1);
        
        double xv1 = geometry.xVoxelCenter;
        double yv1 = geometry.yVoxelCenter;
        double xv2 = xv1 + 1;
        double yv2 = yv1 + 1;

        double C0 = 0; // mM
        
        int dnum = Master.addDiffusant("p", C0, D);

        coordinates.setVoxels(xv1, yv1, 0, xv2, yv2, 0);
        // 2 x 2 voxels, allows cutting geometry to reduce simulation time
        
        double source_spaceVoxels = coordinates.spaceVoxels;
        
        double releaseTime = 0.01;
        
        Master.addSourceImpulse(dnum, coordinates, releaseTime, source_cTotal, "C");

        int di = (int) (distance_between_detectors / project.dx);

        int i = (int) xv1;

        for (int jcnt = 0; jcnt < num_detectors; jcnt++) {
            coordinates.setVoxels(i, yv1, 0, i+1, yv2, 0);
            // 2 x 2 voxels, allows cutting geometry to reduce simulation time
            Master.addDetectorAvg(dnum, coordinates);
            i -= di;
        }

        //Master.addDetectorAvg(dnum); // total concentration should be step function
        
        double cTotal_final = source_cTotal * source_spaceVoxels / geometry.spaceVoxels;
        
        Master.log("final C = " + cTotal_final);

        geometry.resizeWithSpace(width/2, width/2, -1);
        // because of symmetry, we can reduce simulation space
        
        return false;

    }
    
    private boolean initGlomerulus_Nielsen2004() {
        // Thomas A. Nielsen, David A. DiGregorio, and R. Angus Silver
        // Modulation of Glutamate Mobility Reveals the Mechanism Underlying Slow-Rising AMPAR EPSCs
        // Neuron, Vol. 42, 757–771, June 10, 2004
        // Figure 1B and 8G
        
        project.name = "MFT Glutamate Release";
        project.directory = directory;
        project.folder = folder;

        project.simTime = 8;
        project.dx = 0.02;
        project.finiteDifference.stability = 0.4;
        project.set("printDT", 1);
        project.set("saveDT", 0.01);
        //project.set("saveDT", 0.001);
        
        int nDenX = 12; // number of GC dendrites in x,y direction
        int nDenY = 12;
        
        double denWidth = 0.62; // um
        double denHeight = 2; // um
        double cleft = 0.02; // um
        
        double detectorWidth = 0.2; // um // post-synaptic AMPAR
        int idetectorWidthHalf = (int) Math.ceil(0.5 * detectorWidth / project.dx);
        
        int iDenWidth = (int) (denWidth / project.dx); // 31
        int iDenWidthHalf = (int) Math.ceil(0.5 * iDenWidth);
        int iDenHeight = (int) (denHeight / project.dx); // 100
        int iCleft = (int) (cleft / project.dx); // 1
        int kCleft = 0;
        
        double C0 = 0; // mM
        
        int dnum = Master.addDiffusant("glutamate", C0, Utility.D_GLUTAMATE_37C);

        GeometryGlomerulus.create(geometry, nDenX, nDenY, iDenWidth, iDenHeight, iCleft, iCleft);
        
        int molecules = 4000; // glutamate per vesicle
        double trelease = 0.01; // ms
        
        int iDenSource = 8;
        int i0 = (iDenSource - 1) * (iDenWidth + iCleft) + iDenWidthHalf;
        int j0 = i0;
        int k0 = kCleft;
        
        coordinates.xyVoxelSquare(i0, j0, k0, 1);
        
        double ctotal = Utility.N2C(molecules, project.dx, coordinates.spaceVoxels);
        
        Master.addSourceImpulse(dnum, coordinates, trelease, ctotal, "C");
        
        int i1 = i0 - idetectorWidthHalf;
        int j1 = i1;
        int i2 = i0 + idetectorWidthHalf;
        int j2 = i2;

        coordinates.xVoxel1 = i1;
        coordinates.yVoxel1 = j1;
        coordinates.zVoxel1 = kCleft;

        coordinates.xVoxel2 = i2;
        coordinates.yVoxel2 = j2;
        coordinates.zVoxel2 = kCleft;

        Master.addDetectorAvg(dnum, coordinates);
        
        for (int i = 1; i <= 3; i++) {
            coordinates.xVoxel1 = i1 - i * (iDenWidth + iCleft);
            coordinates.xVoxel2 = i2 - i * (iDenWidth + iCleft);
            Master.addDetectorAvg(dnum, coordinates);
        }
        
        coordinates.yVoxel1 = i1 - 1 * (iDenWidth + iCleft);
        coordinates.yVoxel2 = i2 - 1 * (iDenWidth + iCleft);
        
        for (int i = 1; i <= 3; i++) {
            coordinates.xVoxel1 = i1 - i * (iDenWidth + iCleft);
            coordinates.xVoxel2 = i2 - i * (iDenWidth + iCleft);
            Master.addDetectorAvg(dnum, coordinates);
        }
        
        coordinates.yVoxel1 = i1 - 2 * (iDenWidth + iCleft);
        coordinates.yVoxel2 = i2 - 2 * (iDenWidth + iCleft);
        
        for (int i = 2; i <= 3; i++) {
            coordinates.xVoxel1 = i1 - i * (iDenWidth + iCleft);
            coordinates.xVoxel2 = i2 - i * (iDenWidth + iCleft);
            Master.addDetectorAvg(dnum, coordinates);
        }
        
        coordinates.yVoxel1 = i1 - 3 * (iDenWidth + iCleft);
        coordinates.yVoxel2 = i2 - 3 * (iDenWidth + iCleft);

        coordinates.xVoxel1 = i1 - 3 * (iDenWidth + iCleft);
        coordinates.xVoxel2 = i2 - 3 * (iDenWidth + iCleft);
        
        Master.addDetectorAvg(dnum, coordinates);

        return false;

    }
    
    private boolean initFrapAxelrod() {
        // Biophys J. 1976 Sep;16(9):1055-69.
        // Mobility measurement by analysis of fluorescence photobleaching recovery kinetics.
        // Axelrod D, Koppel DE, Schlessinger J, Elson E, Webb WW.
        
        // Parameters from Rothman et al 2016, but with no small probe pulses
        // Elife. 2016 Aug 19;5. pii: e15133. doi: 10.7554/eLife.15133.
        // Physical determinants of vesicle mobility and supply at a central synapse.
        // Rothman JS, Kocsis L, Herzog E Nusser Z, Silver RA.
        // See Figure 2C
        
        double FRAP_onset = 100; // ms
        double FRAP_bleach_length = 0.5; // ms
        
        project.name = "Finite Difference FRAP";
        project.directory = directory;
        project.folder = folder;

        project.simTime = FRAP_onset + 10000;
        project.dx = 0.05;
        project.finiteDifference.stability = 0.0017;
        project.set("printDT", 500); // ms
        project.set("saveDT", 0.1); // ms
        
        double cubeWidth = 3; // um
        
        double psf_fwhm_xy = 0.27; // measured PSF
        double FRAP_kPhoto = 2.0; // produces f = 0.62, or K = 1 in Axelrod Figure 2
        
        //boolean eighth_geometry = false;
        boolean eighth_geometry = true;
        
        geometry.forceEvenVoxels = true;
        geometry.resizeWithSpace(cubeWidth, cubeWidth, cubeWidth);
        
        double psf_stdv_xy = Utility.gaussFWHM2STDV(psf_fwhm_xy);
        double psf_stdv_z = -1; // column in z-axis

        PSFgauss iPSF = new PSFgauss(project, null, psf_stdv_xy, psf_stdv_xy, psf_stdv_z);

        if (eighth_geometry) {
            iPSF.arbitraryVoxelCenter = true;
            iPSF.xySymmetric = false;
        } else {
            iPSF.arbitraryVoxelCenter = false;
            iPSF.xySymmetric = true;
        }

        iPSF.set("xVoxelCenter", geometry.xVoxelCenter);
        iPSF.set("yVoxelCenter", geometry.yVoxelCenter);
        iPSF.set("zVoxelCenter", geometry.zVoxelCenter);
        iPSF.name = "Illumination PSF";
        
        PulseTimer timer = new PulseTimer(project, FRAP_onset, FRAP_bleach_length, FRAP_kPhoto, "1/ms");
        // one bleaching pulse, no small probe pulses
        
        double C0 = 1.0; // initial fluorescence (which is concentration)
        double Dlong = 0.025e-3; // um^2/ms // effective diffusion constant of vesicles
        int dnum = -1; // not used

        DiffusantPhoto vesicles = new DiffusantPhoto(project, "vesicles", C0, Dlong, null, dnum, timer, iPSF);
        //vesicles.save.save2BinaryFile = true; // for saving timer/kPhoto

        int iVesicles = project.addDiffusant(vesicles);
        
        PSFgauss dPSF = new PSFgauss(project, null, psf_stdv_xy, psf_stdv_xy, psf_stdv_z);
        
        if (eighth_geometry) {
            dPSF.arbitraryVoxelCenter = true;
            dPSF.xySymmetric = false;
        } else {
            dPSF.arbitraryVoxelCenter = false;
            dPSF.xySymmetric = true;
        }
        
        dPSF.set("xVoxelCenter", geometry.xVoxelCenter);
        dPSF.set("yVoxelCenter", geometry.yVoxelCenter);
        dPSF.set("zVoxelCenter", geometry.zVoxelCenter);
        dPSF.name = "Detection PSF";
        
        DetectorPSF d = Master.addDetectorPSF(iVesicles, dPSF);
        d.save.save2BinaryFile = true;
        d.save.ydim = "Norm Fluorescence";
        d.save.autoDimensions = false;
        
        if (eighth_geometry) {
            geometry.resizeWithSpace(cubeWidth / 2, cubeWidth / 2, cubeWidth / 2);
            // because of symmetry, we can reduce simulation space
        }
        
        double w = 2 * psf_stdv_xy; // Gauss radius at exp(-2) 
        double tauD = w * w / (4 * Dlong); // characteristic diffusion time
        double gammaD = 1.1; // Axelrod et al 1976, Figure 7
        double tauHalf = tauD * gammaD;
        
        Master.log("FRAP recovery tauD = " + tauD + "ms");
        Master.log("FRAP recovery tauHalf = " + tauHalf + "ms");

        return false;

    }
    
    private boolean init_MFT_Scott_Rusakov() {
        // Main Determinants of Presynaptic Ca2+ Dynamics at Individual Mossy Fiber–CA3 Pyramidal Cell Synapses
        // Ricardo Scott and Dmitri A. Rusakov
        // Journal of Neuroscience 28 June 2006, 26 (26) 7071-7081
        // DOI: https://doi.org/10.1523/JNEUROSCI.0946-06.2006
        //
        // Simulations do not consider F due to Alexa Fluor 594 (Rct).
        // Hence, simulation dF/F0 values are much larger than reported.
        // Concentrations agree.

        double MFBvolume = 10; // um^3 // giant MFT without axon // changing this does not change results
        double dx = Math.pow(MFBvolume, 1.0 / 3.0); // single voxel for entire geometry
        
        String fig = "D"; // Figure 4D
        //String fig = "E";
        //String fig = "F";

        project.directory = directory;
        project.folder = folder;
            
        project.simTime = 1000;
        project.dx = dx;
        project.finiteDifference.stability = 0.0001;
        project.set("printDT", 50); // ms
        project.set("saveDT", 0.01); // ms
        
        double Ca_C0 = 110e-6; // mM
        double Ca_D = 0.223; // um^2/s // but there is no diffusion in simlation
        
        int iCa = Master.addDiffusant("Ca", Ca_C0, Ca_D);
        
        String dye_name;
        double dye_Ctotal, dye_D, dye_Kon, dye_Koff, dye_R;
        
        if (fig.equalsIgnoreCase("D")) {
            dye_name = "Fluo4";
            dye_Ctotal = 0.2;
            dye_D = 0.02; // but no difffusion
            dye_Kon = 600;
            dye_Koff = 0.21;
            dye_R = 110; // ?????
        } else if (fig.equalsIgnoreCase("E")) {
            dye_name = "Fluo4";
            dye_Ctotal = 0.05;
            dye_D = 0.02;
            dye_Kon = 600;
            dye_Koff = 0.21;
            dye_R = 110; // ?????
        } else if (fig.equalsIgnoreCase("F")) {
            dye_name = "Fluo5";
            dye_Ctotal = 0.2;
            dye_D = 0.02;
            dye_Kon = 300;
            dye_Koff = 0.3;
            dye_R = 110; // ?????
        } else {
            error("init_MFT_Scott_Rusakov", "fig", "unknown figure");
            return true;
        }
        
        int iDye = Master.addDiffusantDye(dye_name, dye_Ctotal, dye_D, iCa, dye_Kon, dye_Koff, dye_R);

        double CB_Ctotal = 0.16; // mM // 4 * 0.04 mM (4 binding sites) // Muller et al 2005
        double CB_D = 0.02; // but no diffusion
        double CB_Kon = 30;
        double CB_Koff = 0.009;
            
        int iCB = Master.addDiffusantReactant("Calbindin", CB_Ctotal, CB_D, iCa, CB_Kon, CB_Koff);

        geometry.resizeWithSpace(1, 1, 1); // single voxel
        geometry.update();
        
        double source_Ctotal;
        
        if (fig.equalsIgnoreCase("D")) {
            source_Ctotal = 0.054; // mM
        } else if (fig.equalsIgnoreCase("E")) {
            source_Ctotal = 0.043; // mM
        } else if (fig.equalsIgnoreCase("F")) {
            source_Ctotal = 0.050; // mM
        } else {
            return true;
        }

        double source_tfwhm = 1.1774; // ms
        double source_tstdv = Utility.gaussFWHM2STDV(source_tfwhm); // STDV = 0.5 ms
        // this seems to be error in paper, since they say FWHM = 0.5, but then say sigma (STDV) = 0.5
        // i.e. FWHM and STDV were mixed up.
        // best match to their simulations is when FWHM = 1.1774, STDV = 0.5
        // but authors probably wanted FWHM = 0.5, STDV = 0.21233
        // the difference to simulation is small and does not effect long time scales
        
        double source_tpeak = 150; // ms
        int source_numPulses = 5;
        double source_pulseInterval = 50; // ms // 20 Hz
        
        PulseTimer pt = new PulseTimer(project, source_tpeak, source_tstdv, source_Ctotal, "mM", source_numPulses, source_pulseInterval);

        SourceGauss s = Master.addSourceGauss(iCa, geometry, pt);
        s.saveSelect = "pA";
        s.save.save2TextFile = true;
        
        double uptake_rate;
        
        if (fig.equalsIgnoreCase("D")) {
            uptake_rate = -0.36; // 1/ms
        } else if (fig.equalsIgnoreCase("E")) {
            uptake_rate = -0.43; // 1/ms
        } else if (fig.equalsIgnoreCase("F")) {
            uptake_rate = -0.33; // 1/ms
        } else {
            return true;
        }
        
        SourceUptake su = Master.addSourceUptake(iCa, geometry, uptake_rate);
        
        su.saveSelect = "pA";
        su.save.save2BinaryFile = true;
        
        DetectorDye dd = new DetectorDye(project, Master.detectorName(), iDye, null, null);
        project.addDetector(dd);
        
        Master.addDetectorAvg(iCa);
        Master.addDetectorAvg(iCB);

        project.init();

        return false;

    }
    
}
