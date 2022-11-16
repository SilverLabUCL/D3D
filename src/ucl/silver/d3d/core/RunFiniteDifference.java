package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;
import ucl.silver.d3d.gui.*;

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
public class RunFiniteDifference
        extends ParamVector implements Runnable {
    
    public double stability = 0.4; // stability factor for dt computation ( see below )
    // Note, stability is only for diffusion, not reactions
    // For simulations with reactions, stability factor may have to be reduced
    // Set stability = Double.NaN to stop dt computation
    
    public double maxD = 0; // max diffusion coefficient of all diffusants ( see maxD() )
    public double maxC0 = 0; // maximum concentration (mM) of diffusants ( see maxC0() )

    private transient Geometry geometry;
    public transient Grid grid; // Panel2D voxelGrid

    public transient double[][] diffus = null; // array of diffusant particles
    public transient double[][] diffnext = null; // array of diffusant particles next time step
    public transient double[][] psf = null; // psf of diffusant particles

    private transient double[][] temp; // temp array for switching arrays
    private transient byte[] sum; // sum of adjacent space voxels
    
    private int[] im, ip, jm, jp, km, kp;

    public boolean useCompactArrays = true; // use compact arrays containing only space voxels
    public boolean fastArrayAccess = true; // use fast array access, but this requires creating more arrays
    public boolean singleCompartment = false; // compute simulation using single compartment (i.e. no diffusion)

    public transient int thisPnt;
    private transient int xVoxels, yVoxels, zVoxels;

    public transient int it, itmax;
    public transient double time;

    private transient NativeC ntv = null;
    private transient Thread thrd;

    private transient boolean sourcesExist = false;
    private transient boolean reactionsExist = false;
    private transient boolean detectorsExist = false;
    private transient boolean batchesExist = false;

    private transient boolean loop = false;
    private transient boolean initialized = false;
    private transient boolean preview = false;
    public transient boolean finished = false;

    private transient StopWatch timer1 = null;
    private transient StopWatch timer2 = null;
    
    private transient final int indexLimitRFD = 10000;
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("maxD")) {
            return project.spaceUnits + "^2/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("maxC0")) {
            return project.concUnits;
        }
        return super.units(name);
    }
    
    @Override
    public String help(String name) {
        if (name.equalsIgnoreCase("stabilityFD")) {
            return "used to compute dt, dt=(stabilityFD*dx*dx/(3*maxD);";
        }
        if (name.equalsIgnoreCase("maxD")) {
            return "maximum diffusion constant, used to compute dt";
        }
        if (name.equalsIgnoreCase("maxC0")) {
            return "maximum initial concentration";
        }
        return super.help(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("maxD")) {
            return false;
        }
        if (name.equalsIgnoreCase("maxC0")) {
            return false;
        }
        if (name.equalsIgnoreCase("thisPnt")) {
            return false;
        }
        if (name.equalsIgnoreCase("time")) {
            return false;
        }
        if (name.equalsIgnoreCase("reactionsExist")) {
            return false;
        }
        if (name.equalsIgnoreCase("error")) {
            return false;
        }
        if (name.equalsIgnoreCase("finished")) {
            return false;
        }
        return true;
    }

    public RunFiniteDifference(Project p) {
        super(p);
    }

    public boolean init() {

        if (project.numBatches() > 0) {
            batchesExist = true;
        }
        
        maxC0();
        maxD();
        
        if ((maxD > 0) && Double.isFinite(stability) && (stability > 0)){
            project.dt = (stability * project.dx * project.dx) / (3 * maxD); // "Mathematics of Diffusion", J. Crank, p. 151, Eq. 8.50
        } else {
            project.dt = Double.NaN;
        }

        return false;

    }
    
    public double maxC0() {
        
        maxC0 = 0;

        if (project.diffusants == null) {
            return Double.NaN;
        }

        for (Diffusant d : project.diffusants) {
            if (!(d instanceof DiffusantParticles)) {
                maxC0 = Math.max(maxC0, d.C0);
            }
        }

        return maxC0;

    }

    public double maxD() {

        maxD = 0;

        if (project.diffusants == null) {
            return Double.NaN;
        }

        for (Diffusant d : project.diffusants) {
            if (!(d instanceof DiffusantParticles)) {
                maxD = Math.max(maxD, d.D);
            }
        }

        return maxD;

    }

    public void startSimulation(boolean PREVIEW, boolean startThread) {

        preview = PREVIEW;
        loop = true;

        if (startThread) {
            thrd = new Thread(this);
            thrd.start(); // this calls run() below
        }

    }

    public void pauseSimulation() {
        loop = false;
        Master.log("paused Finite Difference simulation.");
    }

    public void stopSimulation() {
        loop = false;
        it = 10 * itmax; // this stops simulation
    }

    public boolean initSimulation() {
        
        geometry = project.geometry;
        grid = Master.grid();

        timer1 = new StopWatch();
        timer2 = new StopWatch();

        timer1.start();

        reactionsExist = project.reactionExists();
        sourcesExist = project.sourceExists();
        detectorsExist = project.detectorExists();

        xVoxels = geometry.xVoxels;
        yVoxels = geometry.yVoxels;
        zVoxels = geometry.zVoxels;

        it = 0;
        time = 0;
        itmax = project.simPoints(); // total number of time steps

        //if (sourcesExist) {
            
            //for (Source s : project.sources) {
                //if (s == null) {
                    //continue;
                //}
                //if (project.sources[i].array == null) {
                    //continue;
                //}

                //pnts = itmax - project.sources[i].array.length;

                //if (pnts > 0) {
                    //Master.log("WARNING: SOURCE ARRAY IS MISSING TIME POINTS (n=" + pnts + ")");
                //}

            //}
        //}

        if (singleCompartment) {
            createArraysSingleCompartment();
        } else {
            createArrays();
        }
        
        if (preview && (grid != null)) {
            grid.preview = true;
            grid.diffusant = diffus;
            Master.mainframe.panel2D.displayModeSet("Diffusant.0"); // set display to "Diffusant"
        }

        finished = false;
        initialized = true;

        project.updateVectors();

        Master.log("initialized Finite Difference simulation: " + it);

        return false;

    }

    // compute sum, diffusant and psf integration arrays
    public boolean createArrays() {

        int count, vmax, I = 0;
        byte isum;
        boolean isSpace;

        Master.log("creating Finite Difference simulation arrays...");
        
        if (project.diffusants == null) {
            error("createArrays", "diffusants", "no diffusants");
            return true; // no FD diffusants
        }
        
        Diffusant[] diffusants = project.diffusants;

        if (useCompactArrays) {
            vmax = geometry.spaceVoxels;
        } else {
            vmax = xVoxels * yVoxels * zVoxels;
        }

        sum = new byte[vmax];
        diffus = new double[diffusants.length][vmax];
        diffnext = new double[diffusants.length][vmax];
        psf = new double[diffusants.length][vmax];

        if (fastArrayAccess) {
            im = new int[vmax];
            ip = new int[vmax];
            jm = new int[vmax];
            jp = new int[vmax];
            km = new int[vmax];
            kp = new int[vmax];
        }

        for (int k = 0; k < zVoxels; k++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int i = 0; i < xVoxels; i++) {

                    isSpace = geometry.isSpace(i, j, k);

                    if (!useCompactArrays || (useCompactArrays && isSpace)) {

                        if (!useCompactArrays) {
                            I = getI(i, j, k);
                        }

                        if (isSpace) {
                            geometry.setSpace(i, j, k, I); // set space value as pointer to new arrays
                        }

                        for (int ii = 0; ii < diffusants.length; ii++) {

                            diffus[ii][I] = diffusants[ii].getConcentration(i, j, k);

                            if (diffusants[ii].psf != null) {

                                psf[ii][I] = diffusants[ii].psf.getValue(i, j, k);

                                if ((i == xVoxels - 1) && (j == yVoxels - 1) &&
                                        (k == zVoxels - 1)) {
                                    diffusants[ii].psf.getValue(-1, -1, -1); // closes external file if they have been opened
                                }

                            }

                        }

                        isum = 0;

                        if (geometry.getSpace(i - 1, j, k) >= 0) {
                            isum++;
                        }

                        if (geometry.getSpace(i + 1, j, k) >= 0) {
                            isum++;
                        }

                        if (geometry.getSpace(i, j - 1, k) >= 0) {
                            isum++;
                        }

                        if (geometry.getSpace(i, j + 1, k) >= 0) {
                            isum++;
                        }

                        if (geometry.getSpace(i, j, k - 1) >= 0) {
                            isum++;
                        }

                        if (geometry.getSpace(i, j, k + 1) >= 0) {
                            isum++;
                        }

                        sum[I] = isum;

                        I++;

                    }

                }
            }
        }

        if (fastArrayAccess) {
            for (int k = 0; k < zVoxels; k++) {
                for (int j = 0; j < yVoxels; j++) {
                    for (int i = 0; i < xVoxels; i++) {
                        if (geometry.simpleCuboid || geometry.isSpace(i, j, k)) {
                            thisPnt = geometry.getSpace(i, j, k);
                            im[thisPnt] = geometry.getSpace(i - 1, j, k);
                            ip[thisPnt] = geometry.getSpace(i + 1, j, k);
                            jm[thisPnt] = geometry.getSpace(i, j - 1, k);
                            jp[thisPnt] = geometry.getSpace(i, j + 1, k);
                            km[thisPnt] = geometry.getSpace(i, j, k - 1);
                            kp[thisPnt] = geometry.getSpace(i, j, k + 1);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < diffusants.length; i++) {
            if (diffusants[i].psf != null) {
                Master.log("diffusant #" + i + ", computed PSF voxels: " +
                        diffusants[i].psf.computedValues);
            }
        }

        //if (project.sources != null) {
            //for (int i = 0; i < project.sources.length; i++) {
                //count = project.sources[i].coordinates().updateIndex();
                //if (count == 0) {
                //    Master.log("warning: 0 space voxels encountered in Source #" + i);
                //}
            //}
        //}
        
        sourcesUpdateIndex();

        if (project.detectors != null) {
            for (int i = 0; i < project.detectors.length; i++) {
                
                if (project.detectors[i].coordinates().spaceVoxels >= indexLimitRFD) {
                    continue;
                }
                
                if (project.detectors[i].psf != null) {
                    continue; // do not use indexRFD
                }
                
                count = project.detectors[i].coordinates().updateIndexRFD();
                
                if (count == 0) {
                    Master.log("warning: 0 space voxels encountered in Detector #" + i);
                }
                
            }
        }

        return false;

    }
    
    public void sourcesUpdateIndex() {

        int count;
        
        boolean noOverlaps = false;
        //boolean noOverlaps = true;

        if (project.sources == null) {
            return;
        }

        for (int is = 0; is < project.sources.length; is++) {

            if (project.sources[is] == null) {
                continue;
            }

            project.sources[is].indexRFD = null;
            
            if (project.sources[is].coordinates == null) {
                continue;
            }
            
            count = 0;

            for (CoordinatesVoxels c : project.sources[is].coordinates) {

                if (c == null) {
                    continue;
                }

                for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
                    for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                        for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                            if (geometry.isSpace(i, j, k) && project.sources[is].isInside(i, j, k)) {
                                if (noOverlaps) {
                                    if (!project.insideAnySource(i, j, k, is)) {
                                        count++;
                                    }
                                } else {
                                    count++;
                                }
                            }
                        }
                    }
                }

            }

            if (count == 0) {
                Master.log("warning: zero space voxels encountered in Source #" + is);
                continue;
            } else {
                //Master.log("FD source #" + is + " array init: space voxels = " + count);
            }
            
            // repeat again and save values to index array

            project.sources[is].indexRFD = new int[count];

            count = 0;
            
            for (CoordinatesVoxels c : project.sources[is].coordinates) {

                if (c == null) {
                    continue;
                }

                for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
                    for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                        for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                            if (geometry.isSpace(i, j, k) && project.sources[is].isInside(i, j, k)) {
                                if (noOverlaps) {
                                    if (!project.insideAnySource(i, j, k, is)) {
                                        project.sources[is].indexRFD[count] = geometry.space[i][j][k];
                                        count++;
                                    }
                                } else {
                                    project.sources[is].indexRFD[count] = geometry.space[i][j][k];
                                    count++;
                                }
                            }
                        }
                    }
                }

            }

        }

    }

    public boolean createArraysSingleCompartment() {

        int count = 0;
        
        if (project.diffusants == null) {
            error("createArrays", "diffusants", "no diffusants");
            return true; // no FD diffusants
        }
        
        Diffusant[] diffusants = project.diffusants;

        for (int i = 0; i < count; i++) {
            if (project.diffusants[i] instanceof DiffusantParticles) {
                continue; // skip MC
            }
            diffusants[i] = project.diffusants[i];
        }

        fastArrayAccess = false;

        diffus = new double[count][1];
        diffnext = new double[count][1];
        psf = new double[count][1];

        thisPnt = 0;

        Master.log("creating Single Compartment simulation arrays...");

        for (int ii = 0; ii < diffusants.length; ii++) {

            count = 0;

            for (int k = 0; k < zVoxels; k++) {
                for (int j = 0; j < yVoxels; j++) {
                    for (int i = 0; i < xVoxels; i++) {
                        if (geometry.simpleCuboid || geometry.isSpace(i, j, k)) {
                            diffus[ii][thisPnt] += diffusants[ii].getConcentration(i, j, k);
                            count++;
                        }
                    }
                }
            }

            if (count > 1) {
                diffus[ii][thisPnt] /= count;
            }

            if (diffusants[ii].psf != null) {

                count = 0;

                for (int k = 0; k < zVoxels; k++) {
                    for (int j = 0; j < yVoxels; j++) {
                        for (int i = 0; i < xVoxels; i++) {

                            if (geometry.simpleCuboid || geometry.isSpace(i, j, k)) {
                                psf[ii][thisPnt] += diffusants[ii].psf.getValue(i, j, k);
                                count++;
                            }

                            if ((i == xVoxels - 1) && (j == yVoxels - 1) && (k == zVoxels - 1)) {
                                diffusants[ii].psf.getValue(-1, -1, -1); // closes external file if they have been opened
                            }

                        }
                    }
                }

                if (count > 1) {
                    psf[ii][thisPnt] /= count;
                }

            }

        }

        for (int ii = 0; ii < diffusants.length; ii++) {
            if (diffusants[ii].psf != null) {
                Master.log("diffusant #" + ii + ", computed PSF voxels: " +
                        diffusants[ii].psf.computedValues);
            }
        }

        return false;

    }

    public int getI(int i, int j, int k) {
        return i + j * xVoxels + k * xVoxels * yVoxels; // 3D array in 1D
    }

    @Override
    public void run() {

        while (loop) { // loop thru batches

            if (!initialized) {

                if (project.simulationInit(preview)) {
                    stopSimulation(); // ERROR
                } else if (initSimulation()) {
                    stopSimulation(); // ERROR
                }

            }

            if (it < itmax) {
                Master.log("starting Finite Difference simulation at time: " + (time));
            }

            while (loop && (it < itmax)) {
                runDT();
            }

            if (it >= itmax) {
                finish();
            }

        }

    }
    
    public void runDT() {
        
        for (Diffusant d : project.diffusants) {
            if (!(d instanceof DiffusantParticles)) {
                d.save();
            }
        }

        if (sourcesExist) {
            for (Source s : project.sources) {
                s.release(this, geometry);
            }
        }

        if (detectorsExist) {
            for (Detector d : project.detectors) {
                d.detect(this, geometry);
            }
        }

        if (preview && (grid != null)) {
            grid.repaint();
        }

        if (singleCompartment) {
            diffuseReactSingleCompartment();
        } else if (fastArrayAccess) {
            diffuseReactFast();
        } else {
            diffuseReact();
        }
        
        it++; // next time step
        time = it * project.dt;
        timer2.timer(time);
                
    }

    public void diffuseReact() {

        double d;
        int otherPnt;

        for (int k = 0; k < zVoxels; k++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int i = 0; i < xVoxels; i++) {

                    if (geometry.simpleCuboid || geometry.isSpace(i, j, k)) {

                        thisPnt = geometry.space[i][j][k];

                        for (int ii = 0; ii < project.diffusants.length; ii++) { // diffusion
                            
                            if (project.diffusants[ii] instanceof DiffusantParticles) {
                                continue; // skip MC
                            }

                            if (project.diffusants[ii].h <= 0) {

                                diffnext[ii][thisPnt] = diffus[ii][thisPnt];

                            } else {
                                
                                d = 0;
                                otherPnt = geometry.getSpace(i - 1, j, k);

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = geometry.getSpace(i + 1, j, k);

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = geometry.getSpace(i, j - 1, k);

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = geometry.getSpace(i, j + 1, k);

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = geometry.getSpace(i, j, k - 1);

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = geometry.getSpace(i, j, k + 1);

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                diffnext[ii][thisPnt] = diffus[ii][thisPnt] +
                                        project.diffusants[ii].h *
                                        (d - sum[thisPnt] * diffus[ii][thisPnt]);

                            }
                        }

                        if (reactionsExist) {
                            for (int ii = 0; ii < diffus.length; ii++) {
                                project.diffusants[ii].react(this, ii);
                            }
                        }

                    }
                }
            }
        }

        // swap matrixes for next time step

        temp = diffus;
        diffus = diffnext;
        diffnext = temp;

    }

    public void diffuseReactFast() {

        double d;
        int otherPnt;

        for (int k = 0; k < zVoxels; k++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int i = 0; i < xVoxels; i++) {

                    if (geometry.simpleCuboid || geometry.isSpace(i, j, k)) {

                        thisPnt = geometry.space[i][j][k];

                        for (int ii = 0; ii < project.diffusants.length; ii++) { // diffusion
                            
                            if (project.diffusants[ii] instanceof DiffusantParticles) {
                                continue; // skip MC
                            }

                            if (project.diffusants[ii].h <= 0) {

                                diffnext[ii][thisPnt] = diffus[ii][thisPnt];

                            } else {

                                d = 0;
                                otherPnt = im[thisPnt];

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = ip[thisPnt];

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = jm[thisPnt];

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = jp[thisPnt];

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = km[thisPnt];

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                otherPnt = kp[thisPnt];

                                if (otherPnt >= 0) {
                                    d += diffus[ii][otherPnt];
                                }

                                diffnext[ii][thisPnt] = diffus[ii][thisPnt] +
                                        project.diffusants[ii].h *
                                        (d - sum[thisPnt] * diffus[ii][thisPnt]);

                            }
                        }

                        if (reactionsExist) {
                            for (int ii = 0; ii < diffus.length; ii++) {
                                project.diffusants[ii].react(this, ii);
                            }
                        }

                    }
                }
            }
        }

        // swap matrixes for next time step

        temp = diffus;
        diffus = diffnext;
        diffnext = temp;

    }

    public void diffuseReactSingleCompartment() {

        thisPnt = 0;

        for (int ii = 0; ii < diffus.length; ii++) {
            diffnext[ii][thisPnt] = diffus[ii][thisPnt]; // no diffusion
        }

        if (reactionsExist) {
            for (int ii = 0; ii < project.diffusants.length; ii++) {
                if (project.diffusants[ii] instanceof DiffusantParticles) {
                    continue; // skip MC
                }
                project.diffusants[ii].react(this, ii);
            }
        }

        temp = diffus;
        diffus = diffnext;
        diffnext = temp;

    }

    public void finish() {

        if (!initialized) {
            return;
        }

        timer1.stop();

        Master.log("finished Finite Difference simulation, time = " + timer1.toString() + " min");

        timer2.stop();

        project.simulationFinish();

        if (project.sources != null) {
            for (Source s : project.sources) {
                s.indexRFD = null;
            }
        }

        if (project.detectors != null) {
            for (Detector d : project.detectors) {
                d.coordinates().indexRFD = null;
            }
        }

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] instanceof DiffusantParticles) {
                continue; // skip MC
            }
            project.diffusants[i].saveConcentration(this, geometry, i); // save final concentration values
            project.diffusants[i].coordinates().indexRFD = null;
        }

        geometry.checkSpace(); // reset array values to -1 and +1

        diffus = null;
        diffnext = null;
        temp = null;
        sum = null;
        psf = null;

        im = null;
        ip = null;
        jm = null;
        jp = null;
        km = null;
        kp = null;

        if (grid != null) {
            grid.diffusant = null;
            grid.repaint();
        }

        finished = true;

        if (batchesExist) {
            initialized = false;
        } else {
            loop = false;
        }

        if (!loop && (grid != null)) {
            grid.preview = false;
        }

    }

// native C++ code
    public void basicFastNative() { // NO LONGER WORKING!!!

        boolean toFile = false;

        if (ntv == null) {
            ntv = new NativeC();
        }
        /*
        if (proj.detectFile.length() > 0) {
        toFile = true;
        }
        else {
        toFile = initOuts(); // create output arrays
        }
         */
        Master.log("executing NativeC basic integration...");

    //ntv.integrateBasic(proj, c, toFile);

    }

// native C++ code
    public void basicUncageNative() { // NO LONGER WORKING!!!

        boolean toFile = false;

        if (ntv == null) {
            ntv = new NativeC();
        }
        /*
        if (proj.detectFile.length() > 0) {
        toFile = true;
        }
        else {
        toFile = initOuts(); // create output arrays
        }
         */
        Master.log("executing NativeC uncaging integration...");

    //ntv.integrateUncage(proj, c, toFile);

    }

     // do not use this function directly, use ParamVector.set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof RunFiniteDifference)) {
            return false;
        }
        
        if (n.equalsIgnoreCase("stability")) {
            if (v <= 0) {
                return false;
            }
            stability = v;
            init();
            return true;
        }

        if (n.equalsIgnoreCase("useCompactArrays")) {
            useCompactArrays = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("fastArrayAccess")) {
            fastArrayAccess = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("singleCompartment")) {
            singleCompartment = (v == 1);
            return true;
        }

        return super.setMyParams(o, v);

    }

}
