package ucl.silver.d3d.core;

import ucl.silver.d3d.gui.*;
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
public class RunMonteCarlo
        extends ParamVector implements Runnable {
    
    public double maxD = 0; // max diffusion coefficient of all diffusants ( see maxD() )

    public double minParticleRadius; // um
    public double maxParticleRadius; // um
    public double minAllowedDistanceBetweenParticles; // um
    public double minDistanceBetweenParticles; // um

    public int numParticles; // total number of particles
    public double particleDensity; // particles / um^3
    public double particleVolumeFraction;
    public double totalParticleVolume; // um^3

    public boolean initParticleRandom = true;
    public boolean particleLattice = false;

    public double Dcyto = 0; // um^2/ms

    public double voxelVolume = 0; // um^3

    public int mobileParticles = 0;
    public int immobileParticles = 0;
    public double mobilePercent = 0;
    public double immobilePercent = 0;
    public double mobileVolumeFraction = 0;
    public double immobileVolumeFraction = 0;
    public double clusterImmobilePercent = 0;
    public boolean clusterImmobileParticles = false;

    public boolean connectParticles = false;
    public boolean connectorsBinomial = false;
    public int maxNumConnectors = 5;
    public double meanNumConnectors = 1.5;
    public double connectorLength = 0.01; // um
    public double connectRate = 0.7 / 0.1; // P/ms
    public double unconnectRate = 0.175 / 0.1; // P/ms
    public double avgConnectorLifeTime = 0;
    public int connectorLifeTimeCounter = 0;
    public double connectorAllOffTime = 0; // ms

    public double minParticleStep; // um
    public double removeParticleOverlapStep3 = 0.001 / Math.sqrt(3); // um

    //public boolean checkForOverlapsDuringSimulation = false; // set true for testing only!
    
    public boolean PBC = false; // periodic boundary conditions
    public boolean freeDiffusion = false; // particles are allowed to overlap
    public boolean removeParticleOverlap = true;

    public boolean driftOn = false;
    public double driftRateX = 0; // um/ms
    public double driftRateY = 0; // um/ms
    public double driftRateZ = 0; // um/ms
    public double driftOnset = 0; // ms
    private double driftx, drifty, driftz; // um

    // simulation variables
    public transient Grid grid = null; // Panel2D voxel grid
    public transient DiffusantParticles[] diffusants = null;
    
    public transient Thread thread = null;
    public transient Geometry geometry;

    public int itime;
    public double time;
    //public boolean timer = true;
    public double stepx, stepy, stepz;

    public boolean batchesExist = false;
    public boolean autoInit = true;
    public boolean runSimulation = false;
    public boolean cancel = false;
    public boolean initialized = false;
    public boolean preview = false;
    public boolean removingParticleOverlap = false;

    private boolean clusterImmobileOn = false;
    private long clusterImmobileCounter = 0;

    public boolean hydrodynamicsLocalD = false; // hydrodynamic interactions local density computation
    public boolean hydrodynamicsDscale = false; // hydrodynamic interactions scale factor
    public boolean hydroWallZ = false;
    public boolean hydrodynamicsLocalDVoxels = false;

    public boolean sortParticlesByRadius = false;
    public long initParticleRandomTrialLimit = (long) 1e7;
    public long removeParticleOverlapTrialLimit = (long) 1e7;

    public transient StopWatch timer1 = new StopWatch();
    public transient StopWatch timer2 = new StopWatch();

    public boolean saveResidenceTime = false;
    public double residenceTimeStart = 0; // time to start computing mean residence time
    public CoordinatesVoxels[] residenceTimeCoordinates = null;

    public boolean saveMSDspatial = false;
    private Save[] MSDspatial = null;
    public double MSDspatialTbgn = -1;
    private double MSDspatial_t1 = -1, MSDspatial_t2 = -1;
    public double MSDspatial_win = 100; // ms
    public int MSDspatial_numBins = 9;
    public double MSDspatial_binWidth = 0.05; // um

    private final DiffusantParticle testParticle = new DiffusantParticle(project, "test", 0, 0, 0, 0, 0);

    public double particleVolume = Double.NaN; // used for local density computation, but assumes all particles have same volume

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("maxD")) {
            return project.spaceUnits + "^2/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("minParticleRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxParticleRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minAllowedDistanceBetweenParticles")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minDistanceBetweenParticles")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minParticleStep")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("totalParticleVolume")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("particleDensity")) {
            return "particles/" + project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("Dcyto")) {
            return project.spaceUnits + "^2/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("driftRateX")) {
            return project.spaceUnits + "/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("driftRateY")) {
            return project.spaceUnits + "/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("driftRateZ")) {
            return project.spaceUnits + "/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("driftOnset")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("time")) {
            return project.timeUnits;
        }
        return super.units(name);
    }
    
    @Override
    public String help(String name) {
        if (name.equalsIgnoreCase("maxD")) {
            return "maximum diffusion constant, used to compute dt";
        }
        return super.help(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("maxD")) {
            return false;
        }
        if (name.equalsIgnoreCase("minParticleRadius")) {
            return false;
        }
        if (name.equalsIgnoreCase("maxParticleRadius")) {
            return false;
        }
        if (name.equalsIgnoreCase("minAllowedDistanceBetweenParticles")) {
            return false;
        }
        if (name.equalsIgnoreCase("minDistanceBetweenParticles")) {
            return false;
        }
        if (name.equalsIgnoreCase("totalParticleVolume")) {
            return false;
        }
        if (name.equalsIgnoreCase("particleDensity")) {
            return false;
        }
        if (name.equalsIgnoreCase("particleVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("numParticles")) {
            return false;
        }
        if (name.equalsIgnoreCase("mobileParticles")) {
            return false;
        }
        if (name.equalsIgnoreCase("immobileParticles")) {
            return false;
        }
        if (name.equalsIgnoreCase("mobilePercent")) {
            return false;
        }
        if (name.equalsIgnoreCase("immobilePercent")) {
            return false;
        }
        if (name.equalsIgnoreCase("mobileVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("immobileVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("time")) {
            return false;
        }
        if (name.equalsIgnoreCase("batchesExist")) {
            return false;
        }
        if (name.equalsIgnoreCase("runSimulation")) {
            return false;
        }
        if (name.equalsIgnoreCase("cancel")) {
            return false;
        }
        if (name.equalsIgnoreCase("initialized")) {
            return false;
        }
        if (name.equalsIgnoreCase("preview")) {
            return false;
        }
        return super.canEdit(name);
    }

    public RunMonteCarlo(Project p) {

        super(p);
        
        createVector(true);

    }

    public void init() {
        if (project.numBatches() > 0) {
            batchesExist = true;
        }
    }

    public void setOutputRate(double newRate) {
    }

    public void startSimulation(boolean PREVIEW, boolean startThread) {

        preview = PREVIEW;

        grid = Master.grid();

        if ((grid != null) && preview) {
            grid.preview = true;
            //voxelGrid.diffusant = diffus;
            Master.mainframe.panel2D.displayModeSet("Diffusant.0");
        }

        runSimulation = true;
        cancel = false;
        
        if (startThread) {
            thread = new Thread(this);
            thread.start(); // this calls run() below
        }

    }

    public void pauseSimulation() {

        runSimulation = false;

        timer1.stop();
        timer2.stop();

        Master.log("paused Monte Carlo simulation at time " + time + " " + project.timeUnits);

    }

    public void cancelSimulation() {
        runSimulation = false;
        cancel = true;
        Master.log("cancelled Monte Carlo simulation at time " + time + " " + project.timeUnits);
    }

    public boolean initSimulation() {

        if (checkVariables()) {
            return true;
        }

        if (autoInit) {
            if (initAll()) {
                return true;
            }
        }

        itime = 0;
        time = 0;

        if (driftOn) {
            driftx = driftRateX * project.dt;
            drifty = driftRateY * project.dt;
            driftz = driftRateZ * project.dt;
        }

        timer1.start();

        initialized = true;

        printParameters();

        initSave();

        Master.log("initialized Monte Carlo simulation");

        return false;

    }

    public void printParameters() {

        double voxels = geometry.voxels;
        double spaceVoxels = geometry.spaceVoxels;

        Master.log("Geometry xWidth = " + geometry.xWidth);
        Master.log("Geometry yWidth = " + geometry.yWidth);
        Master.log("Geometry zWidth = " + geometry.zWidth);

        Master.log("Geometry % nonspace = " + (100 * (voxels - spaceVoxels) / voxels));

        int j = 0;

        for (DiffusantParticles d : diffusants) {
            Master.log("DiffusantParticles #" + j + ", avg D = " + d.D);
            Master.log("DiffusantParticles #" + j + ", estimated meanStep = " + Math.sqrt(6.0 * d.D * project.dt));
            Master.log("DiffusantParticles #" + j + ", avg meanStep = " + (d.meanStep3 * Math.sqrt(3)));
            Master.log("DiffusantParticles #" + j + ", volume fraction = " + d.volumeFraction);
            Master.log("DiffusantParticles #" + j + ", immobilePercent = " + d.immobilePercent);
            j++;
        }

    }

    public void finishSimulation() {

        finishSave();
        particleStats();

        if (connectParticles) {
            avgConnectorLifeTime /= connectorLifeTimeCounter;
            Master.log("average connector life time (ms): " + avgConnectorLifeTime + " (n=" + Integer.toString(connectorLifeTimeCounter) + ")");
        }

    }

    @Override
    public void run() {

        double simTime = project.simTime;

        while (runSimulation) { // runSimulation thru batches

            if (!initialized) {

                if (project.simulationInit(preview)) {
                    cancelSimulation(); // ERROR
                }

                saveParticlePositions(0);

            }

            if (time < simTime) {
                Master.log("starting Monte Carlo simulation at time " + time + " " + project.timeUnits);
            }

            while (runSimulation && !cancel && (time < simTime)) {

                for (DiffusantParticles d : diffusants) {
                    d.save(); // save diffusant variables
                }

                runMSDspatial();

                //if (project.detectors != null) { // right now there are no MC detectors
                //    for (int d = 0; d < project.detectors.length; d++) {
                //        if (project.detectors[d] instanceof DetectorSnapshot) {
                //            project.detectors[d].detect(this, geometry); // compute detector averages
                //        }
                //    }
                //}
                
                react();

                moveParticlesCichocki(false);
                //testConnectorSeperation();

                if (hydrodynamicsLocalD) {
                    localDensityAll(false);
                }

                if (connectParticles) {
                    connectParticles();
                    unconnectParticles();
                }

                if (clusterImmobileOn) {
                    clusterImmobileParticles();
                }
                
                if (driftOn && (time > driftOnset)) {
                    drift();
                }

                //if (!freeDiffusion && checkForOverlapsDuringSimulation && testParticleOverlap()) {
                //    Master.log("Monte Carlo Simulation has overlapping particles: " + minDistanceBetweenParticles);
                //    runSimulation = false;
                //}
                
                if (saveResidenceTime && (time >= residenceTimeStart)) {
                    updateResidentTimes();
                }

                if ((grid != null) && preview) {
                    grid.repaint();
                }

                itime += 1;
                time += project.dt;
                timer2.timer(time);

            }

            if (cancel || (time >= simTime)) {
                finish();
            }

        }

    }
    
    public void react() {
        // put any reactions here
    }

    public void finish() {

        if (!initialized) {
            return;
        }

        project.simulationFinish();

        saveParticlePositions((int) project.simTime);

        //avgFirstCollision();
        if (saveResidenceTime) {
            saveVoxelResidenceTime();
        }

        //for (int d = 0; d < diffusants.length; d++) {
        //    diffusants[d].saveConc(this, d); // save final concentration values
        //}
        //voxelGrid.diffusant = null;
        if (grid != null) {
            grid.repaint();
        }

        if (batchesExist) {
            initialized = false;
        } else {
            runSimulation = false;
        }

        if ((grid != null) && (!runSimulation)) {
            grid.preview = false;
        }

        timer1.stop();
        timer2.stop();
        Master.log("finished Monte Carlo simulation. time = " + timer1.toString());

    }
    
    public boolean initGeometry() {
        return false; // generic method to initialize geometry
    }

    public boolean initAll() {

        if (checkVariables()) {
            return true;
        }

        if (initDiffusantParticlesArray()) {
            return true;
        }

        if (initDetectors()) {
            return true;
        }

        if (initVoxels()) {
            return true;
        }

        if (initGeometry()) {
            return true;
        }

        if (initParticles()) {
            return true;
        }

        if (checkDX()) {
            return true;
        }

        if (initDT()) {
            return true;
        }

        if (particleLattice) {
            if (initVParticlesLattice()) {
                return true;
            }
        } else if (initParticleRandom) {
            if (initParticlesRandom()) {
                return true;
            }
        }

        if (initParticlesImmobile()) {
            return true;
        }

        if (removeParticleOverlap && removeParticleOverlap()) {
            return true;
        }

        initParticleStartLocations();

        initConnectors();

        particleStats();

        if (clusterImmobileParticles) {
            initClusterImmobileCounter();
        }

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        if (maxNumConnectors == 0) {
            connectParticles = false;
        }
        
        Master.updatePanel2D();

        return false;

    }

    public boolean checkVariables() {

        //if (PBC && !shape.simpleCuboid) {
        //    Master.log("init error: PBC allowed only with simple cuboid geometries");
        //    error = true;
        //     return true;
        //}
        geometry = project.geometry;

        return false;

    }

    public boolean initDiffusantParticlesArray() {

        int j = 0, count = 0;

        if (project.diffusants == null) {
            error("No Diffusant Particles!");
            return true;
        }

        for (Diffusant d : project.diffusants) {
            if (d instanceof DiffusantParticles) {
                count++;
            }
        }

        if (count == 0) {
            error("initDiffusantParticlesArray", "diffusants", "no diffusant particles");
            return true;
        }

        diffusants = new DiffusantParticles[count];

        for (Diffusant d : project.diffusants) {

            if (d instanceof DiffusantParticles) {
                diffusants[j++] = (DiffusantParticles) d;
            }
        }
        
        Master.log("initialized diffusant particle arrays, n = " + count);

        return false;

    }

    public boolean initDetectors() {
        return false;
    }

    public boolean initVoxels() {
        int iPSF = -1; // none
        int dPSF = -1; // none
        
        if (PBC) {
            geometry.initVoxelsPBC(iPSF, dPSF);
        } else {
            geometry.initVoxels(iPSF, dPSF);
        }
        
        voxelVolume = project.dx * project.dx * project.dx;
        
        return false;
        
    }

    public boolean initParticles() {

        double hydroDsDff;
        
        boolean shuffle = true;

        if (diffusants == null) {
            return true;
        }

        if (hydrodynamicsDscale && hydrodynamicsLocalD) {
            Master.exit("hydrodynamicDscale error: hydrodynamicInteractionsLocal is on");
        }

        if (hydrodynamicsDscale) {

            if (Dcyto == 0) {
                Master.exit("hydrodynamicsDscale error: Dcyto = 0");
            }

            if (diffusants.length > 1) {
                Master.exit("hydrodynamicsDscale error: more than one diffusant");
            }

            for (DiffusantParticles d : diffusants) {
                if (d != null) {
                    Master.log("hydrodynamic scaling old D = " + Dcyto);
                    d.D = Dcyto * DiffusantParticle.Dratio_shortNew(d.mobileVolumeFraction, d.immobileVolumeFraction);
                    Master.log("hydrodynamic scaling new D = " + d.D);
                }
            }

        }

        for (DiffusantParticles d : diffusants) {

            if (d != null) {

                d.initParticles();

                if ((d.xyz != null) && (d.particles.length == d.xyz.length)) {
                    for (int j = 0; j < d.particles.length; j++) {
                        setParticleLocation(d.particles[j], d.xyz[j][0], d.xyz[j][1], d.xyz[j][2], true);
                        addToVoxelList(d.particles[j]);
                    }
                }
                
                if (shuffle) {
                    d.shuffleParticles();
                }

            }

        }

        if (hydroWallZ) {

            for (DiffusantParticles d : diffusants) {

                if ((d == null) || (d.particles == null)) {
                    continue;
                }

                hydroDsDff = DiffusantParticle.Dratio_shortNew(d.mobileVolumeFraction, d.immobileVolumeFraction);
                hydroDsDff /= DiffusantParticle.Dff_short_Banchio(d.setVolumeFraction * (1 - d.setImmobilePercent));

                for (DiffusantParticle p : d.particles) {
                    p.DsDff = hydroDsDff;
                }

            }

        }

        return false;

    }

    public boolean initParticlesImmobile() {

        if (diffusants == null) {
            return true;
        }

        for (DiffusantParticles d : diffusants) {

            if (d == null) {
                continue;
            }

            if (d.initImmobileParticles(null)) {
                return true;
            }

        }

        return false;

    }

    public boolean saveParticleDensity(String saveTag) {
        return false;
    }

    public boolean saveParticlePositions(int msec) {

        if (preview || (diffusants == null)) {
            return true;
        }

        for (DiffusantParticles d : diffusants) {

            if (d == null) {
                continue;
            }

            if (d.saveXYZ) {
                d.saveXYZ(msec);
            }

        }

        return false;

    }

    public double spaceVolume() {
        return project.geometry.spaceVolume;
    }

    public void particleStats() {

        double distance;
        double spaceVolume = spaceVolume();

        minParticleRadius = Double.POSITIVE_INFINITY;
        maxParticleRadius = 0;
        minAllowedDistanceBetweenParticles = Double.POSITIVE_INFINITY;

        totalParticleVolume = 0;

        numParticles = 0;
        mobileParticles = 0;
        immobileParticles = 0;

        particleDensity = 0;
        particleVolumeFraction = 0;

        if (diffusants == null) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if (d == null) {
                continue;
            }

            //diffusants[k].init();
            d.particleStats();

            minParticleRadius = Math.min(minParticleRadius, d.radiusMin);
            maxParticleRadius = Math.max(maxParticleRadius, d.radiusMax);

            totalParticleVolume += d.volumeTotal;

            numParticles += d.numParticles;
            mobileParticles += d.mobileParticles;
            immobileParticles += d.immobileParticles;

        }

        mobilePercent = 1.0 * mobileParticles / (immobileParticles + mobileParticles);
        immobilePercent = 1.0 * immobileParticles / (immobileParticles + mobileParticles);

        for (DiffusantParticles d1 : diffusants) {

            if (d1 == null) {
                continue;
            }

            for (DiffusantParticles d2 : diffusants) {

                if (d2 == null) {
                    continue;
                }

                distance = d1.radiusMin + d2.radiusMin;
                minAllowedDistanceBetweenParticles = Math.min(distance, minAllowedDistanceBetweenParticles);

            }
        }

        particleVolume = 4 * Math.PI * minParticleRadius * minParticleRadius * minParticleRadius / 3.0;

        particleDensity = (numParticles * 1.0) / spaceVolume;
        particleVolumeFraction = totalParticleVolume / spaceVolume;

        mobileVolumeFraction = mobilePercent * particleVolumeFraction;
        immobileVolumeFraction = immobilePercent * particleVolumeFraction;

        setParamObject("minParticleRadius", minParticleRadius);
        setParamObject("maxParticleRadius", maxParticleRadius);
        setParamObject("minAllowedDistanceBetweenParticles", minAllowedDistanceBetweenParticles);

        setParamObject("totalParticleVolume", totalParticleVolume);

        setParamObject("numParticles", numParticles);
        setParamObject("mobileParticles", mobileParticles);
        setParamObject("immobileParticles", immobileParticles);
        setParamObject("mobilePercent", mobilePercent);
        setParamObject("immobilePercent", immobilePercent);
        setParamObject("mobileVolumeFraction", mobileVolumeFraction);
        setParamObject("immobileVolumeFraction", immobileVolumeFraction);

        setParamObject("particleDensity", particleDensity);
        setParamObject("particleVolumeFraction", particleVolumeFraction);

        minDistance();

    }
    
    public double maxD() {

        maxD = 0;

        if (project.diffusants == null) {
            return Double.NaN;
        }

        for (Diffusant d : project.diffusants) {
            if (d instanceof DiffusantParticles) {
                maxD = Math.max(maxD, d.D);
            }
        }

        return maxD;

    }

    public boolean checkDX() {

        double maxRadius = 0;

        if (diffusants == null) {
            return true;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (p == null) {
                    continue;
                }

                maxRadius = Math.max(maxRadius, p.radius);

            }

        }

        if (maxRadius <= 0) {
            Master.exit("MonteCarlo: checkDX: bad maxRadius: " + maxRadius);
        }

        if ((int) (100 * project.dx) < (int) (100 * 2 * maxRadius)) {
            //Master.exit("MonteCarlo: checkDX: voxel size is smaller than maximum particle diameter: " + (2 * maxRadius));
        }

        return false;

    }

    public boolean initDT() {

        if (project.diffusants == null) {
            return true;
        }

        maxD();

        if (maxD > 0) {
            project.dt = minParticleStep * minParticleStep / (6.0 * maxD);
        }
        
        //Master.log("maxD=" + maxD);
        //Master.log("dt=" + project.dt);

        return false;

    }

    public static double hydroWall_ll(double radius, double z) {

        double x = radius / z;

        // z = shortest distance from sphere center to wall
        x = Math.min(x, 1);

        double b1 = 1 - (9 * x / 16) + (x * x * x / 8) - (x * x * x * x * 45 / 256) - (x * x * x * x * x / 16);

        return b1;

    }

    public static double hydroWall_T(double radius, double h) {

        double b2;

        //double h = z - radius; // distance between wall and edge of sphere
        h = Math.max(h, 0);
        b2 = 6 * h * h + 2 * radius * h;
        b2 /= 6 * h * h + 9 * radius * h + 2 * radius * radius;

        return b2;

    }

    public boolean setParticleLocation(DiffusantParticle p, double x, double y, double z, boolean initStartLocation) {

        int xVoxel, yVoxel, zVoxel;

        if (PBC) {

            //if (geometry.voxelSpacePBC == null) {
            //    return false;
            //}
            p.x = x;
            p.y = y;
            p.z = z;

            if (initStartLocation) {
                p.x0 = x;
                p.y0 = y;
                p.z0 = z;
            }

            xVoxel = (int) geometry.computeVoxelX(x);
            yVoxel = (int) geometry.computeVoxelY(y);
            zVoxel = (int) geometry.computeVoxelZ(z);

            if ((xVoxel < 0) || (xVoxel >= geometry.voxelSpacePBC.length)) {
                //return false;
                return true;
            }

            if ((yVoxel < 0) || (yVoxel >= geometry.voxelSpacePBC[0].length)) {
                //return false;
                return true;
            }

            if ((zVoxel < 0) || (zVoxel >= geometry.voxelSpacePBC[0][0].length)) {
                //return false;
                return true;
            }

            p.voxel = geometry.voxelSpacePBC[xVoxel][yVoxel][zVoxel];

            return true;

        } else {

            //if (geometry.voxelSpace == null) {
            //    return false;
            //}
            xVoxel = (int) geometry.computeVoxelX(x);
            yVoxel = (int) geometry.computeVoxelY(y);
            zVoxel = (int) geometry.computeVoxelZ(z);

            if ((xVoxel < 0) || (xVoxel >= geometry.voxelSpace.length)) {
                return false;
            }

            if ((yVoxel < 0) || (yVoxel >= geometry.voxelSpace[0].length)) {
                return false;
            }

            if ((zVoxel < 0) || (zVoxel >= geometry.voxelSpace[0][0].length)) {
                return false;
            }

            p.x = x;
            p.y = y;
            p.z = z;

            if (initStartLocation) {
                p.x0 = x;
                p.y0 = y;
                p.z0 = z;
            }

            p.voxel = geometry.voxelSpace[xVoxel][yVoxel][zVoxel];

            return true;

        }

    }

    public boolean outOfBounds(DiffusantParticle p) {

        if (geometry.simpleCuboid) {
            return beyondLimits(p); // simple computation
        } else {
            if (beyondLimits(p)) {
                return true;
            }
            return insideVoxelNonSpace(p); // check nonspace voxels
        }

    }

    public boolean beyondLimits(DiffusantParticle p) {

        double extra;

        if (PBC || freeDiffusion) {
            extra = 0;
        } else {
            extra = p.radius;
        }

        //if (ellipsoid) {
        //    test = x * x / (lx1 * lx1) + y * y / (ly1 * ly1) + z * z / (lz1 * lz1);
        //    if (test > 1) {
        //        return true;
        //    }
        //}
        if ((p.x < geometry.x1 + extra) || (p.x > geometry.x2 - extra)) {
            return true;
        }

        if ((p.y < geometry.y1 + extra) || (p.y > geometry.y2 - extra)) {
            return true;
        }

        if ((p.z < geometry.z1 + extra) || (p.z > geometry.z2 - extra)) {
            return true;
        }

        return false;

    }

    public boolean insideVoxelNonSpace(DiffusantParticle p) {

        double dlimit;
        Voxel ivoxel;

        if (p.voxel == null) {
            return false;
        }

        if (!p.voxel.isSpace) {
            return true;
        }

        if (freeDiffusion) {
            return false;
        }

        dlimit = project.dx * 0.5 + p.radius; // um

        for (int i = 0; i < p.voxel.numNonSpaceNeighbors; i++) {

            ivoxel = p.voxel.nonSpaceNeighbors[i];

            if ((p.x > ivoxel.x - dlimit) && (p.x < ivoxel.x + dlimit)) {
                if ((p.y > ivoxel.y - dlimit) && (p.y < ivoxel.y + dlimit)) {
                    if ((p.z > ivoxel.z - dlimit) && (p.z < ivoxel.z + dlimit)) {
                        //Master.log("nonspace overlap " + ivoxel.x + "," + ivoxel.y + "," + ivoxel.z);
                        return true;
                    }
                }
            }

        }

        return false;

    }

    public void initParticleStartLocations() {

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            d.initStartLocation();

        }

    }

    public boolean initParticleRandom(DiffusantParticle p, Coordinates c, boolean allowOverlaps, long trialLimit) {

        long trial = 0;
        double x, y, z;
        boolean overlap;

        if ((p == null) || (c == null)) {
            return true;
        }

        while (true) {

            x = (Master.mt.nextDouble() * (c.x2 - c.x1)) + c.x1;
            y = (Master.mt.nextDouble() * (c.y2 - c.y1)) + c.y1;
            z = (Master.mt.nextDouble() * (c.z2 - c.z1)) + c.z1;

            if (!setParticleLocation(p, x, y, z, true)) {
                continue;
            }

            if (outOfBounds(p)) {
                continue;
            }

            if (testOverlap(p)) {
                continue;
            }

            if (freeDiffusion) {
                overlap = false;
            } else {
                overlap = (testParticleOverlap(p, null) != null);
            }

            if (!overlap) {
                return (!addToVoxelList(p));
            }
            
            trial++;
            
            if (trial > trialLimit) {
                
                if (allowOverlaps) {
                    return (!addToVoxelList(p));
                } else {
                    Master.log("initParticleRandom warning: reached trial limit");
                    return true;
                }
                
            }

        }

    }
    
    public boolean testOverlap(DiffusantParticle p) {
        return false; // generic function to test particle overlap
    }

    public boolean initParticlesRandom() {

        int nParticles, counter;
        double counter_percent = 1/10.0;

        boolean allowOverlaps = true;

        Coordinates c;

        if (diffusants == null) {
            return true;
        }

        int k;

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            if (d.coordinates == null) {
                c = new Coordinates(project, geometry);
            } else {
                c = new Coordinates(project, d.coordinates);
            }
            
            if (sortParticlesByRadius) {
                d.sortParticlesByRadius();
                //d.sortParticlesByRadiusSlow();
            }

            nParticles = d.particles.length;
            counter = (int) (nParticles * counter_percent);

            Master.log("randomly placing particles for diffusant '" + d.name + "' (n = " + nParticles + ")");

            k = 0;

            for (DiffusantParticle p : d.particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                if (initParticleRandom(p, c, allowOverlaps, initParticleRandomTrialLimit)) {
                    return true;
                }

                if ((k > 0) && (Math.IEEEremainder(k, counter) == 0)) {
                    Master.log("placed particle " + k + " / " + nParticles);
                }

                k++;

            }

        }

        return false;

    }

    public boolean initVParticlesLattice() {

        double newParticleVolume, r, jj;
        double x, y, z;
        double xlimit, ylimit, zlimit;
        int nParticles, iv;
        int i, j, k, ii;

        double maxPacking = Math.PI / (3 * Math.sqrt(2));
        double zstep = Math.sqrt(6) * 2.0 / 3.0;
        double sqrt3 = Math.sqrt(3);

        boolean overlap, finished = false;

        if (diffusants == null) {
            return true;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            nParticles = d.particles.length;
            //newParticleVolume = maxPacking * 1.0 / nParticles;
            //r = Math.pow((3 * newParticleVolume / (4 * Math.PI)), (1 / 3.0));
            r = d.setRadiusMean + 0.00000001;

            //Master.log("new radius " + r);
            if (PBC) {
                xlimit = geometry.x2;
                ylimit = geometry.y2;
                zlimit = geometry.z2;
            } else {
                xlimit = geometry.x2 - d.setRadiusMean;
                ylimit = geometry.y2 - d.setRadiusMean;
                zlimit = geometry.z2 - d.setRadiusMean;
            }

            iv = 0;
            k = 0;

            while (!finished) {
                
                if (iv >= d.particles.length) {
                    finished = true;
                    break;
                }

                z = geometry.z1 + r + k * zstep * r;

                if (z > zlimit) {
                    break;
                }

                j = 0;
                jj = (k % 2) / 3.0;

                while (!finished) {
                    
                    if (iv >= d.particles.length) {
                        finished = true;
                        break;
                    }

                    y = geometry.y1 + r + sqrt3 * (j + jj) * r;

                    if (y > ylimit) {
                        break;
                    }

                    i = 0;
                    ii = (j + k) % 2;

                    while (!finished) {
                        
                        if (iv >= d.particles.length) {
                            finished = true;
                            break;
                        }

                        x = geometry.x1 + r + (2 * i + ii) * r;

                        if (x > xlimit) {
                            break;
                        }

                        if (!setParticleLocation(d.particles[iv], x, y, z, true)) {
                            return true;
                        }

                        overlap = testParticleOverlap(d.particles[iv], null) != null;

                        if (!overlap) {

                            //if ((k % 3 == 1) && (j % 3 == 1) && (k % 3 == 1)) {
                            //    diffusant[id].particle[iv].mobile = true;
                            //} else {
                            //    diffusant[id].particle[iv].mobile = false;
                            //}
                            if (!addToVoxelList(d.particles[iv])) {
                                return true;
                            }

                            //if (iv == d.particles.length - 1) {
                            //    finished = true;
                            //    break;
                           // }

                            iv++;

                        }

                        i++;

                    }

                    j++;

                }

                k++;

            }

            Master.log("created particle lattice (" + iv + "/" + nParticles + ")");

            for (int iiv = iv; iiv < d.particles.length; iiv++) {
                d.particles[iiv].insideGeometry = false;
                d.particles[iiv].mobile = false;
            }

        }

        return false;

    }

    public boolean removeParticleOverlap() {

        int i, j;
        int radiusNM, lastRadiusNM = 0, maxNumParticles = 0;
        long count = 0;
        double halfDistance = 0;
        
        double[][] saveRadius;
        boolean[][] finished;

        if (freeDiffusion) {
            return false; // OK
        }

        if (!particleOverlapExists()) {
            return false; // OK
        }

        if (diffusants == null) {
            return true;
        }

        particleStats(); // will update particle step size

        removingParticleOverlap = true;

        Master.log("elimination of particle overlap (max radius=" + (maxParticleRadius * 1000) + " nm)...");
        
        for (DiffusantParticles d : diffusants) {
            maxNumParticles = Math.max(maxNumParticles, d.particles.length);
        }
        
        if (maxNumParticles == 0) {
            return true;
        }
        
        saveRadius = new double[diffusants.length][maxNumParticles];
        finished = new boolean[diffusants.length][maxNumParticles];
        
        for (double[] dd : saveRadius) {
            for (double d : dd ) {
                d = Double.NaN;
            }
        }
        
        for (boolean[] dd : finished) {
            for (boolean d : dd ) {
                d = false;
            }
        }
        
        i = 0;
        for (DiffusantParticles d : diffusants) {
            j = 0;
            for (DiffusantParticle p : d.particles) {
                saveRadius[i][j] = p.radius;
                j++;
            }
            i++;
        }
        
        while (halfDistance < maxParticleRadius) {

            halfDistance = 0.5 * minDistance();

            radiusNM = (int) (1000 * halfDistance);

            if (radiusNM > lastRadiusNM) {
                Master.log("current radius: " + radiusNM + " nm");
                lastRadiusNM = radiusNM;
                count = 0;
            }
            
            i = 0;

            // set current diameter to minimal distance between particles
            for (DiffusantParticles d : diffusants) {

                //if ((diffusant[k] == null) || (diffusant[k].particle == null)) {
                //    continue;
                //}
                
                j = 0;
                
                for (DiffusantParticle p : d.particles) {
                    if (p.insideGeometry && !finished[i][j]) {
                        //v.radius = Math.min(halfDistance, diffusant[k].meanRadius);
                        //v.radius = Math.min(halfDistance, saveRadius[i][j]);
                        if (halfDistance >= saveRadius[i][j]) {
                            p.radius = saveRadius[i][j];
                            finished[i][j] = true;
                        } else {
                            p.radius = halfDistance;
                        }
                        
                    }
                    j++;
                }
                
                i++;

            }

            moveParticlesCichocki(true);

            count++;

            if (count > removeParticleOverlapTrialLimit) {
                error("removeParticleOverlap: aborted after reaching trial Limit");
                break;
            }

        }

        //i = 0;
        
        for (DiffusantParticles d : diffusants) {

            //if ((d == null) || (d.particles == null)) {
                //continue;
            //}

            //radius = diffusant[k].meanRadius;
            //j = 0;
            
            for (DiffusantParticle p : d.particles) {
                if (p.insideGeometry) {
                    //v.radius = d.meanRadius;
                    //v.radius = saveRadius[i][j];
                    p.x0 = p.x;
                    p.y0 = p.y;
                    p.z0 = p.z;
                }
                //j++;
            }
            
            //i++;

        }

        removingParticleOverlap = false;

        if (particleOverlapExists()) {
            error("removeParticleOverlap: failed to remove particle overlap.");
            return true;
        }

        return initParticlesImmobile();

    }

    public double minDistanceZstack() {

        double sqrDistance, minSqrDist, minSqrDist2;
        double DX, DY, DZ;

        minSqrDist = Double.POSITIVE_INFINITY;
        minDistanceBetweenParticles = Double.NaN;

        boolean found = false;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p1 : d.particles) {

                if (!p1.insideGeometry) {
                    continue;
                }

                minSqrDist2 = Double.POSITIVE_INFINITY;

                for (DiffusantParticle p2 : d.particles) {

                    if (p1 == p2) {
                        continue;
                    }

                    if (!p2.insideGeometry) {
                        continue;
                    }

                    if (p1.z != p2.z) {
                        continue;
                    }

                    DX = p1.x - p2.x;
                    DY = p1.y - p2.y;
                    DZ = p1.z - p2.z;

                    //sqrDistance = DX * DX + DY * DY + DZ * DZ;
                    sqrDistance = DX * DX + DY * DY;

                    if (sqrDistance == 0) {
                        Master.log("warning: zero sqrDistance " + p1 + ", " + p2);
                    }

                    //if (sqrDistance < 0.02) {
                    //Master.log("" + DX + ", " + DY + ", " + DZ);
                    //}
                    if (sqrDistance < minSqrDist) {
                        minSqrDist = sqrDistance;
                        found = true;
                    }

                    if (sqrDistance < minSqrDist2) {
                        minSqrDist2 = sqrDistance;
                    }

                }

                Master.log("" + Math.sqrt(minSqrDist2));

            }

        }

        if (found) {
            minDistanceBetweenParticles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenParticles = Double.NaN;
        }

        return minDistanceBetweenParticles;

    }

    public double minDistanceSlow() {

        int countOverlaps = 0;
        double sqrDistance, minSqrDist, minSqrDist2;
        double DX, DY, DZ;

        minSqrDist = Double.POSITIVE_INFINITY;
        minDistanceBetweenParticles = Double.NaN;

        boolean found = false;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (int j = 0; j < d.particles.length; j++) {

                if (!d.particles[j].insideGeometry) {
                    continue;
                }

                minSqrDist2 = Double.POSITIVE_INFINITY;

                for (int k = j + 1; k < d.particles.length; k++) {

                    //if (j == k) {
                    //    continue;
                    //}
                    if (!d.particles[k].insideGeometry) {
                        continue;
                    }

                    DX = d.particles[j].x - d.particles[k].x;
                    DY = d.particles[j].y - d.particles[k].y;
                    DZ = d.particles[j].z - d.particles[k].z;

                    sqrDistance = DX * DX + DY * DY + DZ * DZ;

                    if (sqrDistance == 0) {
                        Master.log("warning: zero sqrDistance " + j + ", " + k);
                    }

                    if (Math.sqrt(sqrDistance) < 0.04) {
                        if ((Math.abs(DZ) <= 0.03) && (Math.abs(DZ) > 0)) {
                            if ((Math.abs(DX) < 0.01) && (Math.abs(DY) < 0.01)) {
                                //Master.log("" + DX + ", " + DY + ", " + DZ);
                                //Master.log("" + (Math.sqrt(DX * DX + DY * DY)));
                                d.particles[k].insideGeometry = false;
                                d.particles[k].x = Double.NaN;
                                d.particles[k].y = Double.NaN;
                                d.particles[k].z = Double.NaN;
                                countOverlaps++;
                            }
                        }
                    }

                    if (sqrDistance < minSqrDist) {
                        minSqrDist = sqrDistance;
                        found = true;
                    }

                    if (sqrDistance < minSqrDist2) {
                        minSqrDist2 = sqrDistance;
                    }

                }

                //Master.log("" + Math.sqrt(minSqrDist2));
            }

        }

        if (countOverlaps > 0) {
            Master.log("removed overlaps = " + countOverlaps);
        }

        if (found) {
            minDistanceBetweenParticles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenParticles = Double.NaN;
        }

        return minDistanceBetweenParticles;

    }

    public double minDistance() {

        double min, minSqrDist = Double.MAX_VALUE;

        boolean found = false;

        minDistanceBetweenParticles = Double.NaN;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (p == null) {
                    continue;
                }

                if (p.insideGeometry) {

                    min = minSqrDistance(p, p);

                    if (min >= 0) {

                        if (min < minSqrDist) {
                            minSqrDist = min;
                            found = true;
                        }

                    } else {
                        //Master.log("failed to find min distance for particles " + j);
                    }

                }

            }

        }

        if (found) {
            minDistanceBetweenParticles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenParticles = minDistanceSlow();
        }

        return minDistanceBetweenParticles;

    }

    public void testMinDistance() {

        int counter = 0, isteps = 2000;
        double sec, min;

        if (diffusants == null) {
            return;
        }

        timer1.start();

        for (int istep = 0; istep < isteps; istep++) {

            for (DiffusantParticles d : diffusants) {

                if ((d == null) || (d.particles == null)) {
                    continue;
                }

                for (DiffusantParticle p : d.particles) { // runSimulation thru each particles

                    if (p == null) {
                        continue;
                    }

                    if (p.insideGeometry) {
                        min = minSqrDistance(p, p);
                        counter++;
                    }

                }

            }

        }

        timer1.stop();
        sec = (double) (timer1.milliseconds / 1000.0);
        Master.log("Test Time (" + counter + ") = " + sec + " sec");

    }

    public double minSqrDistance(DiffusantParticle p, DiffusantParticle ignoreParticle) {

        double dx, dy, dz, sqrDistance, min = Double.MAX_VALUE;
        int ii, jj, kk;
        boolean found = false;

        DiffusantParticle iParticle;

        Voxel voxel = p.voxel;
        VoxelPBC voxelPBC = null;
        Voxel ivoxel;

        if (voxel == null) {
            return -1;
        }

        if (PBC && (voxel instanceof VoxelPBC)) {
            voxelPBC = (VoxelPBC) voxel;
        }

        for (int i = 0; i < voxel.numNeighbors; i++) {

            ivoxel = voxel.neighbors[i];
            iParticle = ivoxel.firstParticleInChain;

            while (iParticle != null) {

                if ((iParticle != p) && (iParticle != ignoreParticle) && (!Double.isNaN(iParticle.x)) && iParticle.insideGeometry) {

                    dx = p.x - iParticle.x;
                    dy = p.y - iParticle.y;
                    dz = p.z - iParticle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    if (sqrDistance < min) {
                        min = sqrDistance;
                        found = true;
                    }

                    if (sqrDistance <= 0) {
                        Master.log("warning: zero sqrDistance " + p + " and " + iParticle);
                    }

                }

                iParticle = iParticle.nextParticleInChain;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];
                iParticle = ivoxel.firstParticleInChain;

                while (iParticle != null) {

                    if ((iParticle != p) && (iParticle != ignoreParticle) && (!Double.isNaN(iParticle.x)) && iParticle.insideGeometry) {

                        ii = voxelPBC.PBCi[i];
                        jj = voxelPBC.PBCj[i];
                        kk = voxelPBC.PBCk[i];

                        dx = p.x - (ii * 2 * geometry.x2 + iParticle.x);
                        dy = p.y - (jj * 2 * geometry.y2 + iParticle.y);
                        dz = p.z - (kk * 2 * geometry.z2 + iParticle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        if (sqrDistance < min) {
                            min = sqrDistance;
                            found = true;
                        }

                    }

                    iParticle = iParticle.nextParticleInChain;

                }

            }

        }

        if (found) {
            return min;
        } else {
            return -1;
        }

    }

    public boolean particleOverlapExists() {

        if (diffusants == null) {
            return false;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (p == null) {
                    continue;
                }

                if (!p.insideGeometry) {
                    continue;
                }

                if (testParticleOverlap(p, p) != null) {
                    return true;
                }

            }
        }

        return false;

    }

    public DiffusantParticle testParticleOverlap(DiffusantParticle p, DiffusantParticle ignoreParticle) {

        double dx, dy, dz, sqrDistance, minDBV;
        int ii, jj, kk;

        DiffusantParticle iParticle;

        Voxel voxel = p.voxel;
        VoxelPBC voxelPBC = null;
        Voxel ivoxel;

        if (voxel == null) {
            return null;
        }

        if (PBC && (voxel instanceof VoxelPBC)) {
            voxelPBC = (VoxelPBC) voxel;
        }

        for (int i = 0; i < voxel.numNeighbors; i++) {

            ivoxel = voxel.neighbors[i];
            iParticle = ivoxel.firstParticleInChain;

            while (iParticle != null) {

                if ((iParticle != p) && (iParticle != ignoreParticle) && (!Double.isNaN(iParticle.x)) && iParticle.insideGeometry && (iParticle.radius > 0)) {

                    dx = p.x - iParticle.x;
                    dy = p.y - iParticle.y;
                    dz = p.z - iParticle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    minDBV = p.radius + iParticle.radius;

                    if (sqrDistance < minDBV * minDBV) {
                        return iParticle;
                    }

                }

                iParticle = iParticle.nextParticleInChain;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];
                iParticle = ivoxel.firstParticleInChain;

                while (iParticle != null) {

                    if ((iParticle != p) && (iParticle != ignoreParticle) && (!Double.isNaN(iParticle.x)) && iParticle.insideGeometry) {

                        ii = voxelPBC.PBCi[i];
                        jj = voxelPBC.PBCj[i];
                        kk = voxelPBC.PBCk[i];

                        dx = p.x - (ii * 2 * geometry.x2 + iParticle.x);
                        dy = p.y - (jj * 2 * geometry.y2 + iParticle.y);
                        dz = p.z - (kk * 2 * geometry.z2 + iParticle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        minDBV = p.radius + iParticle.radius;

                        if (sqrDistance < minDBV * minDBV) {
                            return iParticle;
                        }

                    }

                    iParticle = iParticle.nextParticleInChain;

                }

            }

        }

        return null;

    }

    public double avgMinSqrDistance() {

        double minSD;
        double avgMinSD = 0, count = 0;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (p == null) {
                    continue;
                }

                minSD = minSqrDistance(p, p);

                if (minSD > 0) {
                    avgMinSD += minSD;
                    count += 1;
                }

            }

        }

        return (avgMinSD / count);

    }

    public void runMSDspatial() {

        double svalue;

        if (!saveMSDspatial || (MSDspatialTbgn < 0)) {
            return;
        }

        if ((MSDspatial == null) && (time >= MSDspatialTbgn)) {
            initMSDspatial(MSDspatialTbgn);
            initD20();
            MSDspatial_t1 = MSDspatialTbgn;
            MSDspatial_t2 = MSDspatial_t1 + MSDspatial_win;
        }

        if ((MSDspatial_t1 < 0) || (MSDspatial_t2 < 0)) {
            return;
        }

        if ((MSDspatial != null) && (time > MSDspatial_t2)) {

            for (Save s : MSDspatial) {
                s.finish("Particles", project.geometry, -1);
            }

            MSDspatial = null;
            MSDspatial_t1 = -1;
            MSDspatial_t2 = -1;
            MSDspatialTbgn = -1;

            return;

        }

        if ((time >= MSDspatial_t1) && (time <= MSDspatial_t2)) {
            for (int i = 0; i < MSDspatial.length; i++) {
                if (MSDspatial[i].skipCounter == 0) {
                    svalue = MSDspatial(i);
                } else {
                    svalue = -1; // does not matter
                }
                MSDspatial[i].saveData(svalue);
            }
        }

    }

    public void initMSDspatial(double time) {

        int bin, it;
        int dataPoints = 1;
        String fname;

        MSDspatial = new Save[MSDspatial_numBins];

        for (int i = 0; i < MSDspatial.length; i++) {
            MSDspatial[i] = new Save(project);
            MSDspatial[i].saveWhileComputing = false;
            MSDspatial[i].save2BinaryFile = true;
            MSDspatial[i].save2TextFile = false;
            MSDspatial[i].init();
            bin = (int) ((i + 1) * MSDspatial_binWidth * 1000);
            it = (int) time;
            fname = "MSD_t" + Integer.toString(it) + "_d" + Integer.toString(bin);
            Master.log("init save " + fname);
            MSDspatial[i].fileName(fname, "");
            MSDspatial[i].xdim = "ms";
            MSDspatial[i].ydim = "Particles" + " (" + project.spaceUnits + "^2)";
            MSDspatial[i].updateVectors();
            MSDspatial[i].init("Particles", project.geometry, -1, dataPoints);
        }

    }

    public void initD20() {
        double dx, dy, dz;

        if (diffusants == null) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if ((p == null) || !p.mobile || !p.insideGeometry) {
                    continue;
                }

                p.x0 = p.x;
                p.y0 = p.y;
                p.z0 = p.z;

                dx = p.x - 0;
                dy = p.y - 0;
                dz = p.z - 0;

                p.d20 = Math.sqrt(dx * dx + dy * dy + dz * dz);

            }

        }

    }

    public double MSDspatial(int ibin) {

        double sumSD = 0.0, count = 0.0;
        double x1 = ibin * MSDspatial_binWidth;
        double x2 = x1 + MSDspatial_binWidth;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if ((p == null) || !p.mobile || !p.insideGeometry) {
                    continue;
                }

                if ((p.d20 >= x1) && (p.d20 < x2)) {
                    sumSD += p.squareDisplacement();
                    count += 1.0;
                    //Master.log("" + dv + " " + dv.d2AZ0);
                }

            }

        }

        return (sumSD / count);

    }

    public boolean initConnectors() {

        if (diffusants == null) {
            return false;
        }

        if (!connectParticles) {
            return false;
        }

        if (connectorsBinomial) {
            return initConnectorBinomialDistribution();
        }

        Master.log("init " + maxNumConnectors + " connectors / particle");

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (!(p instanceof DiffusantParticle)) {
                    continue;
                }

                p.connectTo = new DiffusantParticle[maxNumConnectors];
                p.connectorOffTime = new double[maxNumConnectors];

            }

        }

        return false;

    }

    private boolean initConnectorBinomialDistribution() {

        double avg = 0, count = 0;
        int nConnectors, numConnected = 0;
        //double meanNumConnectors = 0.9;//1.5; // 1.0
        double probability = meanNumConnectors / (1.0 * maxNumConnectors); // Fernandez-Busnadiego 2013

        int[] numConnectorsHisto = new int[maxNumConnectors + 1];

        if (diffusants == null) {
            return false;
        }

        Master.log("max number of connectors = " + maxNumConnectors);
        Master.log("connector probability = " + probability);

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (!(p instanceof DiffusantParticle)) {
                    continue;
                }

                nConnectors = 0;

                for (int k = 0; k < maxNumConnectors; k++) {
                    if (Master.mt.nextDouble() < probability) {
                        nConnectors++;
                    }
                }

                if (nConnectors > 0) {
                    p.connectTo = new DiffusantParticle[nConnectors]; // need twice as many for tracking connectors
                    p.connectorOffTime = new double[nConnectors]; // need twice as many for tracking connectors
                }

                numConnectorsHisto[nConnectors]++;

            }

        }

        for (int i = 0; i < numConnectorsHisto.length; i++) {
            Master.log("" + i + " connectors: " + numConnectorsHisto[i]);
            avg += i * numConnectorsHisto[i];
            count += numConnectorsHisto[i];
            if (i > 0) {
                numConnected += numConnectorsHisto[i];
            }
        }

        avg /= count;

        Master.log("average = " + avg + " connectors / particle");
        Master.log("fraction connected = " + (numConnected / count));

        return false;

    }

    public void connectParticles() {

        int ii, jj, kk;
        double dx, dy, dz, sqrDistance, minDBV;
        double probability_per_time_step = connectRate * project.dt;

        DiffusantParticle kParticle;

        VoxelPBC voxelPBC;

        for (int i = 0; i < diffusants.length; i++) {

            if ((diffusants[i] == null) || (diffusants[i].particles == null)) {
                continue;
            }

            for (DiffusantParticle p : diffusants[i].particles) {

                for (int k = 0; k < p.voxel.numNeighbors; k++) {

                    kParticle = (DiffusantParticle) p.voxel.neighbors[k].firstParticleInChain;

                    while (kParticle != null) {

                        if (p.noMoreConnectors()) {
                            break;
                        }

                        if ((kParticle != p) && !p.isConnectedTo(kParticle) && !kParticle.noMoreConnectors()) {

                            if (p.overlap(kParticle, connectorLength)) {

                                if (Master.mt.nextDouble() < probability_per_time_step) {

                                    if (kParticle.connectToNew(p, unconnectRate, time)) {

                                        avgConnectorLifeTime += p.connectorLifeTime;
                                        connectorLifeTimeCounter += 1;

                                        //Master.log("connected particles " + dv.connectorLifeTime);
                                        if (!p.connectTo(kParticle)) {
                                            Master.exit("connection failure");
                                        }

                                    } else {
                                        Master.exit("connection failure");
                                    }

                                }

                            }

                        }

                        kParticle = kParticle.nextParticleInChain;

                    }

                }

                if (p.noMoreConnectors()) {
                    continue;
                }

                if (!PBC) {
                    continue;
                }

                if (p.voxel instanceof VoxelPBC) {
                    voxelPBC = (VoxelPBC) p.voxel;
                } else {
                    continue;
                }

                if (voxelPBC == null) {
                    continue;
                }

                if (true) {
                    continue;
                }

                for (int k = 0; k < voxelPBC.numPBCneighbors; k++) {

                    kParticle = (DiffusantParticle) voxelPBC.PBCneighbors[k].firstParticleInChain;

                    while (kParticle != null) {

                        if (p.noMoreConnectors()) {
                            break;
                        }

                        if ((kParticle != p) && !p.isConnectedTo(kParticle) && !kParticle.noMoreConnectors()) {

                            ii = voxelPBC.PBCi[i];
                            jj = voxelPBC.PBCj[i];
                            kk = voxelPBC.PBCk[i];

                            dx = p.x - (ii * 2 * geometry.x2 + kParticle.x);
                            dy = p.y - (jj * 2 * geometry.y2 + kParticle.y);
                            dz = p.z - (kk * 2 * geometry.z2 + kParticle.z);

                            sqrDistance = dx * dx + dy * dy + dz * dz;

                            minDBV = p.radius + kParticle.radius + connectorLength;

                            if (sqrDistance < minDBV * minDBV) {

                                if (Master.mt.nextDouble() < connectRate * project.dt) {

                                    if (kParticle.connectToNew(p, unconnectRate, time)) {

                                        avgConnectorLifeTime += p.connectorLifeTime;
                                        connectorLifeTimeCounter += 1;

                                        //Master.log("connected particles PBC " + dv.connectorLifeTime);
                                        if (!p.connectTo(kParticle)) {
                                            Master.exit("connection failure");
                                        }

                                    } else {
                                        Master.exit("connection failure");
                                    }

                                }

                            }

                        }

                        kParticle = kParticle.nextParticleInChain;

                    }

                }

            }

        }

    }

    public void unconnectParticles() {

        DiffusantParticle p2;

        if (diffusants == null) {
            return;
        }

        if (unconnectRate <= 0) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p1 : d.particles) {

                if (p1.connectTo == null) {
                    continue;
                }

                for (int k = 0; k < p1.connectTo.length; k++) {

                    if (p1.connectTo[k] == null) {
                        continue;
                    }

                    if ((connectorAllOffTime > 0) && (time >= connectorAllOffTime)) {
                        p2 = p1.connectTo[k];
                        p2.unconnectFrom(p1);
                        p1.connectTo[k] = null;
                        p1.connectorOffTime[k] = 0;
                    }

                    if ((p1.connectorOffTime[k] > 0) && (time >= p1.connectorOffTime[k])) {
                        //Master.log("unconnected particles at " + time + " ms, offtime = " + (dv.connectorAllOffTime[1][k] * project.dt));
                        //Master.log("unconnected particles at " + time + " ms");
                        p2 = p1.connectTo[k];
                        p2.unconnectFrom(p1);
                        p1.connectTo[k] = null;
                        p1.connectorOffTime[k] = 0;
                    }

                }

            }

        }

        if ((connectorAllOffTime > 0) && (time >= connectorAllOffTime)) {
            connectParticles = false;
        }

    }

    public double[] computeNumberOfConnectors() {

        double avg = 0, connected = 0;
        int n, ntotal = 0;

        if (diffusants == null) {
            return null;
        }

        double[] results = new double[2];

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                n = p.numberOfConnections(true);
                avg += n;
                ntotal++;

                if (n > 0) {
                    connected++;
                }

            }

        }

        if (ntotal > 0) {
            avg /= ntotal;
            connected /= ntotal;
        } else {
            avg = 0;
            connected = 0;
        }

        //Master.log("" + avg + " connectors / particle");
        //Master.log("" + connected + " fraction connected");
        results[0] = avg;
        results[1] = connected;

        return results;

    }

    private void testConnectorSeperation() {

        if (diffusants == null) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (!(p instanceof DiffusantParticle)) {
                    continue;
                }

                if (p.testConnectorSeperation(connectorLength)) {
                    Master.log("connectors seperated " + p);
                }

            }
        }

    }

    public double avgFirstCollision() {

        double avg = 0, count = 0, sqrd = 0;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if ((p == null) || !p.mobile) {
                    continue;
                }

                if (p.firstCollision > 0) {
                    avg += p.firstCollision;
                    sqrd += p.sqrDisplacement;
                    count++;
                    //Master.log("" + vs.particle[j].firstCollision);
                }

            }

        }

        if (count == 0) {
            Master.log("no collisions");
        } else {

            avg /= count;
            sqrd /= count;

            Master.log("avg first collision = " + avg + " ms (n=" + count + ")");
            Master.log("square displacement = " + sqrd + " um^2");
            Master.log("D = " + (sqrd / (6 * avg)) + " um^2/ms");

        }

        return avg;

    }

    public void updateResidentTimes() {

        if (diffusants == null) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (p == null) {
                    continue;
                }

                if (p.mobile && (p.residenceTime >= 0)) {
                    p.residenceTime += project.dt;
                }

            }

        }

    }

    public DiffusantParticle findParticleOverlap(DiffusantParticle p) {

        double dx, dy, dz, sqrDistance, minDBV;
        int ii, jj, kk;

        DiffusantParticle iParticle;

        Voxel voxel = p.voxel;
        VoxelPBC voxelPBC = null;
        Voxel ivoxel;

        if (voxel == null) {
            return null;
        }

        if (PBC && (voxel instanceof VoxelPBC)) {
            voxelPBC = (VoxelPBC) voxel;
        }

        for (int i = 0; i < voxel.numNeighbors; i++) {

            ivoxel = voxel.neighbors[i];
            iParticle = ivoxel.firstParticleInChain;

            while (iParticle != null) {

                if ((iParticle != p) && (!Double.isNaN(iParticle.x)) && iParticle.insideGeometry) {

                    dx = p.x - iParticle.x;
                    dy = p.y - iParticle.y;
                    dz = p.z - iParticle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    minDBV = p.radius + iParticle.radius;

                    if (sqrDistance < minDBV * minDBV) {
                        return iParticle;
                    }

                }

                iParticle = iParticle.nextParticleInChain;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];

                ii = voxelPBC.PBCi[i];
                jj = voxelPBC.PBCj[i];
                kk = voxelPBC.PBCk[i];

                iParticle = ivoxel.firstParticleInChain;

                while (iParticle != null) {

                    if ((iParticle != p) && (!Double.isNaN(iParticle.x)) && iParticle.insideGeometry) {

                        dx = p.x - (ii * 2 * geometry.x2 + iParticle.x);
                        dy = p.y - (jj * 2 * geometry.y2 + iParticle.y);
                        dz = p.z - (kk * 2 * geometry.z2 + iParticle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        minDBV = p.radius + iParticle.radius;

                        if (sqrDistance < minDBV * minDBV) {
                            //return iParticle;
                            return null;
                        }

                    }

                    iParticle = iParticle.nextParticleInChain;

                }

            }

        }

        return null;

    }

    public double localDensityFastHydroWallz(DiffusantParticle p) {

        double h, volumeCap, volumeTotal, localVolumeFraction, localD;
        double nRadiiExt = 4; // Rext

        double dz = geometry.z2 - p.z;
        double radius4 = p.radius * nRadiiExt;

        if (dz >= radius4) {
            return p.step3;
        }

        h = radius4 - dz;
        volumeCap = (Math.PI * h * h / 3) * (3 * radius4 - h);
        volumeTotal = 4 * Math.PI * radius4 * radius4 * radius4 / 3;

        //if (dz < dv.radius) {
        //    Master.log("" + dz);
        //}
        localVolumeFraction = (volumeTotal - volumeCap) * particleVolumeFraction / volumeTotal;

        //localD = Dcyto * DiffusantParticle.Dratio_shortOLD(localVolumeFraction, immobileParticlePercent);
        return Double.NaN; // DiffusantParticle.step3(localD, project.dt);

        //Master.log("" + particleVolumeFraction);
        //return (volumeTotal - volumeCap)/volumeTotal;
    }

    public double localDensityAll(boolean print) {

        double density, avg = 0, avgD = 0, avgDD0 = 0, avgStep = 0, count = 0;
        double min = 99999, max = 0;

        //DiffusantParticle p;
        //if (immobileParticleFraction > 0) {
        //    Master.exit("aborted MC simulation: immobile fraction not allowed with local density");
        //}
        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (p == null) {
                    continue;
                }

                if (p.insideGeometry) {

                    if (PBC) {
                        density = localDensityPBC(p);
                    } else if (hydrodynamicsLocalDVoxels) {
                        density = localDensityVoxels(p);
                    } else {
                        density = localDensity(p);
                    }

                    if (print && (density > 0)) {
                        avg += density;
                        avgStep += p.step3Local;
                        count++;
                        min = Math.min(min, density);
                        max = Math.max(max, density);
                    }

                }

            }

        }

        avg /= count;
        avgStep /= count;

        if (print) {
            Master.log("local density avg = " + avg);
            Master.log("local density min = " + min);
            Master.log("local density max = " + max);
            //Master.log("local D/D0 avg = " + avgDD0 / count);
            Master.log("local D avg = " + DiffusantParticle.D3(avgStep, project.dt));
            Master.log("local step avg = " + avgStep);
        }

        return avg;

    }

    public double localDensity(DiffusantParticle p) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt, nRadiiInt;
        int x0, y0, z0, x1, y1, z1;
        double dx, dy, dz;
        double d, h, Rext, r, totalV, sumV, volumeCap = 0;
        double localMobileVolumeFraction, localD, Ds;

        Voxel voxel;
        DiffusantParticle iParticle;

        if ((p == null) || (geometry.voxelSpace == null)) {
            return Double.NaN;
        }

        //nRadiiInt = 2; // Rint
        //Rint = nRadiiInt * dv.radius;
        dz = geometry.z2 - p.z;

        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * p.radius;
        h = Rext - dz;

        if (h > 0) {
            volumeCap = (Math.PI * h * h / 3) * (3 * Rext - h);
        }

        totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        //sumV = 4.0 * Math.PI * dv.radius * dv.radius * dv.radius / 3.0;
        sumV = particleVolume; // all particles with same radius

        xVoxel = (int) geometry.computeVoxelX(p.x);
        yVoxel = (int) geometry.computeVoxelY(p.y);
        zVoxel = (int) geometry.computeVoxelZ(p.z);

        nVoxels = (int) Math.ceil(Rext / project.dx);

        x0 = xVoxel - nVoxels;
        y0 = yVoxel - nVoxels;
        z0 = zVoxel - nVoxels;

        x1 = xVoxel + nVoxels;
        y1 = yVoxel + nVoxels;
        z1 = zVoxel + nVoxels;

        x0 = Math.max(x0, 0);
        y0 = Math.max(y0, 0);
        z0 = Math.max(z0, 0);

        x1 = Math.min(x1, geometry.voxelSpace.length - 1);
        y1 = Math.min(y1, geometry.voxelSpace[0].length - 1);
        z1 = Math.min(z1, geometry.voxelSpace[0][0].length - 1);

        p.step3Local = 0;

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    voxel = geometry.voxelSpace[i][j][k];
                    iParticle = voxel.firstParticleInChain;

                    while (iParticle != null) {

                        if ((iParticle != p) && iParticle.insideGeometry && iParticle.mobile) {

                            dx = p.x - iParticle.x;
                            dy = p.y - iParticle.y;
                            dz = p.z - iParticle.z;

                            d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            r = iParticle.radius;

                            if (d < Rext - r) {
                                sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            } else if (d < Rext + r) {
                                sumV += Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                            }

                        }

                        iParticle = iParticle.nextParticleInChain;

                    }

                }
            }
        }

        localMobileVolumeFraction = sumV / (totalV - volumeCap);
        //Ds = DiffusantParticle.Dratio_shortOLD(localMobileVolumeFraction, immobileParticlePercent);
        Ds = DiffusantParticle.Dratio_shortNew(localMobileVolumeFraction, immobileVolumeFraction);
        //localD = dv.D * Ds;
        localD = Dcyto * Ds;
        p.step3Local = DiffusantParticle.step3(localD, project.dt);
        p.DsDff = Ds / DiffusantParticle.Dff_short_Banchio(localMobileVolumeFraction);

        return localMobileVolumeFraction;

    }

    public double localDensityVoxels(DiffusantParticle p) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt, nRadiiInt;
        int x0, y0, z0, x1, y1, z1;
        double dx, dy, dz;
        double d, h, Rext, r, totalV = 0, sumV_mobile, sumV_immobile = 0;
        double localMobileVolumeFraction, localImmobileVolumeFraction, localD, Ds;

        Voxel voxel;
        DiffusantParticle iParticle;

        if ((p == null) || (geometry.voxelSpace == null)) {
            return Double.NaN;
        }

        //nRadiiInt = 2; // Rint
        //Rint = nRadiiInt * dv.radius;
        //dz = geometry.z2 - dv.z;
        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * p.radius;
        //h = Rext - dz;

        //if (h > 0) {
        //    volumeCap = (Math.PI * h * h / 3) * (3 * Rext - h);
        //}
        //totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        //sumV = 4.0 * Math.PI * dv.radius * dv.radius * dv.radius / 3.0;
        sumV_mobile = particleVolume; // all particles with same radius

        xVoxel = (int) geometry.computeVoxelX(p.x);
        yVoxel = (int) geometry.computeVoxelY(p.y);
        zVoxel = (int) geometry.computeVoxelZ(p.z);

        nVoxels = (int) Math.ceil(Rext / project.dx);

        x0 = xVoxel - nVoxels;
        y0 = yVoxel - nVoxels;
        z0 = zVoxel - nVoxels;

        x1 = xVoxel + nVoxels;
        y1 = yVoxel + nVoxels;
        z1 = zVoxel + nVoxels;

        x0 = Math.max(x0, 0);
        y0 = Math.max(y0, 0);
        z0 = Math.max(z0, 0);

        x1 = Math.min(x1, geometry.voxelSpace.length - 1);
        y1 = Math.min(y1, geometry.voxelSpace[0].length - 1);
        z1 = Math.min(z1, geometry.voxelSpace[0][0].length - 1);

        p.step3Local = 0;

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    voxel = geometry.voxelSpace[i][j][k];

                    dx = p.x - voxel.x;
                    dy = p.y - voxel.y;
                    dz = p.z - voxel.z;

                    d = Math.sqrt(dx * dx + dy * dy + dz * dz);

                    if (d > Rext) {
                        continue;
                    }

                    if (voxel.isSpace) {
                        totalV += voxelVolume;
                    }

                    iParticle = voxel.firstParticleInChain;

                    while (iParticle != null) {

                        if ((iParticle != p) && iParticle.insideGeometry) {

                            //dx = dv.x - iParticle.x;
                            //dy = dv.y - iParticle.y;
                            //dz = dv.z - iParticle.z;
                            //d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            //r = iParticle.radius;
                            //sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            if (iParticle.mobile) {
                                sumV_mobile += particleVolume;
                            } else {
                                sumV_immobile += particleVolume;
                            }

                            //if (d < Rext - r) {
                            //    sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            //} else if (d < Rext + r) {
                            //    sumV += Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                            //}
                        }

                        iParticle = iParticle.nextParticleInChain;

                    }

                }
            }
        }

        //localDensity =  sumV / (totalV - volumeCap);
        localMobileVolumeFraction = sumV_mobile / totalV;
        localImmobileVolumeFraction = sumV_immobile / totalV;

        //Ds = DiffusantParticle.Dratio_shortOLD(localMobileVolumeFraction, immobileParticlePercent);
        Ds = DiffusantParticle.Dratio_shortNew(localMobileVolumeFraction, immobileVolumeFraction);
        //Ds = DiffusantParticle.Dratio_shortNew(localMobileVolumeFraction, localImmobileVolumeFraction);
        //localD = dv.D * Ds;
        localD = Dcyto * Ds;
        //Master.log("localD = " + localD);
        p.step3Local = DiffusantParticle.step3(localD, project.dt);
        p.DsDff = Ds / DiffusantParticle.Dff_short_Banchio(localMobileVolumeFraction);

        if (Double.isNaN(p.step3Local)) {
            p.step3Local = p.step3;
            Master.log("NaN localSteps: " + localMobileVolumeFraction + "," + p.step3Local);
        }

        return localMobileVolumeFraction;

    }

    public double localDensityPBC(DiffusantParticle p) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt, nRadiiInt;
        int x0, y0, z0, x1, y1, z1, ii, jj, kk;
        double dx, dy, dz;
        double d, Rext, r, totalV, sumV;
        double xdir, ydir, zdir;
        double localMobileVolumeFraction, localD, Ds;

        Voxel voxel;
        DiffusantParticle iParticle;

        if ((p == null) || (geometry.voxelSpacePBC == null)) {
            return Double.NaN;
        }

        //nRadiiInt = 2; // Rint
        //Rint = nRadiiInt * dv.radius;
        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * p.radius;
        totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        //sumV = 4.0 * Math.PI * dv.radius * dv.radius * dv.radius / 3.0;
        sumV = particleVolume; // all particles with same radius

        xVoxel = (int) geometry.computeVoxelX(p.x);
        yVoxel = (int) geometry.computeVoxelY(p.y);
        zVoxel = (int) geometry.computeVoxelZ(p.z);

        nVoxels = (int) Math.ceil(Rext / project.dx);

        x0 = xVoxel - nVoxels;
        y0 = yVoxel - nVoxels;
        z0 = zVoxel - nVoxels;

        x1 = xVoxel + nVoxels;
        y1 = yVoxel + nVoxels;
        z1 = zVoxel + nVoxels;

        p.step3Local = 0;

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    if (i < 0) {
                        ii = geometry.voxelSpacePBC.length + i;
                        xdir = -1;
                    } else if (i >= geometry.voxelSpacePBC.length) {
                        ii = i - geometry.voxelSpacePBC.length;
                        xdir = 1;
                    } else {
                        ii = i;
                        xdir = 0;
                    }

                    if (j < 0) {
                        jj = geometry.voxelSpacePBC[0].length + j;
                        ydir = -1;
                    } else if (j >= geometry.voxelSpacePBC[0].length) {
                        jj = j - geometry.voxelSpacePBC[0].length;
                        ydir = 1;
                    } else {
                        jj = j;
                        ydir = 0;
                    }

                    if (k < 0) {
                        kk = geometry.voxelSpacePBC[0][0].length + k;
                        zdir = -1;
                    } else if (k >= geometry.voxelSpacePBC[0][0].length) {
                        kk = k - geometry.voxelSpacePBC[0][0].length;
                        zdir = 1;
                    } else {
                        kk = k;
                        zdir = 0;
                    }

                    voxel = geometry.voxelSpacePBC[ii][jj][kk];
                    iParticle = voxel.firstParticleInChain;

                    while (iParticle != null) {

                        if ((iParticle != p) && iParticle.insideGeometry && iParticle.mobile) {

                            dx = p.x - (xdir * 2 * geometry.x2 + iParticle.x);
                            dy = p.y - (ydir * 2 * geometry.y2 + iParticle.y);
                            dz = p.z - (zdir * 2 * geometry.z2 + iParticle.z);

                            d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            r = iParticle.radius;

                            if (d < Rext - r) {
                                sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            } else if (d < Rext + r) {
                                sumV += Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                            }

                        }

                        iParticle = iParticle.nextParticleInChain;

                    }

                }
            }
        }

        localMobileVolumeFraction = sumV / totalV;
        Ds = DiffusantParticle.Dratio_shortNew(localMobileVolumeFraction, immobileVolumeFraction);
        //localD = dv.D * Ds;
        localD = Dcyto * Ds;
        p.step3Local = DiffusantParticle.step3(localD, project.dt);
        p.DsDff = Ds / DiffusantParticle.Dff_short_Banchio(localMobileVolumeFraction);

        return localMobileVolumeFraction;

    }

    public void initClusterImmobileCounter() {

        clusterImmobileCounter = Math.round(clusterImmobilePercent * numParticles) - immobileParticles;

        if (clusterImmobileCounter > 0) {
            clusterImmobileOn = true;
        }

    }

    public void clusterImmobileParticles() {

        double dx, dy, dz, sqrDistance, minDBV;
        double tetherDistance = 0.00005;

        DiffusantParticle iParticle;

        Voxel voxel, ivoxel;

        if (diffusants == null) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            d.shuffleParticles();

            for (DiffusantParticle p : d.particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                if (p.mobile) {
                    continue;
                }

                voxel = p.voxel;

                if (voxel == null) {
                    return;
                }

                for (int k = 0; k < voxel.numNeighbors; k++) {

                    ivoxel = voxel.neighbors[k];
                    iParticle = ivoxel.firstParticleInChain;

                    while (iParticle != null) {

                        if ((iParticle != p) && (iParticle.mobile)) {

                            dx = p.x - iParticle.x;
                            dy = p.y - iParticle.y;
                            dz = p.z - iParticle.z;

                            sqrDistance = dx * dx + dy * dy + dz * dz;
                            minDBV = p.radius + iParticle.radius + tetherDistance;

                            if (sqrDistance < minDBV * minDBV) {
                                iParticle.mobile = false;
                                clusterImmobileCounter--;
                                if (clusterImmobileCounter <= 0) {
                                    clusterImmobileOn = false;
                                    return;
                                }
                            }

                        }

                        iParticle = iParticle.nextParticleInChain;

                    }

                }

            }

        }

    }

    public double zMSD() {

        double sumMSD = 0, count = 0;

        DiffusantParticle p;

        for (int i = 0; i <= geometry.xVoxels - 1; i++) {
            for (int j = 0; j <= geometry.yVoxels - 1; j++) {
                for (int k = geometry.zVoxels - 2; k <= geometry.zVoxels - 1; k++) {

                    if (geometry.voxelSpacePBC != null) {
                        p = geometry.voxelSpacePBC[i][j][k].firstParticleInChain;
                    } else if (geometry.voxelSpace != null) {
                        p = geometry.voxelSpace[i][j][k].firstParticleInChain;
                    } else {
                        return Double.NaN;
                    }

                    while (p != null) {
                        sumMSD += p.squareDisplacement();
                        count++;
                        p = p.nextParticleInChain;
                    }

                }
            }
        }

        return (sumMSD / count);

    }

    public void saveVoxelResidenceTime() {

        if (residenceTimeCoordinates != null) {
            for (CoordinatesVoxels c : residenceTimeCoordinates) {
                voxelResidenceTime(c);
            }
        } else {
            saveVoxelResidenceTimeZavg();
        }

    }

    public double voxelResidenceTime(CoordinatesVoxels c) {

        double sum = 0, count = 0;
        double t, n;

        DiffusantParticle iParticle;

        if (c == null) {
            return Double.NaN;
        }

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

                    if (!geometry.voxelSpace[i][j][k].isSpace) {
                        continue;
                    }

                    t = geometry.voxelSpace[i][j][k].sumResidenceTime;
                    n = geometry.voxelSpace[i][j][k].residenceCounter;

                    iParticle = geometry.voxelSpace[i][j][k].firstParticleInChain;

                    while (iParticle != null) {
                        if (iParticle.mobile) {
                            t += iParticle.residenceTime;
                            n++;
                            iParticle = iParticle.nextParticleInChain;
                        }
                    }

                    if ((t > 0) && (n > 0)) {
                        sum += t / n;
                        count++;
                    }

                }
            }
        }

        Master.log("Average voxel residence time = " + (sum / count) + " ms");

        return (sum / count);

    }

    public void saveVoxelResidenceTimeZavg() {

        double sum, count, sumTotal = 0, countTotal = 0;
        double t, n;

        //int xyoffset = 2;
        int xyoffset = geometry.yVoxels / 2 - 2;

        Master.writeToFiles = true;

        Save save_ResidenceTime = new Save(project);

        save_ResidenceTime.saveWhileComputing = true;
        save_ResidenceTime.save2BinaryFile = false;
        save_ResidenceTime.save2TextFile = true;

        save_ResidenceTime.fileName("ResidenceTime", "voxels");
        save_ResidenceTime.xdim = project.spaceUnits;
        save_ResidenceTime.ydim = project.timeUnits;

        save_ResidenceTime.samples2save = geometry.zVoxels;

        if (!save_ResidenceTime.init("ResidenceTime", geometry, -1, 1)) {
            return;
        }

        CoordinatesVoxels c = new CoordinatesVoxels(project);

        c.xVoxel1 = 0 + xyoffset;
        c.xVoxel2 = geometry.xVoxels - 1 - xyoffset;
        c.yVoxel1 = 0 + xyoffset;
        c.yVoxel2 = geometry.yVoxels - 1 - xyoffset;
        c.zVoxel1 = 0;
        c.zVoxel2 = geometry.zVoxels - 1;

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {

            sum = 0;
            count = 0;

            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

                    if (geometry.voxelSpacePBC != null) {
                        t = geometry.voxelSpacePBC[i][j][k].sumResidenceTime;
                        n = geometry.voxelSpacePBC[i][j][k].residenceCounter;
                    } else if (geometry.voxelSpace != null) {
                        t = geometry.voxelSpace[i][j][k].sumResidenceTime;
                        n = geometry.voxelSpace[i][j][k].residenceCounter;
                    } else {
                        return;
                    }

                    if ((t > 0) && (n > 0)) {
                        sum += t / n;
                        count++;
                        sumTotal += t / n;
                        countTotal++;
                    }

                }
            }

            save_ResidenceTime.writeData(sum / count);

        }

        save_ResidenceTime.finish("voxels", geometry, -1);

        Master.log("Average voxel residence time = " + (sumTotal / countTotal) + " ms");

    }

    public void saveArray(double[] data, String fxn, String saveTag) {

        Master.writeToFiles = true;

        Save save = new Save(project);

        save.saveWhileComputing = true;
        save.save2BinaryFile = false;
        save.save2TextFile = true;

        save.fileName(fxn, saveTag);
        save.xdim = project.spaceUnits;
        save.ydim = project.timeUnits;

        save.samples2save = data.length;

        if (!save.init(fxn, geometry, -1, 1)) {
            return;
        }

        for (int i = 0; i < data.length; i++) {
            save.writeData(data[i]);
        }

        save.finish(fxn, geometry, -1);

    }

    public boolean addToVoxelList(DiffusantParticle p) {

        if (p.voxel == null) {
            return false;
        }

        p.nextParticleInChain = p.voxel.firstParticleInChain; // save existing particles
        p.voxel.firstParticleInChain = p; // replace with new particles, creating chain

        return true;

    }

    public boolean removeFromVoxelList(DiffusantParticle p, Voxel voxel) {

        if (voxel == null) {
            return false;
        }

        boolean found = false;

        DiffusantParticle p_test = voxel.firstParticleInChain;
        DiffusantParticle p_next;
        DiffusantParticle p_hold = null;

        while (p_test != null) {

            p_next = p_test.nextParticleInChain;

            if (p_test == p) {
                if (p_hold == null) {
                    voxel.firstParticleInChain = p_next; // dv was first in list
                } else {
                    p_hold.nextParticleInChain = p_next;
                }
                found = true;
            }

            p_hold = p_test;
            p_test = p_next;

        }

        if (found) {
            if (p.residenceTime > 0) {
                voxel.residenceCounter++;
                voxel.sumResidenceTime += p.residenceTime;
            }
            p.residenceTime = 0;
        }

        return found;

    }

    public boolean moveParticleGauss(DiffusantParticle p) {

        double step3;

        if (removingParticleOverlap) {
            step3 = removeParticleOverlapStep3;
        } else if (hydrodynamicsLocalD) {
            step3 = p.step3Local;
        } else {
            step3 = p.step3;
        }

        //if (step <= 0) {
        //    Master.exit("moveParticleGauss Error: bad step: " + step);
        //}
        stepx = Master.randomGauss() * step3;
        stepy = Master.randomGauss() * step3;
        stepz = Master.randomGauss() * step3;

        return setParticleLocation(p, p.x + stepx, p.y + stepy, p.z + stepz, false);

    }

    public boolean moveParticleGaussHydroWallz(DiffusantParticle p) {

        double step3;
        double z, b_ll = 1, b_T = 1;

        if (hydrodynamicsLocalD) {
            step3 = p.step3Local;
        } else if (hydrodynamicsDscale) {
            //step3 = localDensityFastHydroWallz(dv); // using Michailidou model instead
            step3 = p.step3;
        } else {
            step3 = p.step3;
        }

        if (removingParticleOverlap) {

            step3 = removeParticleOverlapStep3;

        } else {

            z = geometry.z2 - p.z;

            //b_ll = Math.sqrt(hydroWall_ll(dv.radius, z));
            //b_T = Math.sqrt(hydroWall_T(dv.radius, z - dv.radius));
            b_ll = hydroWall_ll(p.radius, z);
            b_T = hydroWall_T(p.radius, Math.abs(z - p.radius));

            b_ll = 1 / (1 + p.DsDff * ((1 / b_ll) - 1)); // Michailidou et al. 2009
            b_T = 1 / (1 + p.DsDff * ((1 / b_T) - 1)); // Michailidou et al. 2009

            b_ll = Math.sqrt(b_ll);
            b_T = Math.sqrt(b_T);

        }

        stepx = Master.randomGauss() * step3 * b_ll;
        stepy = Master.randomGauss() * step3 * b_ll;
        stepz = Master.randomGauss() * step3 * b_T;

        //stepx = step3 * b_ll;
        //stepy = step3 * b_ll;
        //stepz = step3 * b_T;
        //dv.localD = (stepx * stepx + stepy * stepy + stepz * stepz) / (6 * project.dt);
        //Master.log("" + (dv.localD/1.397657228105222E-4));
        //dv.localDw = (3 * step3 * step3) / (6 * project.dt);
        return setParticleLocation(p, p.x + stepx, p.y + stepy, p.z + stepz, false);

    }

    public boolean moveParticle(DiffusantParticle p) {

        double step3;

        if (hydrodynamicsLocalD) {
            step3 = p.step3Local;
        } else {
            step3 = p.step3;
        }

        //if (step < 0) {
        //    Master.exit("moveParticle Error: bad step: " + step);
        //}
        if (removingParticleOverlap) {
            step3 = removeParticleOverlapStep3;
        }

        if (Master.mt.nextDouble() < 0.5) {
            stepx = step3;
        } else {
            stepx = -step3;
        }

        if (Master.mt.nextDouble() < 0.5) {
            stepy = step3;
        } else {
            stepy = -step3;
        }

        if (Master.mt.nextDouble() < 0.5) {
            stepz = step3;
        } else {
            stepz = -step3;
        }

        return setParticleLocation(p, p.x + stepx, p.y + stepy, p.z + stepz, false);

    }

    public void moveParticlesCichocki(boolean moveAll) {

        boolean outOfBounds, overlap, sameVoxel, moveConnected, ok;

        DiffusantParticle odv;

        if (diffusants == null) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }
            
            //if (!removeParticleOverlap) {
            d.shuffleParticles();
            //}

            for (DiffusantParticle p : d.particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                if (!moveAll && !p.mobile) {
                    continue;
                }

                testParticle.copy(p);

                if (hydroWallZ) {
                    if (!moveParticleGaussHydroWallz(testParticle)) {
                        continue; // dont move
                    }
                } else {
                    if (!moveParticleGauss(testParticle)) {
                        continue; // dont move
                    }
                }

                outOfBounds = outOfBounds(testParticle);

                if (PBC && outOfBounds) {

                    outOfBounds = wrapAtBorder(testParticle);

                    //if (outOfBounds) {
                    //
                    //Master.log("moveParticlesCichocki error: PBC wrap at border error");
                    //return -1;
                    //} else {
                    //Master.log("wrapped particles");
                    //}
                }

                if (outOfBounds) {
                    continue; // dont move
                }

                if (testOverlap(testParticle)) {
                    if ((time > 0) && (p.firstCollision == 0)) {
                        p.firstCollision = time;
                        p.sqrDisplacement = p.squareDisplacement();
                    }
                    continue; // dont move
                }

                if (!freeDiffusion) {

                    odv = testParticleOverlap(testParticle, p);

                    overlap = odv != null;

                    if (overlap) {
                        //Master.log("square displacement = " + diffusant[k].particle[j].squareDisplacement());
                        //Master.log("time = " + time);
                        //cancel = true;
                        //if (!removingParticleOverlap && diffusant[k].saveDisplacementAfterCollisions && (diffusant[k].particle[j].lastParticleCollision != odv)) {
                        if (!removingParticleOverlap && d.saveDisplacementAfterCollisions) {
                            p.sqrDisplacement = p.squareDisplacement();
                            p.lastParticleCollision = odv;
                            p.x0 = p.x;
                            p.y0 = p.y;
                            p.z0 = p.z;
                        }
                        if ((time > 0) && (p.firstCollision == 0)) {
                            p.firstCollision = time;
                            p.sqrDisplacement = p.squareDisplacement();
                        }
                        continue;
                    }

                }

                if (connectParticles && (p.connectTo != null)) {

                    moveConnected = true;

                    for (DiffusantParticle p2 : p.connectTo) {

                        if (p2 == null) {
                            continue;
                        }

                        if (!testParticle.overlap(p2, connectorLength)) {
                            moveConnected = false;
                            break;
                        }

                    }

                    if (!moveConnected) {
                        continue;
                    }

                }

                sameVoxel = false;

                if (testParticle.voxel == p.voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(p, p.voxel);
                }

                p.copy(testParticle);

                if (!sameVoxel) {
                    addToVoxelList(p);
                }

            }

        }

    }

    public void drift() {
        boolean outOfBounds, sameVoxel, remove;

        if (diffusants == null) {
            return;
        }

        //if (!PBC) {
            //Master.exit("drift error: PBC is not on");
        //}
        
        driftGeometry(driftx, drifty, driftz);

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                testParticle.copy(p);

                setParticleLocation(testParticle, testParticle.x + driftx, testParticle.y + drifty, testParticle.z + driftz, false);

                outOfBounds = outOfBounds(testParticle);
                
                remove = false;

                if (outOfBounds) {
                    if (PBC) {
                        outOfBounds = wrapAtBorder(testParticle);
                        //testParticle.fluorescence = 1.0; // reset F to simulate non-frapped particles moving into PSF
                    } else {
                        remove = true;
                        outOfBounds = false;
                    }
                }

                if (outOfBounds) {
                    Master.exit("drift error: wrap at border error");
                    //continue;
                }

                sameVoxel = false;

                if (testParticle.voxel == p.voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(p, p.voxel);
                }

                p.copy(testParticle);
                
                if (remove) {
                    p.insideGeometry = false;
                }

                if (!sameVoxel) {
                    addToVoxelList(p);
                }

            }

        }

    }
    
    public void driftGeometry(double dx, double dy, double dz) {
        // generic function for creating drift of geometric structures
    }

    public boolean wrapAtBorder(DiffusantParticle p) {

        double x = p.x;
        double y = p.y;
        double z = p.z;

        if (x > geometry.x2) {
            x -= 2 * geometry.x2;
            p.x0 -= 2 * geometry.x2;
        } else if (x < geometry.x1) {
            x -= 2 * geometry.x1;
            p.x0 -= 2 * geometry.x1;
        }

        if (y > geometry.y2) {
            y -= 2 * geometry.y2;
            p.y0 -= 2 * geometry.y2;
        } else if (y < geometry.y1) {
            y -= 2 * geometry.y1;
            p.y0 -= 2 * geometry.y1;
        }

        if (z > geometry.z2) {
            z -= 2 * geometry.z2;
            p.z0 -= 2 * geometry.z2;
        } else if (z < geometry.z1) {
            z -= 2 * geometry.z1;
            p.z0 -= 2 * geometry.z1;
        }

        setParticleLocation(p, x, y, z, false);

        return outOfBounds(p);

    }

    public boolean initSave() {
        project.checkDirectory();
        return true;
    }

    public boolean finishSave() {
        return true;
    }

    @Override
    public boolean addUser(ParamVector pv) {

        if (pv == null) {
            return false;
        }

        super.addUser(pv);

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {
        super.updateVector(v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof RunMonteCarlo)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("minParticleStep")) {
            if (v <= 0) {
                return false;
            }
            minParticleStep = v;
            return true;
        }
        if (n.equalsIgnoreCase("freeDiffusion")) {
            freeDiffusion = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("PBC")) {
            PBC = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("connectParticles")) {
            connectParticles = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("maxNumConnectors")) {
            if (v < 0) {
                return false;
            }
            maxNumConnectors = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("connectRate")) {
            if (v < 0) {
                return false;
            }
            connectRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("unconnectRate")) {
            if (v < 0) {
                return false;
            }
            unconnectRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("hydrodynamicsLocalD")) {
            hydrodynamicsLocalD = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("hydroWallZ")) {
            hydroWallZ = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("MSDspatialTbgn")) {
            MSDspatialTbgn = v;
            return true;
        }
        return super.setMyParams(o, v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof RunMonteCarlo)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        return super.setMyParams(o, s);
    }
}
