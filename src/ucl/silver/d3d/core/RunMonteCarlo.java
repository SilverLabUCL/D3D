package ucl.silver.d3d.core;

import ucl.silver.d3d.gui.*;
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
public class RunMonteCarlo
        extends ParamVector implements Runnable {

    public double minVesicleRadius; // um
    public double maxVesicleRadius; // um
    public double minAllowedDistanceBetweenVesicles; // um
    public double minDistanceBetweenVesicles; // um

    public int numVesicles; // total number of vesicles
    public double vesicleDensity; // vesicles / um^3
    public double vesicleVolumeFraction;
    public double totalVesicleVolume; // um^3

    public boolean initVesicleRandom = true;
    public boolean vesicleLattice = false;

    public double Dcyto = 0; // um^2/ms

    public double voxelVolume = 0; // um^3

    public int mobileVesicles = 0;
    public int immobileVesicles = 0;
    public double mobilePercent = 0;
    public double immobilePercent = 0;
    public double mobileVolumeFraction = 0;
    public double immobileVolumeFraction = 0;
    public double clusterImmobilePercent = 0;
    public boolean clusterImmobileVesicles = false;

    public boolean connectVesicles = false;
    public boolean connectorsBinomial = false;
    public int maxNumConnectors = 5;
    public double meanNumConnectors = 1.5;
    public double connectorLength = 0.01; // um
    public double connectRate = 0.7 / 0.1; // P/ms
    public double unconnectRate = 0.175 / 0.1; // P/ms
    public double avgConnectorLifeTime = 0;
    public int connectorLifeTimeCounter = 0;
    public double connectorAllOffTime = 0; // ms

    public double minVesicleStep; // um
    public double removeVesicleOverlapStep3 = 0.001 / Math.sqrt(3); // um

    //public boolean checkForOverlapsDuringSimulation = false; // set true for testing only!
    
    public boolean PBC = false; // periodic boundary conditions
    public boolean freeDiffusion = false; // vesicles are allowed to overlap
    public boolean removeVesicleOverlap = true;

    public boolean driftOn = false;
    public double driftRateX = 0.0046 / 1000.0; // um/ms
    public double driftRateY = 0.0046 / 1000.0; // um/ms
    public double driftRateZ = 0.0116 / 1000.0; // um/ms (measured May 2012)
    public double driftOnset = 0; // ms
    private double driftx, drifty, driftz; // um

    // simulation variables
    public transient Grid grid = null; // Panel2D voxel grid
    public transient DiffusantVesicles[] diffusants = null;
    
    public transient Thread thread = null;
    public transient Geometry geometry;

    public int itime;
    public double time;
    public boolean timer = true;
    public double stepx, stepy, stepz;

    public boolean batchesExist = false;
    public boolean autoInit = true;
    public boolean runSimulation = false;
    public boolean cancel = false;
    public boolean initialized = false;
    public boolean preview = false;
    public boolean removingVesicleOverlap = false;

    private boolean clusterImmobileOn = false;
    private long clusterImmobileCounter = 0;

    public boolean hydrodynamicsLocalD = false; // hydrodynamic interactions local density computation
    public boolean hydrodynamicsDscale = false; // hydrodynamic interactions scale factor
    public boolean hydroWallZ = false;
    public boolean hydrodynamicsLocalDVoxels = false;

    public int overlapTrialLimit = (int) 1e4;
    public int absLimit = (int) 1e6;

    public transient StopWatch timer1 = new StopWatch();
    public transient StopWatch timer2 = new StopWatch();

    public transient MersenneTwisterFast mt;
    public long seed;

    private double saveRanGauss = 999999; // for Gaussian random number generator

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

    private final DiffusantVesicle testVesicle = new DiffusantVesicle(project, "ready", 0, 0, 0, 0, 0);

    public double vesicleVolume = Double.NaN; // used for local density computation, but assumes all vesicles have same volume

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("minVesicleRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxVesicleRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minAllowedDistanceBetweenVesicles")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minDistanceBetweenVesicles")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("minVesicleStep")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("totalVesicleVolume")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("vesicleDensity")) {
            return "vesicles/" + project.spaceUnits + "^3";
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
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("minVesicleRadius")) {
            return false;
        }
        if (name.equalsIgnoreCase("maxVesicleRadius")) {
            return false;
        }
        if (name.equalsIgnoreCase("minAllowedDistanceBetweenVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("minDistanceBetweenVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("totalVesicleVolume")) {
            return false;
        }
        if (name.equalsIgnoreCase("VesicleDensity")) {
            return false;
        }
        if (name.equalsIgnoreCase("vesicleVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("numVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("mobileVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("immobileVesicles")) {
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

        seed = System.currentTimeMillis();
        seed += Math.random() * 10000000;
        mt = new MersenneTwisterFast(seed); // init random number generator

        Master.log("Mersenne Twister seed = " + seed);

    }

    public void init() {
        if (project.numBatches() > 0) {
            batchesExist = true;
        }
    }

    public void setOutputRate(double newRate) {
    }

    public void startSimulation(boolean PREVIEW) {

        preview = PREVIEW;

        grid = Master.grid();

        if ((grid != null) && preview) {
            grid.preview = true;
            //voxelGrid.diffusant = diffus;
            Master.mainframe.panel2D.displayModeSet("Diffusant.0");
        }

        runSimulation = true;
        cancel = false;
        thread = new Thread(this);
        thread.start(); // this calls run() below

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

        for (DiffusantVesicles d : diffusants) {
            Master.log("DiffusantVesicle #" + j + ", avg D = " + d.D);
            Master.log("DiffusantVesicle #" + j + ", estimated meanStep = " + Math.sqrt(6.0 * d.D * project.dt));
            Master.log("DiffusantVesicle #" + j + ", avg meanStep = " + (d.meanStep3 * Math.sqrt(3)));
            Master.log("DiffusantVesicle #" + j + ", volume fraction = " + d.volumeFraction);
            Master.log("DiffusantVesicle #" + j + ", immobilePercent = " + d.immobilePercent);
            j++;
        }

    }

    public void finishSimulation() {

        finishSave();
        vesicleStats();

        if (connectVesicles) {
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

                saveVesiclePositions(0);

            }

            if (time < simTime) {
                Master.log("starting Monte Carlo simulation at time " + time + " " + project.timeUnits);
            }

            while (runSimulation && !cancel && (time < simTime)) {

                for (DiffusantVesicles d : diffusants) {
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

                moveVesiclesCichocki(false);
                //testConnectorSeperation();

                if (hydrodynamicsLocalD) {
                    localDensityAll(false);
                }

                if (connectVesicles) {
                    connectVesicles();
                    unconnectVesicles();
                }

                if (clusterImmobileOn) {
                    clusterImmobileVesicles();
                }
                
                if (driftOn && (time > driftOnset)) {
                    drift();
                }

                //if (!freeDiffusion && checkForOverlapsDuringSimulation && testVesicleOverlap()) {
                //    Master.log("Monte Carlo Simulation has overlapping vesicles: " + minDistanceBetweenVesicles);
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

                if (timer) {
                    timer2.timer(time);
                }

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

        saveVesiclePositions((int) project.simTime);

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

        Master.log("finished Monte Carlo simulation. time = " + timer1.toString());

        if (timer) {
            timer2.stop();
        }

    }
    
    public boolean initGeometry() {
        return false; // generic method to initialize geometry
    }

    public boolean initAll() {

        if (checkVariables()) {
            return true;
        }

        if (initDiffusantVesiclesArray()) {
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

        if (initVesicles()) {
            return true;
        }

        if (checkDX()) {
            return true;
        }

        if (initDT()) {
            return true;
        }

        if (vesicleLattice) {
            if (initVesiclesLattice()) {
                return true;
            }
        } else if (initVesicleRandom) {
            if (initVesiclesRandom()) {
                return true;
            }
        }

        if (initVesiclesImmobile()) {
            return true;
        }

        if (removeVesicleOverlap && removeVesicleOverlap()) {
            return true;
        }

        initVesicleStartLocations();

        initConnectors();

        vesicleStats();

        if (clusterImmobileVesicles) {
            initClusterImmobileCounter();
        }

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        if (maxNumConnectors == 0) {
            connectVesicles = false;
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

    public boolean initDiffusantVesiclesArray() {

        int count = 0;

        if (project.diffusants == null) {
            error("No Diffusant Vesicles!");
            return true;
        }

        for (Diffusant d : project.diffusants) {
            if (d instanceof DiffusantVesicles) {
                count++;
            }
        }

        if (count == 0) {
            error("initDiffusantVesiclesArray", "diffusants", "no diffusant vesicles");
            return true;
        }

        diffusants = new DiffusantVesicles[project.diffusants.length];

        for (int i = 0; i < project.diffusants.length; i++) {

            if (project.diffusants[i] instanceof DiffusantVesicles) {
                diffusants[i] = (DiffusantVesicles) project.diffusants[i];
            } else {
                diffusants[i] = null;
            }
        }

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

    public boolean initVesicles() {

        double hydroDsDff;

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

            for (DiffusantVesicles d : diffusants) {
                if (d != null) {
                    Master.log("hydrodynamic scaling old D = " + Dcyto);
                    d.D = Dcyto * DiffusantVesicle.Dratio_shortNew(d.mobileVolumeFraction, d.immobileVolumeFraction);
                    Master.log("hydrodynamic scaling new D = " + d.D);
                }
            }

        }

        for (DiffusantVesicles d : diffusants) {

            if (d != null) {

                d.initVesicles();

                if ((d.xyz != null) && (d.vesicles.length == d.xyz.length)) {
                    for (int j = 0; j < d.vesicles.length; j++) {
                        setVesicleLocation(d.vesicles[j], d.xyz[j][0], d.xyz[j][1], d.xyz[j][2], true);
                        addToVoxelList(d.vesicles[j]);
                    }
                }

            }

        }

        if (hydroWallZ) {

            for (DiffusantVesicles d : diffusants) {

                if ((d == null) || (d.vesicles == null)) {
                    continue;
                }

                hydroDsDff = DiffusantVesicle.Dratio_shortNew(d.mobileVolumeFraction, d.immobileVolumeFraction);
                hydroDsDff /= DiffusantVesicle.Dff_short_Banchio(d.setVolumeFraction * (1 - d.setImmobilePercent));

                for (DiffusantVesicle v : d.vesicles) {
                    v.DsDff = hydroDsDff;
                }

            }

        }

        return false;

    }

    public boolean initVesiclesImmobile() {

        if (diffusants == null) {
            return true;
        }

        for (DiffusantVesicles d : diffusants) {

            if (d == null) {
                continue;
            }

            if (d.initImmobileVesicles(null)) {
                return true;
            }

        }

        return false;

    }

    public boolean saveVesicleDensity(String saveTag) {
        return false;
    }

    public boolean saveVesiclePositions(int msec) {

        if (preview || (diffusants == null)) {
            return true;
        }

        for (DiffusantVesicles d : diffusants) {

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

    public void vesicleStats() {

        double distance;
        double spaceVolume = spaceVolume();

        minVesicleRadius = Double.POSITIVE_INFINITY;
        maxVesicleRadius = 0;
        minAllowedDistanceBetweenVesicles = Double.POSITIVE_INFINITY;

        totalVesicleVolume = 0;

        numVesicles = 0;
        mobileVesicles = 0;
        immobileVesicles = 0;

        vesicleDensity = 0;
        vesicleVolumeFraction = 0;

        if (diffusants == null) {
            return;
        }

        for (DiffusantVesicles d : diffusants) {

            if (d == null) {
                continue;
            }

            //diffusants[k].init();
            d.vesicleStats();

            minVesicleRadius = Math.min(minVesicleRadius, d.minRadius);
            maxVesicleRadius = Math.max(maxVesicleRadius, d.maxRadius);

            totalVesicleVolume += d.totalVolume;

            numVesicles += d.numVesicles;
            mobileVesicles += d.mobileVesicles;
            immobileVesicles += d.immobileVesicles;

        }

        mobilePercent = 1.0 * mobileVesicles / (immobileVesicles + mobileVesicles);
        immobilePercent = 1.0 * immobileVesicles / (immobileVesicles + mobileVesicles);

        for (DiffusantVesicles d1 : diffusants) {

            if (d1 == null) {
                continue;
            }

            for (DiffusantVesicles d2 : diffusants) {

                if (d2 == null) {
                    continue;
                }

                distance = d1.minRadius + d2.minRadius;
                minAllowedDistanceBetweenVesicles = Math.min(distance, minAllowedDistanceBetweenVesicles);

            }
        }

        vesicleVolume = 4 * Math.PI * minVesicleRadius * minVesicleRadius * minVesicleRadius / 3.0;

        vesicleDensity = (numVesicles * 1.0) / spaceVolume;
        vesicleVolumeFraction = totalVesicleVolume / spaceVolume;

        mobileVolumeFraction = mobilePercent * vesicleVolumeFraction;
        immobileVolumeFraction = immobilePercent * vesicleVolumeFraction;

        setParamObject("minVesicleRadius", minVesicleRadius);
        setParamObject("maxVesicleRadius", maxVesicleRadius);
        setParamObject("minAllowedDistanceBetweenVesicles", minAllowedDistanceBetweenVesicles);

        setParamObject("vesicleVolume", totalVesicleVolume);

        setParamObject("numVesicles", numVesicles);
        setParamObject("mobileVesicles", mobileVesicles);
        setParamObject("immobileVesicles", immobileVesicles);
        setParamObject("mobilePercent", mobilePercent);
        setParamObject("immobilePercent", immobilePercent);
        setParamObject("mobileVolumeFraction", mobileVolumeFraction);
        setParamObject("immobileVolumeFraction", immobileVolumeFraction);

        setParamObject("vesicleDensity", vesicleDensity);
        setParamObject("vesicleVolumeFraction", vesicleVolumeFraction);

        minDistance();

    }

    public boolean checkDX() {

        double maxRadius = 0;

        if (diffusants == null) {
            return true;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (v == null) {
                    continue;
                }

                maxRadius = Math.max(maxRadius, v.radius);

            }

        }

        if (maxRadius <= 0) {
            Master.exit("MonteCarlo: checkDX: bad maxRadius: " + maxRadius);
        }

        if ((int) (100 * project.dx) < (int) (100 * 2 * maxRadius)) {
            //Master.exit("MonteCarlo: checkDX: voxel size is smaller than maximum vesicle diameter: " + (2 * maxRadius));
        }

        return false;

    }

    public boolean initDT() {

        double dt, maxD = 0;

        if (diffusants == null) {
            return true;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {
                maxD = Math.max(maxD, v.D);
            }

        }

        if (maxD <= 0) {
            Master.exit("initDT: negative maxD: " + maxD);
            return true;
        }

        dt = minVesicleStep * minVesicleStep / (6.0 * maxD);
        project.dt = dt;
        project.stability = dt * 3 * maxD / (project.dx * project.dx);

        // can only compute stability since maxD, dt and dx are set elsewhere
        //Master.log("maxD = " + maxD);
        //Master.log("Monte Carlo dt = " + dt);
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

    public boolean setVesicleLocation(DiffusantVesicle v, double x, double y, double z, boolean initStartLocation) {

        int xVoxel, yVoxel, zVoxel;

        if (PBC) {

            //if (geometry.voxelSpacePBC == null) {
            //    return false;
            //}
            v.x = x;
            v.y = y;
            v.z = z;

            if (initStartLocation) {
                v.x0 = x;
                v.y0 = y;
                v.z0 = z;
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

            v.voxel = geometry.voxelSpacePBC[xVoxel][yVoxel][zVoxel];

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

            v.x = x;
            v.y = y;
            v.z = z;

            if (initStartLocation) {
                v.x0 = x;
                v.y0 = y;
                v.z0 = z;
            }

            v.voxel = geometry.voxelSpace[xVoxel][yVoxel][zVoxel];

            return true;

        }

    }

    public boolean outOfBounds(DiffusantVesicle v) {

        if (geometry.simpleCuboid) {
            return beyondLimits(v); // simple computation
        } else {
            if (beyondLimits(v)) {
                return true;
            }
            return insideVoxelNonSpace(v); // check nonspace voxels
        }

    }

    public boolean beyondLimits(DiffusantVesicle v) {

        double extra;

        if (PBC || freeDiffusion) {
            extra = 0;
        } else {
            extra = v.radius;
        }

        //if (ellipsoid) {
        //    test = x * x / (lx1 * lx1) + y * y / (ly1 * ly1) + z * z / (lz1 * lz1);
        //    if (test > 1) {
        //        return true;
        //    }
        //}
        if ((v.x < geometry.x1 + extra) || (v.x > geometry.x2 - extra)) {
            return true;
        }

        if ((v.y < geometry.y1 + extra) || (v.y > geometry.y2 - extra)) {
            return true;
        }

        if ((v.z < geometry.z1 + extra) || (v.z > geometry.z2 - extra)) {
            return true;
        }

        return false;

    }

    public boolean insideVoxelNonSpace(DiffusantVesicle v) {

        double dlimit;
        Voxel ivoxel;

        if (v.voxel == null) {
            return false;
        }

        if (!v.voxel.isSpace) {
            return true;
        }

        if (freeDiffusion) {
            return false;
        }

        dlimit = project.dx * 0.5 + v.radius; // um

        for (int i = 0; i < v.voxel.numNonSpaceNeighbors; i++) {

            ivoxel = v.voxel.nonSpaceNeighbors[i];

            if ((v.x > ivoxel.x - dlimit) && (v.x < ivoxel.x + dlimit)) {
                if ((v.y > ivoxel.y - dlimit) && (v.y < ivoxel.y + dlimit)) {
                    if ((v.z > ivoxel.z - dlimit) && (v.z < ivoxel.z + dlimit)) {
                        //Master.log("nonspace overlap " + ivoxel.x + "," + ivoxel.y + "," + ivoxel.z);
                        return true;
                    }
                }
            }

        }

        return false;

    }

    public void initVesicleStartLocations() {

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            d.initStartLocation();

        }

    }

    public boolean initVesicleRandom(DiffusantVesicle v, Coordinates c, boolean allowOverlaps, int overlapTrialLimit, int absLimit) {

        int trial = 0;
        double x, y, z;
        boolean overlap;

        if ((v == null) || (c == null)) {
            return true;
        }

        while (true) {

            trial++;

            if (trial >= absLimit) {
                return true;
            }

            x = (mt.nextDouble() * (c.x2 - c.x1)) + c.x1;
            y = (mt.nextDouble() * (c.y2 - c.y1)) + c.y1;
            z = (mt.nextDouble() * (c.z2 - c.z1)) + c.z1;

            if (!setVesicleLocation(v, x, y, z, true)) {
                continue;
            }

            if (outOfBounds(v)) {
                continue;
            }

            if (testOverlap(v)) {
                continue;
            }

            overlap = false;

            if (!freeDiffusion) {
                overlap = (testVesicleOverlap(v, null) != null);
            }

            if (!overlap || (allowOverlaps && (trial > overlapTrialLimit))) {
                return (!addToVoxelList(v));
            }

        }

    }
    
    public boolean testOverlap(DiffusantVesicle v) {
        return false; // generic function to test vesicle overlap
    }

    public boolean initVesiclesRandom() {

        int nVesicles;

        boolean allowOverlaps = true;

        Coordinates c;

        if (diffusants == null) {
            return true;
        }

        int k;

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            if (d.coordinates == null) {
                c = new Coordinates(project, geometry);
            } else {
                c = new Coordinates(project, d.coordinates);
            }

            nVesicles = d.vesicles.length;

            Master.log("randomly placing ready vesicles (" + nVesicles + ")");

            k = 0;

            for (DiffusantVesicle v : d.vesicles) {

                if (!v.insideGeometry) {
                    continue;
                }

                if (initVesicleRandom(v, c, allowOverlaps, overlapTrialLimit, absLimit)) {
                    return true;
                }

                if ((k > 0) && (Math.IEEEremainder(k, 10000) == 0)) {
                    Master.log("placed vesicle " + k + " / " + nVesicles);
                }

                k++;

            }

        }

        return false;

    }

    public boolean initVesiclesLattice() {

        double newVesicleVolume, r, jj;
        double x, y, z;
        double xlimit, ylimit, zlimit;
        int nVesicles, iv;
        int i, j, k, ii;

        double maxPacking = Math.PI / (3 * Math.sqrt(2));
        double zstep = Math.sqrt(6) * 2.0 / 3.0;
        double sqrt3 = Math.sqrt(3);

        boolean overlap, finished = false;

        if (diffusants == null) {
            return true;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            nVesicles = d.vesicles.length;
            newVesicleVolume = maxPacking * 1.0 / nVesicles;
            r = Math.pow((3 * newVesicleVolume / (4 * Math.PI)), (1 / 3.0));

            //Master.log("new radius " + r);
            if (PBC) {
                xlimit = geometry.x2;
                ylimit = geometry.y2;
                zlimit = geometry.z2;
            } else {
                xlimit = geometry.x2 - d.setMeanRadius;
                ylimit = geometry.y2 - d.setMeanRadius;
                zlimit = geometry.z2 - d.setMeanRadius;
            }

            iv = 0;
            k = 0;

            while (!finished) {

                z = geometry.z1 + r + k * zstep * r;

                if (z > zlimit) {
                    break;
                }

                j = 0;
                jj = (k % 2) / 3.0;

                while (!finished) {

                    y = geometry.y1 + r + sqrt3 * (j + jj) * r;

                    if (y > ylimit) {
                        break;
                    }

                    i = 0;
                    ii = (j + k) % 2;

                    while (!finished) {

                        x = geometry.x1 + r + (2 * i + ii) * r;

                        if (x > xlimit) {
                            break;
                        }

                        if (!setVesicleLocation(d.vesicles[iv], x, y, z, true)) {
                            return true;
                        }

                        overlap = testVesicleOverlap(d.vesicles[iv], null) != null;

                        if (!overlap) {

                            //if ((k % 3 == 1) && (j % 3 == 1) && (k % 3 == 1)) {
                            //    diffusant[id].vesicle[iv].mobile = true;
                            //} else {
                            //    diffusant[id].vesicle[iv].mobile = false;
                            //}
                            if (!addToVoxelList(d.vesicles[iv])) {
                                return true;
                            }

                            if (iv == d.vesicles.length - 1) {
                                finished = true;
                                break;
                            }

                            iv++;

                        }

                        i++;

                    }

                    j++;

                }

                k++;

            }

            Master.log("created vesicle lattice (" + iv + "/" + nVesicles + ")");

            for (int iiv = iv; iiv < d.vesicles.length; iiv++) {
                d.vesicles[iiv].insideGeometry = false;
                d.vesicles[iiv].mobile = false;
            }

        }

        return false;

    }

    public boolean removeVesicleOverlap() {

        int radiusNM, lastRadiusNM = 0, count = 0;
        double halfDistance = 0;

        if (freeDiffusion) {
            return false; // OK
        }

        if (!testVesicleOverlap()) {
            return false; // OK
        }

        if (diffusants == null) {
            return true;
        }

        vesicleStats(); // will update vesicles step size

        removingVesicleOverlap = true;

        Master.log("elimination of vesicle overlap...");

        while (halfDistance < minVesicleRadius) {

            halfDistance = 0.5 * minDistance();

            radiusNM = (int) (1000 * halfDistance);

            if (radiusNM > lastRadiusNM) {
                Master.log("current radius: " + radiusNM + " nm");
                lastRadiusNM = radiusNM;
                count = 0;
            }

            // set current diameter to minimal distance between vesicles
            for (DiffusantVesicles d : diffusants) {

                //if ((diffusant[k] == null) || (diffusant[k].vesicle == null)) {
                //    continue;
                //}
                for (DiffusantVesicle v : d.vesicles) {
                    if (v.insideGeometry) {
                        //v.radius = Math.min(halfDistance, diffusant[k].meanRadius);
                        v.radius = halfDistance;
                    }
                }

            }

            moveVesiclesCichocki(true);

            count++;

            if (count > absLimit) {
                break;
            }

        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            //radius = diffusant[k].meanRadius;
            for (DiffusantVesicle v : d.vesicles) {
                if (v.insideGeometry) {
                    v.radius = d.meanRadius;
                    v.x0 = v.x;
                    v.y0 = v.y;
                    v.z0 = v.z;
                }
            }

        }

        removingVesicleOverlap = false;

        if (testVesicleOverlap()) {
            error("removeVesicleOverlap: failed to remove vesicle overlap.");
            return true;
        }

        return initVesiclesImmobile();

    }

    public double minDistanceZstack() {

        double sqrDistance, minSqrDist, minSqrDist2;
        double DX, DY, DZ;

        minSqrDist = Double.POSITIVE_INFINITY;
        minDistanceBetweenVesicles = Double.NaN;

        boolean found = false;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v1 : d.vesicles) {

                if (!v1.insideGeometry) {
                    continue;
                }

                minSqrDist2 = Double.POSITIVE_INFINITY;

                for (DiffusantVesicle v2 : d.vesicles) {

                    if (v1 == v2) {
                        continue;
                    }

                    if (!v2.insideGeometry) {
                        continue;
                    }

                    if (v1.z != v2.z) {
                        continue;
                    }

                    DX = v1.x - v2.x;
                    DY = v1.y - v2.y;
                    DZ = v1.z - v2.z;

                    //sqrDistance = DX * DX + DY * DY + DZ * DZ;
                    sqrDistance = DX * DX + DY * DY;

                    if (sqrDistance == 0) {
                        Master.log("warning: zero sqrDistance " + v1 + ", " + v2);
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
            minDistanceBetweenVesicles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenVesicles = Double.NaN;
        }

        return minDistanceBetweenVesicles;

    }

    public double minDistanceSlow() {

        int countOverlaps = 0;
        double sqrDistance, minSqrDist, minSqrDist2;
        double DX, DY, DZ;

        minSqrDist = Double.POSITIVE_INFINITY;
        minDistanceBetweenVesicles = Double.NaN;

        boolean found = false;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (int j = 0; j < d.vesicles.length; j++) {

                if (!d.vesicles[j].insideGeometry) {
                    continue;
                }

                minSqrDist2 = Double.POSITIVE_INFINITY;

                for (int k = j + 1; k < d.vesicles.length; k++) {

                    //if (j == k) {
                    //    continue;
                    //}
                    if (!d.vesicles[k].insideGeometry) {
                        continue;
                    }

                    DX = d.vesicles[j].x - d.vesicles[k].x;
                    DY = d.vesicles[j].y - d.vesicles[k].y;
                    DZ = d.vesicles[j].z - d.vesicles[k].z;

                    sqrDistance = DX * DX + DY * DY + DZ * DZ;

                    if (sqrDistance == 0) {
                        Master.log("warning: zero sqrDistance " + j + ", " + k);
                    }

                    if (Math.sqrt(sqrDistance) < 0.04) {
                        if ((Math.abs(DZ) <= 0.03) && (Math.abs(DZ) > 0)) {
                            if ((Math.abs(DX) < 0.01) && (Math.abs(DY) < 0.01)) {
                                //Master.log("" + DX + ", " + DY + ", " + DZ);
                                //Master.log("" + (Math.sqrt(DX * DX + DY * DY)));
                                d.vesicles[k].insideGeometry = false;
                                d.vesicles[k].x = Double.NaN;
                                d.vesicles[k].y = Double.NaN;
                                d.vesicles[k].z = Double.NaN;
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
            minDistanceBetweenVesicles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenVesicles = Double.NaN;
        }

        return minDistanceBetweenVesicles;

    }

    public double minDistance() {

        double min, minSqrDist = Double.MAX_VALUE;

        boolean found = false;

        minDistanceBetweenVesicles = Double.NaN;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (v == null) {
                    continue;
                }

                if (v.insideGeometry) {

                    min = minSqrDistance(v, v);

                    if (min >= 0) {

                        if (min < minSqrDist) {
                            minSqrDist = min;
                            found = true;
                        }

                    } else {
                        //Master.log("failed to find min distance for vesicles " + j);
                    }

                }

            }

        }

        if (found) {
            minDistanceBetweenVesicles = Math.sqrt(minSqrDist);
        } else {
            minDistanceBetweenVesicles = minDistanceSlow();
        }

        return minDistanceBetweenVesicles;

    }

    public void testMinDistance() {

        int counter = 0, isteps = 2000;
        double sec, min;

        if (diffusants == null) {
            return;
        }

        timer1.start();

        for (int istep = 0; istep < isteps; istep++) {

            for (DiffusantVesicles d : diffusants) {

                if ((d == null) || (d.vesicles == null)) {
                    continue;
                }

                for (DiffusantVesicle v : d.vesicles) { // runSimulation thru each vesicles

                    if (v == null) {
                        continue;
                    }

                    if (v.insideGeometry) {
                        min = minSqrDistance(v, v);
                        counter++;
                    }

                }

            }

        }

        timer1.stop();
        sec = (double) (timer1.milliseconds / 1000.0);
        Master.log("Test Time (" + counter + ") = " + sec + " sec");

    }

    public double minSqrDistance(DiffusantVesicle v, DiffusantVesicle ignoreVesicle) {

        double dx, dy, dz, sqrDistance, min = Double.MAX_VALUE;
        int ii, jj, kk;
        boolean found = false;

        DiffusantVesicle ivesicle;

        Voxel voxel = v.voxel;
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
            ivesicle = ivoxel.firstReady;

            while (ivesicle != null) {

                if ((ivesicle != v) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry) {

                    dx = v.x - ivesicle.x;
                    dy = v.y - ivesicle.y;
                    dz = v.z - ivesicle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    if (sqrDistance < min) {
                        min = sqrDistance;
                        found = true;
                    }

                    if (sqrDistance <= 0) {
                        Master.log("warning: zero sqrDistance " + v + " and " + ivesicle);
                    }

                }

                ivesicle = ivesicle.nextReady;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];
                ivesicle = ivoxel.firstReady;

                while (ivesicle != null) {

                    if ((ivesicle != v) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry) {

                        ii = voxelPBC.PBCi[i];
                        jj = voxelPBC.PBCj[i];
                        kk = voxelPBC.PBCk[i];

                        dx = v.x - (ii * 2 * geometry.x2 + ivesicle.x);
                        dy = v.y - (jj * 2 * geometry.y2 + ivesicle.y);
                        dz = v.z - (kk * 2 * geometry.z2 + ivesicle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        if (sqrDistance < min) {
                            min = sqrDistance;
                            found = true;
                        }

                    }

                    ivesicle = ivesicle.nextReady;

                }

            }

        }

        if (found) {
            return min;
        } else {
            return -1;
        }

    }

    public boolean testVesicleOverlap() {

        if (diffusants == null) {
            return false;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (v == null) {
                    continue;
                }

                if (!v.insideGeometry) {
                    continue;
                }

                if (testVesicleOverlap(v, v) != null) {
                    return true;
                }

            }
        }

        return false;

    }

    public DiffusantVesicle testVesicleOverlap(DiffusantVesicle v, DiffusantVesicle ignoreVesicle) {

        double dx, dy, dz, sqrDistance, minDBV;
        int ii, jj, kk;

        DiffusantVesicle ivesicle;

        Voxel voxel = v.voxel;
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
            ivesicle = ivoxel.firstReady;

            while (ivesicle != null) {

                if ((ivesicle != v) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry && (ivesicle.radius > 0)) {

                    dx = v.x - ivesicle.x;
                    dy = v.y - ivesicle.y;
                    dz = v.z - ivesicle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    minDBV = v.radius + ivesicle.radius;

                    if (sqrDistance < minDBV * minDBV) {
                        return ivesicle;
                    }

                }

                ivesicle = ivesicle.nextReady;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];
                ivesicle = ivoxel.firstReady;

                while (ivesicle != null) {

                    if ((ivesicle != v) && (ivesicle != ignoreVesicle) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry) {

                        ii = voxelPBC.PBCi[i];
                        jj = voxelPBC.PBCj[i];
                        kk = voxelPBC.PBCk[i];

                        dx = v.x - (ii * 2 * geometry.x2 + ivesicle.x);
                        dy = v.y - (jj * 2 * geometry.y2 + ivesicle.y);
                        dz = v.z - (kk * 2 * geometry.z2 + ivesicle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        minDBV = v.radius + ivesicle.radius;

                        if (sqrDistance < minDBV * minDBV) {
                            return ivesicle;
                        }

                    }

                    ivesicle = ivesicle.nextReady;

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

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (v == null) {
                    continue;
                }

                minSD = minSqrDistance(v, v);

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
                s.finish("Vesicles", project.geometry, -1);
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
            MSDspatial[i].ydim = "Vesicles" + " (" + project.spaceUnits + "^2)";
            MSDspatial[i].updateVectors();
            MSDspatial[i].init("Vesicles", project.geometry, -1, dataPoints);
        }

    }

    public void initD20() {
        double dx, dy, dz;

        if (diffusants == null) {
            return;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if ((v == null) || !v.mobile || !v.insideGeometry) {
                    continue;
                }

                v.x0 = v.x;
                v.y0 = v.y;
                v.z0 = v.z;

                dx = v.x - 0;
                dy = v.y - 0;
                dz = v.z - 0;

                v.d20 = Math.sqrt(dx * dx + dy * dy + dz * dz);

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

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if ((v == null) || !v.mobile || !v.insideGeometry) {
                    continue;
                }

                if ((v.d20 >= x1) && (v.d20 < x2)) {
                    sumSD += v.squareDisplacement();
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

        if (!connectVesicles) {
            return false;
        }

        if (connectorsBinomial) {
            return initConnectorBinomialDistribution();
        }

        Master.log("init " + maxNumConnectors + " connectors / vesicle");

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (!(v instanceof DiffusantVesicle)) {
                    continue;
                }

                v.connectTo = new DiffusantVesicle[maxNumConnectors];
                v.connectorOffTime = new double[maxNumConnectors];

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

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (!(v instanceof DiffusantVesicle)) {
                    continue;
                }

                nConnectors = 0;

                for (int k = 0; k < maxNumConnectors; k++) {
                    if (mt.nextDouble() < probability) {
                        nConnectors++;
                    }
                }

                if (nConnectors > 0) {
                    v.connectTo = new DiffusantVesicle[nConnectors]; // need twice as many for tracking connectors
                    v.connectorOffTime = new double[nConnectors]; // need twice as many for tracking connectors
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

        Master.log("average = " + avg + " connectors / vesicle");
        Master.log("fraction connected = " + (numConnected / count));

        return false;

    }

    public void connectVesicles() {

        int ii, jj, kk;
        double dx, dy, dz, sqrDistance, minDBV;

        DiffusantVesicle kvesicle;

        VoxelPBC voxelPBC;

        for (int i = 0; i < diffusants.length; i++) {

            if ((diffusants[i] == null) || (diffusants[i].vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : diffusants[i].vesicles) {

                for (int k = 0; k < v.voxel.numNeighbors; k++) {

                    kvesicle = (DiffusantVesicle) v.voxel.neighbors[k].firstReady;

                    while (kvesicle != null) {

                        if (v.noMoreConnectors()) {
                            break;
                        }

                        if ((kvesicle != v) && !v.isConnectedTo(kvesicle) && !kvesicle.noMoreConnectors()) {

                            if (v.overlap(kvesicle, connectorLength)) {

                                if (mt.nextDouble() < connectRate * project.dt) {

                                    if (kvesicle.connectToNew(v, mt, unconnectRate, time)) {

                                        avgConnectorLifeTime += v.connectorLifeTime;
                                        connectorLifeTimeCounter += 1;

                                        //Master.log("connected vesicles " + dv.connectorLifeTime);
                                        if (!v.connectTo(kvesicle)) {
                                            Master.exit("connection failure");
                                        }

                                    } else {
                                        Master.exit("connection failure");
                                    }

                                }

                            }

                        }

                        kvesicle = kvesicle.nextReady;

                    }

                }

                if (v.noMoreConnectors()) {
                    continue;
                }

                if (!PBC) {
                    continue;
                }

                if (v.voxel instanceof VoxelPBC) {
                    voxelPBC = (VoxelPBC) v.voxel;
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

                    kvesicle = (DiffusantVesicle) voxelPBC.PBCneighbors[k].firstReady;

                    while (kvesicle != null) {

                        if (v.noMoreConnectors()) {
                            break;
                        }

                        if ((kvesicle != v) && !v.isConnectedTo(kvesicle) && !kvesicle.noMoreConnectors()) {

                            ii = voxelPBC.PBCi[i];
                            jj = voxelPBC.PBCj[i];
                            kk = voxelPBC.PBCk[i];

                            dx = v.x - (ii * 2 * geometry.x2 + kvesicle.x);
                            dy = v.y - (jj * 2 * geometry.y2 + kvesicle.y);
                            dz = v.z - (kk * 2 * geometry.z2 + kvesicle.z);

                            sqrDistance = dx * dx + dy * dy + dz * dz;

                            minDBV = v.radius + kvesicle.radius + connectorLength;

                            if (sqrDistance < minDBV * minDBV) {

                                if (mt.nextDouble() < connectRate * project.dt) {

                                    if (kvesicle.connectToNew(v, mt, unconnectRate, time)) {

                                        avgConnectorLifeTime += v.connectorLifeTime;
                                        connectorLifeTimeCounter += 1;

                                        //Master.log("connected vesicles PBC " + dv.connectorLifeTime);
                                        if (!v.connectTo(kvesicle)) {
                                            Master.exit("connection failure");
                                        }

                                    } else {
                                        Master.exit("connection failure");
                                    }

                                }

                            }

                        }

                        kvesicle = kvesicle.nextReady;

                    }

                }

            }

        }

    }

    public void unconnectVesicles() {

        DiffusantVesicle v2;

        if (diffusants == null) {
            return;
        }

        if (unconnectRate <= 0) {
            return;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v1 : d.vesicles) {

                if (v1.connectTo == null) {
                    continue;
                }

                for (int k = 0; k < v1.connectTo.length; k++) {

                    if (v1.connectTo[k] == null) {
                        continue;
                    }

                    if ((connectorAllOffTime > 0) && (time >= connectorAllOffTime)) {
                        v2 = v1.connectTo[k];
                        v2.unconnectFrom(v1);
                        v1.connectTo[k] = null;
                        v1.connectorOffTime[k] = 0;
                    }

                    if ((v1.connectorOffTime[k] > 0) && (time >= v1.connectorOffTime[k])) {
                        //Master.log("unconnected vesicles at " + time + " ms, offtime = " + (dv.connectorAllOffTime[1][k] * project.dt));
                        //Master.log("unconnected vesicles at " + time + " ms");
                        v2 = v1.connectTo[k];
                        v2.unconnectFrom(v1);
                        v1.connectTo[k] = null;
                        v1.connectorOffTime[k] = 0;
                    }

                }

            }

        }

        if ((connectorAllOffTime > 0) && (time >= connectorAllOffTime)) {
            connectVesicles = false;
        }

    }

    public double[] computeNumberOfConnectors() {

        double avg = 0, connected = 0;
        int n, ntotal = 0;

        if (diffusants == null) {
            return null;
        }

        double[] results = new double[2];

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                n = v.numberOfConnections(true);
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

        //Master.log("" + avg + " connectors / vesicle");
        //Master.log("" + connected + " fraction connected");
        results[0] = avg;
        results[1] = connected;

        return results;

    }

    private void testConnectorSeperation() {

        if (diffusants == null) {
            return;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (!(v instanceof DiffusantVesicle)) {
                    continue;
                }

                if (v.testConnectorSeperation(connectorLength)) {
                    Master.log("connectors seperated " + v);
                }

            }
        }

    }

    public double avgFirstCollision() {

        double avg = 0, count = 0, sqrd = 0;

        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if ((v == null) || !v.mobile) {
                    continue;
                }

                if (v.firstCollision > 0) {
                    avg += v.firstCollision;
                    sqrd += v.sqrDisplacement;
                    count++;
                    //Master.log("" + vs.vesicle[j].firstCollision);
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

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (v == null) {
                    continue;
                }

                if (v.mobile && (v.residenceTime >= 0)) {
                    v.residenceTime += project.dt;
                }

            }

        }

    }

    public DiffusantVesicle findVesicleOverlap(DiffusantVesicle v) {

        double dx, dy, dz, sqrDistance, minDBV;
        int ii, jj, kk;

        DiffusantVesicle ivesicle;

        Voxel voxel = v.voxel;
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
            ivesicle = ivoxel.firstReady;

            while (ivesicle != null) {

                if ((ivesicle != v) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry) {

                    dx = v.x - ivesicle.x;
                    dy = v.y - ivesicle.y;
                    dz = v.z - ivesicle.z;

                    sqrDistance = dx * dx + dy * dy + dz * dz;

                    minDBV = v.radius + ivesicle.radius;

                    if (sqrDistance < minDBV * minDBV) {
                        return ivesicle;
                    }

                }

                ivesicle = ivesicle.nextReady;

            }

        }

        if (PBC && (voxelPBC != null)) {

            for (int i = 0; i < voxelPBC.numPBCneighbors; i++) {

                ivoxel = voxelPBC.PBCneighbors[i];

                ii = voxelPBC.PBCi[i];
                jj = voxelPBC.PBCj[i];
                kk = voxelPBC.PBCk[i];

                ivesicle = ivoxel.firstReady;

                while (ivesicle != null) {

                    if ((ivesicle != v) && (!Double.isNaN(ivesicle.x)) && ivesicle.insideGeometry) {

                        dx = v.x - (ii * 2 * geometry.x2 + ivesicle.x);
                        dy = v.y - (jj * 2 * geometry.y2 + ivesicle.y);
                        dz = v.z - (kk * 2 * geometry.z2 + ivesicle.z);

                        sqrDistance = dx * dx + dy * dy + dz * dz;

                        minDBV = v.radius + ivesicle.radius;

                        if (sqrDistance < minDBV * minDBV) {
                            //return ivesicle;
                            return null;
                        }

                    }

                    ivesicle = ivesicle.nextReady;

                }

            }

        }

        return null;

    }

    public double localDensityFastHydroWallz(DiffusantVesicle v) {

        double h, volumeCap, volumeTotal, localVolumeFraction, localD;
        double nRadiiExt = 4; // Rext

        double dz = geometry.z2 - v.z;
        double radius4 = v.radius * nRadiiExt;

        if (dz >= radius4) {
            return v.step3;
        }

        h = radius4 - dz;
        volumeCap = (Math.PI * h * h / 3) * (3 * radius4 - h);
        volumeTotal = 4 * Math.PI * radius4 * radius4 * radius4 / 3;

        //if (dz < dv.radius) {
        //    Master.log("" + dz);
        //}
        localVolumeFraction = (volumeTotal - volumeCap) * vesicleVolumeFraction / volumeTotal;

        //localD = Dcyto * DiffusantVesicle.Dratio_shortOLD(localVolumeFraction, immobileVesiclePercent);
        return Double.NaN; // DiffusantVesicle.step3(localD, project.dt);

        //Master.log("" + vesicleVolumeFraction);
        //return (volumeTotal - volumeCap)/volumeTotal;
    }

    public double localDensityAll(boolean print) {

        double density, avg = 0, avgD = 0, avgDD0 = 0, avgStep = 0, count = 0;
        double min = 99999, max = 0;

        //DiffusantVesicle v;
        //if (immobileVesicleFraction > 0) {
        //    Master.exit("aborted MC simulation: immobile fraction not allowed with local density");
        //}
        if (diffusants == null) {
            return Double.NaN;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (v == null) {
                    continue;
                }

                if (v.insideGeometry) {

                    if (PBC) {
                        density = localDensityPBC(v);
                    } else if (hydrodynamicsLocalDVoxels) {
                        density = localDensityVoxels(v);
                    } else {
                        density = localDensity(v);
                    }

                    if (print && (density > 0)) {
                        avg += density;
                        avgStep += v.localStep3;
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
            Master.log("local D avg = " + DiffusantVesicle.D3(avgStep, project.dt));
            Master.log("local step avg = " + avgStep);
        }

        return avg;

    }

    public double localDensity(DiffusantVesicle v) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt, nRadiiInt;
        int x0, y0, z0, x1, y1, z1;
        double dx, dy, dz;
        double d, h, Rext, r, totalV, sumV, volumeCap = 0;
        double localMobileVolumeFraction, localD, Ds;

        Voxel voxel;
        DiffusantVesicle ivesicle;

        if ((v == null) || (geometry.voxelSpace == null)) {
            return Double.NaN;
        }

        //nRadiiInt = 2; // Rint
        //Rint = nRadiiInt * dv.radius;
        dz = geometry.z2 - v.z;

        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * v.radius;
        h = Rext - dz;

        if (h > 0) {
            volumeCap = (Math.PI * h * h / 3) * (3 * Rext - h);
        }

        totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        //sumV = 4.0 * Math.PI * dv.radius * dv.radius * dv.radius / 3.0;
        sumV = vesicleVolume; // all vesicles with same radius

        xVoxel = (int) geometry.computeVoxelX(v.x);
        yVoxel = (int) geometry.computeVoxelY(v.y);
        zVoxel = (int) geometry.computeVoxelZ(v.z);

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

        v.localStep3 = 0;

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    voxel = geometry.voxelSpace[i][j][k];
                    ivesicle = voxel.firstReady;

                    while (ivesicle != null) {

                        if ((ivesicle != v) && ivesicle.insideGeometry && ivesicle.mobile) {

                            dx = v.x - ivesicle.x;
                            dy = v.y - ivesicle.y;
                            dz = v.z - ivesicle.z;

                            d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            r = ivesicle.radius;

                            if (d < Rext - r) {
                                sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            } else if (d < Rext + r) {
                                sumV += Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                            }

                        }

                        ivesicle = ivesicle.nextReady;

                    }

                }
            }
        }

        localMobileVolumeFraction = sumV / (totalV - volumeCap);
        //Ds = DiffusantVesicle.Dratio_shortOLD(localMobileVolumeFraction, immobileVesiclePercent);
        Ds = DiffusantVesicle.Dratio_shortNew(localMobileVolumeFraction, immobileVolumeFraction);
        //localD = dv.D * Ds;
        localD = Dcyto * Ds;
        v.localStep3 = DiffusantVesicle.step3(localD, project.dt);
        v.DsDff = Ds / DiffusantVesicle.Dff_short_Banchio(localMobileVolumeFraction);

        return localMobileVolumeFraction;

    }

    public double localDensityVoxels(DiffusantVesicle v) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt, nRadiiInt;
        int x0, y0, z0, x1, y1, z1;
        double dx, dy, dz;
        double d, h, Rext, r, totalV = 0, sumV_mobile, sumV_immobile = 0;
        double localMobileVolumeFraction, localImmobileVolumeFraction, localD, Ds;

        Voxel voxel;
        DiffusantVesicle ivesicle;

        if ((v == null) || (geometry.voxelSpace == null)) {
            return Double.NaN;
        }

        //nRadiiInt = 2; // Rint
        //Rint = nRadiiInt * dv.radius;
        //dz = geometry.z2 - dv.z;
        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * v.radius;
        //h = Rext - dz;

        //if (h > 0) {
        //    volumeCap = (Math.PI * h * h / 3) * (3 * Rext - h);
        //}
        //totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        //sumV = 4.0 * Math.PI * dv.radius * dv.radius * dv.radius / 3.0;
        sumV_mobile = vesicleVolume; // all vesicles with same radius

        xVoxel = (int) geometry.computeVoxelX(v.x);
        yVoxel = (int) geometry.computeVoxelY(v.y);
        zVoxel = (int) geometry.computeVoxelZ(v.z);

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

        v.localStep3 = 0;

        for (int i = x0; i <= x1; i++) {
            for (int j = y0; j <= y1; j++) {
                for (int k = z0; k <= z1; k++) {

                    voxel = geometry.voxelSpace[i][j][k];

                    dx = v.x - voxel.x;
                    dy = v.y - voxel.y;
                    dz = v.z - voxel.z;

                    d = Math.sqrt(dx * dx + dy * dy + dz * dz);

                    if (d > Rext) {
                        continue;
                    }

                    if (voxel.isSpace) {
                        totalV += voxelVolume;
                    }

                    ivesicle = voxel.firstReady;

                    while (ivesicle != null) {

                        if ((ivesicle != v) && ivesicle.insideGeometry) {

                            //dx = dv.x - ivesicle.x;
                            //dy = dv.y - ivesicle.y;
                            //dz = dv.z - ivesicle.z;
                            //d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            //r = ivesicle.radius;
                            //sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            if (ivesicle.mobile) {
                                sumV_mobile += vesicleVolume;
                            } else {
                                sumV_immobile += vesicleVolume;
                            }

                            //if (d < Rext - r) {
                            //    sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            //} else if (d < Rext + r) {
                            //    sumV += Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                            //}
                        }

                        ivesicle = ivesicle.nextReady;

                    }

                }
            }
        }

        //localDensity =  sumV / (totalV - volumeCap);
        localMobileVolumeFraction = sumV_mobile / totalV;
        localImmobileVolumeFraction = sumV_immobile / totalV;

        //Ds = DiffusantVesicle.Dratio_shortOLD(localMobileVolumeFraction, immobileVesiclePercent);
        Ds = DiffusantVesicle.Dratio_shortNew(localMobileVolumeFraction, immobileVolumeFraction);
        //Ds = DiffusantVesicle.Dratio_shortNew(localMobileVolumeFraction, localImmobileVolumeFraction);
        //localD = dv.D * Ds;
        localD = Dcyto * Ds;
        //Master.log("localD = " + localD);
        v.localStep3 = DiffusantVesicle.step3(localD, project.dt);
        v.DsDff = Ds / DiffusantVesicle.Dff_short_Banchio(localMobileVolumeFraction);

        if (Double.isNaN(v.localStep3)) {
            v.localStep3 = v.step3;
            Master.log("NaN localSteps: " + localMobileVolumeFraction + "," + v.localStep3);
        }

        return localMobileVolumeFraction;

    }

    public double localDensityPBC(DiffusantVesicle v) {

        int xVoxel, yVoxel, zVoxel, nVoxels, nRadiiExt, nRadiiInt;
        int x0, y0, z0, x1, y1, z1, ii, jj, kk;
        double dx, dy, dz;
        double d, Rext, r, totalV, sumV;
        double xdir, ydir, zdir;
        double localMobileVolumeFraction, localD, Ds;

        Voxel voxel;
        DiffusantVesicle ivesicle;

        if ((v == null) || (geometry.voxelSpacePBC == null)) {
            return Double.NaN;
        }

        //nRadiiInt = 2; // Rint
        //Rint = nRadiiInt * dv.radius;
        nRadiiExt = 4; // Rext
        Rext = nRadiiExt * v.radius;
        totalV = 4.0 * Math.PI * Rext * Rext * Rext / 3.0;
        //sumV = 4.0 * Math.PI * dv.radius * dv.radius * dv.radius / 3.0;
        sumV = vesicleVolume; // all vesicles with same radius

        xVoxel = (int) geometry.computeVoxelX(v.x);
        yVoxel = (int) geometry.computeVoxelY(v.y);
        zVoxel = (int) geometry.computeVoxelZ(v.z);

        nVoxels = (int) Math.ceil(Rext / project.dx);

        x0 = xVoxel - nVoxels;
        y0 = yVoxel - nVoxels;
        z0 = zVoxel - nVoxels;

        x1 = xVoxel + nVoxels;
        y1 = yVoxel + nVoxels;
        z1 = zVoxel + nVoxels;

        v.localStep3 = 0;

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
                    ivesicle = voxel.firstReady;

                    while (ivesicle != null) {

                        if ((ivesicle != v) && ivesicle.insideGeometry && ivesicle.mobile) {

                            dx = v.x - (xdir * 2 * geometry.x2 + ivesicle.x);
                            dy = v.y - (ydir * 2 * geometry.y2 + ivesicle.y);
                            dz = v.z - (zdir * 2 * geometry.z2 + ivesicle.z);

                            d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                            r = ivesicle.radius;

                            if (d < Rext - r) {
                                sumV += 4.0 * Math.PI * r * r * r / 3.0;
                            } else if (d < Rext + r) {
                                sumV += Math.PI * (Rext + r - d) * (Rext + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * Rext + 6 * r * Rext - 3 * Rext * Rext) / (12 * d);
                            }

                        }

                        ivesicle = ivesicle.nextReady;

                    }

                }
            }
        }

        localMobileVolumeFraction = sumV / totalV;
        Ds = DiffusantVesicle.Dratio_shortNew(localMobileVolumeFraction, immobileVolumeFraction);
        //localD = dv.D * Ds;
        localD = Dcyto * Ds;
        v.localStep3 = DiffusantVesicle.step3(localD, project.dt);
        v.DsDff = Ds / DiffusantVesicle.Dff_short_Banchio(localMobileVolumeFraction);

        return localMobileVolumeFraction;

    }

    public void initClusterImmobileCounter() {

        clusterImmobileCounter = Math.round(clusterImmobilePercent * numVesicles) - immobileVesicles;

        if (clusterImmobileCounter > 0) {
            clusterImmobileOn = true;
        }

    }

    public void clusterImmobileVesicles() {

        double dx, dy, dz, sqrDistance, minDBV;
        double tetherDistance = 0.00005;

        DiffusantVesicle ivesicle;

        Voxel voxel, ivoxel;

        if (diffusants == null) {
            return;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            d.shuffleVesicles();

            for (DiffusantVesicle v : d.vesicles) {

                if (!v.insideGeometry) {
                    continue;
                }

                if (v.mobile) {
                    continue;
                }

                voxel = v.voxel;

                if (voxel == null) {
                    return;
                }

                for (int k = 0; k < voxel.numNeighbors; k++) {

                    ivoxel = voxel.neighbors[k];
                    ivesicle = ivoxel.firstReady;

                    while (ivesicle != null) {

                        if ((ivesicle != v) && (ivesicle.mobile)) {

                            dx = v.x - ivesicle.x;
                            dy = v.y - ivesicle.y;
                            dz = v.z - ivesicle.z;

                            sqrDistance = dx * dx + dy * dy + dz * dz;
                            minDBV = v.radius + ivesicle.radius + tetherDistance;

                            if (sqrDistance < minDBV * minDBV) {
                                ivesicle.mobile = false;
                                clusterImmobileCounter--;
                                if (clusterImmobileCounter <= 0) {
                                    clusterImmobileOn = false;
                                    return;
                                }
                            }

                        }

                        ivesicle = ivesicle.nextReady;

                    }

                }

            }

        }

    }

    public double zMSD() {

        double sumMSD = 0, count = 0;

        DiffusantVesicle v;

        for (int i = 0; i <= geometry.xVoxels - 1; i++) {
            for (int j = 0; j <= geometry.yVoxels - 1; j++) {
                for (int k = geometry.zVoxels - 2; k <= geometry.zVoxels - 1; k++) {

                    if (geometry.voxelSpacePBC != null) {
                        v = geometry.voxelSpacePBC[i][j][k].firstReady;
                    } else if (geometry.voxelSpace != null) {
                        v = geometry.voxelSpace[i][j][k].firstReady;
                    } else {
                        return Double.NaN;
                    }

                    while (v != null) {
                        sumMSD += v.squareDisplacement();
                        count++;
                        v = v.nextReady;
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

        DiffusantVesicle ivesicle;

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

                    ivesicle = geometry.voxelSpace[i][j][k].firstReady;

                    while (ivesicle != null) {
                        if (ivesicle.mobile) {
                            t += ivesicle.residenceTime;
                            n++;
                            ivesicle = ivesicle.nextReady;
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

    public boolean addToVoxelList(DiffusantVesicle v) {

        if (v.voxel == null) {
            return false;
        }

        v.nextReady = v.voxel.firstReady; // save existing vesicles
        v.voxel.firstReady = v; // replace with new vesicles, creating chain

        return true;

    }

    public boolean removeFromVoxelList(DiffusantVesicle v, Voxel voxel) {

        if (voxel == null) {
            return false;
        }

        boolean found = false;

        DiffusantVesicle vtest = voxel.firstReady;
        DiffusantVesicle vnext;
        DiffusantVesicle vhold = null;

        while (vtest != null) {

            vnext = vtest.nextReady;

            if (vtest == v) {
                if (vhold == null) {
                    voxel.firstReady = vnext; // dv was first in list
                } else {
                    vhold.nextReady = vnext;
                }
                found = true;
            }

            vhold = vtest;
            vtest = vnext;

        }

        if (found) {
            if (v.residenceTime > 0) {
                voxel.residenceCounter++;
                voxel.sumResidenceTime += v.residenceTime;
            }
            v.residenceTime = 0;
        }

        return found;

    }

    double ranGauss() { // 0.0 mean, 1.0 stdv

        double v1, v2, w, step;

        if (saveRanGauss != 999999) {
            step = saveRanGauss;
            saveRanGauss = 999999;
            return step;
        }

        do {
            v1 = (2.0 * mt.nextDouble() - 1.0);
            v2 = (2.0 * mt.nextDouble() - 1.0);
            w = v1 * v1 + v2 * v2;
        } while (w >= 1.0);

        w = Math.sqrt((-2.0 * Math.log(w)) / w);

        step = v1 * w; // first value
        saveRanGauss = v2 * w; // save second value

        return step;

    }

    public boolean moveVesicleGauss(DiffusantVesicle v) {

        double step3;

        if (removingVesicleOverlap) {
            step3 = removeVesicleOverlapStep3;
        } else if (hydrodynamicsLocalD) {
            step3 = v.localStep3;
        } else {
            step3 = v.step3;
        }

        //if (step <= 0) {
        //    Master.exit("moveVesicle Error: bad step: " + step);
        //}
        stepx = ranGauss() * step3;
        stepy = ranGauss() * step3;
        stepz = ranGauss() * step3;

        return setVesicleLocation(v, v.x + stepx, v.y + stepy, v.z + stepz, false);

    }

    public boolean moveVesicleGaussHydroWallz(DiffusantVesicle v) {

        double step3;
        double z, b_ll = 1, b_T = 1;

        if (hydrodynamicsLocalD) {
            step3 = v.localStep3;
        } else if (hydrodynamicsDscale) {
            //step3 = localDensityFastHydroWallz(dv); // using Michailidou model instead
            step3 = v.step3;
        } else {
            step3 = v.step3;
        }

        if (removingVesicleOverlap) {

            step3 = removeVesicleOverlapStep3;

        } else {

            z = geometry.z2 - v.z;

            //b_ll = Math.sqrt(hydroWall_ll(dv.radius, z));
            //b_T = Math.sqrt(hydroWall_T(dv.radius, z - dv.radius));
            b_ll = hydroWall_ll(v.radius, z);
            b_T = hydroWall_T(v.radius, Math.abs(z - v.radius));

            b_ll = 1 / (1 + v.DsDff * ((1 / b_ll) - 1)); // Michailidou et al. 2009
            b_T = 1 / (1 + v.DsDff * ((1 / b_T) - 1)); // Michailidou et al. 2009

            b_ll = Math.sqrt(b_ll);
            b_T = Math.sqrt(b_T);

        }

        stepx = ranGauss() * step3 * b_ll;
        stepy = ranGauss() * step3 * b_ll;
        stepz = ranGauss() * step3 * b_T;

        //stepx = step3 * b_ll;
        //stepy = step3 * b_ll;
        //stepz = step3 * b_T;
        //dv.localD = (stepx * stepx + stepy * stepy + stepz * stepz) / (6 * project.dt);
        //Master.log("" + (dv.localD/1.397657228105222E-4));
        //dv.localDw = (3 * step3 * step3) / (6 * project.dt);
        return setVesicleLocation(v, v.x + stepx, v.y + stepy, v.z + stepz, false);

    }

    public boolean moveVesicle(DiffusantVesicle v) {

        double step3;

        if (hydrodynamicsLocalD) {
            step3 = v.localStep3;
        } else {
            step3 = v.step3;
        }

        //if (step < 0) {
        //    Master.exit("moveVesicle Error: bad step: " + step);
        //}
        if (removingVesicleOverlap) {
            step3 = removeVesicleOverlapStep3;
        }

        if (mt.nextDouble() < 0.5) {
            stepx = step3;
        } else {
            stepx = -step3;
        }

        if (mt.nextDouble() < 0.5) {
            stepy = step3;
        } else {
            stepy = -step3;
        }

        if (mt.nextDouble() < 0.5) {
            stepz = step3;
        } else {
            stepz = -step3;
        }

        return setVesicleLocation(v, v.x + stepx, v.y + stepy, v.z + stepz, false);

    }

    public void moveVesiclesCichocki(boolean moveAll) {

        boolean outOfBounds, overlap, sameVoxel, moveConnected, ok;

        DiffusantVesicle odv;

        if (diffusants == null) {
            return;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            d.shuffleVesicles();

            for (DiffusantVesicle v : d.vesicles) {

                if (!v.insideGeometry) {
                    continue;
                }

                if (!moveAll && !v.mobile) {
                    continue;
                }

                testVesicle.copy(v);

                if (hydroWallZ) {
                    if (!moveVesicleGaussHydroWallz(testVesicle)) {
                        continue; // dont move
                    }
                } else {
                    if (!moveVesicleGauss(testVesicle)) {
                        continue; // dont move
                    }
                }

                outOfBounds = outOfBounds(testVesicle);

                if (PBC && outOfBounds) {

                    outOfBounds = wrapAtBorder(testVesicle);

                    //if (outOfBounds) {
                    //
                    //Master.log("moveVesicles error: PBC wrap at border error");
                    //return -1;
                    //} else {
                    //Master.log("wrapped vesicles");
                    //}
                }

                if (outOfBounds) {
                    continue; // dont move
                }

                if (testOverlap(testVesicle)) {
                    if ((time > 0) && (v.firstCollision == 0)) {
                        v.firstCollision = time;
                        v.sqrDisplacement = v.squareDisplacement();
                    }
                    continue; // dont move
                }

                if (!freeDiffusion) {

                    odv = testVesicleOverlap(testVesicle, v);

                    overlap = odv != null;

                    if (overlap) {
                        //Master.log("square displacement = " + diffusant[k].vesicle[j].squareDisplacement());
                        //Master.log("time = " + time);
                        //cancel = true;
                        //if (!removingVesicleOverlap && diffusant[k].saveDisplacementAfterCollisions && (diffusant[k].vesicle[j].lastVesicleCollision != odv)) {
                        if (!removingVesicleOverlap && d.saveDisplacementAfterCollisions) {
                            v.sqrDisplacement = v.squareDisplacement();
                            v.lastVesicleCollision = odv;
                            v.x0 = v.x;
                            v.y0 = v.y;
                            v.z0 = v.z;
                        }
                        if ((time > 0) && (v.firstCollision == 0)) {
                            v.firstCollision = time;
                            v.sqrDisplacement = v.squareDisplacement();
                        }
                        continue;
                    }

                }

                if (connectVesicles && (v.connectTo != null)) {

                    moveConnected = true;

                    for (DiffusantVesicle v2 : v.connectTo) {

                        if (v2 == null) {
                            continue;
                        }

                        if (!testVesicle.overlap(v2, connectorLength)) {
                            moveConnected = false;
                            break;
                        }

                    }

                    if (!moveConnected) {
                        continue;
                    }

                }

                sameVoxel = false;

                if (testVesicle.voxel == v.voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(v, v.voxel);
                }

                v.copy(testVesicle);

                if (!sameVoxel) {
                    addToVoxelList(v);
                }

            }

        }

    }

    public void drift() {
        boolean outOfBounds, sameVoxel;

        if (diffusants == null) {
            return;
        }

        if (!PBC) {
            Master.exit("drift error: PBC is not on");
        }
        
        driftGeometry(driftx, drifty, driftz);

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (!v.insideGeometry) {
                    continue;
                }

                testVesicle.copy(v);

                setVesicleLocation(testVesicle, testVesicle.x + driftx, testVesicle.y + drifty, testVesicle.z + driftz, false);

                outOfBounds = outOfBounds(testVesicle);

                if (outOfBounds) {
                    outOfBounds = wrapAtBorder(testVesicle);
                    //testVesicle.fluorescence = 1.0; // reset F to simulate non-frapped vesicles moving into PSF
                }

                if (outOfBounds) {
                    Master.exit("drift error: wrap at border error");
                    //continue;
                }

                sameVoxel = false;

                if (testVesicle.voxel == v.voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(v, v.voxel);
                }

                v.copy(testVesicle);

                if (!sameVoxel) {
                    addToVoxelList(v);
                }

            }

        }

    }
    
    public void driftGeometry(double dx, double dy, double dz) {
        // generic function for creating drift of geometric structures
    }

    public boolean wrapAtBorder(DiffusantVesicle v) {

        double x = v.x;
        double y = v.y;
        double z = v.z;

        if (x > geometry.x2) {
            x -= 2 * geometry.x2;
            v.x0 -= 2 * geometry.x2;
        } else if (x < geometry.x1) {
            x -= 2 * geometry.x1;
            v.x0 -= 2 * geometry.x1;
        }

        if (y > geometry.y2) {
            y -= 2 * geometry.y2;
            v.y0 -= 2 * geometry.y2;
        } else if (y < geometry.y1) {
            y -= 2 * geometry.y1;
            v.y0 -= 2 * geometry.y1;
        }

        if (z > geometry.z2) {
            z -= 2 * geometry.z2;
            v.z0 -= 2 * geometry.z2;
        } else if (z < geometry.z1) {
            z -= 2 * geometry.z1;
            v.z0 -= 2 * geometry.z1;
        }

        setVesicleLocation(v, x, y, z, false);

        return outOfBounds(v);

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

        if (n.equalsIgnoreCase("minVesicleStep")) {
            if (v <= 0) {
                return false;
            }
            minVesicleStep = v;
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
        if (n.equalsIgnoreCase("connectVesicles")) {
            connectVesicles = (v == 1);
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
