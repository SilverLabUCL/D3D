package ucl.silver.d3d.core;

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
public class RunMonteCarloAZ
        extends RunMonteCarlo {

    public boolean dockingOn = false;
    public boolean dockingOffAfterInit = false;
    public double dockingOnTime = 0; // ms

    public int dockedVesicles = 0; // number of docked vesicles at active zone
    public int dockedVesiclesMax = -1; // use -1 for no limit
    private DiffusantVesicleAZ firstDocked = null;

    public int reserveVesicles = 0; // number of reserve vesicles
    public int reserveVesiclesMax = -1; // use -1 for no limit
    private DiffusantVesicleAZ firstReserve = null;

    public Coordinates activeZone = null; // active zone (az) coordinates

    public double[][][] azXYZ = null;

    public double azWidth = 0.16; // active zone xy-width
    //public double azHeight = 0.025; // active zone z-height
    public double azHeight = 0.05; // active zone z-height
    public boolean azHemisphere = false; // active zone is semi-circle

    public double azx1[] = null;
    public double azx2[] = null;
    public double azy1[] = null;
    public double azy2[] = null;
    public double azz1[] = null;
    public double azz2[] = null;
    public String azPlane = "xy";

    public double azKeepClearExtra = 0; // 0.05; // um (keeps AZ clear of immobile vesicles)
    public boolean azExcludeReadyinit = true; // keep AZ clear when initiating ready vesicles

    public double releaseRate = 0; // release docked, move to reserve (kHz)
    public double releaseStartTime = 0; // ms, min time to allow release
    public int releaseCounter = 0;
    public int releaseLimit = 0; // use value > 0 to end simulation
    public double releaseTimeOutLimit = 0; // abort simulation if no releases occur before this time
    public double releaseProb = 1.0; // release probability

    public double replenishRate = 0; // move reserve to ready (kHz)
    public double replenishStartTime = 0; // ms, min time to allow replacement
    private int replenishCounter = 0;

    public Coordinates replenishCoordinates = null; // where reserve vesicles are placed
    public Coordinates replenishCoordinatesExclude = null;
    public double replenishFromAZdistance = -1; // number of voxels from AZ for setting replenish coordinates

    public boolean tetherVesicles = false;
    public double azPermitConnectorRadius = 0.2; // um
    public int azNumTethers = 1;
    public double azTetherLength = 0.008; // um
    public Coordinates activeZoneTethering = null;
    public DiffusantParticle azDV = null;

    public boolean saveRelease = false;
    public boolean saveNumDocked = false;
    public boolean saveConnectorStats = false;
    public boolean saveVesicleDensityAZ = false;
    public boolean saveDNearAZ = false;

    public PulseTimer saveXYZincrement = null;

    public int[] azVesicleDensityZoltan = {6356, 3596, 2233, 1468, 938}; // cumulative density from AZ

    private long iRelease, iReleaseMax;
    private long iReplenish, iReplenishMax;

    Save save_Release = null; // saving data to file and/or internal array
    Save save_NumDocked = null; // saving data to file and/or internal array
    Save save_ConnectorStats = null;
    Save save_VesicleDensityAZ = null;

    private double[] connectorStats;
    private double[] vesicleDensityAZ;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("azWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("releaseRate")) {
            return project.freqUnits;
        }
        if (name.equalsIgnoreCase("replenishRate")) {
            return project.freqUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("dockedVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("tetheredVesicles")) {
            return false;
        }
        if (name.equalsIgnoreCase("reserveVesicles")) {
            return false;
        }
        return super.canEdit(name);
    }

    public RunMonteCarloAZ(Project p) {

        super(p);

        save_Release = new Save(p);
        save_NumDocked = new Save(p);
        save_ConnectorStats = new Save(p);
        save_VesicleDensityAZ = new Save(p);

        save_Release.name = "Release";
        save_NumDocked.name = "NumDockedVesicles";
        save_ConnectorStats.name = "Connector Stats";
        save_VesicleDensityAZ.name = "AZ Vesicle Density";

        activeZone = new Coordinates(project, 0, 0, 0, 0, 0, 0);
        activeZone.name = "Active Zone";
        activeZone.setShape("cylinderz");
        activeZone.color.name = "AZcolor";
        activeZone.color.setColor("[r=204,g=0,b=0]");

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        if (saveRelease) {
            save_Release.init();
            save_Release.save2TextFile = true;
            save_Release.fileName("Release", "");
            save_Release.ydim = project.timeUnits;
            save_Release.saveWhileComputing = true;
            save_Release.saveWhileComputingAppend = true;
        }

        if (saveNumDocked) {
            save_NumDocked.init();
            save_NumDocked.save2TextFile = true;
            save_NumDocked.fileName("NumDocked", "");
        }

        if (saveConnectorStats) {
            save_ConnectorStats.init();
            save_ConnectorStats.save2TextFile = false;
            save_ConnectorStats.save2BinaryFile = true;
            save_ConnectorStats.saveWhileComputing = true;
            save_ConnectorStats.fileName("ConnectorStats", "");
        }

        if (saveVesicleDensityAZ) {
            save_VesicleDensityAZ.init();
            save_VesicleDensityAZ.save2TextFile = false;
            save_VesicleDensityAZ.save2BinaryFile = true;
            save_VesicleDensityAZ.saveWhileComputing = true;
            save_VesicleDensityAZ.fileName("VesicleDensityAZ", "");
        }

        if (PBC) {
            //Master.exit("RunMonteCarloAZ: init: PBC not allowed in active zone simulations");
        }

    }

    @Override
    public void setOutputRate(double newRate) {

        super.setOutputRate(newRate);

        if (save_Release != null) {
            save_Release.setOutputRate(newRate);
        }

        if (save_NumDocked != null) {
            save_NumDocked.setOutputRate(newRate);
        }

        if (save_ConnectorStats != null) {
            save_ConnectorStats.setOutputRate(newRate);
        }

        if (save_VesicleDensityAZ != null) {
            save_VesicleDensityAZ.setOutputRate(newRate);
        }

    }

    @Override
    public boolean initSimulation() {

        if (checkVariables()) {
            return true;
        }

        if (autoInit) {
            if (initAll()) {
                return true;
            }
        }

        iRelease = 0;
        iReplenish = 0;

        if (releaseRate > 0) {
            iReleaseMax = (int) ((1.0 / releaseRate) / project.dt);
            iReleaseMax = Math.max(iReleaseMax, 1); // at least one time rstep
            //Master.log("iReleaseMax: " + iReleaseMax);
            iRelease = iReleaseMax; // start at time = 0
        }

        if (replenishRate > 0) {
            iReplenishMax = (int) (1.0 / (replenishRate * project.dt));
            iReplenishMax = Math.max(iReplenishMax, 1); // at least one time rstep
            replenishCounter = 0;
        }

        time = 0;

        timer1.start();

        //saveXYZincrement.impulse = true;
        //saveXYZincrement.name = "saveXYZincrement";
        //saveXYZincrement.initTimer();

        initialized = true;

        printParameters();

        Master.log("initialized Monte Carlo simulation");

        return false;

    }

    @Override
    public void finishSimulation(){

        finishSave();

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        particleStats();

        if (connectParticles) {
            avgConnectorLifeTime /= connectorLifeTimeCounter;
            Master.log("average connector life time (ms): " + avgConnectorLifeTime + " (n=" + Integer.toString(connectorLifeTimeCounter) + ")");
        }

        Master.log("final docked vesicles: " + dockedVesicles);
        Master.log("final reserve vesicles: " + reserveVesicles);
        Master.log("final vesicle volume fraction: " + particleVolumeFraction);

        if (replenishRate > 0) {
            Master.log("replenished vesicles: " + replenishCounter);
        }

        //computeVesicleDensityAZ();

    }

    @Override
    public double spaceVolume() {
        return project.geometry.spaceVolume - activeZone.geometryVolume();
    }

    @Override
    public void run() {

        int icount, numReturned;

        //double simTime = project.simTime;

        DiffusantVesicleAZ vaz;

        while (runSimulation) { // runSimulation thru batches

            if (!initialized) {

                if (project.simulationInit(preview)) {
                    cancelSimulation(); // ERROR
                }

                saveParticlePositions(0);
                saveParticleDensity("start");

            }

            if (time < project.simTime) {
                Master.log("starting Monte Carlo simulation at time " + time + " " + project.timeUnits);
            }

            while (runSimulation && !cancel && (time < project.simTime)) {

                for (DiffusantParticles d : diffusants) {
                    d.save(); // save diffusant variables
                }

                runMSDspatial();

                //if ((saveXYZincrement != null) && (saveXYZincrement.timer != null)) {

                    //if (saveXYZincrement.timer[this.it] > 0) {
                        //saveVesiclePositions((int) time);
                    //}

                //}

                //if (project.detectors != null) {
                    //for (int d = 0; d < project.detectors.length; d++) {
                        //if (project.detectors[d] instanceof DetectorSnapshot) {
                        //    project.detectors[d].detect(this, geometry); // compute detector averages
                        //}
                    //}
                //}

                //if (!preview && saveNumDocked) {
                //    saveNumDockedData();
                //}

                //if (!preview && saveConnectorStats) {
                //    saveConnectorStats();
                //}

                if (!preview && saveVesicleDensityAZ) {
                    saveVesicleDensityAZ();
                }

                moveVesiclesActiveZone();

                if (hydrodynamicsLocalD) {
                    localDensityAll(false);
                }

                //if (hydroWallZ) {
                //    hydroWallZ2();
                //}

                //if (!freeDiffusion && checkForOverlapsDuringSimulation && testVesicleOverlap()) {
                //    Master.log("Monte Carlo Simulation has overlapping vesicles: " + minDistanceBetweenVesicles);
                //    runSimulation = false;
                //}

                if (connectParticles) {
                    connectParticles();
                    unconnectParticles();
                }

                if (saveResidenceTime) {
                    updateResidentTimes();
                }

                if ((releaseRate > 0) && (iRelease == iReleaseMax)) {

                    //if (releaseRate == Double.POSITIVE_INFINITY) {
                    //    icount = 9999;
                    //} else {
                    //    icount = dockedVesiclesMax;
                    //}

                    icount = dockedVesicles;

                    if (time < releaseStartTime) {
                        icount = 0;
                    }

                    for (int i = 0; i < icount; i++) {
                        
                        vaz = lastDockedVesicle();
                        
                        if (vaz == null) {
                            break;
                        }

                        releaseDockedVesicle(vaz, time);

                    }

                    iRelease = 0; // reset counter

                }

                if ((replenishRate > 0) && (iReplenish == iReplenishMax)) {

                    if (replenishRate == Double.POSITIVE_INFINITY) {
                        icount = 999;
                    } else {
                        icount = 1;
                    }

                    if (time < replenishStartTime) {
                        icount = 0;
                    }

                    numReturned = moveReserveToReady(icount);

                    //if (numReturned > 0) {
                    //    Master.log("added reserve vesicles: " + numReturned);
                    //}

                    replenishCounter += numReturned;

                    iReplenish = 0; // reset counter

                }

                if ((grid != null) && preview) {
                    grid.repaint();
                }

                iRelease++;
                iReplenish++;

                itime += 1;
                time += project.dt;
                timer2.timer(time);

                if ((releaseLimit > 0) && (releaseCounter >= releaseLimit)) {
                    cancel = true;
                }

                if ((releaseTimeOutLimit > 0) && (time >= releaseTimeOutLimit) && (releaseCounter == 0)) {
                    Master.log("cancelled MC AZ simulation: no release events within " + releaseTimeOutLimit + " ms");
                    cancel = true;
                }

            }

            if (cancel || (time >= project.simTime)) {
                finish();
            }

        }

    }

    @Override
    public boolean initAll(){

        releaseCounter = 0;

        if (checkVariables()) {
            return true;
        }

        if (initDiffusantParticlesArray()) {
            return true;
        }

        if (initVoxels()) {
            return true;
        }

        if (initActiveZone()) {
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

        if (initVesiclesDocked()){
            return true;
        }

        if (initParticleRandom && initParticlesRandom()) {
            return true;
        }

        //minDistanceSlow();

        if (initDockedList()) {
            return true;
        }

        if (initReserveList()) {
            return true;
        }

        if (removeParticleOverlap && removeParticleOverlap()) {
            return true;
        }

        if (initParticlesImmobile()) {
            return true;
        }

        initConnectors();

        particleStats();

        replenishCoordinates = new Coordinates(project, geometry);

        if (!freeDiffusion && (particleVolumeFraction > 0.44)) {
            replenishRate = 0; // TAKES TOO LONG
            Master.log("TURNED OFF REPLENISHING: replenishRate = 0");
        }

        //initAZ_patch();

        //Master.updatePanel2D();

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        if (dockingOffAfterInit) {
            dockingOn = false; // TEMPORARY removal of docking
        }

        return false;

    }

    public void initAZ_patch() {

        DiffusantVesicleAZ vaz;

        computeVesicleDensityAZ2(true);

        if (true) {
            return;
        }

        azHemisphere = true;
        azWidth = 0.3;
        azHeight = 0.3;
        initActiveZone();

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle v : d.particles) {

                vaz = (DiffusantVesicleAZ) v;

                if (activeZone.isInside(vaz.x, vaz.y, vaz.z, vaz.radius)) {
                    //dv.isReserve = true;
                    //dv.isReady = false;
                    //dv.isDocked = false;
                    //dv.x = Double.NaN;
                    //dv.y = Double.NaN;
                    //dv.z = Double.NaN;
                    vaz.mobile = true;
                    //Master.log("inside AZ " + diffusant[j].vesicle[j]);
                }

            }

        }

    }

    @Override
    public boolean checkVariables(){

        if (PBC) {
           //Master.exit("checkVariables: PBC not allowed with active zone");
        }

        return super.checkVariables();

    }

    public boolean initActiveZone() {

        double radius, height, SA, diameter;
        double x0, x1, y0, y1, z0, z1;

        if ((geometry == null) || (activeZone == null)) {
            return true;
        }

        radius = azWidth / 2.0;
        height = azHeight / 2.0;

        if (azPlane.equalsIgnoreCase("xy")) {
            x0 = 0 - radius;
            y0 = 0 - radius;
            z0 = geometry.z2 - height;
            x1 = 0 + radius;
            y1 = 0 + radius;
            z1 = geometry.z2 + height;
            activeZone.setShape("cylinderz");
        } else if (azPlane.equalsIgnoreCase("yz")) {
            x0 = geometry.x2 - height;
            y0 = 0 - radius;
            z0 = 0 - radius;
            x1 = geometry.x2 + height;
            y1 = 0 + radius;
            z1 = 0 + radius;
            activeZone.setShape("cylinderx");
        } else if (azPlane.equalsIgnoreCase("zx")) {
            x0 = 0 - radius;
            y0 = geometry.y2 - height;
            z0 = 0 - radius;
            x1 = 0 + radius;
            y1 = geometry.y2 + height;
            z1 = 0 + radius;
            activeZone.setShape("cylindery");
        } else {
            return true;
        }

        activeZone.setCoordinates(x0, y0, z0, x1, y1, z1);

        if (azHemisphere) {
            azHeight = azWidth;
            activeZone.setShape("ellipsoid");
        }

        

        if (tetherVesicles) {

            radius += azTetherLength;
            height += azTetherLength;

            activeZoneTethering = new Coordinates(project, 0, 0, 0, 0, 0, 0);
            activeZoneTethering.setCoordinates(x0, y0, z0, x1, y1, z1);
            activeZoneTethering.setShape(activeZone.shape);

            if (azNumTethers > 0) {
                azDV = new DiffusantParticle(project, "AZ", 0, 0, 0, 0, 0);
                azDV.connectTo = new DiffusantParticle[azNumTethers];
                azDV.connectorOffTime = new double[azNumTethers];
            }

        }

        return false;

    }

    @Override
    public boolean initParticlesImmobile() {

        double xx1, xx2, yy1, yy2, zz1, zz2;

        Coordinates keepClear = null;

        if (diffusants == null) {
            return true;
        }

        if (azKeepClearExtra >= 0) {

            xx1 = activeZone.x1 - 1 * azKeepClearExtra;
            xx2 = activeZone.x2 + 1 * azKeepClearExtra;
            yy1 = activeZone.y1 - 1 * azKeepClearExtra;
            yy2 = activeZone.y2 + 1 * azKeepClearExtra;
            zz1 = activeZone.z1 - 1 * azKeepClearExtra;
            zz2 = activeZone.z2 + 1 * azKeepClearExtra;

            keepClear = new Coordinates(project, xx1, yy1, zz1, xx2, yy2, zz2);

        }

        for (DiffusantParticles d : diffusants) {

            if (d == null) {
                continue;
            }

            //Master.log(keepClear.printDimensions());

            if (d.initImmobileParticles(keepClear)){
                return true;
            }

        }

        return false;

    }

    public boolean initVesicleRandom(DiffusantParticle dv, Coordinates c, Coordinates exclude, boolean allowOverlaps, boolean excludeAZ, double distanceFromAZ, long trialLimit, double az_z) {

        long trial = 0;
        double x, y, z;
        double dx, dy, dz, d2AZ;
        boolean overlap;

        if (dv == null) {
            return true;
        }

        if (c == null) {
            return true;
        }

        while (true) {

            x = (Master.mt.nextDouble() * (c.x2 - c.x1)) + c.x1;
            y = (Master.mt.nextDouble() * (c.y2 - c.y1)) + c.y1;

            if (Double.isNaN(az_z)) {
                z = (Master.mt.nextDouble() * (c.z2 - c.z1)) + c.z1;
            } else {
                z = az_z;
            }

            if (!c.isInside(x, y, z, dv.radius)) {
                continue;
            }

            if ((exclude != null) && exclude.isInside(x, y, z, dv.radius)) {
                continue;
            }

            if (distanceFromAZ > 0) {

                dx = x - activeZone.xCenter;
                dy = y - activeZone.yCenter;
                dz = z - activeZone.zCenter;

                d2AZ = Math.sqrt(dx * dx + dy * dy + dz * dz);

                if (d2AZ < distanceFromAZ) {
                    continue;
                }

            }

            if (excludeAZ && (activeZone != null)) {
                if (azHemisphere) {
                    if (!activeZone.isOutsideEllipsoid(x, y, z, dv.radius)) {
                        continue;
                    }
                } else {
                    if (activeZone.isInside(x, y, z)) {
                        continue;
                    }
                }
            }

            if (!setParticleLocation(dv, x, y, z, true)){
                continue;
            }

            if (outOfBounds(dv)) {
                continue;
            }

            if (testOverlap(dv)) {
                continue;
            }

            overlap = false;

            if (!freeDiffusion) {
                overlap = testParticleOverlap(dv, null) != null;
            }
            
            if (!overlap) {
                return (!addToVoxelList(dv));
            }
            
            trial++;
            
            if (trial > trialLimit) {
                if (allowOverlaps) {
                    return (!addToVoxelList(dv));
                } else {
                    Master.log("initVesicleRandom warning: reached trial limit");
                    return true;
                }
                
            }

        }

    }

    @Override
    public boolean setParticleLocation(DiffusantParticle dv, double x, double y, double z, boolean initStartLocation) {

        int xVoxel, yVoxel, zVoxel;

        if (PBC) {

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

            dv.x = x;
            dv.y = y;
            dv.z = z;

            if (initStartLocation) {
                dv.x0 = x;
                dv.y0 = y;
                dv.z0 = z;
            }

            dv.voxel = geometry.voxelSpacePBC[xVoxel][yVoxel][zVoxel];

            return true;

        } else {

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

            dv.x = x;
            dv.y = y;
            dv.z = z;

            if (initStartLocation) {
                dv.x0 = x;
                dv.y0 = y;
                dv.z0 = z;
            }

            dv.voxel = geometry.voxelSpace[xVoxel][yVoxel][zVoxel];

            return true;

        }

    }

    @Override
    public boolean initParticlesRandom() {

        boolean allowOverlaps = true;

        Coordinates c;
        DiffusantVesicleAZ vaz;

        if (diffusants == null) {
            return true;
        }

        Master.log("randomly placing ready vesicles...");

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            if (d.coordinates == null) {
                c = new Coordinates(project, geometry);
            } else {
                c = new Coordinates(project, d.coordinates);
            }

            for (DiffusantParticle v : d.particles) {

                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (!vaz.isReady) {
                    continue;
                }

                if (initVesicleRandom(vaz, c, null, allowOverlaps, azExcludeReadyinit, -1, initParticleRandomTrialLimit, Double.NaN)) {
                    return true;
                }

            }

        }

        return false;

    }

    @Override
    public boolean removeParticleOverlap() {

        int radiusNM, lastRadiusNM = 0, count = 0;
        double halfDistance = 0;
        boolean saveDockingOn = dockingOn;

        if (freeDiffusion) {
            return false; // OK
        }

        if (!particleOverlapExists()){
            return false; // OK
        }

        if (diffusants == null) {
            return true;
        }

        particleStats(); // will update vesicles step size

        removingParticleOverlap = true;
        dockingOn = false;

        Master.log("elimination of vesicle overlap...");

        while (halfDistance < minParticleRadius) {

            halfDistance = 0.5 * minDistance();

            radiusNM = (int) (1000 * halfDistance);

            if (radiusNM > lastRadiusNM) {
                Master.log("current radius: " + radiusNM + " nm");
                lastRadiusNM = radiusNM;
                count = 0;
            }

            // set current diameter to minimal distance between vesicles

            for (DiffusantParticles d : diffusants) {

                //if ((diffusant[j] == null) || (diffusant[j].vesicle == null)) {
                //    continue;
                //}

                for (DiffusantParticle v : d.particles) {
                    //if (v.insideGeometry) {
                        //v.radius = Math.min(halfDistance, diffusant[j].meanRadius);
                        v.radius = halfDistance;
                    //}
                }

            }

            moveVesiclesActiveZone();

            count++;

            if (count > removeParticleOverlapTrialLimit) {
                break;
            }

        }

        dockingOn = saveDockingOn;
        removingParticleOverlap = false;

        if (particleOverlapExists()) {

            error("removeVesicleOverlap: failed to remove vesicle overlap.");

            for (DiffusantParticles d : diffusants) {

                //if ((d == null) || (d.vesicles == null)) {
                    //continue;
                //}

                for (DiffusantParticle v : d.particles) {
                    v.radius = d.radiusMean;
                    v.x = v.x0;
                    v.y = v.y0;
                    v.z = v.z0;
                }

            }

            return true;

        }

        for (DiffusantParticles d : diffusants) {

            //if ((d == null) || (d.vesicles == null)) {
                //continue;
            //}

            for (DiffusantParticle v : d.particles) {
                v.radius = d.radiusMean;
                v.x0 = v.x;
                v.y0 = v.y;
                v.z0 = v.z;
            }

        }

        return false;

    }

    public boolean initVesiclesDocked() {

        double az_z;
        double zwidth;
        double vesVolume;

        int count = 0;
        int initNum, j;

        boolean allowOverlaps = true;
        boolean excludeAZ = false;

        DiffusantVesiclesAZ daz;
        DiffusantVesicleAZ vaz;

        if (!dockingOn) {
            return false;
        }

        if (diffusants == null) {
            return true;
        }

        zwidth = 0.5 * azHeight;

        //Master.log("randomly placing docked vesicles...");

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            if (!(d instanceof DiffusantVesiclesAZ)) {
                continue;
            }

            daz = (DiffusantVesiclesAZ) d;

            if (daz.dockedVesiclesInit > 0) {
                initNum = daz.dockedVesiclesInit;
            } else if (daz.dockedVesiclesInitFraction > 0) {
                vesVolume = 4.0 * Math.PI * Math.pow(daz.setRadiusMean, 3) / 3.0;
                initNum = (int) (daz.dockedVesiclesInitFraction * activeZone.geometryVolume() / vesVolume) ;
            } else {
                continue;
            }

            //initVesiclesLattice(activeZone, diffusant[j], initNum);

            //if (true) {
                //continue;
            //}

            //Master.log("placing docked vesicles n = " + initNum);

            count = 0;
            j = 0;

            for (DiffusantParticle v : d.particles) {

                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (!vaz.isReady) {
                    continue;
                }

                if (azHemisphere) {
                    az_z = Double.NaN;
                } else {
                    az_z = geometry.z2 - zwidth;
                }

                if (initVesicleRandom(vaz, activeZone, null, allowOverlaps, excludeAZ, -1, initParticleRandomTrialLimit, az_z)) {
                    Master.log("failed to place docked vesicle #" + j++);
                    return true;
                }

                vaz.setType("docked");

                count++;

                if (count >= initNum) {
                    break; // finished
                }

            }

        }

        Master.log("placed docked vesicles: " + count);

        return false;

    }

    public int initVesiclesLattice(Coordinates c, DiffusantParticles dvs, int numVesicles) {

        double r, jj;
        double x, y, z;
        int iv;
        int i, j, k, ii;

        double maxPacking = Math.PI / (3 * Math.sqrt(2));
        double zstep = Math.sqrt(6) * 2.0 / 3.0;
        double sqrt3 = Math.sqrt(3);

        double xlimit = c.x2 - 0 * dvs.setRadiusMean;
        double ylimit = c.y2 - 0 * dvs.setRadiusMean;
        double zlimit = c.z2 - 0 * dvs.setRadiusMean;

        boolean finished = false;

        DiffusantVesicleAZ vaz;

        //nVesicles = dv.length;
        //newVesicleVolume = maxPacking * 1.0 / nVesicles;
        //r = Math.pow((3 * newVesicleVolume / (4 * Math.PI)), (1 / 3.0));
        r = dvs.setRadiusMean + 0.00000001;

        iv = 0;
        k = 0;

        while (!finished) {

            z = geometry.z2 - r - k * zstep * r;

            if (z < c.z1) {
                break;
            }

            j = 0;
            jj = (k % 2) / 3.0;

            while (!finished) {

                y = c.y1 + r + sqrt3 * (j + jj) * r;

                if (y > ylimit) {
                    break;
                }

                i = 0;
                ii = (j + k) % 2;

                while (!finished) {

                    x = c.x1 + r + (2 * i + ii) * r;

                    if (x > xlimit) {
                        break;
                    }

                    if (geometry.isInside(x, y, z) && c.isInside(x, y, z)) {

                        vaz = (DiffusantVesicleAZ) dvs.particles[iv];

                        if (!setParticleLocation(vaz, x, y, z, true)) {
                            return iv;
                        }

                        //Master.log("z = " + z);

                        if (!addToVoxelList(vaz)) {
                            return iv;
                        }

                        vaz.setType("docked");
                        vaz.insideGeometry = true;

                        iv++;

                        if ((iv == numVesicles) || (iv == dvs.particles.length)) {
                            finished = true;
                            break;
                        }

                    }

                    i++;

                }

                j++;

            }

            k++;

        }

        Master.log("placed docked vesicles : " + iv + " / " + numVesicles);

        return iv;

    }

    public boolean initDockedList() {

        DiffusantVesicleAZ vaz;

        dockedVesicles = 0;

        firstDocked = null;

        if (diffusants == null) {
            return true;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle v : d.particles) {

                if (v == null) {
                    continue;
                }

                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (vaz.isDocked) {
                    addToDockedList(vaz, -99999);
                }

            }

        }

        Master.log("initial docked vesicles: " + dockedVesicles);

        return false;

    }

    public boolean addToDockedList(DiffusantVesicleAZ vaz, double time) {

        if ((dockedVesiclesMax >= 0) && (dockedVesicles >= dockedVesiclesMax)) {
            return false;
        }

        vaz.setType("docked");
        vaz.dockTime = time;
        vaz.nextDocked = firstDocked;
        firstDocked = vaz;
        dockedVesicles += 1;

        //System.out.println("docked vesicle " + dv.dockTime + ", type " + dv.name);

        return true;

    }

    public boolean removeFromDockedList(DiffusantVesicleAZ vaz) {

        DiffusantVesicleAZ idv = firstDocked;
        DiffusantVesicleAZ idvnext;
        DiffusantVesicleAZ idvhold = null;

        boolean found = false;

        while (idv != null) {

            idvnext = idv.nextDocked;

            if (idv == vaz) {
                if (idvhold == null) {
                    firstDocked = idvnext; // dv was first list
                } else {
                    idvhold.nextDocked = idvnext;
                }
                found = true;
            }

            idvhold = idv;
            idv = idvnext;

        }

        if (found) {
            dockedVesicles -= 1;
            vaz.dockTime = -1;
            vaz.name = "";
        }

        return found;

    }

    public DiffusantVesicleAZ lastDockedVesicle() {

        DiffusantVesicleAZ idv = firstDocked;

        while (idv != null) {

            if (idv.nextDocked == null) {
                return idv;
            }

            idv = idv.nextDocked;

        }

        return null;

    }

    public boolean releaseDockedVesicle(DiffusantVesicleAZ vaz, double time) {

        double d;
        double dlimit;

        if (vaz == null) {
            //Master.log("failed to release " + time);
            return false;
        }

        //if (dockedVesicles <= 0) {
            //Master.log("release error, no docked vesicles " + dv);
            //return false;
        //}

        if (vaz.dockRefractoryPeriod > 0) {
            if (time - vaz.dockTime < vaz.dockRefractoryPeriod) {
                return false; // still in refractory period
            }
        }

        if ((releaseProb < 1.0) && (Master.mt.nextDouble() > releaseProb)) {
            //Master.log("release failure at time " + time);
            return false;
        }

        d = vaz.squareDisplacement();
        d = Math.sqrt(d);
        //dlimit = Math.sqrt(6 * time * dv.D);

        //if ( (d > 0.5) && (d > dlimit)) {
        //    Master.log("unusually fast vesicle at t = " + time);
        //    Master.log("" + dv);
        //    Master.log("d=" + d + ", dlimit=" + dlimit);
        //    Master.log("x0=" + dv.x0 + ", y0=" + dv.y0 + ",z0=" + dv.z0);
        //    Master.log("x=" + dv.x + ", y=" + dv.y + ",z=" + dv.z);
        //    Master.log("dx=" + Math.abs(dv.x - dv.x0) + ", dy=" + Math.abs(dv.y - dv.y0) + ",dz=" + Math.abs(dv.z - dv.z0));
        //}

        if (connectParticles) {

            vaz.unconnectAll();

            if (vaz.isTetheredToAZ) {
                azDV.unconnectFrom(vaz);
                vaz.isTetheredToAZ = false;
            }

        }

        //dv.z = geometry.z2 + project.dx * 0.5; // move out of the way
        vaz.x = geometry.x2 + project.dx * 1.0; // move out of the way
        vaz.y = geometry.y1 + releaseCounter * project.dx;

        vaz.residenceTime = -1; // do not include released vesicles in residence time

        //Master.log("released docked vesicle " + releaseCounter + " at time " + time + ", docking time " + vaz.dockTime);

        removeFromVoxelList(vaz, vaz.voxel);
        removeFromDockedList(vaz);
        addToReserveList(vaz);
        releaseCounter++;

        //Toolkit.getDefaultToolkit().beep();

        //d2az = this.activeZone.distanceToCenter(d, d, d);

        
        //Master.log("released docked vesicle " + releaseCounter + ", docking time " + (time - dv.dockTime));
        //Master.log("" + d + ", " + dv.d2AZ);

        if (!preview && saveRelease) {
            if (vaz.d2AZ == 0) {
                saveReleaseData(time, 0); // reserve vesicle
            } else {
                saveReleaseData(time, d);
            }
        }

        return true;

    }

    private boolean releaseDockedVesicleCooperative(DiffusantVesicleAZ vaz) {

        double dx, dy, dz, dr, ux, uy, uz;
        double dvx, dvy, dvz;
        double rstep = 0.0001; // um
        int zsteps;
        int jselect = 0;

        boolean sameVoxel, move;

        double[] moved;

        DiffusantVesicleAZ testVesicle = new DiffusantVesicleAZ(project, "ready", 0, 0, 0, 0, 0);
        DiffusantParticle dv2;
        DiffusantParticle dv3, odv;

        int numConnections = vaz.numberOfConnections(false);

        dz = geometry.z2 + vaz.radius - vaz.z;
        zsteps = (int) (dz / rstep);

        dvx = vaz.x;
        dvy = vaz.y;
        dvz = vaz.z;

        vaz.z = geometry.z2 + vaz.radius;

        if (numConnections == 0) {
            return true;
        }

        moved = new double[vaz.connectTo.length];

        Master.log("releasing vesicle...");
        Master.log("connections: " + numConnections);
        //Master.log("dz (nm): " + (1000 * dz)); 50 nm
        //Master.log("zsteps: " + zsteps);

        for (int j = 0; j < vaz.connectTo.length; j++) {
            if (vaz.connectTo[j] == null) {
                moved[j] = Double.NaN;
            } else {
                moved[j] = 0;
                jselect = j;
            }
        }

        for (int i = 0; i < zsteps; i++) {

            dvz += rstep;

            //for (int j = 0; j < dv.connectTo.length; j++) {
            for (int j = jselect; j <= jselect; j++) {

                if (vaz.connectTo[j] == null) {
                    continue;
                }

                dv2 = vaz.connectTo[j];

                move = true;

                testVesicle.copy(dv2);

                dx = dvx - dv2.x;
                dy = dvy - dv2.y;
                dz = dvz - dv2.z;

                dr = Math.sqrt(dx * dx + dy * dy + dz * dz);

                ux = dx / dr; // unit vector towards dv
                uy = dy / dr; // unit vector towards dv
                uz = dz / dr; // unit vector towards dv

                dx = ux * rstep;
                dy = uy * rstep;
                dz = uz * rstep;

                //dx = 0;
                //dy = 0;
                //dz = rstep;

                if (!setParticleLocation(testVesicle, testVesicle.x + dx, testVesicle.y + dy, testVesicle.z + dz, false)) {
                    move = false;
                    //Master.log("failed to set vesicle location " + (1000 * j * rstep));
                }

                if (move && outOfBounds(testVesicle)) {
                    move = false;
                    //Master.log("out of bounds " + (1000 * j * rstep));
                }

                if (move) {

                    odv = testParticleOverlap(testVesicle, dv2);

                    if (odv != null) {
                        move = false;
                        //Master.log("vesicle overlap " + odv + ", " + (1000 * j * rstep));
                    }

                }

                if (move) {

                    for (int k = 0; k < dv2.connectTo.length; k++) {

                        dv3 = dv2.connectTo[k];

                        if ((dv3 == null) || (dv3 == vaz)) {
                            continue;
                        }

                        if (!testVesicle.overlap(dv3, connectorLength)) {
                            dv3.unconnectFrom(dv2);
                            dv2.connectTo[k] = null;
                            dv2.connectorOffTime[k] = 0;
                            //Master.log("lost connection with " + dv3 + ", " + (1000 * j * rstep));
                        }

                    }

                }

                if (!move) {
                    // break connection
                    dv2.unconnectFrom(vaz);
                    vaz.connectTo[j] = null;
                    vaz.connectorOffTime[j] = 0;
                    continue;
                }

                // OK to move vesicle
                
                sameVoxel = (testVesicle.voxel == dv2.voxel);

                if (!sameVoxel) {
                    removeFromVoxelList(dv2, dv2.voxel);
                }

                dv2.copy(testVesicle);

                if (!sameVoxel) {
                    addToVoxelList(dv2);
                }

                moved[j] += rstep;

                //Master.log("moved vesicle " + dv2);

            }

        }

        //Master.log("released vesicle " + dv);

        return true;

    }

    @Override
    public void connectParticles() {

        int neighbors;
        double sqrDistance;
        boolean inside;

        DiffusantVesicleAZ vaz, kvesicle;

        if (tetherVesicles && (azDV == null)) {
            Master.exit("azDV is null");
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle v : d.particles) {

                vaz = (DiffusantVesicleAZ) v;

                if (vaz.isReserve) {
                    continue;
                }

                if (tetherVesicles && !vaz.isTetheredToAZ && !azDV.noMoreConnectors()) {

                    inside = isInsideActiveZoneTethering(vaz);

                    if (inside) {
                        if (azDV.connectToNew(vaz, Double.POSITIVE_INFINITY, time)) {
                            vaz.isTetheredToAZ = true;
                            //Master.log("tethered to AZ " + dv);
                        } else {
                            Master.exit("tether to AZ failure ");
                        }
                    }

                }

                if (vaz.noMoreConnectors()) {
                    continue;
                }

                sqrDistance = vaz.squareDistanceFromActiveZone(activeZone);

                if (sqrDistance >= azPermitConnectorRadius * azPermitConnectorRadius) {
                    continue;
                }

                neighbors = vaz.voxel.numNeighbors;

                for (int k = 0; k < neighbors; k++) {

                    kvesicle = (DiffusantVesicleAZ) vaz.voxel.neighbors[k].firstParticleInChain;

                    while (kvesicle != null) {

                        if (vaz.noMoreConnectors()) {
                            break;
                        }

                        if (kvesicle.isReady && (kvesicle != vaz) && !vaz.isConnectedTo(kvesicle) && !kvesicle.noMoreConnectors()) {

                            if (vaz.overlap(kvesicle, connectorLength)) {

                                sqrDistance = kvesicle.squareDistanceFromActiveZone(activeZone);

                                if (sqrDistance < azPermitConnectorRadius * azPermitConnectorRadius) {

                                    if (Master.mt.nextDouble() < connectRate * project.dt) {

                                        if (kvesicle.connectToNew(vaz, unconnectRate, time)) {

                                            avgConnectorLifeTime += vaz.connectorLifeTime;
                                            connectorLifeTimeCounter += 1;

                                            if (!vaz.connectTo(kvesicle)) {
                                                Master.exit("connection failure");
                                            } else {
                                                //Master.log("connected vesicles");
                                            }

                                        } else {
                                            Master.exit("connection failure");
                                        }

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

    public double computeDshort() {

        double sqrDistance, avgD = 0, count = 0;

        DiffusantVesicleAZ vaz;

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle v : d.particles) {

                vaz = (DiffusantVesicleAZ) v;

                if (vaz.isReserve) {
                    continue;
                }

                sqrDistance = vaz.squareDistanceFromActiveZone(activeZone);

                if (sqrDistance >= azPermitConnectorRadius * azPermitConnectorRadius) {
                    continue;
                }

                avgD += vaz.Dlocal;
                count++;

            }

        }

        if (count <= 0) {
            return 0;
        }

        avgD /= count;

        Master.log("Dshort = " + avgD + " um^2/ms");

        return avgD;

    }

    @Override
    public double[] computeNumberOfConnectors() {

        double sqrDistance, avg = 0, connected = 0;
        int n, ntotal = 0;

        double[] results = new double[2];

        DiffusantVesicleAZ vaz;

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle v : d.particles) {

                vaz = (DiffusantVesicleAZ) v;

                if (vaz.isReserve) {
                    continue;
                }

                sqrDistance = vaz.squareDistanceFromActiveZone(activeZone);

                if (sqrDistance >= azPermitConnectorRadius * azPermitConnectorRadius) {
                    continue;
                }

                n = vaz.numberOfConnections(true);
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

    public double[] computeVesicleDensityAZ() {

        double dmin;
        double zwin, zstep, zlimit = 0.5;
        double AZradius = azWidth / 2.0;
        double vesVolume = 0;

        DiffusantVesicleAZ vaz;

        zstep = 0.03;

        int zsteps = (int) (zlimit / zstep);
        int iz;

        double[] count = new double[zsteps];
        double[] volume = new double[zsteps];
        double[] density = new double[zsteps];

        for (int zcnt = 0; zcnt < zsteps; zcnt++) {
            zwin = (zcnt + 1) * zstep;
            volume[zcnt] = ((AZradius * AZradius) * Math.PI * zwin) + (2.0 * (zwin * zwin * zwin) * Math.PI / 3.0) + (0.5 * AZradius * (zwin * zwin) * (Math.PI * Math.PI));
            count[zcnt] = 0;
        }

        for (int zcnt = zsteps - 1; zcnt >= 1; zcnt--) {
            volume[zcnt] = volume[zcnt] - volume[zcnt - 1];
        }

        //for (int zcnt = 0; zcnt < zsteps; zcnt++) {
        //    Master.log("" + zcnt + " volume = " + volume[zcnt]);
        //}

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle v : d.particles) {

                vaz = (DiffusantVesicleAZ) v;

                if (vaz.isReserve) {
                    continue;
                }

                if (vesVolume == 0) {
                    vesVolume = DiffusantParticle.volume(vaz.radius); // all the same
                }

                dmin = minDistanceToAZ(vaz.x, vaz.y, vaz.z);

                iz = (int) (dmin / zstep);

                if ((iz >= 0) && (iz < count.length)) {
                    count[iz]++;
                }

            }

        }

        for (int zcnt = 0; zcnt < zsteps; zcnt++) {
            density[zcnt] = count[zcnt] / volume[zcnt];
            density[zcnt] *= vesVolume;
            //Master.log("" + count[zcnt]);
            Master.log("" + density[zcnt]);
        }

        return density;

    }

    public double minDistanceToAZ(double x, double y, double z) {

        double d, dx2, dy2;

        double AZx0 = 0;
        double AZy0 = 0;
        double AZradius = azWidth * 0.5;

        double dx = AZx0 - x;
        double dy = AZy0 - y;
        double dz = geometry.z2 - z;
        double dxy = Math.sqrt(dx * dx + dy * dy);

        if (dxy < AZradius) {

            d = dz;

        } else {

            dx2 = AZradius * dx / dxy;
            dy2 = AZradius * dy / dxy;
            dx2 = dx - dx2;
            dy2 = dy - dy2;

            d = Math.sqrt(dx2 * dx2 + dy2 * dy2 + dz * dz);

        }

        return d;

    }

    public double[] computeVesicleDensityAZ2(boolean print) {

        int numVes;
        double d, zstep, zlimit = 0.5;

        double vesVolume = 0;

        DiffusantParticle dv;

        zstep = 0.03;

        int zsteps = (int) (zlimit / zstep);
        int ibin;

        double[] count = new double[zsteps];
        double[] volume = new double[zsteps];
        double[] density = new double[zsteps];

        for (int zcnt = 0; zcnt < zsteps; zcnt++) {
            count[zcnt] = 0;
            volume[zcnt] = 0;
        }
        
        for (Voxel[][] i : geometry.voxelSpace) {
            for (Voxel[] j : i) {
                for (Voxel v : j) {

                    if (!v.isSpace) {
                        continue;
                    }

                    numVes = 0;

                    dv = v.firstParticleInChain;

                    if ((dv != null) && (vesVolume == 0)) {
                        vesVolume = DiffusantParticle.volume(dv.radius); // all the same
                    }

                    while (dv != null) {
                        numVes++;
                        dv = dv.nextParticleInChain;
                    }

                    d = minDistanceToAZ(v.x, v.y, v.z);

                    if (Double.isNaN(d) || Double.isInfinite(d)) {
                        continue;
                    }

                    ibin = (int) (d / zstep);

                    if ((ibin >= 0) && (ibin < volume.length)) {
                        volume[ibin] += voxelVolume;
                        count[ibin] += numVes;
                    }

                }
            }
        }

        for (int zcnt = 0; zcnt < zsteps; zcnt++) {
            density[zcnt] = count[zcnt] / volume[zcnt];
            density[zcnt] *= vesVolume;
            //Master.log("" + count[zcnt]);
            if (print) {
                Master.log("" + density[zcnt]);
            }
        }

        return density;

    }

    private double distanceToAZ(Coordinates c) {

        double dx, dy, dz;
        double d = Double.NaN;

        if (c == null) {
            return Double.NaN;
        }

        switch(c.shapeNum) {
                case 2:
                    dy = c.yCenter - activeZone.yCenter;
                    dz = c.zCenter - activeZone.zCenter;
                    d = Math.sqrt(dy * dy + dz * dz);
                    d -= c.yWidth * 0.5;
                    break;
                case 3:
                    dx = c.xCenter - activeZone.xCenter;
                    dz = c.zCenter - activeZone.zCenter;
                    d = Math.sqrt(dx * dx + dz * dz);
                    d -= c.xWidth * 0.5;
                    break;
                case 4:
                    dx = c.xCenter - activeZone.xCenter;
                    dy = c.yCenter - activeZone.yCenter;
                    d = Math.sqrt(dx * dx + dy * dy);
                    d -= c.xWidth * 0.5;
                    break;
            }

        return d;

    }

    public boolean initReserveList() {

        DiffusantVesicleAZ vaz;

        reserveVesicles = 0;

        firstReserve = null;

        if (diffusants == null) {
            return true;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle v : d.particles) {

                if (v == null) {
                    continue;
                }

                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (vaz.isReserve) {
                    addToReserveList(vaz);
                }

            }

        }

        Master.log("initial reserve vesicles: " + reserveVesicles);

        return false;

    }

    public boolean addToReserveList(DiffusantVesicleAZ vaz) {

        if ((reserveVesiclesMax >= 0) && (reserveVesicles >= reserveVesiclesMax)) {
            return false;
        }

        vaz.setType("reserve");
        vaz.nextReserve = firstReserve;
        firstReserve = vaz;
        reserveVesicles += 1;

        return true;

    }

    public boolean removeFromReserveList(DiffusantVesicleAZ vaz) {

        DiffusantVesicleAZ idv = firstReserve;
        DiffusantVesicleAZ idvnext;
        DiffusantVesicleAZ idvhold = null;

        boolean found = false;

        while (idv != null) {

            idvnext = idv.nextReserve;

            if (idv == vaz) {
                if (idvhold == null) {
                    firstReserve = idvnext; // dv was first list
                } else {
                    idvhold.nextReserve = idvnext;
                }
                found = true;
            }

            idvhold = idv;
            idv = idvnext;

        }

        if (found) {
            reserveVesicles -= 1;
        }

        if (vaz.isReserve){
            vaz.name = "";
        }

        return found;

    }

     private int moveReserveToReady(int nvesicles) {

        int count = 0;
        double x, y, z, dx, dy, dz;
        double d2AZ;
        boolean error;

        boolean allowOverlaps = false;
        boolean excludeAZ = true;

        DiffusantVesicleAZ vaz;

        if (nvesicles <= 0) {
            return 0;
        }

        for (int i = 0; i < nvesicles; i++) {

            vaz = firstReserve;

            if (vaz == null) {
                return count;
            }

            x = vaz.x;
            y = vaz.y;
            z = vaz.z;

            error = initVesicleRandom(vaz, replenishCoordinates, replenishCoordinatesExclude, allowOverlaps, excludeAZ, replenishFromAZdistance, initParticleRandomTrialLimit, Double.NaN);
            
            if (error) {

                vaz.x = x;
                vaz.y = y;
                vaz.z = z;

                return count;

            } else {

                //dx = dv.x - activeZone.xCenter;
                //dy = dv.y - activeZone.yCenter;
                //dz = dv.z - activeZone.zCenter;
                //d2AZ = Math.sqrt(dx * dx + dy * dy + dz * dz);
                //Master.log("added reserve vesicle: " + dv + ", distance to AZ = " + d2AZ);

                vaz.setType("ready");
                removeFromReserveList(vaz);
                vaz.isTetheredToAZ = false;
                vaz.isDocked = false;
                vaz.isReserve = false;
                
                count++;

            }

        }

        return count;

    }

    //public boolean setParticleLocation(DiffusantVesicleAZ vaz, double x, double y, double z, boolean initStartLocation) {
        //return super.setParticleLocation(vaz, x, y, z, initStartLocation);
    //}

    public boolean moveVesicle(DiffusantVesicleAZ vaz) {

        double step3;
        //double stepx, stepy, stepz;

        if (hydrodynamicsLocalD || hydroWallZ) {
            step3 = vaz.step3Local;
        } else {
            step3 = vaz.step3;
        }

        //if (rstep < 0) {
        //    Master.exit("moveParticle: bad rstep: " + rstep);
        //    return false;
        //}

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

        if (vaz.isDocked) {
            stepz = 0;
        }

        return setParticleLocation(vaz, vaz.x + stepx, vaz.y + stepy, vaz.z + stepz, false);

    }

    public boolean moveVesicleGauss(DiffusantVesicleAZ vaz) {

        double step3;
        //double stepx, stepy, stepz;

        if (hydrodynamicsLocalD) {
            step3 = vaz.step3Local;
        } else {
            step3 = vaz.step3;
        }

        if (removingParticleOverlap) {
            step3 = removeParticleOverlapStep3;
        }

        //if (rstep <= 0) {
        //    Master.exit("moveParticleGauss: bad rstep: " + rstep);
        //}

        stepx = Master.randomGauss() * step3;
        stepy = Master.randomGauss() * step3;
        stepz = Master.randomGauss() * step3;

        if (vaz.isDocked) {
            stepz = 0;
        }

        return setParticleLocation(vaz, vaz.x + stepx, vaz.y + stepy, vaz.z + stepz, false);

    }

    public boolean moveVesicleGaussHydroWallz(DiffusantVesicleAZ vaz) {

        double step3;
        double dz, b_ll, b_T;

        if (hydrodynamicsLocalD) {
            step3 = vaz.step3Local;
        } else if (hydrodynamicsDscale) {
            //step3 = localDensityFastHydroWallz(dv); // using Michailidou model instead
            step3 = vaz.step3;
        } else {
            step3 = vaz.step3;
        }

        if (removingParticleOverlap) {

            step3 = removeParticleOverlapStep3;

            stepx = Master.randomGauss() * step3;
            stepy = Master.randomGauss() * step3;
            stepz = Master.randomGauss() * step3;

        } else {

            dz = geometry.z2 - vaz.z;

            //b_ll = Math.sqrt(hydroWall_ll(dv.radius, z));
            //b_T = Math.sqrt(hydroWall_T(dv.radius, z - dv.radius));

            b_ll = hydroWall_ll(vaz.radius, dz);
            b_T = hydroWall_T(vaz.radius, Math.abs(dz - vaz.radius));

            b_ll = 1 / (1 + vaz.DsDff * ((1 / b_ll) - 1)); // Michailidou et al. 2009
            b_T = 1 / (1 + vaz.DsDff * ((1 / b_T) - 1)); // Michailidou et al. 2009

            b_ll = Math.sqrt(b_ll);
            b_T = Math.sqrt(b_T);

            stepx = Master.randomGauss() * step3 * b_ll;
            stepy = Master.randomGauss() * step3 * b_ll;
            stepz = Master.randomGauss() * step3 * b_T;

        }

        //stepx = step3 * b_ll;
        //stepy = step3 * b_ll;
        //stepz = step3 * b_T;

        //dv.localD = (stepx * stepx + stepy * stepy + stepz * stepz) / (6 * project.dt);
        //Master.log("" + (dv.localD/1.397657228105222E-4));
        //dv.localDw = (3 * step3 * step3) / (6 * project.dt);

        return setParticleLocation(vaz, vaz.x + stepx, vaz.y + stepy, vaz.z + stepz, false);

    }

    public int moveVesiclesActiveZone() {

        int moved = 0;
        double distance;
        boolean outOfBounds, overlap, sameVoxel, insideAZ = false, outsideAZ = true, moveTethered;
        boolean ok;

        DiffusantVesicleAZ vaz;
        DiffusantVesicleAZ testVesicle = new DiffusantVesicleAZ(project, "ready", 0, 0, 0, 0, 0);

        if (diffusants == null) {
            return 0;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            d.shuffleParticles();

            for (DiffusantParticle v : d.particles) {

                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (!vaz.mobile || vaz.isDocked || vaz.isReserve || vaz.isBound) {
                    continue;
                }

                testVesicle.copy(vaz);

                if (hydroWallZ) {
                    ok = moveVesicleGaussHydroWallz(testVesicle);
                    if (!ok) {
                        //    continue;
                    }
                } else {
                    ok = moveVesicleGauss(testVesicle);
                    if (!ok) {
                    //    continue;
                    }
                }

                if (ok && !removingParticleOverlap) {

                    if (azHemisphere) {
                        outsideAZ = activeZone.isOutsideEllipsoid(testVesicle.x, testVesicle.y, testVesicle.z, vaz.radius);
                    }

                    insideAZ = isInsideActiveZone(testVesicle);

                    if (dockingOn && vaz.isDocked && !insideAZ) {
                        ok = false;
                        //continue; // dont let docked vesicles move off AZ
                    }

                }

                if (ok && tetherVesicles && vaz.isTetheredToAZ) {
                    if (!isInsideActiveZoneTethering(testVesicle)) {
                        ok = false;
                        //continue; // dont let tethered vesicles move outside AZ tethering region
                    }
                }

                if (ok) {
                    outOfBounds = outOfBounds(testVesicle);
                } else {
                    outOfBounds = false;
                }

                if (ok && outOfBounds) {
                    if (dockingOn && insideAZ) {

                        testVesicle.z = activeZone.z1;

                    } else if (PBC) {

                        outOfBounds = wrapAtBorder(testVesicle);

                        //if (azPerimeter) {
                        //    insideAZ = activeZone.isInsideCylinderZ_Perimeter(testVesicle.x, testVesicle.y, testVesicle.z, azPerimeterRatio);
                        //} else {
                        //    insideAZ = activeZone.isInsideCylinderZ(testVesicle.x, testVesicle.y, testVesicle.z);
                        //}

                        //if (insideAZ) {
                        //    continue;
                        //}

                        

                    } else {
                        ok = false;
                        //continue;
                    }
                }

                if (outOfBounds) {
                    ok = false;
                    //continue;
                }

                if (ok) {
                    if (testOverlap(testVesicle)) {
                        ok = false;
                        //continue;
                    }
                }

                if (ok && !freeDiffusion) {

                    overlap = testParticleOverlap(testVesicle, vaz) != null;

                    if (overlap) {
                        ok = false;
                        //continue;
                    }

                }

                if (ok && connectParticles && (vaz.connectTo != null)) {

                    moveTethered = true;

                    for (DiffusantParticle vaz2 : vaz.connectTo) {

                        if (vaz2 == null) {
                            continue;
                        }

                        if (!testVesicle.overlap(vaz2, connectorLength)) {
                            moveTethered = false;
                            break;
                        }

                    }

                    if (!moveTethered) {
                        ok = false;
                        //continue;
                    }

                }

                if (!ok) {
                    continue;
                }

                moved++;

                if (dockingOn && (time >= dockingOnTime) && !vaz.isDocked) {

                    if (azHemisphere && !outsideAZ) {
                        if (!addToDockedList(vaz, time)) {
                            continue;
                        }
                    } else if (insideAZ) {
                        if (!addToDockedList(vaz, time)) {
                            continue;
                        }
                    }

                }

                sameVoxel = false;

                if (testVesicle.voxel == vaz.voxel) {
                    sameVoxel = true;
                }

                if (!sameVoxel) {
                    removeFromVoxelList(vaz, vaz.voxel);
                }

                setParticleLocation(vaz, testVesicle.x, testVesicle.y, testVesicle.z, false);
                //dv.localD = testVesicle.localD;

                if (!sameVoxel) {
                    addToVoxelList(vaz);
                }

                //if (time > 0) {
                //    dv.totalDisplacement += Math.sqrt(stepx * stepx + stepy * stepy + stepz * stepz);
                //    Master.log("" + (dv.totalDisplacement * dv.totalDisplacement / (6 * time * itime)));
                //}

            }

        }

        return moved;

    }

    public boolean isInsideActiveZone(DiffusantVesicleAZ vaz) {

        if (azHemisphere) {
            return activeZone.isInsideEllipsoid(vaz.x, vaz.y, vaz.z, vaz.radius);
        } else {
            return activeZone.isInsideCylinderZ(vaz.x, vaz.y, vaz.z);
        }

    }

    public boolean isInsideActiveZoneTethering(DiffusantVesicleAZ vaz) {
        return activeZoneTethering.isInsideCylinderZ(vaz.x, vaz.y, vaz.z);
    }

    @Override
    public boolean initSave() {

        super.initSave();

        int dataPoints = 1;

        if (save_Release != null) {
            save_Release.init("Monte Carlo AZ Release Times", null, -1, dataPoints);
        }

        if (save_NumDocked != null) {
            save_NumDocked.init("Monte Carlo AZ Docked Vesicles", null, -1, dataPoints);
        }

        if (save_ConnectorStats != null) {
            save_ConnectorStats.init("Monte Carlo AZ Connector Stats", null, -1, 2);
        }

        if (save_VesicleDensityAZ != null) {
            save_VesicleDensityAZ.init("Monte Carlo AZ Vesicle Density", null, -1, azVesicleDensityZoltan.length);
        }

        return true;

    }

    @Override
    public boolean finishSave() {

        super.finishSave();

        if (save_Release != null) {
            save_Release.finish("Monte Carlo AZ Release Times", null, -1);
        }

        if (save_NumDocked != null) {
            save_NumDocked.finish("Monte Carlo AZ Docked Vesicles", null, -1);
        }

        if (save_ConnectorStats != null) {
            save_ConnectorStats.finish("Monte Carlo AZ Connector Stats", null, -1);
        }

        if (save_VesicleDensityAZ != null) {
            save_VesicleDensityAZ.finish("Monte Carlo AZ Vesicle Density", null, -1);
        }

        return true;

    }

    public boolean saveReleaseData(double t, double d) {
        if (d == Double.NaN) {
            return save_Release.writeData(t);
        } else {
            return save_Release.writeString(Double.toString(t) + '\t' + Double.toString(d));
        }
    }

    public boolean saveNumDockedData() {
        return save_NumDocked.writeData(dockedVesicles);
    }

    public boolean saveConnectorStats() {

        if (save_ConnectorStats.skipCounter == 0) {
            connectorStats = computeNumberOfConnectors();
            //computeDshort();
        } else {
            for (int i = 0; i < connectorStats.length; i++) {
                connectorStats[i] = 0;
            }
        }

        save_ConnectorStats.saveData(connectorStats);

        return true;

    }

    public boolean saveVesicleDensityAZ() {

        if (save_VesicleDensityAZ.skipCounter == 0) {
            vesicleDensityAZ = computeVesicleDensityAZ2(false);
        } else {
            for (int i = 0; i < vesicleDensityAZ.length; i++) {
                vesicleDensityAZ[i] = 0;
            }
        }

        save_VesicleDensityAZ.saveData(vesicleDensityAZ);

        return true;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        if (pv == null) {
            return false;
        }

        super.addUser(pv);

        if (activeZone != null) {
            activeZone.addUser(pv);
        }

        if (save_Release != null) {
            save_Release.addUser(pv);
        }

        if (save_NumDocked != null) {
            save_NumDocked.addUser(pv);
        }

        if (save_ConnectorStats != null) {
            save_ConnectorStats.addUser(pv);
        }

        if (save_VesicleDensityAZ != null) {
            save_VesicleDensityAZ.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (activeZone != null) {
            addBlankParam();
            activeZone.createVector(true);
            addVector(activeZone.getVector());
            activeZone.addUser(this);
        }

        if (save_Release != null) {
            addBlankParam();
            save_Release.createVector(true);
            addVector(save_Release.getVector());
            save_Release.addUser(this);
        }

        if (save_NumDocked != null) {
            addBlankParam();
            save_NumDocked.createVector(true);
            addVector(save_NumDocked.getVector());
            save_NumDocked.addUser(this);
        }

        if (save_ConnectorStats != null) {
            addBlankParam();
            save_ConnectorStats.createVector(true);
            addVector(save_ConnectorStats.getVector());
            save_ConnectorStats.addUser(this);
        }

        if (save_VesicleDensityAZ != null) {
            addBlankParam();
            save_VesicleDensityAZ.createVector(true);
            addVector(save_VesicleDensityAZ.getVector());
            save_VesicleDensityAZ.addUser(this);
        }

        if (saveXYZincrement != null) {
            addBlankParam();
            saveXYZincrement.createVector(true);
            addVector(saveXYZincrement.getVector());
            saveXYZincrement.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (activeZone != null) {
            activeZone.updateVector(v);
        }

        if (save_Release != null) {
            save_Release.updateVector(v);
        }

        if (save_NumDocked != null) {
            save_NumDocked.updateVector(v);
        }

        if (save_ConnectorStats != null) {
            save_ConnectorStats.updateVector(v);
        }

        if (save_VesicleDensityAZ != null) {
            save_VesicleDensityAZ.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof RunMonteCarloAZ)) {
            return false;
        }

        if (n.equalsIgnoreCase("freeDiffusion")) {
            freeDiffusion = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("releaseRate")) {
            releaseRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("releaseStartTime")) {
            releaseStartTime = v;
            return true;
        }
        if (n.equalsIgnoreCase("releaseLimit")) {
            releaseLimit = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("replenishRate")) {
            replenishRate = v;
            return true;
        }
        if (n.equalsIgnoreCase("azWidth")) {
            azWidth = Math.abs(v);
            return true;
        }
        if (n.equalsIgnoreCase("azHeight")) {
            azHeight = Math.abs(v);
            return true;
        }
        if (n.equalsIgnoreCase("azHemisphere")) {
            azHemisphere = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("tetherVesicles")) {
            tetherVesicles = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("azNumTethers")) {
            azNumTethers = (int) v;
            return true;
        }
        return super.setMyParams(o, v);
    }

}
