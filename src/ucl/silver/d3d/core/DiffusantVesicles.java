package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.MersenneTwisterFast;

import java.awt.Color;

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
public class DiffusantVesicles extends Diffusant {

    public double minRadius; // um
    public double maxRadius; // um
    public double meanRadius; // um
    public double setMeanRadius; // um

    public double meanVolume; // um^3
    public double totalVolume; // um^3

    public double meanStep3; // um

    public double maxX0, maxY0, maxZ0; // um

    public transient DiffusantVesicle[] vesicles; // the array of vesicles
    public transient double xyz[][] = null;

    public double setDensity; // vesicles / um^3
    public double setVolumeFraction;

    public double density; // vesicles / um^3
    public double volumeFraction;

    public int numVesicles; // total number of vesicles in array

    public int mobileVesicles, immobileVesicles;
    public double mobileVolumeFraction, immobileVolumeFraction;
    public double mobilePercent, immobilePercent; // e.g. 0.25
    public double setImmobilePercent = 0; // set immobile fraction

    public ColorD3D colorReady = new ColorD3D( "colorReady", Color.white);
    public ColorD3D colorImmobile = new ColorD3D( "colorImmobile", Color.gray);
    public ColorD3D colorConnected = new ColorD3D( "colorConnected", new Color(102,102,0));

    public boolean PV_includeVesicleArray = false;

    public boolean saveXYZ = false; // save vesicles xyz position
    public boolean saveMSD = false; // save mean square distance
    public boolean saveDisplacementAfterCollisions = false;

    public Save save_XYZ = null;
    public Save save_MSD = null;

    //CoordinatesVoxels initCoordinates = null; // coordinates for initial vesicles placement ( leave null for entire shape )

    public MersenneTwisterFast mt = new MersenneTwisterFast(); // create/init random number generator

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("minRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("setMeanRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanVolume")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("totalVolume")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("meanStep3")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanExtraX")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanExtraY")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("meanExtraZ")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxX0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxY0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("maxZ0")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("density")) {
            return "vesicles/" + project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("setDensity")) {
            return "vesicles/" + project.spaceUnits + "^3";
        }
        return super.units(name);
    }

    public DiffusantVesicles(Project p, String NAME, double InitialConcentration, double DiffusionConstant, CoordinatesVoxels c,
            double radius) {

        super(p, NAME, InitialConcentration, DiffusionConstant, c);

        setMeanRadius = radius;

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        if (vesicles != null) {
            for (DiffusantVesicle v : vesicles) {
                if (v != null) {
                    v.init();
                }
            }
        }

        vesicleStats();

        if (save != null) {
            save.init();
            saveFileName();
            saveDimensions();
            save.updateVectors();
        }

        if (saveXYZ) {

            if (save_XYZ == null) {
                save_XYZ = new Save(project);
                save_XYZ.samples2save = 0;
                save_XYZ.skipSamples = 0;
                save_XYZ.sampleRate = -1;
                save_XYZ.sampleInterval = 0;
                save_XYZ.saveWhileComputing = true;
                save_XYZ.save2BinaryFile = false;
                save_XYZ.save2TextFile = true;
            }

            save_XYZ.init();
            saveXYZFileName(-1);
            saveXYZDimensionsInit();
            save_XYZ.updateVectors();

        }

        if (saveMSD) {

            if (save_MSD == null) {
                save_MSD = new Save(project);
                save_MSD.saveWhileComputing = false;
                save_MSD.save2BinaryFile = true;
                save_MSD.save2TextFile = false;
            }

            save_MSD.init();
            saveMSDFileName();
            saveMSDDimensions();
            save_MSD.updateVectors();

        }

    }

    public void vesicleStats() {

        int count = 0;
        double volume;
        double spaceVolume;

        if (project.monteCarlo == null) {
            return;
        }

        if (vesicles == null) {
            return;
        }

        spaceVolume = project.monteCarlo.spaceVolume();

        minRadius = Double.POSITIVE_INFINITY;
        maxRadius = 0;
        meanRadius = 0;

        meanVolume = 0;
        totalVolume = 0;

        D = 0;
        meanStep3 = 0;

        maxX0 = 0;
        maxY0 = 0;
        maxZ0 = 0;

        mobileVesicles = 0;
        immobileVesicles = 0;
        mobilePercent = 0;
        immobilePercent = 0;

        density = 0;
        volumeFraction = 0;
        mobileVolumeFraction = 0;
        immobileVolumeFraction = 0;

        for (DiffusantVesicle v : vesicles) {

            if (v == null) {
                continue;
            }

            if (!v.insideGeometry) {
                continue;
            }

            count++;

            volume = v.getVolume();

            minRadius = Math.min(v.radius, minRadius);
            maxRadius = Math.max(v.radius, maxRadius);

            meanRadius += v.radius;
            meanVolume += volume;
            totalVolume += volume;

            D += v.D;
            meanStep3 += v.step3;

            maxX0 = Math.max(maxX0, Math.abs(v.x0));
            maxY0 = Math.max(maxY0, Math.abs(v.y0));
            maxZ0 = Math.max(maxZ0, Math.abs(v.z0));

            if (v.mobile) {
                mobileVesicles++;
            } else {
                immobileVesicles++;
            }

        }

        numVesicles = immobileVesicles + mobileVesicles;

        meanRadius /= numVesicles;
        meanVolume /= numVesicles;
        D /= numVesicles;
        meanStep3 /= numVesicles;

        immobilePercent = immobileVesicles / ( 1.0 * numVesicles);
        mobilePercent = mobileVesicles / ( 1.0 * numVesicles);

        density = (1.0 * numVesicles) / spaceVolume;
        volumeFraction = totalVolume / spaceVolume;
        mobileVolumeFraction = mobilePercent * volumeFraction;
        immobileVolumeFraction = immobilePercent * volumeFraction;

        setParamObject("minRadius", minRadius);
        setParamObject("maxRadius", maxRadius);
        setParamObject("meanRadius", meanRadius);

        setParamObject("meanVolume", meanVolume);
        setParamObject("totalVolume", totalVolume);

        setParamObject("D", D);
        setParamObject("meanStep3", meanStep3);

        setParamObject("maxX0", maxX0);
        setParamObject("maxY0", maxY0);
        setParamObject("maxZ0", maxZ0);

        setParamObject("numVesicles", numVesicles);
        setParamObject("mobileVesicles", mobileVesicles);
        setParamObject("immobileVesicles", immobileVesicles);
        setParamObject("mobileVesicleFraction", mobilePercent);
        setParamObject("immobilePercent", immobilePercent);

        setParamObject("density", density);
        setParamObject("volumeFraction", volumeFraction);
        setParamObject("mobileVolumeFraction", mobileVolumeFraction);
        setParamObject("immobileVolumeFraction", immobileVolumeFraction);

    }

    public void computeDensity() {

        double sumVolume = 0.0, count = 0.0;

        double spaceVolume = project.monteCarlo.spaceVolume();

        for (DiffusantVesicle v : vesicles) {

            if ((v != null) && v.insideGeometry) {
                count++;
                sumVolume += v.getVolume();
            }

        }

        density = count / spaceVolume;
        volumeFraction = sumVolume / spaceVolume;

        setParamObject("density", density);
        setParamObject("volumeFraction", volumeFraction);

    }

    public static int numVesiclesPossible(double volumeGeometry, double vesicleRadius, double vFraction) {

        double vesicleVolume = 4.0 * Math.PI * Math.pow(vesicleRadius, 3) / 3.0;

        int nVesicles = (int) (vFraction * volumeGeometry / vesicleVolume);

        return nVesicles;

    }

    //@Override
    //public double displayValue(int i, int j, int k) { // already in Diffusant

        //if (concentration != null) {
            //return concentration[i][j][k];
        //}

        //if (psf == null) {
            //return C0;
        //}

        //psf.checkExists();

        //return psf.getArrayValue(i, j, k);

    //}

    @Override
    public boolean saveInit() {

        int dataPoints = 1;

        if (save != null) {
            if (!save.init(name, coordinates(), -1, dataPoints)) {
                return false;
            }
        }
        
        if (save_MSD != null) {
            if (!save_MSD.init(name, coordinates(), -1, dataPoints)) {
                return false;
            }
        }

        return true;

    }

    @Override
    public boolean saveFinish() {

        if (save != null) {
            save.finish(name, coordinates(), -1);
        }

        if (save_MSD != null) {
            save_MSD.finish(name, coordinates(), -1);
        }

        return true;

    }

    @Override
    public boolean save() {

        double svalue = -1;

        if (save_MSD != null) {
            if (save_MSD.skipCounter == 0) {
                if (saveDisplacementAfterCollisions) {
                    svalue = meanDisplacementAfterCollision();
                } else {
                    svalue = meanSquareDisplacement();
                }
            }
            save_MSD.saveData(svalue);
        }

        return true;

    }

    public boolean saveXYZFileName(int msec) {

        if (save_XYZ == null) {
            return false;
        }

        //save_XYZ.fileName("XYZF", name);
        save_XYZ.fileName("XYZ", name);

        if (msec >= 0) {
            save_XYZ.outputFile += "_" + Integer.toString(msec) + project.timeUnits;
        }

        return true;

    }

    public void saveXYZDimensionsInit() {

        if ((save_XYZ == null) || (!save_XYZ.autoDimensions)) {
            return;
        }

        save_XYZ.xdim = "XYZ positions";

        if ((name == null) || (name.length() == 0)) {
            save_XYZ.ydim = project.spaceUnits;
        } else {
            save_XYZ.ydim = name + " (" + project.spaceUnits + ")";
        }

    }

    public boolean saveXYZ(int msec) {

        double f;
        String ostr;

        int dataPoints = 1;

        if (!saveXYZ || (vesicles == null)) {
            return false;
        }

        saveXYZFileName(msec);

        if (!save_XYZ.init(name, coordinates(), msec, dataPoints)) {
            return false;
        }

        for (DiffusantVesicle v : vesicles) {

            if (v == null) {
                continue;
            }

            if (!v.insideGeometry) {
                continue;
            }

            f = v.fluorescence;
            //f = vesicles[i].fluorescence * vesicles[i].voxel.PSFd;
            //f = vesicles[i].localD;

            //ostr = Double.toString(vesicle[i].x) + '\t' + Double.toString(vesicle[i].y) + '\t' + Double.toString(vesicle[i].z) + '\t' + Double.toString(f);
            ostr = Double.toString(v.x) + '\t' + Double.toString(v.y) + '\t' + Double.toString(v.z);

            save_XYZ.writeString(ostr);

        }

        save_XYZ.finish(name, coordinates(), msec);

        //Master.log("Saved vesicles positions at " + msec + " " + project.timeUnits);

        return false;

    }

    public boolean saveMSDFileName() {

        if (save_MSD == null) {
            return false;
        }

        save_MSD.fileName("MSD", name);

        return true;

    }

    public void saveMSDDimensions() {

        if ((save_MSD == null) || (!save_MSD.autoDimensions)) {
            return;
        }

        save_MSD.xdim = "ms";

        if ((name == null) || (name.length() == 0)) {
            save_MSD.ydim = project.spaceUnits;
        } else {
            save_MSD.ydim = name + " (" + project.spaceUnits + "^2)";
        }

    }

    public boolean checkVesicleNum(int i) {
        return ((vesicles != null) && (i >= 0) && (i < vesicles.length));
    }

    public void initVesicles() {

        double spaceVolume = project.monteCarlo.spaceVolume();

        if ((setVolumeFraction > 0) && (spaceVolume > 0)) {
            numVesicles = numVesiclesPossible(spaceVolume, setMeanRadius, setVolumeFraction);
        } else if ((setDensity > 0) && (spaceVolume > 0)) {
            numVesicles = (int) (setDensity * spaceVolume);
        }

        Master.log("numVesicles = " + numVesicles);

        if (numVesicles <= 0) {
            Master.exit("DiffusantVesicles: initVesicles: bad number of vesicles: " + numVesicles);
        }

        vesicles = new DiffusantVesicle[numVesicles];

        for (int i = 0; i < vesicles.length; i++) {
            vesicles[i] = new DiffusantVesicle(project, "ready", setMeanRadius, D, Double.NaN, Double.NaN, Double.NaN);
        }

    }

    public void shuffleVesicles() {

        int r, size = vesicles.length;

        DiffusantVesicle temp;

        for (int i = 0; i < size; i++) {

            r = (int) (mt.nextDouble() * size);

            if ((i != r) && (r < size)) { // swap
                temp = vesicles[i];
                vesicles[i] = vesicles[r];
                vesicles[r] = temp;
            }

        }

    }

    public void initStartLocation() {

        for (DiffusantVesicle v : vesicles) {
            v.initStartLocation();
        }

    }

    public double meanSquareDisplacement() {

        double sumSD = 0.0, count = 0.0;

        for (DiffusantVesicle v : vesicles) {
            if (v.mobile && v.insideGeometry) {
                sumSD += v.squareDisplacement();
                count += 1.0;
            }
        }

        return sumSD / count;

    }

    public double meanDisplacementAfterCollision() {

        double sumSD = 0.0, count = 0.0;

        for (DiffusantVesicle v : vesicles) {
            if (v.mobile && v.insideGeometry && (v.sqrDisplacement > 0)) {
                sumSD += v.sqrDisplacement;
                count += 1.0;
            }
        }

        return Math.sqrt(sumSD / count);

    }

    public boolean initImmobileVesicles(Coordinates keepClear) {

        long immobile = Math.round(setImmobilePercent * numVesicles);

        mobileVesicles = 0;
        immobileVesicles = 0;

        for (DiffusantVesicle v : vesicles) {

            if (v == null) {
                continue;
            }

            v.mobile = true;

            if (v.insideGeometry) {
                mobileVesicles++;
            }

        }

        if (immobile == 0) {
            return false; // nothing to do
        }

        shuffleVesicles();

        for (DiffusantVesicle v : vesicles) {

            if (v == null) {
                continue;
            }

            if (!v.insideGeometry) {
                continue;
            }

            if ((keepClear != null) && keepClear.isInside(v.x, v.y, v.z)) {
                continue;
            }

            v.mobile = false;
            immobileVesicles++;
            mobileVesicles--;

            if (immobileVesicles >= immobile) {
                break; // finished
            }

        }

        if (immobileVesicles != immobile) {
            Master.exit("DiffusantVesicles: initImmobileVesicles: failed to initialize immobile vesicles: " + (immobile - immobileVesicles));
        }

        immobilePercent = 1.0 * immobileVesicles / (immobileVesicles + mobileVesicles) ;

        return false;

    }

    int countVesicles(String vesicleType) {

        int count = 0;

        if (vesicleType.equalsIgnoreCase("mobile")) {

            for (DiffusantVesicle v : vesicles) {
                if (v.mobile) {
                    count++;
                }
            }

        } else if (vesicleType.equalsIgnoreCase("immobile")) {

            for (DiffusantVesicle v : vesicles) {
                if (!v.mobile) {
                    count++;
                }
            }

        } else { // "ready" or "docked" or "reserve"

            for (DiffusantVesicle v : vesicles) {
                if (v.name.equalsIgnoreCase(vesicleType)) {
                    count++;
                }
            }

        }

        return count;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (save_XYZ != null) {
            save_XYZ.addUser(pv);
        }

        if (save_MSD != null) {
            save_MSD.addUser(pv);
        }

        if (colorReady != null) {
            colorReady.addUser(pv);
        }

        if (colorImmobile != null) {
            colorImmobile.addUser(pv);
        }

        if (colorConnected != null) {
            colorConnected.addUser(pv);
        }

        if (PV_includeVesicleArray && (vesicles != null)) {
            for (DiffusantVesicle v : vesicles) {
                if (v != null) {
                    v.addUser(pv);
                }
            }
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {
 
        if (!super.createVector(false)) {
            return false;
        }

        if (save_XYZ != null) {
            addBlankParam();
            save_XYZ.createVector(true);
            addVector(save_XYZ.getVector());
            save_XYZ.addUser(this);
        }

        if (save_MSD != null) {
            addBlankParam();
            save_MSD.createVector(true);
            addVector(save_MSD.getVector());
            save_MSD.addUser(this);
        }

        if (colorReady != null) {
            addBlankParam();
            colorReady.createVector(true);
            addVector(colorReady.getVector());
            colorReady.addUser(this);
        }

        if (colorImmobile != null) {
            addBlankParam();
            colorImmobile.createVector(true);
            addVector(colorImmobile.getVector());
            colorImmobile.addUser(this);
        }

        if (colorConnected != null) {
            addBlankParam();
            colorConnected.createVector(true);
            addVector(colorConnected.getVector());
            colorConnected.addUser(this);
        }

        if (PV_includeVesicleArray && (vesicles != null)) {
            for (DiffusantVesicle v : vesicles) {
                if (v != null) {
                    addBlankParam();
                    v.createVector(true);
                    addVector(v.getVector());
                    v.addUser(this);
                }
            }
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void clearVector() {

        super.clearVector();

        if (PV_includeVesicleArray && (vesicles != null)) {
            for (DiffusantVesicle v : vesicles) {
                if (v != null) {
                    v.clearVector();
                }
            }
        }

    }

    @Override
    public void updateVector(ParamObject[] ov) {
        
        super.updateVector(ov);

        if (save_XYZ != null) {
            save_XYZ.updateVector(ov);
        }

        if (save_MSD != null) {
            save_MSD.updateVector(ov);
        }

        if (colorReady != null) {
            colorReady.updateVector(ov);
        }

        if (colorImmobile != null) {
            colorImmobile.updateVector(ov);
        }

        if (colorConnected != null) {
            colorConnected.updateVector(ov);
        }

        if (PV_includeVesicleArray && (vesicles != null)) {
            for (DiffusantVesicle v : vesicles) {
                if (v != null) {
                    v.updateVector(ov);
                }
            }
        }

    }

    public boolean setVesicles(String varName, double value) {

        boolean ok, atLeastOne = false;

        if (vesicles == null) {
            return false;
        }

        for (DiffusantVesicle v : vesicles) {
            if (v != null) {
                ok = v.set(varName, value);
                if (ok) {
                    atLeastOne = true;
                }
            }
        }

        vesicleStats();

        return atLeastOne;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantVesicles)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("setMeanRadius")) {
            if (v <= 0) {
                return false;
            }
            return setVesicles("radius",v);
        }
        if (n.equalsIgnoreCase("D")) {
            if (v < 0) {
                return false;
            }
            D = v;
            return setVesicles("D",v);
        }
        if (n.equalsIgnoreCase("meanStep3")) {
            if (v < 0) {
                return false;
            }
            return setVesicles("step3",v);
        }
        if (n.equalsIgnoreCase("meanExtraX")) {
            return setVesicles("extraX",v);
        }
        if (n.equalsIgnoreCase("meanExtraY")) {
            return setVesicles("extraY",v);
        }
        if (n.equalsIgnoreCase("meanExtraZ")) {
            return setVesicles("extraZ",v);
        }
        if (n.equalsIgnoreCase("meanDockRefractoryPeriod")) {
            if (v < 0) {
                return false;
            }
            return setVesicles("dockRefractoryPeriod",v);
        }
        if (n.equalsIgnoreCase("setDensity")) {
            if (v < 0) {
                return false;
            }
            setDensity = v;
            return true;
        }
        if (n.equalsIgnoreCase("setVolumeFraction")) {
            if (v < 0) {
                return false;
            }
            setVolumeFraction = v;
            return true;
        }
        if (n.equalsIgnoreCase("setImmobilePercent")) {
            if (v < 0) {
                return false;
            }
            setImmobilePercent = v;
            return true;
        }
        if (n.equalsIgnoreCase("savePositions")) {
            saveXYZ = (v == 1);
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

        if (!(o.paramVector instanceof DiffusantVesicles)) {
            return false;
        }
        
        String n = o.getName();

        return super.setMyParams(o, s);
    }

}
