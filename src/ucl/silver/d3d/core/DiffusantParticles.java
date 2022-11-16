package ucl.silver.d3d.core;

//import ucl.silver.d3d.utils.MersenneTwisterFast;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import java.awt.Color;
import java.util.Arrays;
import ucl.silver.d3d.utils.Utility;

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
public class DiffusantParticles extends Diffusant {

    public double radiusMin; // um
    public double radiusMax; // um
    public double radiusMean; // um
    public double radiusStdv; // um
    public double setRadiusMean; // um
    public double setRadiusStdv; // um
    public double setRadiusMin = 0; // um
    public double setRadiusMax = Double.POSITIVE_INFINITY; // um
    
    public double setRadiusChiDF = 0;
    public double setRadiusChiBeta = 0;

    public double volumeMean; // um^3
    public double volumeTotal; // um^3

    public double meanStep3; // um

    public double maxX0, maxY0, maxZ0; // um

    public transient DiffusantParticle[] particles; // the array of particles
    public transient double xyz[][] = null;

    public double setDensity; // particles / um^3
    public double setVolumeFraction;

    public double density; // particles / um^3
    public double volumeFraction;

    public int numParticles; // total number of particles in array
    public int setNumParticles;

    public int mobileParticles, immobileParticles;
    public double mobileVolumeFraction, immobileVolumeFraction;
    public double mobilePercent, immobilePercent; // e.g. 0.25
    public double setImmobilePercent = 0; // set immobile fraction

    public ColorD3D colorReady = new ColorD3D("colorReady", Color.white);
    public ColorD3D colorImmobile = new ColorD3D("colorImmobile", Color.gray);
    public ColorD3D colorConnected = new ColorD3D("colorConnected", new Color(102, 102, 0));

    public boolean PV_includeParticleArray = false;

    public boolean saveXYZ = false; // save particle xyz position
    public boolean saveMSD = false; // save mean square distance
    public boolean saveDisplacementAfterCollisions = false;

    public Save save_XYZ = null;
    public Save save_MSD = null;
    
    NumberFormat formatter = new DecimalFormat("0.####");
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("radiusMin")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("radiusMax")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("radiusMean")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("setRadiusMean")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("setRadiusMin")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("setRadiusMax")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("radiusStdv")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("setRadiusStdv")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("setRadiusGammaBeta")) {
            return "/" + project.spaceUnits;
        }
        if (name.equalsIgnoreCase("volumeMean")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("volumeTotal")) {
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
            return "particles/" + project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("setDensity")) {
            return "particles/" + project.spaceUnits + "^3";
        }
        return super.units(name);
    }

    public DiffusantParticles(Project p, String NAME, double InitialConcentration, double DiffusionConstant, CoordinatesVoxels c,
            double radius_mean, double radius_stdv) {

        super(p, NAME, InitialConcentration, DiffusionConstant, c);
        
        if (radius_mean < 0) {
            Master.exit("bad radius_mean: " + radius_mean);
        }
        
        if (radius_stdv < 0) {
            Master.exit("bad radius_stdv: " + radius_stdv);
        }

        setRadiusMean = radius_mean;
        setRadiusStdv = radius_stdv;

        createVector(true);

    }
    
    public boolean setChi(double df, double beta) {

        if (df < 1) {
            Master.exit("bad chi df: " + df);
        }
        
        if (beta < 0) {
            Master.exit("bad chi beta: " + beta);
        }
        
        setRadiusChiDF = df;
        setRadiusChiBeta = beta;

        setRadiusMean = Utility.chiMean(df, beta);
        setRadiusStdv = Utility.chiStdv(df, beta);

        createVector(true);
        
        return false;

    }

    @Override
    public void init() {

        super.init();

        if (particles != null) {
            for (DiffusantParticle p : particles) {
                if (p != null) {
                    p.init();
                }
            }
        }

        particleStats();

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
            save_MSD.updateVectors();

        }

    }

    public void particleStats() {

        double volume;
        double spaceVolume;

        if (project.monteCarlo == null) {
            return;
        }

        if (particles == null) {
            return;
        }

        spaceVolume = project.monteCarlo.spaceVolume();

        radiusMin = Double.POSITIVE_INFINITY;
        radiusMax = 0;
        radiusMean = 0;
        radiusStdv = 0;

        volumeMean = 0;
        volumeTotal = 0;

        D = 0;
        meanStep3 = 0;

        maxX0 = 0;
        maxY0 = 0;
        maxZ0 = 0;

        mobileParticles = 0;
        immobileParticles = 0;
        mobilePercent = 0;
        immobilePercent = 0;

        density = 0;
        volumeFraction = 0;
        mobileVolumeFraction = 0;
        immobileVolumeFraction = 0;

        for (DiffusantParticle p : particles) {

            if (p == null) {
                continue;
            }

            if (!p.insideGeometry) {
                continue;
            }

            volume = p.getVolume();

            radiusMin = Math.min(p.radius, radiusMin);
            radiusMax = Math.max(p.radius, radiusMax);

            radiusMean += p.radius;
            radiusStdv += p.radius * p.radius;
            volumeMean += volume;
            volumeTotal += volume;

            D += p.D;
            meanStep3 += p.step3;

            maxX0 = Math.max(maxX0, Math.abs(p.x0));
            maxY0 = Math.max(maxY0, Math.abs(p.y0));
            maxZ0 = Math.max(maxZ0, Math.abs(p.z0));

            if (p.mobile) {
                mobileParticles++;
            } else {
                immobileParticles++;
            }

        }

        numParticles = immobileParticles + mobileParticles;
        radiusStdv = Math.sqrt((radiusStdv - radiusMean * radiusMean / numParticles) / (numParticles - 1.0));
        radiusMean /= numParticles;
        volumeMean /= numParticles;
        D /= numParticles;
        meanStep3 /= numParticles;

        immobilePercent = immobileParticles / (1.0 * numParticles);
        mobilePercent = mobileParticles / (1.0 * numParticles);

        density = (1.0 * numParticles) / spaceVolume;
        volumeFraction = volumeTotal / spaceVolume;
        mobileVolumeFraction = mobilePercent * volumeFraction;
        immobileVolumeFraction = immobilePercent * volumeFraction;

        setParamObject("radiusMin", radiusMin);
        setParamObject("radiusMax", radiusMax);
        setParamObject("radiusMean", radiusMean);
        setParamObject("radiusStdv", radiusStdv);

        setParamObject("volumeMean", volumeMean);
        setParamObject("volumeTotal", volumeTotal);

        setParamObject("D", D);
        setParamObject("meanStep3", meanStep3);

        setParamObject("maxX0", maxX0);
        setParamObject("maxY0", maxY0);
        setParamObject("maxZ0", maxZ0);

        setParamObject("numParticles", numParticles);
        setParamObject("mobileParticles", mobileParticles);
        setParamObject("immobileParticles", immobileParticles);
        setParamObject("mobilePercent", mobilePercent);
        setParamObject("immobilePercent", immobilePercent);

        setParamObject("density", density);
        setParamObject("volumeFraction", volumeFraction);
        setParamObject("mobileVolumeFraction", mobileVolumeFraction);
        setParamObject("immobileVolumeFraction", immobileVolumeFraction);

    }

    public void computeDensity() {

        double sumVolume = 0.0, count = 0.0;

        double spaceVolume = project.monteCarlo.spaceVolume();

        for (DiffusantParticle p : particles) {

            if ((p != null) && p.insideGeometry) {
                count++;
                sumVolume += p.getVolume();
            }

        }

        density = count / spaceVolume;
        volumeFraction = sumVolume / spaceVolume;

        setParamObject("density", density);
        setParamObject("volumeFraction", volumeFraction);

    }

    public static int numParticlesPossible(double volumeGeometry, double particleRadius, double vFraction) {
        double particleVolume = 4.0 * Math.PI * Math.pow(particleRadius, 3) / 3.0;
        return (int) (vFraction * volumeGeometry / particleVolume); // all particles are the same size
    }

    @Override
    public void saveFileName() {

        if (save_MSD != null) {
            save_MSD.fileName(name, "MSD");
        }

        if (save_XYZ != null) {
            save_XYZ.fileName(name, "XYZ");
        }

        super.saveFileName();

    }

    @Override
    public void saveDimensions() {

        if ((save_MSD != null) && save_MSD.autoDimensions) {

            save_MSD.xdim = "ms";

            if ((name == null) || (name.length() == 0)) {
                save_MSD.ydim = project.spaceUnits;
            } else {
                save_MSD.ydim = name + " (" + project.spaceUnits + "^2)";
            }

        }

        if ((save_XYZ != null) && save_XYZ.autoDimensions) {

            save_XYZ.xdim = "XYZ positions";

            if ((name == null) || (name.length() == 0)) {
                save_XYZ.ydim = project.spaceUnits;
            } else {
                save_XYZ.ydim = name + " (" + project.spaceUnits + ")";
            }

        }

        super.saveDimensions();

    }

    @Override
    public boolean saveInit() {

        int dataPoints = 1;

        if (save_MSD != null) {
            if (!save_MSD.init(name, coordinates(), -1, dataPoints)) {
                return false;
            }
        }

        return super.saveInit();

    }

    @Override
    public boolean saveFinish() {

        if (save_MSD != null) {
            save_MSD.finish(name, coordinates(), -1);
        }

        return super.saveFinish();

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

        return super.save();

    }

    public void saveXYZFileName(int msec) {

        if (save_XYZ != null) {

            //save_XYZ.fileName("XYZF", name);
            save_XYZ.fileName("XYZ", name);

            if (msec >= 0) {
                save_XYZ.outputFile += "_" + Integer.toString(msec) + project.timeUnits;
            }

        }

    }

    public boolean saveXYZ(int msec) {

        double f;
        String ostr;

        int dataPoints = 1;

        if (!saveXYZ || (particles == null)) {
            return false;
        }

        saveXYZFileName(msec);

        if (!save_XYZ.init(name, coordinates(), msec, dataPoints)) {
            return false;
        }

        for (DiffusantParticle p : particles) {

            if (p == null) {
                continue;
            }

            if (!p.insideGeometry) {
                continue;
            }

            f = p.fluorescence;
            //f = particles[j].fluorescence * particles[j].voxel.PSFd;
            //f = particles[j].localD;

            //ostr = Double.toString(v.x) + '\t' + Double.toString(v.y) + '\t' + Double.toString(v.z) + '\t' + Double.toString(f);
            ostr = Double.toString(p.x) + '\t' + Double.toString(p.y) + '\t' + Double.toString(p.z);

            save_XYZ.writeString(ostr);

        }

        save_XYZ.finish(name, coordinates(), msec);

        //Master.log("Saved particle positions at " + msec + " " + project.timeUnits);
        return false;

    }

    public boolean checkParticleNum(int i) {
        return ((particles != null) && (i >= 0) && (i < particles.length));
    }
    
    public double[] getParameter(String select){
        
        if ((particles == null) || (particles.length == 0)) {
            return null;
        }
        
        double[] parray = new double[particles.length];
        
        for (int i = 0; i < particles.length; i++) {
            if (particles[i] == null) {
                parray[i] = Double.NaN;
                continue;
            }
            switch (select) {
                case "radius":
                    parray[i] = particles[i].radius;
                    break;
                case "diameter":
                    parray[i] = particles[i].radius * 2;
                    break;
                case "volume":
                    parray[i] = particles[i].getVolume();
                    break;
                case "D":
                    parray[i] = particles[i].D;
                    break;
                case "x0":
                    parray[i] = particles[i].x0;
                    break; 
                case "y0":
                    parray[i] = particles[i].y0;
                    break;
                case "z0":
                    parray[i] = particles[i].z0;
                    break;
                case "x":
                    parray[i] = particles[i].x;
                    break; 
                case "y":
                    parray[i] = particles[i].y;
                    break;
                case "z":
                    parray[i] = particles[i].z;
                    break;
                case "fluorescence":
                    parray[i] = particles[i].fluorescence;
                    break;
                default:
                    parray[i] = Double.NaN;
            }
            
        }

        return parray;

    }

    public void initParticles() {

        double r, vfraction, v, vsum = 0;
        double[] rtemp = null;
        
        int rtrials = 100; // trials for computing random radius

        double spaceVolume = project.monteCarlo.spaceVolume();

        if ((setVolumeFraction > 0) && (spaceVolume > 0)) {

            numParticles = numParticlesPossible(spaceVolume, setRadiusMean, setVolumeFraction);
            
            if (setRadiusChiDF >= 1) {
                
                rtemp = new double[numParticles * 2];

                for (int i = 0; i < rtemp.length; i++) {
                    r = 0;
                    for (int j = 0; j < rtrials; j++) {
                        r = Master.randomChi(setRadiusChiDF, setRadiusChiBeta);
                        if ((r > setRadiusMin) && (r < setRadiusMax)) {
                            break;
                        }
                    }
                    if ((r > setRadiusMin) && (r < setRadiusMax)) {
                        rtemp[i] = r;
                    } else {
                        Master.exit("bad particle radius: " + r);
                    }
                    v = 4.0 * Math.PI * Math.pow(r, 3) / 3.0;
                    vsum += v;
                    vfraction = vsum / spaceVolume;
                    if (vfraction >= setVolumeFraction) {
                        numParticles = i;
                        break;
                    }
                }

            } else if (setRadiusStdv > 0) {

                rtemp = new double[numParticles * 2];

                for (int i = 0; i < rtemp.length; i++) {
                    r = 0;
                    for (int j = 0; j < rtrials; j++) {
                        r = setRadiusMean + Master.randomGauss() * setRadiusStdv;
                        if ((r > setRadiusMin) && (r < setRadiusMax)) {
                            break;
                        }
                    }
                    if ((r > setRadiusMin) && (r < setRadiusMax)) {
                        rtemp[i] = r;
                    } else {
                        Master.exit("bad particle radius: " + r);
                    }
                    v = 4.0 * Math.PI * Math.pow(r, 3) / 3.0;
                    vsum += v;
                    vfraction = vsum / spaceVolume;
                    if (vfraction >= setVolumeFraction) {
                        numParticles = i;
                        break;
                    }
                }

            }

        } else if ((setDensity > 0) && (spaceVolume > 0)) {

            numParticles = (int) (setDensity * spaceVolume);

        } else if ((setNumParticles > 0) && (spaceVolume > 0)) {

            numParticles = setNumParticles;

        }

        Master.log("numParticles = " + numParticles);

        if (numParticles <= 0) {
            Master.exit("DiffusantParticles: initParticles: bad number of particles: " + numParticles);
        }

        particles = new DiffusantParticle[numParticles];

        r = setRadiusMean;

        for (int i = 0; i < particles.length; i++) {
            
            if (setRadiusChiDF >= 1) {
                if ((rtemp != null) && (i < rtemp.length)) {
                    r = rtemp[i];
                } else {
                    for (int j = 0; j < rtrials; j++) {
                        r = Master.randomChi(setRadiusChiDF, setRadiusChiBeta);
                        if ((r > setRadiusMin) && (r < setRadiusMax)) {
                            break;
                        }
                    }
                }
            } else if (setRadiusStdv > 0) {
                if ((rtemp != null) && (i < rtemp.length)) {
                    r = rtemp[i];
                } else {
                    for (int j = 0; j < rtrials; j++) {
                        r = setRadiusMean + Master.randomGauss() * setRadiusStdv;
                        if ((r > setRadiusMin) && (r < setRadiusMax)) {
                            break;
                        }
                    }
                }
            }
            if ((r > setRadiusMin) && (r < setRadiusMax)) {
                particles[i] = new DiffusantParticle(project, "ready", r, D, Double.NaN, Double.NaN, Double.NaN);
                //Master.log("" + r);
            } else {
                Master.exit("bad particle radius: " + r);
            }
        }

    }

    public void shuffleParticles() {

        int r, size = particles.length;

        DiffusantParticle temp;

        for (int i = 0; i < size; i++) {

            r = (int) (Master.mt.nextDouble() * size);

            if ((i != r) && (r < size)) { // swap
                temp = particles[i];
                particles[i] = particles[r];
                particles[r] = temp;
            }

        }

    }
    
    public boolean sortParticlesByRadius() {

        double[] r_sorted;

        if (particles == null) {
            return true;
        }

        if (particles.length <= 1) {
            return true;
        }

        r_sorted = new double[particles.length];

        for (int i = 0; i < particles.length; i++) {
            r_sorted[i] = particles[i].radius;
        }
        
        //Master.log("before sorting diameter = " + 2000 * Utility.average(r_sorted) + " ± " + 2000 * Utility.stdev(r_sorted));

        Arrays.sort(r_sorted);
        
        //Master.log("after sorting diameter = " + 2000 * Utility.average(r_sorted) + " ± " + 2000 * Utility.stdev(r_sorted));
        
        for (int i = r_sorted.length - 1, j = 0; i >= 0; i--, j++) {
            particles[j].radius = r_sorted[i];
        }
        
        Master.log("sorted particles by size");

        return false;

    }
    
    public boolean sortParticlesByRadiusSlow() {
        
        int jmax;
        DiffusantParticle rmax;

        if (particles == null) {
            return true;
        }

        if (particles.length <= 1) {
            return true;
        }
        
        DiffusantParticle[] v2 = new DiffusantParticle[particles.length];
        
        for (int i = 0; i < v2.length; i++) {
            
            rmax = null;
            jmax = -1;
            
            for (int j = 0; j < particles.length; j++) {
                
                if (particles[j] == null) {
                    continue;
                }
                
                if (rmax == null) {
                    rmax = particles[j];
                    jmax = j;
                } else if (particles[j].radius > rmax.radius) {
                    rmax = particles[j];
                    jmax = j;
                }
                
            }
            
            if (rmax == null) {
                Master.exit("failed to find maximum radius for particle i = " + i);
            }
            
            v2[i] = rmax;
            particles[jmax] = null;

        }
        
        System.arraycopy(v2, 0, particles, 0, particles.length);
        
        Master.log("sorted particles by size");

        return false;

    }

    public void initStartLocation() {

        for (DiffusantParticle p : particles) {
            p.initStartLocation();
        }

    }

    public double meanSquareDisplacement() {

        double sumSD = 0.0, count = 0.0;

        for (DiffusantParticle p : particles) {
            if (p.mobile && p.insideGeometry) {
                sumSD += p.squareDisplacement();
                count += 1.0;
            }
        }

        return sumSD / count;

    }

    public double meanDisplacementAfterCollision() {

        double sumSD = 0.0, count = 0.0;

        for (DiffusantParticle p : particles) {
            if (p.mobile && p.insideGeometry && (p.sqrDisplacement > 0)) {
                sumSD += p.sqrDisplacement;
                count += 1.0;
            }
        }

        return Math.sqrt(sumSD / count);

    }

    public boolean initImmobileParticles(Coordinates keepClear) {

        long immobile = Math.round(setImmobilePercent * numParticles);

        mobileParticles = 0;
        immobileParticles = 0;

        for (DiffusantParticle p : particles) {

            if (p == null) {
                continue;
            }

            p.mobile = true;

            if (p.insideGeometry) {
                mobileParticles++;
            }

        }

        if (immobile == 0) {
            return false; // nothing to do
        }

        shuffleParticles();

        for (DiffusantParticle p : particles) {

            if (p == null) {
                continue;
            }

            if (!p.insideGeometry) {
                continue;
            }

            if ((keepClear != null) && keepClear.isInside(p.x, p.y, p.z)) {
                continue;
            }

            p.mobile = false;
            immobileParticles++;
            mobileParticles--;

            if (immobileParticles >= immobile) {
                break; // finished
            }

        }

        if (immobileParticles != immobile) {
            Master.exit("DiffusantParticles: initImmobileParticles: failed to initialize immobile particles: " + (immobile - immobileParticles));
        }

        immobilePercent = 1.0 * immobileParticles / (immobileParticles + mobileParticles);

        return false;

    }

    int countParticles(String particleType) {

        int count = 0;

        if (particleType.equalsIgnoreCase("mobile")) {

            for (DiffusantParticle p : particles) {
                if (p.mobile) {
                    count++;
                }
            }

        } else if (particleType.equalsIgnoreCase("immobile")) {

            for (DiffusantParticle p : particles) {
                if (!p.mobile) {
                    count++;
                }
            }

        } else { // e.g. "ready" or "docked" or "reserve"

            for (DiffusantParticle p : particles) {
                if (p.name.equalsIgnoreCase(particleType)) {
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

        if (PV_includeParticleArray && (particles != null)) {
            for (DiffusantParticle p : particles) {
                if (p != null) {
                    p.addUser(pv);
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

        if (PV_includeParticleArray && (particles != null)) {
            for (DiffusantParticle p : particles) {
                if (p != null) {
                    addBlankParam();
                    p.createVector(true);
                    addVector(p.getVector());
                    p.addUser(this);
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

        if (PV_includeParticleArray && (particles != null)) {
            for (DiffusantParticle p : particles) {
                if (p != null) {
                    p.clearVector();
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

        if (PV_includeParticleArray && (particles != null)) {
            for (DiffusantParticle p : particles) {
                if (p != null) {
                    p.updateVector(ov);
                }
            }
        }

    }

    public boolean setParticles(String varName, double value) {

        boolean ok, atLeastOne = false;

        if (particles == null) {
            //return false;
            return true; // need true to get batches to run
        }

        for (DiffusantParticle p : particles) {
            if (p != null) {
                ok = p.set(varName, value);
                if (ok) {
                    atLeastOne = true;
                }
            }
        }

        particleStats();

        return atLeastOne;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantParticles)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("D")) {
            if (v < 0) {
                return false;
            }
            D = v;
            return setParticles("D", v);
        }
        if (n.equalsIgnoreCase("meanStep3")) {
            if (v < 0) {
                return false;
            }
            return setParticles("step3", v);
        }
        if (n.equalsIgnoreCase("meanExtraX")) {
            return setParticles("extraX", v);
        }
        if (n.equalsIgnoreCase("meanExtraY")) {
            return setParticles("extraY", v);
        }
        if (n.equalsIgnoreCase("meanExtraZ")) {
            return setParticles("extraZ", v);
        }
        if (n.equalsIgnoreCase("meanDockRefractoryPeriod")) {
            if (v < 0) {
                return false;
            }
            return setParticles("dockRefractoryPeriod", v);
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

        if (!(o.paramVector instanceof DiffusantParticles)) {
            return false;
        }

        String n = o.getName();

        return super.setMyParams(o, s);
    }

}
