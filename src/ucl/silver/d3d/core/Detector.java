package ucl.silver.d3d.core;

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
public class Detector extends ParamVector {

    public String diffusantName = null;
    public int diffusantNum = -1; // diffusant array number in Project

    public boolean useGeometryCoordinates = true;

    private final Color defaultColor = new Color(102, 153, 0);

    public ColorD3D color = null;
    private CoordinatesVoxels coordinates = null;
    public Save save = null;
    public PSF psf = null;
    public PulseTimer pulseTimer = null;

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("useGeometryCoordinates")) {
            return false;
        }
        return true;
    }

    public Detector(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Detector";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;
        save = new Save(p);
        //save.save2TextFile = true;
        save.save2BinaryFile = true;
        //save.saveWhileComputing = true;
        save.saveWhileComputing = false;

        if (c == null) {
            useGeometryCoordinates = true;
            coordinates = null;
            save.dataPoints = p.geometry.voxels;
        } else {
            useGeometryCoordinates = false;
            coordinates = new CoordinatesVoxels(p);
            coordinates.matchVoxels(c);
            save.dataPoints = coordinates.voxels;
        }

        createVector(true);

    }

    public void matchCoordinates(CoordinatesVoxels c) {
        if ((!useGeometryCoordinates) && (coordinates != null)) {
            coordinates.matchVoxels(c);
        }
    }

    public void updateCoordinates() {

        if ((psf == null) || (psf.useGeometryCoordinates)) {
            return;
        }

        psf.coordinates().matchVoxels(coordinates);

    }
    
    public CoordinatesVoxels coordinates() {
        if ((useGeometryCoordinates) || (coordinates == null)) {
            return project.geometry;
        } else  {
            return coordinates;
        }
    }

    public void check() {

        //if (coordinates().voxels != coordinates().spaceVoxels) {
        //    log("Warning: " + name + " is located at non-space voxel(s)");
        //}

        if (psf != null) {
            psf.checkExists();
        }

    }

    public boolean setDiffusantNum(int arrayNum) {

        Diffusant d = project.getDiffusant(arrayNum);

        if (d != null) {
            diffusantNum = arrayNum;
            diffusantName = d.name;
            updateVectors();
            return true;
        } else {
            Master.log("Detector SetDiffusantNum error: failed to find diffusant #" + arrayNum);
            updateVectors();
            return false;
        }

    }

    public boolean setDiffusantName(String dName) {

        int i = project.getDiffusantNum(dName);

        if (i >= 0) {
            diffusantName = dName;
            diffusantNum = i;
            updateVectors();
            return true;
        } else {
            Master.log("Detector SetDiffusantName error: failed to find diffusant " + dName);
            updateVectors();
            return false;
        }

    }

    public void init() {

        if (save == null) {

            Master.exit("Encountered null Save class for detector " + name);

        } else {

            save.init();
            saveFileName();
            saveDimensions();

            if ((save.save2TextFile == false) && (save.save2BinaryFile == false)) {
                Master.log("Warning!!! Saving to external data file is disabled for detector " + name);
            }

        }

        setParamError("diffusantNum", null);

        diffusantName = project.diffusantName(diffusantNum);

        if (!project.checkDiffusantNum(diffusantNum)) {
            error("init", "diffusantNum", "out of range");
        }

        if (coordinates != null) {
            coordinates.update();
        }

        if (psf != null) {
            psf.init();
        }

        updateCoordinates();

    }

    public void initPulseTimer() {
        if (pulseTimer != null) {
            pulseTimer.initTimer();
        }
    }

    public void detect(RunFiniteDifference fd, Geometry geometry) {

        if (fd.diffus[0].length == 1) {
            save.saveData(fd.diffus[diffusantNum][0]);
        } else if (coordinates().indexRFD == null) {
            // code needs to be written
        } else {
            save.saveData(fd.diffus, diffusantNum, coordinates().indexRFD); // save all voxel concentrations
        }

    }

    public double displayValue(int xVoxel, int yVoxel, int zVoxel) {
        
        if (psf != null) {

            psf.checkExists();

            int xVoxel1 = coordinates().xVoxel1;
            int yVoxel1 = coordinates().yVoxel1;
            int zVoxel1 = coordinates().zVoxel1;

            double psfValue = psf.getArrayValue(xVoxel - xVoxel1, yVoxel - yVoxel1, zVoxel - zVoxel1);

            return Math.max(0, psfValue);

        }

        if (coordinates().isInside(xVoxel, yVoxel, zVoxel)) {
            return 0.5; // half color scale
        } else {
            return 0;
        }
        
    }

    public Color displayColor(int xVoxel, int yVoxel, int zVoxel, double min, double max) {

        double psfValue;

        if (color == null) {
            return null;
        }

        if (!coordinates().isInside(xVoxel, yVoxel, zVoxel)) {
            return null;
        }

        if (color.isColorScale() && psf != null) {

            psf.checkExists();

            psfValue = psf.getArrayValue(xVoxel - coordinates().xVoxel1, yVoxel - coordinates().yVoxel1, zVoxel - coordinates().zVoxel1);

            psfValue = Math.max(0, psfValue);

            return color.getColor(psfValue, 0, 1);

        }

        return color.getColor(((max - min) / 2.0), min, max);

    }

    public boolean saveFileName(){

        String dname = "";

        Diffusant[] diffusants = project.diffusants;

        if (save == null) {
            return false;
        }

        if ((diffusantNum >= 0) && (diffusantNum < diffusants.length)) {
            dname = diffusants[diffusantNum].name;
        }

        save.fileName(name, dname);

        return true;

    }

    public void saveDimensions() {

        if ((save == null) || (!save.autoDimensions)) {
            return;
        }

        save.xdim = project.timeUnits;

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            save.ydim = project.concUnits;
        } else {
            save.ydim = diffusantName + " (" + project.concUnits + ")";
        }

    }

    public boolean saveInit() {

        if (save == null) {
            return false;
        }

        int dataPoints = coordinates().voxels;

        return save.init(name, coordinates(), -1, dataPoints);

    }

    public boolean saveFinish() {

        if (save == null) {
            return false;
        }

        return save.finish(name, coordinates(), -1);

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (color != null) {
            color.addUser(pv);
        }

        if (coordinates != null) {
            coordinates.addUser(pv);
        }

        if (save != null) {
            save.addUser(pv);
        }

        if (psf != null) {
            psf.addUser(pv);
        }

        if (pulseTimer != null) {
            pulseTimer.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (color != null) {
            addBlankParam();
            color.createVector(true);
            addVector(color.getVector());
            color.addUser(this);
        }

        if (coordinates != null) {
            addBlankParam();
            coordinates.createVector(true);
            addVector(coordinates.getVector());
            coordinates.addUser(this);
        }

        if (save != null) {
            addBlankParam();
            save.createVector(true);
            addVector(save.getVector());
            save.addUser(this);
        }

        if (psf != null) {
            addBlankParam();
            psf.createVector(true);
            addVector(psf.getVector());
            psf.addUser(this);
        }

        if (pulseTimer != null) {
            addBlankParam();
            pulseTimer.createVector(true);
            addVector(pulseTimer.getVector());
            pulseTimer.addUser(this);
        }

        if (close){
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {
        
        super.updateVector(v);

        if (color != null) {
            color.updateVector(v);
        }

        if (coordinates != null) {
            coordinates.updateVector(v);
        }
        
        if (save != null) {
            save.updateVector(v);
        }

        if (psf != null) {
            psf.updateVector(v);
        }

        if (pulseTimer != null) {
            pulseTimer.updateVector(v);
        }
        
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            if (coordinates.setMyParams(o, v)) {
                updateCoordinates();
                return true;
            }
            return false;
        }

        if (!(o.paramVector instanceof Detector)) {
            return false;
        }
        
        String n = o.getName();
        
        if (n.equalsIgnoreCase("diffusantNum")) {
            if (v < 0) {
                return false;
            }
            return setDiffusantNum((int)v);
        }
        return super.setMyParams(o, v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {
        
        if ((o == null) || (s == null)) {
            return false;
        }

        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            if (coordinates.setMyParams(o, s)) {
                updateCoordinates();
                return true;
            }
            return false;
        }

        if (!(o.paramVector instanceof Detector)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        if (n.equalsIgnoreCase("diffusantName")) {
            return setDiffusantName(s);
        }
        return super.setMyParams(o, s);
    }
}
