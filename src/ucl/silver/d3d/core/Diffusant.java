package ucl.silver.d3d.core;

import java.awt.Color;
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
public class Diffusant extends ParamVector {
    
    public double C0; // default initial concentration (mM)
    public double D; // default diffusion coefficient (um^2/ms)
    public double h; // h = dt * D / dx^2
    public double charge; // for conversion between molar flux and amperes
    
    public boolean reaction = false;
    public boolean saveConc = false;
    
    public double concentration[][][] = null; // where concentration values are saved

    public boolean useGeometryCoordinates = true;

    //public String displaySelect = "C0"; // "C0" or "D"
    public double displayMinC = 0; // mM
    public double displayMaxC = 1; // mM

    public ColorD3D color = null;
    public CoordinatesVoxels coordinates = null; // coordinates of diffusant
    public Save save = null;
    public PSF psf = null; // excitation reaction
    public PulseTimer pulseTimer = null; // excitation reaction
    
    public double saveValue = Double.NaN; // the value that saved to data file during simulation

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("C0")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("D")) {
            return project.spaceUnits + "^2/" + project.timeUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("h")) {
            return false;
        }
        if (name.equalsIgnoreCase("charge")) {
            return false; // not editable since Sources need to be updated if this changes
        }
        if (name.equalsIgnoreCase("reaction")) {
            return false; // parameter is set by constructors
        }
        if (name.equalsIgnoreCase("useGeometryCoordinates")) {
            return false;
        }
        return super.canEdit(name);
    }

    public Diffusant(Project p, String NAME, double InitialConcentration, double DiffusionConstant, CoordinatesVoxels c) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Diffusant";
        }
        
        charge = Utility.getCharge(NAME);

        color = new ColorD3D(name + "_color", "Heat");

        C0 = InitialConcentration;
        D = DiffusionConstant;

        save = new Save(p);
        
        if (c == null) {
            useGeometryCoordinates = true;
            coordinates = null;
        } else {
            useGeometryCoordinates = false;
            coordinates = new CoordinatesVoxels(p);
            coordinates.matchVoxels(c);
        }

        createVector(true);

    }

    public void matchCoordinates(CoordinatesVoxels c) {
        if ((!useGeometryCoordinates) && (coordinates != null)) {
            coordinates.matchVoxels(c);
        }
    }

    public CoordinatesVoxels coordinates() {
        if ((useGeometryCoordinates) || (coordinates == null)) {
            return project.geometry;
        } else  {
            return coordinates;
        }
    }

    public void updateCoordinates() {
        // nothing to do
    }

    public void init() {

        if (save != null) {
            save.init();
            saveFileName();
            saveDimensions();
        }

        h = computeH(D);

        if (psf != null) {
            psf.init();
        }

    }

    public void initPulseTimer() {
        if (pulseTimer != null) {
            pulseTimer.initTimer();
        }
    }

    public double computeH(double d) {
        return project.dt * d / (Math.pow(project.dx, 2.0));
    }

    public boolean checkConcentrationXYZ(int xVoxel, int yVoxel, int zVoxel) {
        if (concentration == null) {
            return false;
        }
        if ((xVoxel < 0) || (xVoxel >= concentration.length)) {
            return false;
        }
        if ((yVoxel < 0) || (yVoxel >= concentration[0].length)) {
            return false;
        }
        if ((zVoxel < 0) || (zVoxel >= concentration[0][0].length)) {
            return false;
        }
        return true; // OK
    }

    public double displayValue(int xVoxel, int yVoxel, int zVoxel) {

        if (concentration != null) {
            if (checkConcentrationXYZ(xVoxel, yVoxel, zVoxel)) {
               return concentration[xVoxel][yVoxel][zVoxel]; // display saved concentration if it exists
            } else {
                return -1;
            }
        }

        if (psf != null) {
            return psf.getArrayValue(xVoxel, yVoxel, zVoxel); // or PSF if it exists
        }

        return C0;

    }

    public Color displayColor(int xVoxel, int yVoxel, int zVoxel, double concentration) {

        double psfValue;

        if (color == null) {
            return null;
        }

        if (!coordinates().isInside(xVoxel, yVoxel, zVoxel)) {
            return null;
        }

        if (psf != null) {

            psf.checkExists();

            psfValue = psf.getArrayValue(xVoxel - coordinates().xVoxel1, yVoxel - coordinates().yVoxel1, zVoxel - coordinates().zVoxel1);
            psfValue = Math.max(0, psfValue);

            return color.getColor(psfValue, 0, 1);

        }

        return color.getColor(concentration, displayMinC, displayMaxC);

    }

    public Color displayColor(int xVoxel, int yVoxel, int zVoxel) {

        double conc, psfValue;

        if (color == null) {
            return null;
        }

        if (!coordinates().isInside(xVoxel, yVoxel, zVoxel)) {
            return null;
        }

        if (psf != null) {

            psf.checkExists();

            psfValue = psf.getArrayValue(xVoxel - coordinates().xVoxel1, yVoxel - coordinates().yVoxel1, zVoxel - coordinates().zVoxel1);
            psfValue = Math.max(0, psfValue);

            return color.getColor(psfValue, 0, 1);

        }

        conc = getConcentration(xVoxel, yVoxel, zVoxel);

        return color.getColor(conc, displayMinC, displayMaxC);

    }

    public double getConcentration(int xVoxel, int yVoxel, int zVoxel) {

        if (concentration != null) {
            if (checkConcentrationXYZ(xVoxel, yVoxel, zVoxel)) {
                return concentration[xVoxel][yVoxel][zVoxel]; // use saved concentration if it exists
            } else {
                return -1;
            }
        }

        return C0;

    }

    public PSF getPSF() {
        return psf;
    }

    public void check() {
        if (psf != null) {
            psf.checkExists();
        }
    }

    public void react(RunFiniteDifference fd, int abnum) {
        // no reaction
    }

    public void saveConcentration(RunFiniteDifference fd, Geometry geometry, int dnum) {

        CoordinatesVoxels c = geometry;

        if (!saveConc || (fd.diffus == null) || (dnum >= fd.diffus.length) ) {
            return;
        }

        concentration = new double[geometry.xVoxels][geometry.yVoxels][geometry.zVoxels];

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                    if (geometry.isSpace(i, j, k)) {
                        concentration[i][j][k] = fd.diffus[dnum][geometry.space[i][j][k]];
                    } else {
                        concentration[i][j][k] = 0;
                    }
                }
            }
        }

    }

    public boolean savedConcExists() {
        return (concentration != null);
    }

    public boolean saveFileName(){

        if (save == null) {
            return false;
        }

        save.fileName(name, "");

        return true;

    }

    public void saveDimensions() {

        if ((save == null) || (!save.autoDimensions)) {
            return;
        }

        save.xdim = project.timeUnits;

        if ((name == null) || (name.length() == 0)) {
            save.ydim = project.concUnits;
        } else {
            save.ydim = name + " (" + project.concUnits + ")";
        }

    }

    public boolean saveInit() {

        int dataPoints = 1;

        if (save == null) {
            return false;
        }

        return save.init(name, coordinates(), -1, dataPoints);

    }

    public boolean save() { // should be called from simulation "run" function
        
        if (save == null) {
            return false;
        }
        
        return save.saveData(saveValue);
        
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
        
        if ((o.paramVector instanceof Pulse) && (pulseTimer != null) && (pulseTimer.pulses != null)) {
            for (Pulse p : pulseTimer.pulses) {
                if (o.paramVector.equals(p)) {
                    return (p.setMyParams(o, v));
                }
            }
            return false;
        }

        if (!(o.paramVector instanceof Diffusant)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("C0")) {
            if (v < 0) {
                return false;
            }
            C0 = v;
            return true;
        }
        if (n.equalsIgnoreCase("D")) {
            if (v < 0) {
                return false;
            }
            D = v;
            return true;
        }
        if (n.equalsIgnoreCase("charge")) {
            charge = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("displayMinC")) {
            if (v < 0) {
                return false;
            }
            displayMinC = v;
            return true;
        }
        if (n.equalsIgnoreCase("displayMaxC")) {
            if (v < 0) {
                return false;
            }
            displayMaxC = v;
            return true;
        }
        if (n.equalsIgnoreCase("reaction")) {
            reaction = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("saveConc")) {
            saveConc = (v == 1);
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

        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            if (coordinates.setMyParams(o, s)) {
                updateCoordinates();
                return true;
            }
            return false;
        }
        
        if ((o.paramVector instanceof Pulse) && (pulseTimer != null) && (pulseTimer.pulses != null)) {
            for (Pulse p : pulseTimer.pulses) {
                if (o.paramVector.equals(p)) {
                    return (p.setMyParams(o, s));
                }
            }
            return false;
        }

        if (!(o.paramVector instanceof Diffusant)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        //if (n.equalsIgnoreCase("displaySelect")) {
            //displaySelect = s;
            //return true;
        //}
        return super.setMyParams(o, s);
    }
}
