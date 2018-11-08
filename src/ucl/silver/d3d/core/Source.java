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
public class Source extends ParamVector {

    public String diffusantName = null;
    public int diffusantNum = -1; // diffusant array number in Project

    //public double Qmax = Double.NaN; // mM/ms
    public double Ctotal = Double.NaN; // mM, M = moles/L
    private double Cvoxel = Double.NaN; // mM, for source with no pulse timer
    public double Ntotal = Double.NaN; // molecules
    public double Nvoxel = Double.NaN; // molecules

    //public String initSelect = "Ctotal";
    public String saveSelect = "mM/ms"; // "mM/ms" or "mM" or "pA"
    public int voxels, spaceVoxels;
    public boolean totalVoxelsForSpaceVoxelsOnly = true;

    public boolean clamp = false; // source "clamps" voxels to Ctotal
    public boolean noOverlapCoordinates = false;

    public double C2I_conversionFactor = Double.NaN;

    private final Color defaultColor = new Color(153, 0, 51);

    public ColorD3D color = null;
    public CoordinatesVoxels[] coordinates = null; // coordinates of sources
    public transient int[] indexRFD = null; // index counter for FD simulations
    public Save save = null;
    public PulseTimer pulseTimer = null; // timer for source

    @Override
    public String units(String name) {

        //if (name.equalsIgnoreCase("Qmax")) {
            //return project.concUnits + "/" + project.timeUnits;
        //}
        if (name.equalsIgnoreCase("Ctotal")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("Ntotal")) {
            return "#";
        }
        if (name.equalsIgnoreCase("Nvoxel")) {
            return "#";
        }
        if (name.equalsIgnoreCase("C2I_conversionFactor")) {
            return project.currentUnits + "/" + project.concUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        //if (name.equalsIgnoreCase("Qmax")) {
            //return false;
        //}
        if (name.equalsIgnoreCase("Ctotal")) {
            if (pulseTimer != null) {
                return false; // edit via pulse timer amplitude
            }
        }
        if (name.equalsIgnoreCase("Ntotal")) {
            return false;
        }
        if (name.equalsIgnoreCase("Nvoxel")) {
            return false;
        }
        if (name.equalsIgnoreCase("voxels")) {
            return false;
        }
        if (name.equalsIgnoreCase("spaceVoxels")) {
            return false;
        }
        if (name.equalsIgnoreCase("C2I_conversionFactor")) {
            return false;
        }
        if (name.equalsIgnoreCase("currentUnits")) {
            return false;
        }
        return super.canEdit(name);
    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double onset, double duration, double ctotal) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;

        pulseTimer = new PulseTimer(project, onset, duration, ctotal);
        save = new Save(p);

        Ctotal = ctotal;

        if (c == null) {
            coordinates = null;
        } else {
            coordinates = new CoordinatesVoxels[1];
            coordinates[0] = new CoordinatesVoxels(p);
            coordinates[0].matchVoxels(c);
        }

        init();

        createVector(true); // sets ParamVector

    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;

        pulseTimer = pt;
        save = new Save(p);

        if (c == null) {
            coordinates = null;
        } else {
            coordinates = new CoordinatesVoxels[1];
            coordinates[0] = new CoordinatesVoxels(p);
            coordinates[0].matchVoxels(c);
        }

        init();

        createVector(true); // sets ParamVector

    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double clampValue) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;

        pulseTimer = null;
        save = new Save(p);

        Ctotal = clampValue;
        clamp = true;

        if (c == null) {
            coordinates = null;
        } else {
            coordinates = new CoordinatesVoxels[1];
            coordinates[0] = new CoordinatesVoxels(p);
            coordinates[0].matchVoxels(c);
        }

        init();

        createVector(true); // sets ParamVector

    }

    public void addCoordinates(CoordinatesVoxels c) {

        int inew = 1 + coordinates.length;

        CoordinatesVoxels[] ctemp = new CoordinatesVoxels[inew];

        System.arraycopy(coordinates, 0, ctemp, 0, coordinates.length);

        ctemp[inew - 1] = new CoordinatesVoxels(project);
        ctemp[inew - 1].matchVoxels(c);

        coordinates = ctemp;

        init();

    }

    public void removeCoordinates(int iremove) {

        int j = 0, inew;

        if ((iremove < 0) || (iremove >= coordinates.length)) {
            return;
        }

        inew = coordinates.length - 1;

        CoordinatesVoxels[] ctemp = new CoordinatesVoxels[inew];

        for (int ic = 0; ic < coordinates.length; ic++) {
            if (ic == iremove) {
                continue;
            }
            ctemp[j] = coordinates[ic];
            j++;
        }

        coordinates = ctemp;

        init();

    }

    public CoordinatesVoxels coordinates(int i) {

        if (coordinates == null) {
            return null;
        }

        if ((i >= 0) && (i < coordinates.length)) {
            return coordinates[i];
        }

        return null;

    }

    public void updateCoordinates() {

        voxels = 0;
        spaceVoxels = 0;

        if (coordinates == null) {
            return;
        }
        
        for (CoordinatesVoxels c :  coordinates) {

            if (c == null) {
                continue;
            }

            c.update();

            for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
                for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                    for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                        if (c.isInside(i, j, k)) {
                            voxels++;
                            if (project.geometry.isSpace(i, j, k)) {
                                spaceVoxels++;
                            }
                        }
                    }
                }
            }

            //voxels += coordinates[ic].voxels;
            //spaceVoxels += coordinates[ic].spaceVoxels;
        }

    }

    public boolean isInside(int xVoxel, int yVoxel, int zVoxel) {

        if (noOverlapCoordinates) {
            return isInsideNoOverlaps(xVoxel, yVoxel, zVoxel);
        }

        if (coordinates == null) {
            return false;
        }

        for (CoordinatesVoxels c :  coordinates) {
            if ((c != null) && c.isInside(xVoxel, yVoxel, zVoxel)) {
                return true;
            }
        }

        return false;

    }

    public boolean isInsideNoOverlaps(int xVoxel, int yVoxel, int zVoxel) {

        boolean inside, inside2;

        if (coordinates == null) {
            return false;
        }

        for (CoordinatesVoxels c1 :  coordinates) {

            if (c1 == null) {
                continue;
            }

            inside = c1.isInside(xVoxel, yVoxel, zVoxel);

            if (!inside) {
                continue;
            }

            for (CoordinatesVoxels c2 :  coordinates) {
                
                if ((c2 == null) || (c2 == c1)) {
                    continue;
                }

                if (c2.shell) {
                    c2.shell = false;
                    inside2 = c2.isInside(xVoxel, yVoxel, zVoxel);
                    c2.shell = true;
                } else {
                    inside2 = c2.isInside(xVoxel, yVoxel, zVoxel);
                }

                if (inside2) {
                    return false; // overlap
                }

            }

            return true;

        }

        return false;

    }

    public boolean SetDiffusantNum(int arrayNum) {

        Diffusant d = project.getDiffusant(arrayNum);

        if (d != null) {
            diffusantNum = arrayNum;
            diffusantName = d.name;
            updateVectors();
            return true;
        } else {
            error("SetDiffusantNum", "diffusantNum", "failed to find diffusant #" + arrayNum);
            updateVectors();
            return false;
        }

    }

    public boolean SetDiffusantName(String dName) {

        int i = project.getDiffusantNum(dName);

        if (i >= 0) {
            diffusantName = dName;
            diffusantNum = i;
            updateVectors();
            return true;
        } else {
            Master.log("Source SetDiffusantName error: failed to find diffusant " + dName);
            updateVectors();
            return false;
        }

    }

    public void init() {

        if (save != null) {
            save.init();
            saveFileName();
            saveDimensions();
        }

        setParamError("diffusantNum", null);

        diffusantName = project.diffusantName(diffusantNum);

        if (!project.checkDiffusantNum(diffusantNum)) {
            error("init", "diffusantNum", "out of range");
        }

        updateCoordinates();
        updateStats();

    }

    public int numVoxels() {
        if (totalVoxelsForSpaceVoxelsOnly) {
            return spaceVoxels;
        } else {
            return voxels;
        }
    }

    public void updateStats() {
        double charge;
        double molesTotal, sumAmp, sumDuration;

        double litersPerVoxel = Utility.litersPerVoxel(project.dx);
        int numVoxels = numVoxels();
        
        String initSelect = "Ctotal"; // this is fixed to Ctotal now

        if (clamp) {
            //Qmax = Double.NaN;
            // Ctotal is used for clamping
            Ntotal = Double.NaN;
            Nvoxel = Double.NaN;
            updateVectors();
            return;
        }

        charge = project.diffusants[diffusantNum].charge;

        C2I_conversionFactor = 1e12 * charge * Utility.FARADAY * litersPerVoxel / project.dt; // convert mM to pA

        if (pulseTimer != null) {

            sumAmp = pulseTimer.sumAmplitude();
            sumDuration = pulseTimer.sumDuration();

            molesTotal = molesTotal(initSelect, sumAmp, sumDuration);
            //Ctotal = 1e3 * molesTotal / (numVoxels * litersPerVoxel);
            Ctotal = sumAmp;
            Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            Nvoxel = Ntotal / numVoxels;
            //Qmax = Ctotal / sumDuration;

        } else if (initSelect.equalsIgnoreCase("Qmax")) { // mM/ms

            //molesTotal = molesTotal(initSelect, Qmax, project.simTime);
            //Ctotal = 1e3 * molesTotal / (numVoxels * litersPerVoxel);
            //Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            //Nvoxel = Ntotal / numVoxels;
            //Qmax = Ctotal / project.simTime;

        } else if (initSelect.equalsIgnoreCase("Ctotal")) { // mM, M = moles/L

            molesTotal = molesTotal(initSelect, Ctotal, project.simTime);
            //Ctotal = 1e3 * molesTotal / (numVoxels * litersPerVoxel);
            Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            Nvoxel = Ntotal / numVoxels;
            //Qmax = Ctotal / project.simTime;

        } else if (initSelect.equalsIgnoreCase("Ntotal")) {

            molesTotal = molesTotal(initSelect, Ntotal, project.simTime);
            Ctotal = 1e3 * molesTotal / (numVoxels * litersPerVoxel);
            //Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            Nvoxel = Ntotal / numVoxels;
            //Qmax = Ctotal / project.simTime;

        } else if (initSelect.equalsIgnoreCase("Nvoxel")) {

            molesTotal = molesTotal(initSelect, Nvoxel, project.simTime);
            Ctotal = 1e3 * molesTotal / (numVoxels * litersPerVoxel);
            Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            //Nvoxel = Ntotal / numVoxels;
            //Qmax = Ctotal / project.simTime;

        }

        updateVectors();

    }

    public double molesTotal(String select, double amp, double time) {

        double litersPerVoxel = Utility.litersPerVoxel(project.dx);
        int numVoxels = numVoxels();

        if (select.equalsIgnoreCase("Qmax")) { // mM/ms
            return (amp * 1e-3) * numVoxels * litersPerVoxel * time;
        } else if (select.equalsIgnoreCase("Ctotal")) { // mM, M = moles/L
            return (amp * 1e-3) * numVoxels * litersPerVoxel;
        } else if (select.equalsIgnoreCase("Ntotal")) { // #
            return amp / Utility.AVOGADRO;
        } else if (select.equalsIgnoreCase("Nvoxel")) { // #
            return (amp / Utility.AVOGADRO) * numVoxels;
        }

        return Double.NaN;

    }

    public void initPulseTimer() {

        if ((pulseTimer == null) || (pulseTimer.pulses == null)) {

            init();

            if (clamp) {
                Cvoxel = Ctotal;
            } else {
                Cvoxel = Ctotal / project.simPoints(); // spread evenly across each dt
            }

            return;
            
        }

        Cvoxel = Double.NaN;

        pulseTimer.initTimer();

        updateVectors();

    }

    public double cvoxel(Pulse p) { // NOT USED ANYMORE
        double litersPerVoxel = Utility.litersPerVoxel(project.dx);
        int numVoxels = numVoxels();
        double molesTotal = molesTotal("Ctotal", p.amplitude, p.duration);
        double ctotal = 1e3 * molesTotal / (numVoxels * litersPerVoxel); // mM

        double cvoxel = ctotal / (p.itend - p.itbgn + 1); // spread evenly across pulse

        return cvoxel;

    }

    public double release(RunFiniteDifference fd, Geometry geometry) {
        double cvoxel, c2, cdelta = 0;

        if (spaceVoxels == 0) {
            return 0;
        }

        if ((fd.diffus == null) || (fd.it >= fd.itmax)) {
            return 0;
        }

        if ((pulseTimer == null) || (pulseTimer.pulses == null)) {
            cvoxel = Cvoxel; // source is on for all time
        } else {
            cvoxel = pulseTimer.high(fd.it);
        }

        if (cvoxel > 0) {

            if (fd.diffus[0].length == 1) { // single compartment

                if (clamp) {
                    c2 = cvoxel - fd.diffus[diffusantNum][0];
                    fd.diffus[diffusantNum][0] = cvoxel;
                    cdelta = c2;
                } else {
                    fd.diffus[diffusantNum][0] += cvoxel;
                    cdelta = cvoxel;
                }

            } else {

                //if (indexRFD == null) {
                    //Master.exit("Source error: coordinates index array has not be initialized.");
                    //return 0;
                //}

                for (int i : indexRFD) {

                    if (clamp) {
                        c2 = cvoxel - fd.diffus[diffusantNum][i];
                        fd.diffus[diffusantNum][i] = cvoxel;
                        cdelta += c2;
                    } else {
                        fd.diffus[diffusantNum][i] += cvoxel;
                        cdelta += cvoxel;
                    }

                }

            }

        }

        saveValue(cdelta);

        return cdelta;

    }

    public double displayValue(int xVoxel, int yVoxel, int zVoxel) {

        if (isInside(xVoxel, yVoxel, zVoxel)) {
            return Ctotal;
        }

        return -1;
    }

    public Color displayColor(int xVoxel, int yVoxel, int zVoxel, double min, double max) {

        if ((color != null) && isInside(xVoxel, yVoxel, zVoxel)) {
            return color.getColor(Ctotal, min, max);
        }

        return null;

    }

    public boolean check() {

        boolean ok = true;

        updateCoordinates();

        if (voxels != spaceVoxels) {
            Master.log("Warning: " + name + " is located at non-space voxel(s)");
        }

        return ok;

    }

    public boolean saveFileName() {

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
            save.ydim = saveSelect;
        } else {
            save.ydim = diffusantName + " (" + saveSelect + ")";
        }

    }

    public boolean saveInit() {

        int dataPoints = 1;

        if (save == null) {
            return false;
        }

        return save.init(name, coordinates(0), -1, dataPoints);

    }

    public boolean saveFinish() {

        if (save == null) {
            return false;
        }

        return save.finish(name, coordinates(0), -1);

    }

    public boolean saveValue(double value) {

        if (saveSelect.equalsIgnoreCase("mM/ms")) {
            return save.saveData(value / project.dt);
        } else if (saveSelect.equalsIgnoreCase("mM")) {
            return save.saveData(value); // this varies depending on dt
        } else if (saveSelect.equalsIgnoreCase("pA")) {
            return save.saveData(value * C2I_conversionFactor);
        }

        return false;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (color != null) {
            color.addUser(pv);
        }

        if (coordinates != null) {
            for (CoordinatesVoxels c :  coordinates) {
                if (c != null) {
                    c.addUser(pv);
                }
            }
        }

        if (save != null) {
            save.addUser(pv);
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
            for (CoordinatesVoxels c :  coordinates) {
                if (c != null) {
                    c.createVector(true);
                    addVector(c.getVector());
                    c.addUser(this);
                }
            }
        }

        if (save != null) {
            addBlankParam();
            save.createVector(true);
            addVector(save.getVector());
            save.addUser(this);
        }

        if (pulseTimer != null) {
            addBlankParam();
            pulseTimer.createVector(true);
            addVector(pulseTimer.getVector());
            pulseTimer.addUser(this);
        }

        if (close) {
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
            for (CoordinatesVoxels c :  coordinates) {
                if (c != null) {
                    c.updateVector(v);
                }
            }
        }

        if (save != null) {
            save.updateVector(v);
        }

        if (pulseTimer != null) {
            pulseTimer.updateVector(v);
        }

    }
    
    public void setCtotal(double ntotal) {
        double molesTotal = ntotal / Utility.AVOGADRO;
        double litersPerVoxel = Utility.litersPerVoxel(project.dx);
        Ctotal = 1e3 * molesTotal / (numVoxels() * litersPerVoxel);
        updateStats();
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }
        
        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            for (CoordinatesVoxels c :  coordinates) {
                if (o.paramVector.equals(c)) {
                    if (c.setMyParams(o, v)) {
                        init();
                        return true;
                    } else {
                        return false;
                    }
                }
            }
            return false;
        }
        
        if ((o.paramVector instanceof Pulse) && (pulseTimer != null) && (pulseTimer.pulses != null)) {

            for (Pulse p : pulseTimer.pulses) {
                if (o.paramVector.equals(p)) {
                    if (p.setMyParams(o, v)) {
                        init();
                        return true;
                    } else {
                        return false;
                    }

                }

            }

            return false;

        }

        if (!(o.paramVector instanceof Source)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("diffusantNum")) {
            if (v < 0) {
                return false;
            }
            return SetDiffusantNum((int) v);
        }
        if (n.equalsIgnoreCase("Ctotal")) {
            if (v < 0) {
                return false;
            }
            Ctotal = v;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("Ntotal")) {
            if (v < 0) {
                return false;
            }
            setCtotal(v);
            init();
            return true;
        }
        if (n.equalsIgnoreCase("Nvoxel")) {
            if (v < 0) {
                return false;
            }
            setCtotal(v * numVoxels());
            init();
            return true;
        }
        if (n.equalsIgnoreCase("clamp")) {
            clamp = (v == 1);
            init();
            return true;
        }
        if (n.equalsIgnoreCase("totalVoxelsForSpaceVoxelsOnly")) {
            totalVoxelsForSpaceVoxelsOnly = (v == 1);
            init();
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

        String n = o.getName();
        
        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            for (CoordinatesVoxels c :  coordinates) {
                if (o.paramVector.equals(c)) {
                    if (c.setMyParams(o, s)) {
                        init();
                        return true;
                    } else {
                        return false;
                    }
                }
            }
            return false;
        }
       
        if (!(o.paramVector instanceof Source)) {
            return false;
        }

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        if (n.equalsIgnoreCase("diffusantName")) {
            return SetDiffusantName(s);
        }
        if (n.equalsIgnoreCase("saveSelect")) {
            if (saveSelect.equalsIgnoreCase("mM/ms")) {
                saveSelect = s;
            } else if (saveSelect.equalsIgnoreCase("mM")) {
                saveSelect = s;
            } else if (saveSelect.equalsIgnoreCase("pA")) {
                saveSelect = s;
            } else {
                return false;
            }
            return true;
        }
        return super.setMyParams(o, s);
    }
}
