package ucl.silver.d3d.core;

import java.awt.Color;
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
public class Source extends ParamVector {
    
    // simple step source
    // on for all time, unless PulseTimer exists
    
    // Pulse time = onset (ms)
    // Pulse duration (ms)
    // Pulse amplitude = Q (mM/ms) or Ctotal (mM) or Ntotal (molecules) (see type)

    public String diffusantName = null;
    public int diffusantNum = -1; // diffusant array number in Project
    
    public String type = "Q";
    // "Q"  (mM/ms)
    // "C" for Ctotal (mM)
    // "Cclamp" for clamping to Cvoxel (mM)
    // "N" for Ntotal (molecules)
    // "K" (1/ms) (e.g. uptake)

    public double Q = Double.NaN; // mM/ms
    public double Ctotal = Double.NaN; // mM, M = moles/L
    public double Cvoxel = Double.NaN; // mM // Ctotal / numVoxels
    private double Cvoxel_dt = Double.NaN; // mM per dt, used in "release" if no PulseTimer
    public double Ntotal = Double.NaN; // molecules
    public double Nvoxel = Double.NaN; // molecules
    public boolean clamp = false; // source clamps voxels to Cvoxel
    
    public int voxels, spaceVoxels;
    public boolean totalVoxelsForSpaceVoxelsOnly = true;
    public boolean noOverlapCoordinates = false;
    private double liters = Double.NaN;
    
    public String saveSelect = "mM/ms"; // "mM/ms" or "mM" or "pA", output file units
    private double C2I = Double.NaN; // conversion scale factor for saving in pA

    private final Color defaultColor = new Color(153, 0, 51);

    public ColorD3D color = null;
    public CoordinatesVoxels[] coordinates = null; // coordinates of sources
    public transient int[] indexRFD = null; // index counter for FD simulations
    public Save save = null;
    public PulseTimer pulseTimer = null; // timer for source

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("Q")) {
            return project.concUnits + "/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("Ctotal")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("Cvoxel")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("Ntotal")) {
            return "#";
        }
        if (name.equalsIgnoreCase("Nvoxel")) {
            return "#";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        boolean timer = (pulseTimer != null);
        if (name.equalsIgnoreCase("type")) {
            return false;
        }
        if (name.equalsIgnoreCase("Q")) {
            return type.equalsIgnoreCase("Q") && !timer;
        }
        if (name.equalsIgnoreCase("Ctotal")) {
            return type.equalsIgnoreCase("C") && !timer && !clamp;
        }
        if (name.equalsIgnoreCase("Cvoxel")) {
            return type.equalsIgnoreCase("C") && !timer && clamp;
        }
        if (name.equalsIgnoreCase("Ntotal")) {
            return type.equalsIgnoreCase("N") && !timer;
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
        return super.canEdit(name);
    }
    
    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double valueQCN, String TYPE) {
        
        super(p);
        
        if (TYPE.equalsIgnoreCase("Q")) {
            Q = valueQCN;
        } else if (TYPE.equalsIgnoreCase("C") || TYPE.equalsIgnoreCase("Ctotal") || TYPE.equalsIgnoreCase("Cclamp")) {
            if (TYPE.equalsIgnoreCase("Cclamp")) {
                Cvoxel = valueQCN;
            } else {
                Ctotal = valueQCN;
            }
        } else if (TYPE.equalsIgnoreCase("N") || TYPE.equalsIgnoreCase("Ntotal")) {
            Ntotal = valueQCN;
        }

        sourceInit(NAME, DiffusantNum, c, null, TYPE);
        
    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double onset, double duration, double valueQCN, String TYPE) {
        super(p);
        pulseTimer = new PulseTimer(project, onset, duration, valueQCN, typeUnits(TYPE));
        sourceInit(NAME, DiffusantNum, c, pulseTimer, TYPE);
    }

    public Source(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt, String TYPE) {
        super(p);
        sourceInit(NAME, DiffusantNum, c, pt, TYPE);
    }
    
    private void sourceInit(String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt, String TYPE) {

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Source";
        }

        color = new ColorD3D(name + "_color", defaultColor);

        diffusantNum = DiffusantNum;
        
        clamp = TYPE.equalsIgnoreCase("Cclamp");
        
        if (TYPE.equalsIgnoreCase("Q")) {
            type = "Q";
        } else if (TYPE.equalsIgnoreCase("C") || TYPE.equalsIgnoreCase("Ctotal") || TYPE.equalsIgnoreCase("Cclamp")) {
            type = "C";
        } else if (TYPE.equalsIgnoreCase("N") || TYPE.equalsIgnoreCase("Ntotal")) {
            type = "N";
        } else if (TYPE.equalsIgnoreCase("K")) {
            type = "K";
        } else {
            error("Source", "type", "bad value" + TYPE);
        }

        pulseTimer = pt;
        pulseTimer.ampUnits = typeUnits(TYPE);
        
        save = new Save(project);

        if (c == null) {
            coordinates = null;
        } else {
            coordinates = new CoordinatesVoxels[1];
            coordinates[0] = new CoordinatesVoxels(project);
            coordinates[0].matchVoxels(c);
        }

        init();
        
        createVector(true); // sets ParamVector
        
    }
    
    public static String typeUnits(String TYPE) {
        if (TYPE.equalsIgnoreCase("Q")) {
            return "mM/ms";
        } else if (TYPE.equalsIgnoreCase("C") || TYPE.equalsIgnoreCase("Ctotal") || TYPE.equalsIgnoreCase("Cclamp")) {
            return "mM";
        } else if (TYPE.equalsIgnoreCase("N") || TYPE.equalsIgnoreCase("Ntotal")) {
            return "#";
        } else if (TYPE.equalsIgnoreCase("K")) {
            return "1/ms";
        } else {
            return "";
        }
    }
    
    public void init() {
        
        double q_dt, charge;

        if (save != null) {
            save.init();
            saveFileName();
            saveDimensions();
            if (clamp) {
                saveSelect = "mM";
            }
        }

        setParamError("diffusantNum", null);

        diffusantName = project.diffusantName(diffusantNum);

        if (!project.checkDiffusantNum(diffusantNum)) {
            error("init", "diffusantNum", "out of range");
        }
        
        liters();
        
        q_dt = 1.0 / project.dt; // mM/ms // unit flux
        charge = project.diffusants[diffusantNum].charge;

        C2I = Utility.Q2I(q_dt, charge, liters); // conversion scale factor

        updateCoordinates();
        scaleQ10();
        updateStats();

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
    
    public void scaleQ10() {
        // see SourceUptake
    }

    public int numVoxels() {
        if (totalVoxelsForSpaceVoxelsOnly) {
            return spaceVoxels;
        } else {
            return voxels;
        }
    }
    
    public double liters() {
        liters = Utility.litersPerVoxel(project.dx) * numVoxels();
        return liters;
    }

    public void updateStats() {
        
        double molesTotal, numVoxels;

        if (clamp) {

            Ctotal = Double.NaN;
            Ntotal = Double.NaN;
            Nvoxel = Double.NaN;
            Q = Double.NaN;

            if (pulseTimer != null) {
                Cvoxel = Double.NaN; // use pulse timer amplitude
            }

            updateVectors();

            return;

        }
        
        if (type.equalsIgnoreCase("K")) {
            
            Ctotal = Double.NaN;
            Cvoxel = Double.NaN;
            Ntotal = Double.NaN;
            Nvoxel = Double.NaN;
            Q = Double.NaN;
            
            updateVectors();

            return;
            
        }
        
        numVoxels = numVoxels();
        
        if (pulseTimer != null) {
            
            Ctotal = 0;
            Ntotal = 0;
            Q = 0;
            
            for (Pulse p : pulseTimer.pulses) {
                molesTotal = molesTotal(type, p.amplitude, p.duration);
                Ctotal += 1e3 * molesTotal / liters; // mM
                Ntotal += molesTotal * Utility.AVOGADRO; // molecules
                Q += (1e3 * molesTotal / liters) / p.duration; // mM/ms
            }
            
        } else if (type.equalsIgnoreCase("Q")) { // mM/ms
            molesTotal = molesTotal(type, Q, project.simTime);
            Ctotal = 1e3 * molesTotal / liters; // mM
            Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            //Q = Ctotal / project.simTime;
        } else if (type.equalsIgnoreCase("C")) { // mM, M = moles/L
            molesTotal = molesTotal(type, Ctotal, project.simTime);
            //Ctotal = 1e3 * molesTotal / liters; // mM
            Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            Q = Ctotal / project.simTime;
        } else if (type.equalsIgnoreCase("N")) {
            molesTotal = molesTotal(type, Ntotal, project.simTime);
            Ctotal = 1e3 * molesTotal / liters; // mM
            //Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            Q = Ctotal / project.simTime;
        } else {
            Ctotal = Double.NaN;
            Ntotal = Double.NaN;
            Q = Double.NaN;
        }
        
        Cvoxel = Ctotal / numVoxels;
        Nvoxel = Ntotal / numVoxels;

        updateVectors();

    }
    
    public double molesTotal(String TYPE, double valueQCN, double time_ms) {
        if (TYPE.equalsIgnoreCase("Q")) { // mM/ms
            return (valueQCN * 1e-3) * liters * time_ms;
        } else if (TYPE.equalsIgnoreCase("C") || TYPE.equalsIgnoreCase("Ctotal")) { // mM, M = moles/L
            return (valueQCN * 1e-3) * liters;
        } else if (TYPE.equalsIgnoreCase("N") || TYPE.equalsIgnoreCase("Ntotal")) { // #
            return valueQCN / Utility.AVOGADRO;
        } else {
            return Double.NaN;
        }
    }

    public void initPulseTimer() {
        
        boolean normalizeAmplitudesPerDT = !clamp;
        
        if (pulseTimer == null) {
            // no timer, source is active for all time

            init();

            if (clamp) {
                Cvoxel_dt = Cvoxel;
            } else {
                Cvoxel_dt = Cvoxel / project.simPoints(); // spread evenly across each dt
            }

            return;

        }
        
        // init PulseTimer...
        
        Cvoxel_dt = Double.NaN;
        
        for (Pulse p : pulseTimer.pulses) {
            p.update_it_params(normalizeAmplitudesPerDT);
            p.it_amp = cVoxel(type, p.it_amp, p.duration); // mM
        }
        
        pulseTimer.initTimer(normalizeAmplitudesPerDT, false);
        
        updateVectors();

    }
    
    public double cVoxel(String TYPE, double valueQCN, double time_ms) {

        double molesTotal;
        double ctotal = Double.NaN; // mM

        if (TYPE.equalsIgnoreCase("Q")) { // mM/ms
            ctotal = valueQCN * time_ms;
        } else if (TYPE.equalsIgnoreCase("C") || TYPE.equalsIgnoreCase("Ctotal")) { // mM, M = moles/L
            ctotal = valueQCN;
        } else if (TYPE.equalsIgnoreCase("N") || TYPE.equalsIgnoreCase("Ntotal")) { // #
            molesTotal = valueQCN / Utility.AVOGADRO;
            ctotal = 1e3 * molesTotal / liters;
        }

        if (clamp) {
            return ctotal;
        } else {
            return ctotal / (1.0 * numVoxels());
        }

    }

    public double release(RunFiniteDifference fd, Geometry geometry) {
        double c, c2, cdelta = 0;

        if (fd.it >= fd.itmax) {
            return 0;
        }
        
        if (pulseTimer == null) {
            c = Cvoxel_dt; // source is on for all time
        } else {
            c = pulseTimer.high(fd.it);
        }

        if (c >= 0) { // include 0 for clamp

            if (fd.diffus[0].length == 1) { // single compartment

                if (clamp) {
                    c2 = c - fd.diffus[diffusantNum][0];
                    fd.diffus[diffusantNum][0] = c;
                    cdelta = c2;
                } else {
                    fd.diffus[diffusantNum][0] += c;
                    cdelta = c;
                }

            } else {

                for (int i : indexRFD) {
                    if (clamp) {
                        c2 = c - fd.diffus[diffusantNum][i];
                        fd.diffus[diffusantNum][i] = c;
                        cdelta += c2;
                    } else {
                        fd.diffus[diffusantNum][i] += c;
                        cdelta += c;
                    }
                }

            }

        }

        saveValue(cdelta);

        return cdelta;

    }

    public double displayValue(int xVoxel, int yVoxel, int zVoxel) {

        if (isInside(xVoxel, yVoxel, zVoxel)) {
            return Cvoxel;
        }

        return -1;
    }

    public Color displayColor(int xVoxel, int yVoxel, int zVoxel, double min, double max) {

        if ((color != null) && isInside(xVoxel, yVoxel, zVoxel)) {
            return color.getColor(Cvoxel, min, max);
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
        
        if (clamp) {
            saveSelect = "mM";
        }

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

    public boolean saveValue(double value_mM) {

        if (saveSelect.equalsIgnoreCase("mM/ms")) {
            return save.saveData(value_mM / project.dt);
        } else if (saveSelect.equalsIgnoreCase("mM")) {
            return save.saveData(value_mM); // this varies depending on dt (except clamping)
        } else if (saveSelect.equalsIgnoreCase("pA")) {
            return save.saveData(C2I * value_mM);
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
        if (n.equalsIgnoreCase("Q")) {
            if (v < 0) {
                return false;
            }
            Q = v;
            type = "Q";
            init();
            return true;
        }
        if (n.equalsIgnoreCase("Ctotal")) {
            if (v < 0) {
                return false;
            }
            Ctotal = v;
            type = "C";
            init();
            return true;
        }
        if (n.equalsIgnoreCase("Cvoxel")) {
            if (v < 0) {
                return false;
            }
            Cvoxel = v;
            type = "C";
            clamp = true;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("Ntotal")) {
            if (v < 0) {
                return false;
            }
            Ntotal = v;
            type = "N";
            init();
            return true;
        }
        if (n.equalsIgnoreCase("clamp")) {
            clamp = (v == 1);
            if (clamp) {
                type = "C";
            }
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
            if (s.equalsIgnoreCase("mM/ms") || s.equalsIgnoreCase("mM") || s.equalsIgnoreCase("pA")) {
                saveSelect = s;
                return true;
            } else {
                return false;
            }
        }
        return super.setMyParams(o, s);
    }
}
