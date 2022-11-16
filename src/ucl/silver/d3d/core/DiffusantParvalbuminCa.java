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
public class DiffusantParvalbuminCa extends Diffusant {

    // [PB-Mg] = [Mg] + [PB] + [Ca] = [PB-Ca]
    //
    // Ctotal = [PB-Mg] + [PB] + [PB-Ca]
    //
    // d[PB]/dt = [PB-Mg]K1off - [PB][Mg]K1on - [PB][Ca]K2on + [PB-Ca]K2off
    // d[PB-Mg]dt = [PB][Mg]K1on - [PB-Mg]K1off
    // d[PB-Ca]/dt = [PB][Ca]K2on - [PB-Ca]K2off
    //
    // this Diffusant represents product [PB-Ca]
    // [PB-Mg] is not saved, but can be computed from: Ctotal = [PB-Mg] + [PB] + [PB-Ca]
    // from Eggermann and Jonas 2012
    // "How the ‘slow’ Ca2+ buffer parvalbumin affects transmitter release in nanodomain-coupling regimes"
    // VOLUME 15 | NUMBER 1 | JANUARY 2012 nature neuroscience
    // Room Temperature 22-24C
    //
    public int iMg = -1; // diffusant number for free Mg
    public int iCa = -1; // diffusant number for free Ca
    public int iPB = -1; // diffusant number for free parvalbumin
    public int iPBCa = -1; // this diffusant
    
    public double Ctotal; // total concentration of parvalbumin (mM)

    public double K1on = 1; // on reaction rate (1/(mM*ms))
    public double K1off = 0.03; // off reaction rate (1/ms)
    public double K2on = 400; // on reaction rate (1/(mM*ms))
    public double K2off = 0.004; // off reaction rate (1/ms)
    public double Kd1 = K1off / K1on;
    public double Kd2 = K2off / K2on;
    
    public Q10 Q10_K1on = null; // Q10 temperature scaling for K1on
    public Q10 Q10_K1off = null; // Q10 temperature scaling for K1off
    public Q10 Q10_K2on = null; // Q10 temperature scaling for K2on
    public Q10 Q10_K2off = null; // Q10 temperature scaling for K2off

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("Ctotal")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("K1on")) {
            return "1/" + project.concUnits + "*" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("K1off")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("K2on")) {
            return "1/" + project.concUnits + "*" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("K2off")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("Kd1")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("Kd2")) {
            return project.concUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("iPBCa")) {
            return false;
        }
        if (name.equalsIgnoreCase("Kd1")) {
            return false;
        }
        if (name.equalsIgnoreCase("Kd2")) {
            return false;
        }
        return super.canEdit(name);
    }

    public DiffusantParvalbuminCa(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int IMG, int ICA, int IPB) { // init reaction particle [ab]

        super(p, NAME, 0, diffusionConstant, c);

        Ctotal = CTOTAL;
        iMg = IMG;
        iCa = ICA;
        iPB = IPB;

        reaction = true;

        createVector(true);

    }

    public DiffusantParvalbuminCa(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int IMG, int ICA, int IPB, double K1ON, double K1OFF, double K2ON, double K2OFF) { // init reaction particle [ab]

        super(p, NAME, 0, diffusionConstant, c);

        Ctotal = CTOTAL;
        iMg = IMG;
        iCa = ICA;
        iPB = IPB;
        K1on = K1ON;
        K1off = K1OFF;
        K2on = K2ON;
        K2off = K2OFF;
        Kd1 = K1off / K1on;
        Kd2 = K2off / K2on;

        reaction = true;

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        reaction = true;
        iPBCa = -1;

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] == this) {
                iPBCa = i;
                break;
            }
        }

        if (!project.checkDiffusantNum(iMg)) {
            error("init", "iMg", "out of range");
        }

        if (!project.checkDiffusantNum(iCa)) {
            error("init", "iCa", "out of range");
        }

        if (!project.checkDiffusantNum(iPB)) {
            error("init", "iPB", "out of range");
        }

        if (iPBCa >= 0) {
            if (iPBCa == iMg) {
                error("init", "iMg", "value conflicts with iPBCa");
            }
            if (iPBCa == iCa) {
                error("init", "iCa", "value conflicts with iPBCa");
            }
            if (iPBCa == iPB) {
                error("init", "iPB", "value conflicts with iPBCa");
            }
        }

        computeC0();

    }
    
    @Override
    public void Q10execute() {
        double k1on, k1off, k2on, k2off;
        super.Q10execute();
        if ((Q10_K1on != null) && Q10_K1on.on) {
            k1on = Q10_K1on.getScaledValue();
            if (Double.isFinite(k1on)) {
                K1on = k1on;
            }
        }
        if ((Q10_K1off != null) && Q10_K1off.on) {
            k1off = Q10_K1off.getScaledValue();
            if (Double.isFinite(k1off)) {
                K1off = k1off;
            }
        }
        if ((Q10_K2on != null) && Q10_K2on.on) {
            k2on = Q10_K2on.getScaledValue();
            if (Double.isFinite(k2on)) {
                K2on = k2on;
            }
        }
        if ((Q10_K2off != null) && Q10_K2off.on) {
            k2off = Q10_K2off.getScaledValue();
            if (Double.isFinite(k2off)) {
                K2off = k2off;
            }
        }
        Kd1 = K1off / K1on;
        Kd2 = K2off / K2on;
    }

    public double computeC0() {

        double r1, r2;

        Diffusant Mg, Ca, PB;

        Kd1 = K1off / K1on;
        Kd2 = K2off / K2on;
        C0 = Double.NaN;

        if ((iMg >= 0) && (iMg < project.diffusants.length)) {
            if ((iCa >= 0) && (iCa < project.diffusants.length)) {
                if ((iPB >= 0) && (iPB < project.diffusants.length)) {

                    Mg = project.diffusants[iMg];
                    Ca = project.diffusants[iCa];
                    PB = project.diffusants[iPB];

                    r1 = Mg.C0 * Kd2 / (Ca.C0 * Kd1);
                    r2 = ((r1 * K1off + K2off) / (Mg.C0 * K1on + Ca.C0 * K2on)) + r1 + 1;
                    C0 = Ctotal / r2;

                    PB.C0 = C0 * K2off / (Ca.C0 * K2on);

                }
            }
        }

        return C0;

    }

    @Override
    public void react(RunFiniteDifference fd, int iPBCA) {

        int ipnt = fd.thisPnt;
        double Mg, Ca, PB, PBMg, PBCa, dMg, dCa, dPB, dPBCa, dPBMg;

        if (iPBCA != iPBCa) {
            return;
        }
        
        if ((pulseTimer != null) && (pulseTimer.high(fd.it) == 0)) {
            return;
        }

        Mg = fd.diffus[iMg][ipnt];
        Ca = fd.diffus[iCa][ipnt];
        PB = fd.diffus[iPB][ipnt];
        PBCa = fd.diffus[iPBCA][ipnt];
        PBMg = (Ctotal - PB - PBCa);

        dPB = (PBMg * K1off - PB * Mg * K1on - PB * Ca * K2on + PBCa * K2off) * project.dt;
        dPBMg = (PB * Mg * K1on - PBMg * K1off) * project.dt;
        dPBCa = (PB * Ca * K2on - PBCa * K2off) * project.dt;

        dMg = -dPBMg;
        dCa = -dPBCa;

        fd.diffnext[iMg][ipnt] += dMg;
        fd.diffnext[iCa][ipnt] += dCa;
        fd.diffnext[iPB][ipnt] += dPB;
        fd.diffnext[iPBCA][ipnt] += dPBCa;

        if (fd.diffnext[iMg][ipnt] < 0) {
            error("react", "fd.diffnext[iMg]", "negative concentration for diffusant #" + iMg);
            fd.stopSimulation();
        }

        if (fd.diffnext[iCa][ipnt] < 0) {
            error("react", "fd.diffnext[iCa]", "negative concentration for diffusant #" + iCa);
            fd.stopSimulation();
        }

        if (fd.diffnext[iPB][ipnt] < 0) {
            error("react", "fd.diffnext[iPB]", "negative concentration for diffusant #" + iPB);
            fd.stopSimulation();
        }

        if (fd.diffnext[iPBCA][ipnt] < 0) {
            error("react", "fd.diffnext[iPBCa]", "negative concentration for diffusant #" + iPBCA);
            fd.stopSimulation();
        }

    }
    
    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);
        
        if (Q10_K1on != null) {
            Q10_K1on.addUser(pv);
        }
        
        if (Q10_K1off != null) {
            Q10_K1off.addUser(pv);
        }
        
        if (Q10_K2on != null) {
            Q10_K2on.addUser(pv);
        }
        
        if (Q10_K2off != null) {
            Q10_K2off.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }
        
        if (Q10_K1on != null) {
            addBlankParam();
            Q10_K1on.createVector(true);
            addVector(Q10_K1on.getVector());
            Q10_K1on.addUser(this);
        }
        
        if (Q10_K1off != null) {
            addBlankParam();
            Q10_K1off.createVector(true);
            addVector(Q10_K1off.getVector());
            Q10_K1off.addUser(this);
        }
        
        if (Q10_K2on != null) {
            addBlankParam();
            Q10_K2on.createVector(true);
            addVector(Q10_K2on.getVector());
            Q10_K2on.addUser(this);
        }
        
        if (Q10_K2off != null) {
            addBlankParam();
            Q10_K2off.createVector(true);
            addVector(Q10_K2off.getVector());
            Q10_K2off.addUser(this);
        }

        if (close){
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);
        
        if (Q10_K1on != null) {
            Q10_K1on.updateVector(v);
        }
        
        if (Q10_K1off != null) {
            Q10_K1off.updateVector(v);
        }
        
        if (Q10_K2on != null) {
            Q10_K2on.updateVector(v);
        }
        
        if (Q10_K2off != null) {
            Q10_K2off.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof Diffusant)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("Ctotal")) {
            if (v < 0) {
                return false;
            }
            Ctotal = v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("iMg")) {
            if (v < 0) {
                return false;
            }
            iMg = (int) v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("iCa")) {
            if (v < 0) {
                return false;
            }
            iCa = (int) v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("iPB")) {
            if (v < 0) {
                return false;
            }
            iPB = (int) v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("K1on")) {
            if (v < 0) {
                return false;
            }
            K1on = v;
            if (Q10_K1on != null) {
                Q10_K1on.value = v;
            }
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("K1off")) {
            if (v < 0) {
                return false;
            }
            K1off = v;
            if (Q10_K1off != null) {
                Q10_K1off.value = v;
            }
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("K2on")) {
            if (v < 0) {
                return false;
            }
            K2on = v;
            if (Q10_K2on != null) {
                Q10_K2on.value = v;
            }
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("K2off")) {
            if (v < 0) {
                return false;
            }
            K2off = v;
            if (Q10_K2off != null) {
                Q10_K2off.value = v;
            }
            computeC0();
            return true;
        }
        return super.setMyParams(o, v);
    }

}
