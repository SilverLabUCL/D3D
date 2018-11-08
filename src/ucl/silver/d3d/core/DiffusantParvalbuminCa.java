package ucl.silver.d3d.core;

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
public class DiffusantParvalbuminCa extends Diffusant {

    // [PB-Mg] = [Mg] + [PB] + [Ca] = [PB-Ca]
    //
    // cTotal = [PB-Mg] + [PB] + [PB-Ca]
    //
    // d[PB]/dt = [PB-Mg]k1off - [PB][Mg]k1on - [PB][Ca]k2on + [PB-Ca]k2off
    // d[PB-Mg]dt = [PB][Mg]k1on - [PB-Mg]k1off
    // d[PB-Ca]/dt = [PB][Ca]k2on - [PB-Ca]k2off
    //
    // this Diffusant represents product [PB-Ca]
    // [PB-Mg] is not saved, but can be computed from: cTotal = [PB-Mg] + [PB] + [PB-Ca]
    // from Eggermann and Jonas 2012
    // "How the ‘slow’ Ca2+ buffer parvalbumin affects transmitter release in nanodomain-coupling regimes"
    // VOLUME 15 | NUMBER 1 | JANUARY 2012 nature neuroscience
    //
    public int iMg = -1; // diffusant number for free Mg
    public int iCa = -1; // diffusant number for free Ca
    public int iPB = -1; // diffusant number for free parvalbumin
    public int iPBCa = -1; // this diffusant
    
    public double cTotal; // total concentration of parvalbumin (mM)

    public double k1on = 1; // on reaction rate (1/(mM*ms))
    public double k1off = 0.03; // off reaction rate (1/ms)
    public double k2on = 400; // on reaction rate (1/(mM*ms))
    public double k2off = 0.004; // off reaction rate (1/ms)
    public double kd1 = k1off / k1on;
    public double kd2 = k2off / k2on;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("cTotal")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("k1on")) {
            return "1/" + project.concUnits + "*" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("k1off")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("k2on")) {
            return "1/" + project.concUnits + "*" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("k2off")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("kd1")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("kd2")) {
            return project.concUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("iPBCa")) {
            return false;
        }
        if (name.equalsIgnoreCase("kd1")) {
            return false;
        }
        if (name.equalsIgnoreCase("kd2")) {
            return false;
        }
        return super.canEdit(name);
    }

    public DiffusantParvalbuminCa(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int IMG, int ICA, int IPB) { // init reaction particle [ab]

        super(p, NAME, 0, diffusionConstant, c);

        cTotal = CTOTAL;
        iMg = IMG;
        iCa = ICA;
        iPB = IPB;

        reaction = true;

        createVector(true);

    }

    public DiffusantParvalbuminCa(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int IMG, int ICA, int IPB, double K1ON, double K1OFF, double K2ON, double K2OFF) { // init reaction particle [ab]

        super(p, NAME, 0, diffusionConstant, c);

        cTotal = CTOTAL;
        iMg = IMG;
        iCa = ICA;
        iPB = IPB;
        k1on = K1ON;
        k1off = K1OFF;
        k2on = K2ON;
        k2off = K2OFF;
        kd1 = k1off / k1on;
        kd2 = k2off / k2on;

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

    public double computeC0() {

        double r1, r2;

        Diffusant Mg, Ca, PB;

        kd1 = k1off / k1on;
        kd2 = k2off / k2on;
        C0 = Double.NaN;

        if ((iMg >= 0) && (iMg < project.diffusants.length)) {
            if ((iCa >= 0) && (iCa < project.diffusants.length)) {
                if ((iPB >= 0) && (iPB < project.diffusants.length)) {

                    Mg = project.diffusants[iMg];
                    Ca = project.diffusants[iCa];
                    PB = project.diffusants[iPB];

                    r1 = Mg.C0 * kd2 / (Ca.C0 * kd1);
                    r2 = ((r1 * k1off + k2off) / (Mg.C0 * k1on + Ca.C0 * k2on)) + r1 + 1;
                    C0 = cTotal / r2;

                    PB.C0 = C0 * k2off / (Ca.C0 * k2on);

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
        PBMg = (cTotal - PB - PBCa);

        dPB = (PBMg * k1off - PB * Mg * k1on - PB * Ca * k2on + PBCa * k2off) * project.dt;
        dPBMg = (PB * Mg * k1on - PBMg * k1off) * project.dt;
        dPBCa = (PB * Ca * k2on - PBCa * k2off) * project.dt;

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

        if (n.equalsIgnoreCase("cTotal")) {
            if (v < 0) {
                return false;
            }
            cTotal = v;
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
        if (n.equalsIgnoreCase("k1on")) {
            if (v < 0) {
                return false;
            }
            k1on = v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("k1off")) {
            if (v < 0) {
                return false;
            }
            k1off = v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("k2on")) {
            if (v < 0) {
                return false;
            }
            k2on = v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("k2off")) {
            if (v < 0) {
                return false;
            }
            k2off = v;
            computeC0();
            return true;
        }
        return super.setMyParams(o, v);
    }

}
