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
public class DiffusantCalretininCa extends Diffusant {

    // [TT] + [Ca] = [TRCa] + [Ca] = [RCaRCa]
    //
    // cTotal = [TT] + [TRCa] + [RCaRCa]
    //
    // d[TRCa]/dt = [TT][Ca]k1on - [TRCa]k1off - [TRCa][Ca]k2on + [RCaRCa]k2off
    // d[TT]dt = [TRCa]k1off - [TT][Ca]k1on
    // d[RCaRCa]/dt = [TRCa][Ca]k2on - [RCaRCa]k2off
    //
    // this Diffusant represents [RCaRCa]
    // [TT] is not saved, but can be computed from: cTotal = [TT] + [TRCa] + [RCaRCa]
    // from Faas et al 2007
    // "Resolving the Fast Kinetics of Cooperative Binding: Ca2+ Buffering by Calretinin"
    // PLoS Biology | www.plosbiology.org
    // November 2007 | Volume 5 | Issue 11 | e311
    //
    public int iCa = -1; // diffusant number for free Ca
    public int iTRCa = -1; // diffusant number for TRCa
    public int iRCaRCa = -1; // this diffusant

    public double cTotal; // total concentration (mM) // [TT] + [TRCa] + [RCaRCa]

    public double k1on = 3.6; // on reaction rate (1/(mM*ms))
    public double k1off = 0.053; // off reaction rate (1/ms)
    public double k2on = 310.0; // on reaction rate (1/(mM*ms))
    public double k2off = 0.040; // off reaction rate (1/ms)
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
        if (name.equalsIgnoreCase("iRCaRCa")) {
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

    public DiffusantCalretininCa(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int ICA, int ITRCA) { // init reaction particle [ab]

        super(p, NAME, 0, diffusionConstant, c);

        cTotal = CTOTAL;

        iCa = ICA;
        iTRCa = ITRCA;

        reaction = true;

        createVector(true);

    }

    public DiffusantCalretininCa(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int ICA, int ITRCA, double K1ON, double K1OFF, double K2ON, double K2OFF) { // init reaction particle [ab]

        super(p, NAME, 0, diffusionConstant, c);

        cTotal = CTOTAL;
        iCa = ICA;
        iTRCa = ITRCA;
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
        iRCaRCa = -1;

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] == this) {
                iRCaRCa = i;
                break;
            }
        }

        if (iRCaRCa >= 0) {
            if (iRCaRCa == iCa) {
                error("init", "iCa", "value conflicts with iRCaRCa");
            }
            if (iRCaRCa == iTRCa) {
                error("init", "iTRCa", "value conflicts with iRCaRCa");
            }
        }

        if (!project.checkDiffusantNum(iCa)) {
            error("init", "iCa", "out of range");
        }

        if (!project.checkDiffusantNum(iTRCa)) {
            error("init", "iTRCa", "out of range");
        }
        
        computeC0();

    }

    public double computeC0() {

        double K1, K2;

        Diffusant Ca, TRCa;

        C0 = Double.NaN;

        if ((iCa >= 0) && (iCa < project.diffusants.length)) {
            if ((iTRCa >= 0) && (iTRCa < project.diffusants.length)) {

                Ca = project.diffusants[iCa];
                TRCa = project.diffusants[iTRCa];

                //r1 = Ca.C0 * k1on - k2off + (k2off / (Ca.C0 * k2on)) * (Ca.C0 * k1on + k1off + Ca.C0 * k2on);
                //C0 = cTotal * Ca.C0 * k1on / r1;
                //TRCa.C0 = C0 * k2off / (Ca.C0 * k2on);
                kd1 = k1off / k1on;
                kd2 = k2off / k2on;

                K1 = k1on / k1off;
                K2 = k2on / k2off;

                TRCa.C0 = K1 * Ca.C0 * cTotal / (1 + K1 * Ca.C0 + K1 * K2 * Ca.C0 * Ca.C0);
                C0 = K1 * K2 * Ca.C0 * Ca.C0 * cTotal / (1 + K1 * Ca.C0 + K1 * K2 * Ca.C0 * Ca.C0);

            }
        }

        return C0;

    }

    @Override
    public void react(RunFiniteDifference fd, int iRCARCA) {

        int ipnt = fd.thisPnt;
        double Ca, dCa, TT, dTT, TRCa, dTRCa, RCaRCa, dRCaRCa;

        if (iRCARCA != iRCaRCa) {
            return;
        }
        
        if ((pulseTimer != null) && (pulseTimer.high(fd.it) == 0)) {
            return;
        }

        Ca = fd.diffus[iCa][ipnt];
        TRCa = fd.diffus[iTRCa][ipnt];
        RCaRCa = fd.diffus[iRCARCA][ipnt];
        TT = (cTotal - TRCa - RCaRCa);

        dTRCa = (TT * Ca * k1on - TRCa * k1off - TRCa * Ca * k2on + RCaRCa * k2off) * project.dt;
        dTT = (TRCa * k1off - TT * Ca * k1on) * project.dt; // not used
        dRCaRCa = (TRCa * Ca * k2on - RCaRCa * k2off) * project.dt;

        dCa = -2 * dRCaRCa - dTRCa; // 2 dRCaRCa for 2 binding sites

        fd.diffnext[iCa][ipnt] += dCa;
        fd.diffnext[iTRCa][ipnt] += dTRCa;
        fd.diffnext[iRCARCA][ipnt] += dRCaRCa;

        if (fd.diffnext[iCa][ipnt] < 0) {
            error("react", "fd.diffnext[iCa]", "negative concentration for diffusant #" + iCa);
            fd.stopSimulation();
        }

        if (fd.diffnext[iTRCa][ipnt] < 0) {
            error("react", "fd.diffnext[iTRCa]", "negative concentration for diffusant #" + iTRCa);
            fd.stopSimulation();
        }

        if (fd.diffnext[iRCARCA][ipnt] < 0) {
            error("react", "fd.diffnext[iRCaRCa]", "negative concentration for diffusant #" + iRCARCA);
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
        if (n.equalsIgnoreCase("iCa")) {
            if (v < 0) {
                return false;
            }
            iCa = (int) v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("iTRCa")) {
            if (v < 0) {
                return false;
            }
            iTRCa = (int) v;
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
