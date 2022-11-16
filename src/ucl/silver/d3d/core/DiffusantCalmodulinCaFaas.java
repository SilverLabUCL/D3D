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
public class DiffusantCalmodulinCaFaas extends Diffusant {
    // T - Tense state
    // R - Relaxed state
    //
    // [TT] + [Ca] = [CaTR] + [Ca] = [CaRCaR]
    //
    // Ctotal = [TT] + [CaTR] + [CaRCaR]
    //
    // d[CaTR]/dt = [TT][Ca]K1on - [CaTR]K1off - [CaTR][Ca]K2on + [CaRCaR]K2off
    // d[TT]dt = [CaTR]K1off - [TT][Ca]K1on
    // d[CaRCaR]/dt = [CaTR][Ca]K2on - [CaRCaR]K2off
    //
    // this Diffusant represents [CaRCaR]
    // [TT] is not saved, but can be computed from: Ctotal = [TT] + [CaTR] + [CaRCaR]
    //
    // "Calmodulin as a direct detector of Ca2+ signals"
    // Faas GC, Raghavachari S, Lisman JE, Mody I
    // Nat Neurosci. 2011 Mar;14(3):301-4. doi: 10.1038/nn.2746. Epub 2011 Jan 23
    // temperature = 35C
    //
    public int iCa = -1; // diffusant number for free Ca
    public int iCaTR = -1; // diffusant number for CaTR
    public int iCaRCaR = -1; // this diffusant

    public double Ctotal; // total concentration (mM) // [TT] + [CaTR] + [CaRCaR]

    public double K1on = Double.NaN; // on reaction rate (1/(mM*ms))
    public double K1off = Double.NaN; // off reaction rate (1/ms)
    public double K2on = Double.NaN; // on reaction rate (1/(mM*ms))
    public double K2off = Double.NaN; // off reaction rate (1/ms)
    public double Kd1 = K1off / K1on;
    public double Kd2 = K2off / K2on;
    
    public double temp = 35; // Faas experiments
    
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
        if (name.equalsIgnoreCase("temp")) {
            return project.tempUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("iCa")) {
            return false;
        }
        if (name.equalsIgnoreCase("iCaTR")) {
            return false;
        }
        if (name.equalsIgnoreCase("iCaRCaR")) {
            return false;
        }
        if (name.equalsIgnoreCase("Kd1")) {
            return false;
        }
        if (name.equalsIgnoreCase("Kd2")) {
            return false;
        }
        if (name.equalsIgnoreCase("temp")) {
            return false;
        }
        return super.canEdit(name);
    }

    public DiffusantCalmodulinCaFaas(Project p, String NAME, double CTOTAL, double diffusionCoefficient, CoordinatesVoxels c,
            int ICA, int ITRCA, String lobeNorC) {

        super(p, NAME, 0, diffusionCoefficient, c);

        iCa = ICA;
        iCaTR = ITRCA;

        reaction = true;
        
        // Ctotal = 0.100 / 2.0; // mM // split between N and C lobes // Faas site:
        // Biber A, Schmid G, Hempel K. 
        // Calmodulin content in specific brain areas. 
        // Exp Brain Res. 1984;56:323–326.
        // 0.03 - 0.12 mM
        //
        // Banay-Schwartz M, Kenessey A, DeGuzman T, Lajtha A, Palkovits M. 
        // Protein content of various regions of rat brain and adult aging human brain.
        // Age. 1992;15:51–54.
        
        // D = 0.05; // um^2/ms // Faas site:
        // Bloodgood BL, Sabatini BL. 
        // Neuronal activity regulates diffusion across the neck of dendritic spines. 
        // Science. 2005;310:866–869. doi: 10.1126/science.1114816. 310/5749/866
        // PAGFP D = 0.037 um^2/s // room temp // 28 kD
        // using Q10 = 1.3 for diffusion, D = 0.052 um^2/s at 35C
        // CaM is 17 kD
        
        if (lobeNorC.equalsIgnoreCase("N") || lobeNorC.equalsIgnoreCase("NL")) {
            name = "CM_NL_CaRCaR";
            K1on = 2 * 770;
            K1off = 160;
            // Kd1 = 0.1039 mM
            K2on = 32000;
            K2off = 2 * 22;
            // Kd2 = 0.00137 mM
        } else if (lobeNorC.equalsIgnoreCase("C") || lobeNorC.equalsIgnoreCase("CL")) {
            name = "CM_CL_CaRCaR";
            K1on = 2 * 84;
            K1off = 2.6;
            // Kd1 = 0.01548 mM
            K2on = 25;
            K2off = 2 * 0.0065;
            // Kd2 = 0.00052 mM
        }

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        reaction = true;
        iCaRCaR = -1;

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] == this) {
                iCaRCaR = i;
                break;
            }
        }

        if (iCaRCaR >= 0) {
            if (iCaRCaR == iCa) {
                error("init", "iCa", "value conflicts with iCaRCaR");
            }
            if (iCaRCaR == iCaTR) {
                error("init", "iCaTR", "value conflicts with iCaRCaR");
            }
        }

        if (!project.checkDiffusantNum(iCa)) {
            error("init", "iCa", "out of range");
        }

        if (!project.checkDiffusantNum(iCaTR)) {
            error("init", "iCaTR", "out of range");
        }
        
        computeC0();

    }
    
    public void Q10init() {
        Q10_D = new Q10(project, "D", D, temp, Q10.Q10_D);
        Q10_K1on = new Q10(project, "K1on", K1on, temp, Q10.Q10_K);
        Q10_K1off = new Q10(project, "K1off", K1off, temp, Q10.Q10_K);
        Q10_K2on = new Q10(project, "K2on", K2on, temp, Q10.Q10_K);
        Q10_K2off = new Q10(project, "K2off", K2off, temp, Q10.Q10_K);
    }
    
    @Override
    public void Q10execute() {
        double k1on, k1off, k2on, k2off;
        super.Q10execute();
        if ((Q10_K1on != null) && Q10_K1on.on) {
            Q10_K1on.temp = temp;
            k1on = Q10_K1on.getScaledValue();
            if (Double.isFinite(k1on)) {
                K1on = k1on;
            }
        }
        if ((Q10_K1off != null) && Q10_K1off.on) {
            Q10_K1off.temp = temp;
            k1off = Q10_K1off.getScaledValue();
            if (Double.isFinite(k1off)) {
                K1off = k1off;
            }
        }
        if ((Q10_K2on != null) && Q10_K2on.on) {
            Q10_K2on.temp = temp;
            k2on = Q10_K2on.getScaledValue();
            if (Double.isFinite(k2on)) {
                K2on = k2on;
            }
        }
        if ((Q10_K2off != null) && Q10_K2off.on) {
            Q10_K2off.temp = temp;
            k2off = Q10_K2off.getScaledValue();
            if (Double.isFinite(k2off)) {
                K2off = k2off;
            }
        }
        Kd1 = K1off / K1on;
        Kd2 = K2off / K2on;
    }

    public double computeC0() {

        double K1, K2;

        Diffusant Ca, CaTR;

        C0 = Double.NaN;

        if ((iCa >= 0) && (iCa < project.diffusants.length)) {
            if ((iCaTR >= 0) && (iCaTR < project.diffusants.length)) {

                Ca = project.diffusants[iCa];
                CaTR = project.diffusants[iCaTR];

                //r1 = Ca.C0 * K1on - K2off + (K2off / (Ca.C0 * K2on)) * (Ca.C0 * K1on + K1off + Ca.C0 * K2on);
                //C0 = Ctotal * Ca.C0 * K1on / r1;
                //CaTR.C0 = C0 * K2off / (Ca.C0 * K2on);
                Kd1 = K1off / K1on;
                Kd2 = K2off / K2on;

                K1 = K1on / K1off;
                K2 = K2on / K2off;

                CaTR.C0 = K1 * Ca.C0 * Ctotal / (1 + K1 * Ca.C0 + K1 * K2 * Ca.C0 * Ca.C0);
                C0 = K1 * K2 * Ca.C0 * Ca.C0 * Ctotal / (1 + K1 * Ca.C0 + K1 * K2 * Ca.C0 * Ca.C0);

            }
        }

        return C0;

    }

    @Override
    public void react(RunFiniteDifference fd, int iRCARCA) {

        int ipnt = fd.thisPnt;
        double Ca, dCa, TT, dTT, CaTR, dCaTR, CaRCaR, dCaRCaR;

        if (iRCARCA != iCaRCaR) {
            return;
        }
        
        if ((pulseTimer != null) && (pulseTimer.high(fd.it) == 0)) {
            return;
        }

        Ca = fd.diffus[iCa][ipnt];
        CaTR = fd.diffus[iCaTR][ipnt];
        CaRCaR = fd.diffus[iRCARCA][ipnt];
        TT = (Ctotal - CaTR - CaRCaR);

        dCaTR = (TT * Ca * K1on - CaTR * K1off - CaTR * Ca * K2on + CaRCaR * K2off) * project.dt;
        dTT = (CaTR * K1off - TT * Ca * K1on) * project.dt; // not used
        dCaRCaR = (CaTR * Ca * K2on - CaRCaR * K2off) * project.dt;

        dCa = -2 * dCaRCaR - dCaTR; // 2 dCaRCaR for 2 binding sites

        fd.diffnext[iCa][ipnt] += dCa;
        fd.diffnext[iCaTR][ipnt] += dCaTR;
        fd.diffnext[iRCARCA][ipnt] += dCaRCaR;

        if (fd.diffnext[iCa][ipnt] < 0) {
            error("react", "fd.diffnext[iCa]", "negative concentration for diffusant #" + iCa);
            fd.stopSimulation();
        }

        if (fd.diffnext[iCaTR][ipnt] < 0) {
            error("react", "fd.diffnext[iCaTR]", "negative concentration for diffusant #" + iCaTR);
            fd.stopSimulation();
        }

        if (fd.diffnext[iRCARCA][ipnt] < 0) {
            error("react", "fd.diffnext[iCaRCaR]", "negative concentration for diffusant #" + iRCARCA);
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
        if (n.equalsIgnoreCase("iCa")) {
            if (v < 0) {
                return false;
            }
            iCa = (int) v;
            computeC0();
            return true;
        }
        if (n.equalsIgnoreCase("iCaTR")) {
            if (v < 0) {
                return false;
            }
            iCaTR = (int) v;
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
