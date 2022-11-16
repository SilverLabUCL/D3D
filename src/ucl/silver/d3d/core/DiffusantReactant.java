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
public class DiffusantReactant extends Diffusant {

    // [a] + [b] = [ab]
    //
    // Ctotal = [a] + [ab]
    //
    // d[ab]/dt = [a][b]kon - [ab]koff
    // d[a]/dt = d[b]/dt = -d[ab]/dt
    //
    // this Diffusant represents [ab]
    // [a] is not saved, but can be computed as [a] = Ctotal - [ab]
    //
    public String bDiffusantName = null; // diffusant name for [b]
    public int bDiffusantNum = -1; // diffusant number for [b]
    public int abDiffusantNum; // this diffusant

    public double Ctotal; // [a] + [ab] (mM)

    public double Kon; // on reaction rate (1/(mM*ms))
    public double Koff; // off reaction rate (1/ms)
    public double Kd; // mM // Kd = Koff / Kon;
    
    public Q10 Q10_Kon = null; // Q10 temperature scaling for Kon
    public Q10 Q10_Koff = null; // Q10 temperature scaling for Koff

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("Ctotal")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("Kon")) {
            return "1/" + project.concUnits + "*" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("Koff")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("Kd")) {
            return project.concUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("abDiffusantNum")) {
            return false;
        }
        return super.canEdit(name);
    }

    public DiffusantReactant(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int DiffusantNumB, double KON, double KOFF) {

        super(p, NAME, 0, diffusionConstant, c);

        Ctotal = CTOTAL;
        bDiffusantNum = DiffusantNumB;
        Kon = KON;
        Koff = KOFF;
        Kd = Koff / Kon;
        reaction = true;

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        reaction = true;
        abDiffusantNum = -1;

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] == this) {
                abDiffusantNum = i;
                break;
            }
        }

        if (!project.checkDiffusantNum(bDiffusantNum)) {
            error("init", "bDiffusantNum", "out of range");
            return;
        }
        
        bDiffusantName = project.diffusantName(bDiffusantNum);

        if (abDiffusantNum >= 0) {
            if (abDiffusantNum == bDiffusantNum) {
                error("init", "bDiffusantNum", "value conflicts with abDiffusantNum");
                return;
            }
        }

        computeC0();

    }
    
    @Override
    public void Q10execute() {
        double kon, koff;
        super.Q10execute();
        if ((Q10_Kon != null) && Q10_Kon.on) {
            kon = Q10_Kon.getScaledValue();
            if (Double.isFinite(kon)) {
                Kon = kon;
            }
        }
        if ((Q10_Koff != null) && Q10_Koff.on) {
            koff = Q10_Koff.getScaledValue();
            if (Double.isFinite(koff)) {
                Koff = koff;
            }
        }
        Kd = Koff / Kon;
    }
    
    public void Q10init(double temp) {
        Q10_D = new Q10(project, "D", D, temp, Q10.Q10_D);
        Q10_Kon = new Q10(project, "Kon", Kon, temp, Q10.Q10_K);
        Q10_Koff = new Q10(project, "Koff", Koff, temp, Q10.Q10_K);
    }

    public double computeC0() {

        Diffusant b;

        C0 = Double.NaN;

        if ((bDiffusantNum >= 0) && (bDiffusantNum < project.diffusants.length)) {

            b = project.diffusants[bDiffusantNum];

            Kd = Koff / Kon;

            C0 = Ctotal * b.C0 / (b.C0 + Kd); // initial value of [ab]

            if (C0 < 0) {
                error("computeC0", "C0", "out of range: " + C0);
            }

        }

        super.setParamObject("C0", C0);

        return C0;

    }

    public boolean SetDiffusantNumB(int diffusantNum) {

        Diffusant d = project.getDiffusant(diffusantNum);

        if (d != null) {
            bDiffusantNum = diffusantNum;
            bDiffusantName = d.name;
            init();
            updateVectors();
            return true;
        } else {
            Master.log("DiffusantReactant SetDiffusantNumB error: failed to find diffusant #" + diffusantNum);
            updateVectors();
            return false;
        }

    }

    public boolean SetDiffusantNameB(String diffusantName) {

        int i = project.getDiffusantNum(diffusantName);

        if (i >= 0) {
            bDiffusantName = diffusantName;
            bDiffusantNum = i;
            init();
            updateVectors();
            return true;
        } else {
            Master.log("DiffusantReactant SetDiffusantNameB error: failed to find diffusant " + diffusantName);
            updateVectors();
            return false;
        }

    }

    @Override
    public void react(RunFiniteDifference fd, int abNum) {

        int ipnt = fd.thisPnt;
        double a, b, ab, d_b, d_ab;

        if (abNum != abDiffusantNum) {
            return;
        }
        
        if ((pulseTimer != null) && (pulseTimer.high(fd.it) == 0)) {
            return;
        }

        ab = fd.diffus[abNum][ipnt];
        b = fd.diffus[bDiffusantNum][ipnt];
        a = Ctotal - ab;

        d_ab = (a * b * Kon - ab * Koff) * project.dt;
        d_b = -d_ab;
        
        fd.diffnext[abNum][ipnt] += d_ab;
        fd.diffnext[bDiffusantNum][ipnt] += d_b;

        if (fd.diffnext[abNum][ipnt] < 0) {
            error("react", "fd.diffnext[abDiffusantNum]", "negative concentration for diffusant #" + abNum);
            fd.stopSimulation();
        }

        if (fd.diffnext[bDiffusantNum][ipnt] < 0) {
            error("react", "fd.diffnext[bDiffusantNum]", "negative concentration for diffusant #" + bDiffusantNum);
            fd.stopSimulation();
        }

    }
    
    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);
        
        if (Q10_Kon != null) {
            Q10_Kon.addUser(pv);
        }
        
        if (Q10_Koff != null) {
            Q10_Koff.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }
        
        if (Q10_Kon != null) {
            addBlankParam();
            Q10_Kon.createVector(true);
            addVector(Q10_Kon.getVector());
            Q10_Kon.addUser(this);
        }
        
        if (Q10_Koff != null) {
            addBlankParam();
            Q10_Koff.createVector(true);
            addVector(Q10_Koff.getVector());
            Q10_Koff.addUser(this);
        }

        if (close){
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);
        
        if (Q10_Kon != null) {
            Q10_Kon.updateVector(v);
        }
        
        if (Q10_Koff != null) {
            Q10_Koff.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantReactant)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("Ctotal")) {
            if (v < 0) {
                return false;
            }
            Ctotal = v;
            init();
            updateVectors();
            return true;
        }
        if (n.equalsIgnoreCase("bDiffusantNum")) {
            if (v < 0) {
                return false;
            }
            return SetDiffusantNumB((int) v);
        }
        if (n.equalsIgnoreCase("Kon")) {
            
            if (v < 0) {
                return false;
            }
            
            Kon = v;
            //Kd = Koff / Kon;
            Koff = Kon * Kd;
            
            if (Q10_Kon != null) {
                Q10_Kon.value = v;
            }
            
            if (Q10_Koff != null) {
                Q10_Koff.value = v;
            }
            
            init();
            updateVectors();
            return true;
            
        }
        if (n.equalsIgnoreCase("Koff")) {
            
            if (v < 0) {
                return false;
            }
            
            Koff = v;
            //Kd = Koff / Kon;
            Kon = Koff / Kd;
            
            if (Q10_Kon != null) {
                Q10_Kon.value = v;
            }
            
            if (Q10_Koff != null) {
                Q10_Koff.value = v;
            }
            
            init();
            updateVectors();
            return true;
            
        }
        if (n.equalsIgnoreCase("Kd")) {
            
            if (v < 0) {
                return false;
            }
            
            Kd = v;
            Kon = Koff / Kd;
            //Koff = Kon * Kd;
            
            if (Q10_Kon != null) {
                Q10_Kon.value = v;
            }
            
            init();
            updateVectors();
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

        if (!(o.paramVector instanceof DiffusantReactant)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("bDiffusantName")) {
            return SetDiffusantNameB(s);
        }
        return super.setMyParams(o, s);
    }

}
