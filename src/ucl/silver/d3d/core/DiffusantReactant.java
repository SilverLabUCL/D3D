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
public class DiffusantReactant extends Diffusant {

    // [a] + [b] = [ab]
    //
    // cTotal = [a] + [ab]
    //
    // d[ab]/dt = [a][b]kon - [ab]koff
    // d[a]/dt = d[b]/dt = -d[ab]/dt
    //
    // this Diffusant represents [ab]
    // [a] is not saved, but can be computed as [a] = cTotal - [ab]
    //
    public String bDiffusantName = null; // diffusant name for [b]
    public int bDiffusantNum = -1; // diffusant number for [b]
    public int abDiffusantNum; // this diffusant

    public double cTotal; // [a] + [ab] (mM)

    public double kon; // on reaction rate (1/(mM*ms))
    public double koff; // off reaction rate (1/ms)
    public double kd; // mM

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("cTotal")) {
            return project.concUnits;
        }
        if (name.equalsIgnoreCase("kon")) {
            return "1/" + project.concUnits + "*" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("koff")) {
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("kd")) {
            return project.concUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("abDiffusantNum")) {
            return false;
        }
        if (name.equalsIgnoreCase("kd")) {
            return false;
        }
        return super.canEdit(name);
    }

    public DiffusantReactant(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int DiffusantNumB, double Kon, double Koff) {

        super(p, NAME, 0, diffusionConstant, c);

        cTotal = CTOTAL;
        bDiffusantNum = DiffusantNumB;
        kon = Kon;
        koff = Koff;
        kd = koff / kon;
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

    public double computeC0() {

        Diffusant b;

        C0 = Double.NaN;

        if ((bDiffusantNum >= 0) && (bDiffusantNum < project.diffusants.length)) {

            b = project.diffusants[bDiffusantNum];

            kd = koff / kon;

            C0 = cTotal * b.C0 / (b.C0 + kd); // initial value of [ab]

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
        a = cTotal - ab;

        d_ab = (a * b * kon - ab * koff) * project.dt;
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

        if (n.equalsIgnoreCase("cTotal")) {
            if (v < 0) {
                return false;
            }
            cTotal = v;
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
        if (n.equalsIgnoreCase("kon")) {
            if (v < 0) {
                return false;
            }
            kon = v;
            kd = koff / kon;
            init();
            updateVectors();
            return true;
        }
        if (n.equalsIgnoreCase("koff")) {
            if (v < 0) {
                return false;
            }
            koff = v;
            kd = koff / kon;
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
