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
public class DiffusantNOR2 extends Diffusant {

    // NO receptor (NOR)
    // [NO] + [R] = [NOR1] = [NOR2]   (receptor can be in state 1 or 2)
    //
    // [a] + [b] = [ab] = [abc]
    //
    // cTotal = [a] + [ab] + [abc]
    //
    // d[ab]/dt = [a][b]k1on - [ab]k1off - [ab]k2on + [abc]k2off
    // d[abc]/dt = [ab]k2on - [abc]k2off
    // d[a]/dt = [ab]k1off - [a][b]k1on
    //
    // this Diffusant represents [abc]
    // [a] is not saved, but can be computed as [a] = cTotal - [ab] - [abc]
    //
    public String bReactantName = null; // diffusant name for reactant [b]
    public int bReactantNum = -1; // diffusant number for reactant [b]

    public String abReactantName = null; // diffusant name for reactant [ab]
    public int abReactantNum = -1; // diffusant number for reactant [ab]

    public double cTotal; // total concentration of a, cTotal = [a] + [ab] + [abc] (mM)

    public double k1on;
    public double k1off;
    public double k2on;
    public double k2off;

    private int thisDiffusantNum;

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
            return "1/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("k2off")) {
            return "1/" + project.timeUnits;
        }
        return super.units(name);
    }

    public DiffusantNOR2(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int BReactantNum, int ABReactantNum, double K1ON, double K1OFF, double K2ON, double K2OFF) {

        super(p, NAME, 0, diffusionConstant, c);

        bReactantNum = BReactantNum;
        abReactantNum = ABReactantNum;
        cTotal = CTOTAL;
        k1on = K1ON;
        k1off = K1OFF;
        k2on = K2ON;
        k2off = K2OFF;
        reaction = true;

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        thisDiffusantNum = -1;

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] == this) {
                thisDiffusantNum = i;
                break;
            }
        }

        if (!project.checkDiffusantNum(bReactantNum)) {
            error("init", "bReactantNum", "out of range");
        }

        if (!project.checkDiffusantNum(abReactantNum)) {
            error("init", "productNum", "out of range");
        }

        if (thisDiffusantNum >= 0) {
            if (thisDiffusantNum == bReactantNum) {
                Master.exit("conflicting diffusant number for reactantNum: " + this);
            }
            if (thisDiffusantNum == abReactantNum) {
                Master.exit("conflicting diffusant number for productNum: " + this);
            }
        }

        bReactantName = project.diffusantName(bReactantNum);
        abReactantName = project.diffusantName(abReactantNum);

        computeC0();

    }

    public double computeC0() {

        double x0, x1, x2 = Double.NaN;

        Diffusant b, ab;

        C0 = Double.NaN;

        if ((bReactantNum >= 0) && (bReactantNum < project.diffusants.length)) {
            if ((abReactantNum >= 0) && (abReactantNum < project.diffusants.length)) {

                b = project.diffusants[bReactantNum];
                ab = project.diffusants[abReactantNum];

                x0 = cTotal / (1 + (k1on * b.C0 / k1off) + (k1on * k2on * b.C0 / (k1off * k2off)));
                x1 = k1on * b.C0 * x0 / k1off;
                x2 = k2on * x1 / k2off;

                ab.C0 = x1;
                C0 = x2;

                ab.setParamObject("C0", x1);

            }
        }

        super.setParamObject("C0", x2);

        return C0;

    }

    @Override
    public void react(RunFiniteDifference fd, int abcDiffusant) {

        int ipnt = fd.thisPnt;
        double a, b, ab, abc, da, db, dab, dabc;

        if (!reaction || (abcDiffusant != thisDiffusantNum)) {
            return;
        }
        
        if ((pulseTimer != null) && (pulseTimer.high(fd.it) == 0)) {
            return;
        }

        b = fd.diffus[bReactantNum][ipnt];
        ab = fd.diffus[abReactantNum][ipnt];
        abc = fd.diffus[abcDiffusant][ipnt];
        a = (cTotal - ab - abc);

        dab = (a * b * k1on - ab * k1off - ab * k2on + abc * k2off) * project.dt;
        da = (ab * k1off - a * b * k1on) * project.dt;
        dabc = (ab * k2on - abc * k2off) * project.dt;
        //db = da;

        fd.diffnext[bReactantNum][ipnt] += da; // db = da
        fd.diffnext[abReactantNum][ipnt] += dab;
        fd.diffnext[abcDiffusant][ipnt] += dabc;

        if (fd.diffnext[bReactantNum][ipnt] < 0) {
            Master.log("reaction error: negative concentration, diffusant #"
                    + bReactantNum);
        }

        if (fd.diffnext[abReactantNum][ipnt] < 0) {
            Master.log("reaction error: negative concentration, diffusant #"
                    + abReactantNum);
        }

        if (fd.diffnext[abcDiffusant][ipnt] < 0) {
            Master.log("reaction error: negative concentration, diffusant #"
                    + abcDiffusant);
        }

    }

    public boolean SetReactantNumB(int diffusantNum) {

        Diffusant d = project.getDiffusant(diffusantNum);

        if (d != null) {
            bReactantNum = diffusantNum;
            bReactantName = d.name;
            updateVectors();
            return true;
        } else {
            Master.log("DiffusantReactant2 SetReactantNumB error: failed to find diffusant #" + diffusantNum);
            updateVectors();
            return false;
        }

    }

    public boolean SetReactantNameB(String diffusantName) {

        int i = project.getDiffusantNum(diffusantName);

        if (i >= 0) {
            bReactantName = diffusantName;
            bReactantNum = i;
            updateVectors();
            return true;
        } else {
            Master.log("DiffusantReactant2 SetReactantNameB error: failed to find diffusant " + diffusantName);
            updateVectors();
            return false;
        }

    }

    public boolean SetProductNum(int diffusantNum) {

        Diffusant d = project.getDiffusant(diffusantNum);

        if (d != null) {
            abReactantNum = diffusantNum;
            abReactantName = d.name;
            updateVectors();
            return true;
        } else {
            Master.log("DiffusantReactant2 SetProductNum error: failed to find diffusant #" + diffusantNum);
            updateVectors();
            return false;
        }

    }

    public boolean SetProductName(String diffusantName) {

        int i = project.getDiffusantNum(diffusantName);

        if (i >= 0) {
            abReactantName = diffusantName;
            abReactantNum = i;
            updateVectors();
            return true;
        } else {
            Master.log("DiffusantReactant2 SetProductName error: failed to find diffusant " + diffusantName);
            updateVectors();
            return false;
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantNOR2)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("cTotal")) {
            if (v < 0) {
                return false;
            }
            cTotal = v;
            return true;
        }
        if (n.equalsIgnoreCase("bReactantNum")) {
            if (v < 0) {
                return false;
            }
            return SetReactantNumB((int) v);
        }
        if (n.equalsIgnoreCase("productNum")) {
            if (v < 0) {
                return false;
            }
            return SetProductNum((int) v);
        }
        if (n.equalsIgnoreCase("k1on")) {
            if (v < 0) {
                return false;
            }
            k1on = v;
            return true;
        }
        if (n.equalsIgnoreCase("k1off")) {
            if (v < 0) {
                return false;
            }
            k1off = v;
            return true;
        }
        if (n.equalsIgnoreCase("k2on")) {
            if (v < 0) {
                return false;
            }
            k2on = v;
            return true;
        }
        if (n.equalsIgnoreCase("k2off")) {
            if (v < 0) {
                return false;
            }
            k2off = v;
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

        if (!(o.paramVector instanceof DiffusantNOR2)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("bReactantName")) {
            return SetReactantNameB(s);
        }
        if (n.equalsIgnoreCase("productName")) {
            return SetProductName(s);
        }
        return super.setMyParams(o, s);
    }

}
