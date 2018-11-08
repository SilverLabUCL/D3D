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
public class DiffusantDye extends DiffusantReactant {
    
    // [dye] + [b] = [dye-b]
    //
    // cTotal = [dye] + [dye-b]
    //
    // d[dye-b]/dt = [dye][b]kon - [dye-b]koff
    // d[dye]/dt = d[b]/dt = -d[dye-b]/dt

    // this Diffusant represents product [dye-b]
    // [dye] is not saved, but can be computed as [dye] = cTotal - [dye-b]

    public double R; // fmax/fmin, for f/f0 computation
    public double f0;
    public double fmax;

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("f0")) {
            return false;
        }
        if (name.equalsIgnoreCase("fmax")) {
            return false;
        }
        return super.canEdit(name);
    }

    public DiffusantDye(Project p, String NAME, double CTOTAL, double diffusionConstant, CoordinatesVoxels c,
            int ReactantNum, double Kon, double Koff, double r) {

        super(p, NAME, CTOTAL, diffusionConstant, c, ReactantNum, Kon, Koff);

        R = r;

        createVector(true);

    }
    
    @Override
    public void init() {
        super.init();
        computeF0();
    }

    public double computeF0() {

        f0 = computeF(C0, "F");
        fmax = cTotal * R; // all dye is bound

        super.setParamObject("f0", f0);
        super.setParamObject("fmax", fmax);

        return f0;

    }
    
    public double computeF(double dye_bound, String norm) {
        
        double f;

        // convert [dye] to F
        // F = Fb + Fu ... (b)ound and (u)nbound
        // Fb = B[dye-b]
        // Fu = U[dye]
        // B/U = Fmax/Fmin = R
        
        double unboundDye = cTotal - dye_bound;

        f = dye_bound * R + unboundDye;
        
        if (norm.equalsIgnoreCase("dF/F0") || norm.equalsIgnoreCase("(F-F0)/F0")) {
            f = (f - f0) / f0;
        } else if (norm.equalsIgnoreCase("F/F0")) {
            f /= f0;
        } else if (norm.equalsIgnoreCase("F/Fmax")) {
            f /= fmax;
        }
        
        return f;

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

        if (n.equalsIgnoreCase("R")) {
            if (v < 0) {
                return false;
            }
            R = v;
            computeF0();
            updateVectors();
            return true;
        }
        return super.setMyParams(o, v);
    }

}
