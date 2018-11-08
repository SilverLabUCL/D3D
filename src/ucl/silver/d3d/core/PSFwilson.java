package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;

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
public class PSFwilson extends PSF {

    public int integrationSteps = 1000; // numerical integration steps

    public double numericalAperture = 1.0;
    public double waveLength = 0.488; // light wavelength
    public double refractiveIndex = 1.338; // refractive index of H20 at 37C

    public double rscale = 1.0;

    private double sinAlpha, sinHalfAlphaSqr, alpha, dTheta, kindex;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("waveLength")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    public PSFwilson(Project p, CoordinatesVoxels c) {
        super(p, c);
        //xySymmetric = true; // should manually set to true
        name = "Wilson PSF";
        createVector(true);
    }

    @Override
    public void init() {

        sinAlpha = numericalAperture / refractiveIndex;
        alpha = Math.asin(numericalAperture / refractiveIndex);
        sinHalfAlphaSqr = Math.pow(Math.sin(alpha / 2.0), 2);
        kindex = 2.0 * Math.PI * refractiveIndex / waveLength;
        dTheta = 1.0 / integrationSteps;

        super.init();

    }

    @Override
    public double computeVoxel(double i, double j, double k) {

        double x, y, z, r, v, uhalf;
        double w, theta, J, C, S, a;
        double Re = 0, Im = 0;

        dx = project.dx;

        i -= xVoxelCenter();
        j -= yVoxelCenter();
        k -= zVoxelCenter();

        rotateAll(i, j, k); // results saved in irotated, jrotated, krotated

        x = iRotated() * dx; // um
        y = jRotated() * dx; // um
        z = kRotated() * dx; // um

        r = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2)) * rscale;

        v = kindex * r * sinAlpha;
        uhalf = 0.5 * 4 * kindex * z * sinHalfAlphaSqr;

        // Simpson's Rule

        for (int itheta = 0; itheta <= integrationSteps; itheta += 1) {

            if ((itheta == 0) || (itheta == integrationSteps)) {
                w = 1;
            } else if (itheta % 2 == 0) {
                w = 2;
            } else {
                w = 4;
            }

            theta = itheta * dTheta;

            J = MoreMath.j0(v * theta);
            a = uhalf * theta * theta;
            C = Math.cos(a);
            S = Math.sin(a);

            Re += J * C * w * theta;
            Im += J * S * w * theta;

        }

        Re *= dTheta / 3.0;
        Im *= dTheta / 3.0;

        return 4.0 * (Im * Im + Re * Re); // || h(u,v) ||^2

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }
        
        if (!(o.paramVector instanceof PSFwilson)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("integrationSteps")) {
            if (v < 5) {
                return false;
            }
            integrationSteps = (int) v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("numericalAperture")) {
            if (v < 0) {
                return false;
            }
            numericalAperture = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("waveLength")) {
            if (v < 0) {
                return false;
            }
            waveLength = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("refractiveIndex")) {
            if (v < 0) {
                return false;
            }
            refractiveIndex = v;
            array = null;
            return true;
        }
        return super.setMyParams(o, v);
    }
}
