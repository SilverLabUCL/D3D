package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;

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
public class PSFgauss extends PSF {

    public double xSTDV, ySTDV, zSTDV; // Gaussian standard deviations
    public double xFWHM, yFWHM, zFWHM; // Gaussian full-width half-max

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("xSTDV")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("ySTDV")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("zSTDV")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("xFWHM")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("yFWHM")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("zFWHM")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    public PSFgauss(Project p, CoordinatesVoxels c, double STDVx, double STDVy, double STDVz) {
        super(p, c);
        xSTDV = STDVx;
        ySTDV = STDVy;
        zSTDV = STDVz;
        xFWHM = Utility.gaussSTDV2FWHM(xSTDV);
        yFWHM = Utility.gaussSTDV2FWHM(ySTDV);
        zFWHM = Utility.gaussSTDV2FWHM(zSTDV);
        //xySymmetric = true; // should manually set to true
        name = "Gaussian PSF";
        createVector(true);
    }

    @Override
    public double computeVoxel(double xVoxel, double yVoxel, double zVoxel) {

        double x, y, z;
        double v1 = 1, v2 = 1, v3 = 1;
        double e1 = 1, e2 = 1, e3 = 1;

        dx = project.dx;
        
        xVoxel -= xVoxelCenter();
        yVoxel -= yVoxelCenter();
        zVoxel -= zVoxelCenter();

        rotateAll(xVoxel, yVoxel, zVoxel); // results saved in irotated, jrotated, krotated

        x = iRotated() * dx; // um
        y = jRotated() * dx; // um
        z = kRotated() * dx; // um

        if (xSTDV > 0) {
            //v1 = Math.pow(x, 2.0);
            //e1 = Math.exp(-v1 / (2.0 * Math.pow(xSTDV, 2.0)));
            e1 = Utility.gaussNormPeak(x, 0, xSTDV);
        }

        if (ySTDV > 0) {
            //v2 = Math.pow(y, 2.0);
            //e2 = Math.exp(-v2 / (2.0 * Math.pow(ySTDV, 2.0)));
            e2 = Utility.gaussNormPeak(y, 0, ySTDV);
        }

        if (zSTDV > 0) {
            //v3 = Math.pow(z, 2.0);
            //e3 = Math.exp(-v3 / (2.0 * Math.pow(zSTDV, 2.0)));
            e3 = Utility.gaussNormPeak(z, 0, zSTDV);
        }

        return e1 * e2 * e3;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }
        
        if (!(o.paramVector instanceof PSFgauss)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("xSTDV")) {
            xSTDV = v;
            xFWHM = Utility.gaussSTDV2FWHM(xSTDV);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ySTDV")) {
            ySTDV = v;
            yFWHM = Utility.gaussSTDV2FWHM(ySTDV);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("zSTDV")) {
            zSTDV = v;
            zFWHM = Utility.gaussSTDV2FWHM(zSTDV);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("xFWHM")) {
            xFWHM = v;
            xSTDV = Utility.gaussFWHM2STDV(xFWHM);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("yFWHM")) {
            yFWHM = v;
            ySTDV = Utility.gaussFWHM2STDV(yFWHM);
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("zFWHM")) {
            zFWHM = v;
            zSTDV = Utility.gaussFWHM2STDV(zFWHM);
            array = null;
            return true;
        }
        return super.setMyParams(o, v);
    }
}
