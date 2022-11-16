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
public class PSFconfocal extends PSFgauss {

    // Angus Silver Lab one-photon confocal rig (Dec 2010)
    // illumination wavelength = 488 nm
    // detection wavelength = 527 nm
    // 0.4 um pinhole (40 um / 100 magnification)
    public final static double NA_MIN = 0.7; // polynomial fit lower limit
    public final static double NA_MAX = 1.0; // polynomial fit upper limit

    public double numericalAperture = 1.0;

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("NA_MIN")) {
            return false;
        }
        if (name.equalsIgnoreCase("NA_MAX")) {
            return false;
        }
        if (name.equalsIgnoreCase("xSTDV")) {
            return false;
        }
        if (name.equalsIgnoreCase("ySTDV")) {
            return false;
        }
        if (name.equalsIgnoreCase("zSTDV")) {
            return false;
        }
        if (name.equalsIgnoreCase("xFWHM")) {
            return false;
        }
        if (name.equalsIgnoreCase("yFWHM")) {
            return false;
        }
        if (name.equalsIgnoreCase("zFWHM")) {
            return false;
        }
        return super.canEdit(name);
    }

    public PSFconfocal(Project p, CoordinatesVoxels c, double NA) {
        super(p, c, 0.1, 0.1, 0.1);
        setNumericalAperture(NA);
        xySymmetric = true;
        name = "Confocal PSF";
        createVector(true);
    }

    public final boolean setNumericalAperture(double NA) {

        if (NA < NA_MIN) {
            Master.log("PSFconfocal error: NA is only defined for values greater than " + NA_MIN);
            return false; // not defined
        }

        if (NA > NA_MAX) {
            Master.log("PSFconfocal error: NA is only defined for values less than " + NA_MAX);
            return false; // not defined
        }

        numericalAperture = NA;

        xFWHM = cPSF_FWHM_XY(NA);
        yFWHM = cPSF_FWHM_XY(NA);
        zFWHM = cPSF_FWHM_Z(NA);

        xSTDV = Utility.gaussFWHM2STDV(xFWHM);
        ySTDV = Utility.gaussFWHM2STDV(yFWHM);
        zSTDV = Utility.gaussFWHM2STDV(zFWHM);

        array = null;

        updateVectors();

        return true;

    }

    public static double cPSF_FWHM_XY(double NA) {

        double p0 = 0.6545; // polynomial fit 0.7 <= NA <= 1.0
        double p1 = -0.74776;
        double p2 = 0.30224;

        // Angus Silver Lab one-photon confocal rig (Dec 2010)
        // illumination wavelength = 488 nm
        // detection wavelength = 527 nm
        // 40/100 um pinhole
        return p0 + p1 * NA + p2 * NA * NA;

    }

    public static double cPSF_FWHM_Z(double NA) {

        double p0 = 7.3019; // polynomial fit 0.7 <= NA <= 1.0
        double p1 = -11.717;
        double p2 = 5.1535;

        // Angus Silver Lab one-photon confocal rig (Dec 2010)
        // illumination wavelength = 488 nm
        // detection wavelength = 527 nm
        // 40/100 um pinhole
        return p0 + p1 * NA + p2 * NA * NA;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }
        
        if (!(o.paramVector instanceof PSFconfocal)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("numericalAperture")) {
            return setNumericalAperture(v);
        }
        return super.setMyParams(o, v);
    }

}
