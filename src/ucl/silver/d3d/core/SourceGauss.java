package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.Utility;

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
public class SourceGauss extends SourceExp {
    //
    // initSelect = "C", using Ctotal (mM)
    //
    // Requires PulseTimer
    // Pulse time = Gauss peak (ms)
    // Pulse duration = Gauss stdv (ms)
    // Pulse amplitude = Ctotal (mM)
    //
    // Q (mM/ms) = Ctotal * Utility.gaussPDF(t, tpeak, stdv)
    // integral Q = Ctotal (mM)
    // Qpeak = Ctotal / (stdv * SQRT2PI)
    //
    
    public SourceGauss(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double tpeak_ms, double stdv_ms, double cTotal) {
        super(p, NAME, DiffusantNum, c, tpeak_ms, stdv_ms, cTotal);
        createVector(true); // sets ParamVector
    }

    public SourceGauss(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer PT_Ctotal) {
        super(p, NAME, DiffusantNum, c, PT_Ctotal);
        createVector(true); // sets ParamVector
    }

    @Override
    public double q(double t, double tpeak, double stdv, double ctotal) {
        if (integral) {
            return ctotal * Utility.gaussCDF(t, tpeak, stdv); // mM
        } else {
            return ctotal * Utility.gaussPDF(t, tpeak, stdv); // mM/ms
        }   
    }

}
