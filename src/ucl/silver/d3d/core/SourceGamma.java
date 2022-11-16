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
public class SourceGamma extends SourceExp {
    
    // initSelect = "C", using Ctotal (mM)
    
    // Gamma distribution (PDF) with parameters alpha and beta
    // https://en.wikipedia.org/wiki/Gamma_distribution

    // Requires PulseTimer
    // Pulse time = onset time (ms)
    // Pulse duration = tau = 1/beta (ms)
    // Pulse amplitude = Ctotal (mM)
 
    // Q (mM/ms) = Ctotal * Utility.gamma(t, tonset, alpha, beta)
    // integral Q = Ctotal (mM)
    // Qpeak = Ctotal * Utility.gamma(timeOfPeak, 0, alpha, beta);
    // timeOfPeak = (alpha - 1) / beta
    
    public int alpha = 2; // alpha > 0
    
    public SourceGamma(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double onset, double tau, double ctotal) {
        super(p, NAME, DiffusantNum, c, onset, tau, ctotal);
        createVector(true); // sets ParamVector
    }
    
    public SourceGamma(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        super(p, NAME, DiffusantNum, c, pt);
        createVector(true); // sets ParamVector
    }
    
    @Override
    public double q(double t, double tonset, double tau, double ctotal) {
        if (integral) {
            return ctotal * Utility.gammaCDF(t, tonset, alpha, 1.0 / tau); // mM
        } else {
            return ctotal * Utility.gammaPDF(t, tonset, alpha, 1.0 / tau); // mM/ms
        }
    }
    
    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }
  
        String n = o.getName();

        if (!(o.paramVector instanceof SourceGamma)) {
            return false;
        }

        if (n.equalsIgnoreCase("alpha")) {
            if (v > 0) {
                alpha = (int) v;
                return true;
            }
        }
        return super.setMyParams(o, v);
    }
    
}
