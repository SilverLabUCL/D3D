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
public class SourceExp extends Source {
    //
    // initSelect = "C", using Ctotal (mM)
    //
    // Requires PulseTimer
    // Pulse time = onset (ms)
    // Pulse duration = exp decay tau (ms)
    // Pulse amplitude = Ctotal (mM)
    //
    // Q (mM/ms) = ( Ctotal / tau ) * Math.exp(-(t - tonset) / tau ))
    // integral Q = Ctotal (mM)
    // Qpeak = Ctotal / tau (mM/ms)
    //
    
    public boolean integral = false; // source is integral of Exp function, can use when clamp = true
    
    public SourceExp(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double tpeak, double tau, double ctotal) {
        super(p, NAME, DiffusantNum, c, tpeak, tau, ctotal, "C");
        createVector(true); // sets ParamVector
    }

    public SourceExp(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        super(p, NAME, DiffusantNum, c, pt, "C");
        createVector(true); // sets ParamVector
    }

    @Override
    public void updateStats() {

        double molesTotal, numVoxels;
        double liters = liters();

        Q = Double.NaN;

        if (clamp) {

            Ctotal = Double.NaN;
            Cvoxel = Double.NaN;
            Ntotal = Double.NaN;
            Nvoxel = Double.NaN;

            updateVectors();

            return;

        }

        numVoxels = numVoxels();

        Ctotal = 0;
        Ntotal = 0;

        for (Pulse p : pulseTimer.pulses) {
            molesTotal = (p.amplitude * 1e-3) * liters;
            Ctotal += p.amplitude; // mM
            Ntotal += molesTotal * Utility.AVOGADRO; // molecules
        }

        Cvoxel = Ctotal / numVoxels;
        Nvoxel = Ntotal / numVoxels;

        updateVectors();

    }

    @Override
    public void initPulseTimer() {

        double q, cvoxel_dt;
        int itmax;

        double nvoxels = numVoxels();

        if (pulseTimer == null) {
            error("initPulseTimer", "pulseTimer", "no timer");
            return;
        }

        pulseTimer.timer = null;

        if (pulseTimer.pulses == null) {
            return;
        }

        if (!pulseTimer.useTimerArray) {
            return;
        }

        itmax = project.simPoints();

        if (itmax <= 0) {
            return;
        }

        pulseTimer.timer = new double[itmax];

        Master.log("computing Exponential pulse timer...");

        for (Pulse p : pulseTimer.pulses) {
            for (int it = 0; it < pulseTimer.timer.length; it++) {
                
                q = q(it * project.dt, p.time, p.duration, p.amplitude); // mM/ms
                
                if (clamp) {
                    cvoxel_dt = q; // mM // see integral parameter
                } else {
                    cvoxel_dt = q * project.dt / nvoxels; // mM
                }
                
                pulseTimer.timer[it] += cvoxel_dt;
                
            }
        }

        Q = 0;

        for (int it = 0; it < pulseTimer.timer.length; it++) {
            q = pulseTimer.timer[it] * nvoxels / project.dt;
            Q = Math.max(Q, q);
        }

    }

    public double q(double t, double tonset, double tau, double ctotal) {
        if (integral) {
            return ctotal * (1 - Math.exp(-(t - tonset) / tau)); // mM
        } else {
            return ctotal * Utility.expNormArea(t, tonset, tau); // mM/ms
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

        if (n.equalsIgnoreCase("integral")) {
            integral = (v == 1);
            return true;
        }
        return super.setMyParams(o, v);
    }

}
