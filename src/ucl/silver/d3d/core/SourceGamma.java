package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.Utility;

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
public class SourceGamma extends Source {
    
    // Gamma distribution with parameters alpha and beta
    // https://en.wikipedia.org/wiki/Gamma_distribution
    // beta = 1 / tau
    // tau is set by pulse duration
    
    public int alpha = 2;
    
    public SourceGamma(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double onset, double tau, double ctotal) {
        super(p, NAME, DiffusantNum, c, onset, tau, ctotal);
        createVector(true); // sets ParamVector
    }
    
    public SourceGamma(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        super(p, NAME, DiffusantNum, c, pt);
        createVector(true); // sets ParamVector
    }
    
    public static double Ctotal2Ipeak(double cTotal_mM, double tau_ms, double charge, double volume_um3, int alpha) {
        double liters = Utility.litersPerUM3(volume_um3);
        double t = 1.0 / Utility.gamma((tau_ms * 1e-3), 0, alpha, 1.0 / (tau_ms * 1e-3)); // seconds // for Gamma
        double qmax = cTotal_mM * 1e-3 / t; // M/s
        double molespersecond = qmax * liters; // moles/s
        double amps = molespersecond * charge * Utility.FARADAY; // amps
        return amps * 1e12; // pA
    }

    public static double Ipeak2Ctotal(double iPeak_pA, double tau_ms, double charge, double volume_um3, int alpha) {
        double liters = Utility.litersPerUM3(volume_um3);
        double t = 1.0 / Utility.gamma((tau_ms * 1e-3), 0, alpha, 1.0 / (tau_ms * 1e-3)); // seconds // for Gamma
        double molespersecond = ((iPeak_pA * 1e-12) / (charge * Utility.FARADAY)); // moles/s
        double qmax = molespersecond / liters; // M/s
        double ctotal = qmax * 1e3 * t; // mM
        return ctotal;
    }
    
    @Override
    public void updateStats() {
        double charge, molesTotal, sumDuration;

        double litersPerVoxel = Utility.litersPerVoxel(project.dx);
        int numVoxels = numVoxels();

        charge = project.diffusants[diffusantNum].charge;

        C2I_conversionFactor = 1e12 * charge * Utility.FARADAY * litersPerVoxel / project.dt; // convert mM to pA

        if (pulseTimer == null) {
            error("updateSourceStats", "pulseTimer", "SourceGamma requires a pulse timer");
        } else {
            Ctotal = pulseTimer.sumAmplitude();
            //sumDuration = pulseTimer.sumDuration();
            molesTotal = (Ctotal * 1e-3) * numVoxels * litersPerVoxel;
            Ntotal = molesTotal * Utility.AVOGADRO; // molecules
            Nvoxel = Ntotal / numVoxels;
            //Qmax = Ctotal / (sumDuration * Math.sqrt(2 * Math.PI));
        }

    }
    
    @Override
    public void initPulseTimer() {
        
        double t, tonset, tau, ctotal, qmax;
        int itmax, ip = 0;

        if (pulseTimer == null) {
            error("initPulseTimer", "pulseTimer", "SourceGamma requires a pulse timer");
            return;
        }
        
        itmax = project.simPoints();

        if (itmax <= 0) {
            pulseTimer.timer = null;
            return;
        }

        pulseTimer.timer = new double[itmax];
        
        Master.log("computing Gamma pulse timer...");
        
        for (Pulse p : pulseTimer.pulses) {

            tonset = p.time;
            tau = p.duration;
            ctotal = p.amplitude;
            qmax = ctotal * Utility.gamma(tonset + tau, tonset, alpha, 1.0 / tau);
            Master.log("pulse " + ip++ + ", Ctotal = " + ctotal + " mM, Qmax = " + qmax + " mM/ms");
            
            for (int it = 0; it < itmax; it++) {
                t = it * project.dt;
                pulseTimer.timer[it] += ctotal * Utility.gamma(t, tonset, alpha, 1.0 / tau) * project.dt;
            }

        }
       
        updateVectors();

    }
    
    @Override
    public double release(RunFiniteDifference fd, Geometry geometry) {

        double cvoxel, cdelta = 0;

        if (spaceVoxels == 0) {
            return 0;
        }

        if (fd.it >= pulseTimer.timer.length) {
            return 0;
        }

        cvoxel = pulseTimer.timer[fd.it];

        if (fd.diffus[0].length == 1) {

            fd.diffus[diffusantNum][0] += cvoxel;
            cdelta += cvoxel;

        } else {

            //if (indexRFD == null) {
            //Master.exit("SourceGuass error: coordinates index array has not be initialized.");
            //}

            for (int i : indexRFD) {
                fd.diffus[diffusantNum][i] += cvoxel;
                cdelta += cvoxel;
            }

        }

        saveValue(cdelta);

        return cdelta;

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
