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
public class SourceGauss extends Source {

    // pulse tpeak is time of Gauss peak (ms)
    // pulse duration = Gauss stdv (ms)
    // pulse amplitude is Ctotal (mM)
    // Q = [ Ctotal / (stdv * SQRT2PI) ] * Math.exp(-Math.pow((t - tpeak) / (SQRT2 * stdv), 2)) (mM/ms)
    // integral Q = Ctotal (mM)
    public SourceGauss(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double tpeak, double stdv, double cTotal) {
        super(p, NAME, DiffusantNum, c, tpeak, stdv, cTotal);
        createVector(true); // sets ParamVector
    }

    public SourceGauss(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        super(p, NAME, DiffusantNum, c, pt);
        createVector(true); // sets ParamVector
    }
    
    public static double Ctotal2Ipeak(double cTotal_mM, double stdv_ms, double charge, double volume_um3) {
        double liters = Utility.litersPerUM3(volume_um3);
        double t = (stdv_ms * 1e-3) * Math.sqrt(2 * Math.PI); // seconds // for Gaussian
        double qmax = cTotal_mM * 1e-3 / t; // M/s
        double molespersecond = qmax * liters; // moles/s
        double amps = molespersecond * charge * Utility.FARADAY; // amps
        return amps * 1e12; // pA
    }

    public static double Ipeak2Ctotal(double iPeak_pA, double stdv_ms, double charge, double volume_um3) {
        double liters = Utility.litersPerUM3(volume_um3);
        double t = (stdv_ms * 1e-3) * Math.sqrt(2 * Math.PI); // seconds // for Gaussian
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
            error("updateSourceStats", "pulseTimer", "SourceGauss requires a pulse timer");
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

        double t, tpeak, stdv, ctotal, qmax;
        int itmax, ip = 0;

        if (pulseTimer == null) {
            error("initPulseTimer", "pulseTimer", "SourceGauss requires a pulse timer");
            return;
        }

        itmax = project.simPoints();

        if (itmax <= 0) {
            pulseTimer.timer = null;
            return;
        }

        pulseTimer.timer = new double[itmax];

        Master.log("computing Gaussian pulse timer...");
        
        for (Pulse p : pulseTimer.pulses) {
            
            tpeak = p.time;
            stdv = p.duration;
            ctotal = p.amplitude;
            qmax = ctotal / (stdv * Math.sqrt(2 * Math.PI));
            Master.log("pulse " + ip++ + ", Ctotal = " + ctotal + " mM, Qmax = " + qmax + " mM/ms");
            
            for (int it = 0; it < pulseTimer.timer.length; it++) {
                t = it * project.dt;
                pulseTimer.timer[it] += ctotal * Utility.gaussNormArea(t, tpeak, stdv) * project.dt; // mM
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
}
