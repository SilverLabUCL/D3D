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
public class SourceStepWithTail extends Source {

    // tail current is 2-exp decay

    public double cStepScaleFactor = 1.0; // amplitude at end of step = C0 * cStepScaleFactor
    public double a1Fraction = 0.5; // a1Fraction = a1 / ( a1 + a2 ), a2 = 1 - a1
    public double a2Fraction = 0.5;
    public double tau1 = 1.0, tau2 = 10.0; // time constants

    private int stepCounter = 0;
    private double tStepOff;

    public SourceStepWithTail(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        super(p, NAME, DiffusantNum, c, pt, "C");
        createVector(true); // sets ParamVector
    }

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("tau1")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("tau2")) {
            return project.timeUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("a2Fraction")) {
            return false;
        }
        return true;
    }

    @Override
    public void init() {

        String thisObj = "SourceStepWithTail";

        if ((pulseTimer == null) || (pulseTimer.pulses == null) || (pulseTimer.pulses.length == 0)) {
            error("init", "pulseTimer", "not defined");
        }

        if (clamp == true) {
            error("init", "clamp", "not appropriate");
        }

        if ((a1Fraction < 0) || (a1Fraction > 1)) {
            error("init", "a1Fraction", "must be fractional value between 0 and 1");
        }

        if (tau1 < 0) {
            error("init", "tau1", "bad time constant");
        }

        if (tau2 < 0) {
            error("init", "tau2", "bad time constant");
        }

        stepCounter = 0;
        a2Fraction = 1.0 - a1Fraction;
        tStepOff = pulseTimer.pulses[0].time + pulseTimer.pulses[0].duration;

        super.init();

    }

    @Override
    public double release(RunFiniteDifference fd, Geometry geometry) {

        int indexNum, numVoxels;
        double scale, a1, a2, cValue, avgC = 0;
        boolean timerHigh;

        if (spaceVoxels == 0) {
            return 0;
        }

        if ((fd.diffus == null) || (fd.it >= fd.itmax)) {
            return 0;
        }
        
        if (pulseTimer != null) {
            timerHigh = pulseTimer.high(fd.it) > 0;
        } else {
            timerHigh = true; // on for all time
        }

        if (timerHigh) {
            stepCounter = 1; // step on
        } else if (stepCounter == 1) {
            stepCounter = 2; // step off, tail current on
        }

        if (stepCounter > 0) {

            if (fd.time >= tStepOff) {
                a1 = a1Fraction * Math.exp(-(fd.time - tStepOff) / tau1);
                a2 = a2Fraction * Math.exp(-(fd.time - tStepOff) / tau2);
                scale = cStepScaleFactor * (a1 + a2);
                cValue = Ctotal * scale;
                //iValue = Ivoxel * scale;
            } else {
                cValue = Ctotal;
                //iValue = Ivoxel;
            }

            if (fd.diffus[0].length == 1) { // single compartment

                fd.diffus[diffusantNum][0] += cValue;
                avgC = cValue;
                //avgI = iValue;

            } else {

                if (indexRFD == null) {
                    //Master.exit("Source error: coordinates index array has not be initialized.");
                    return 0;
                }

                numVoxels = indexRFD.length;

                for (int i = 0; i < numVoxels; i++) {

                    indexNum = indexRFD[i];

                    fd.diffus[diffusantNum][indexNum] += cValue;
                    avgC += cValue;
                    //avgI += iValue;

                }

                if (numVoxels > 1) {
                    avgC /= 1.0 * numVoxels;
                    //avgI /= 1.0 * numVoxels;
                }

            }

        }

        saveValue(avgC);

        return avgC;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("cStepScaleFactor")) {
            cStepScaleFactor = v;
            return true;
        }
        if (n.equalsIgnoreCase("a1Fraction")) {
            if ((v < 0) && (v > 1)) {
                return false;
            }
            a1Fraction = v;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("tau1")) {
            if (v <= 0) {
                return false;
            }
            tau1 = v;
            return true;
        }
        if (n.equalsIgnoreCase("tau2")) {
            if (v <= 0) {
                return false;
            }
            tau2 = v;
            return true;
        }
        return super.setMyParams(o, v);
    }

}
