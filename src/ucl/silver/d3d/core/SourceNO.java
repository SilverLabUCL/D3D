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
public class SourceNO extends Source {
    
    // y = M1*(1-exp(-x/t1))^P * exp(-x/t2) + M3*(1-exp(-(invTau3*x)))
    
    public double M1=16.2/1000.0; // molecules/ms
    public double C1=0; // mM/ms/voxel
    public double tau1=0.0107*1000.0; // ms
    public double Power=4.931;
    public double tau2=0.1134*1000.0; // ms
    public double M3=1.552/1000.0; // molecules/ms
    public double C3=0; // mM/ms/voxel
    public double invTau3=9.056/1000.0; // kHz
    
    public double Cscale = 1.0; // extra scaling factor

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("M1")) {
            return "molecules/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("M3")) {
            return "molecules/" + project.timeUnits;
        }
        if (name.equalsIgnoreCase("C1")) {
            return project.concUnits + "/" + project.timeUnits + "/voxel";
        }
        if (name.equalsIgnoreCase("C3")) {
            return project.concUnits + "/" + project.timeUnits + "/voxel";
        }
        if (name.equalsIgnoreCase("tau1")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("tau2")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("invTau3")) {
            return "1/" + project.timeUnits;
        }
        return super.units(name);
    }

    public SourceNO(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        super(p, NAME, DiffusantNum, c, pt);
        createVector(true); // sets ParamVector
    }

    @Override
    public void init() {
        super.init();
        updateC();
    }

    public void updateC() {

        double moles;
        double AN = 6.0221415E23; // Avogadro's number
        double dx = project.dx;
        double litersPerVoxel = Math.pow(dx, 3) * 1e-18 * 1000; // liters per voxel

        moles = M1 / AN;
        C1 = 1000 * moles / litersPerVoxel; // mM / ms
        C1 /= voxels; // spread evenly between voxels
        C1 *= Cscale; // additional scale factor

        moles = M3 / AN;
        C3 = 1000 * moles / litersPerVoxel; // mM / ms
        C3 /= voxels; // spread evenly between voxels
        C3 *= Cscale; // additional scale factor

        updateVectors();

        // NOTE, C0 is not used in this class

    }

    @Override
    public double release(RunFiniteDifference fd, Geometry geometry) {

        int indexNum, numVoxels;
        double trel, e1, e2, e3;
        double conc = 0, avgC = 0;
        double dt = project.dt;

        if (spaceVoxels == 0) {
            return 0;
        }

        if ((fd.diffus == null) || (fd.it >= fd.itmax)) {
            return 0;
        }

        if ((pulseTimer == null) || (pulseTimer.pulses == null)) {
            return 0;
        }

        for (int i = 0; i < pulseTimer.pulses.length; i += 1) { // loop thru pulses

            trel = fd.time - pulseTimer.pulses[i].time;

            e1 = C1 * Math.pow(1.0 - Math.exp(-trel / tau1), Power);
            e2 = Math.exp(-trel / tau2);
            e3 = C3 * (1.0 - Math.exp(-trel * invTau3));

            conc += e1 * e2 + e3; // mM/ms

        }

        conc *= dt; // mM

        if (fd.diffus[0].length == 1) {

            fd.diffus[diffusantNum][0] += conc;
            avgC = conc;

        } else {

            if (indexRFD == null) {
                Master.exit("SourceNO error: coordinates index array has not be initialized.");
            }

            numVoxels = indexRFD.length;

            for (int i = 0; i < numVoxels; i++) {
                indexNum = indexRFD[i];
                fd.diffus[diffusantNum][indexNum] += conc;
                avgC += conc;
            }

            if (numVoxels > 1) {
                avgC /= 1.0 * numVoxels;
            }

        }

        save.saveData(avgC);

        return avgC;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {
        
        if (o == null) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("M1")) {
            if (v < 0) {
                return false;
            }
            M1 = v;
            updateC();
            return true;
        }
        if (n.equalsIgnoreCase("tau1")) {
            if (v <= 0) {
                return false;
            }
            tau1 = v;
            return true;
        }
        if (n.equalsIgnoreCase("Power")) {
            if (v < 0) {
                return false;
            }
            Power = v;
            return true;
        }
        if (n.equalsIgnoreCase("tau2")) {
            if (v <= 0) {
                return false;
            }
            tau2 = v;
            return true;
        }
        if (n.equalsIgnoreCase("M3")) {
            if (v < 0) {
                return false;
            }
            M3 = v;
            updateC();
            return true;
        }
        if (n.equalsIgnoreCase("invTau3")) {
            if (v < 0) {
                return false;
            }
            invTau3 = v;
            return true;
        }
        if (n.equalsIgnoreCase("Cscale")) {
            if (v <= 0) {
                return false;
            }
            Cscale = v;
            return true;
        }
        return super.setMyParams(o, v);
    }
}
