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
public class PulseTimer extends ParamVector {

    public int numPulses = 0;

    public boolean useTimerArray = true; // create timer array to be used during simulations
    // this is faster than looping thru each pulse on each time step
    // but creates an array of data equal to number of simulation points
    // which could be large depending on time step and simulation length
    
    public String ampUnits = ""; // amplitude units for all pulses

    public Pulse[] pulses = null;

    public double[] timer = null; // pulse array used during simulations if useTimerArray = true

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("numPulses")) {
            return false;
        }
        return true;
    }

    public PulseTimer(Project p) {
        super(p);
    }

    public PulseTimer(Project p, double TIME) {
        super(p);
        add(TIME);
    }

    public PulseTimer(Project p, double TIME, double DURATION) {
        super(p);
        add(TIME, DURATION);
    }
    
    public PulseTimer(Project p, double TIME, double DURATION, double AMPLITUDE, String AMPUNITS) {
        super(p);
        ampUnits = AMPUNITS;
        add(TIME, DURATION, AMPLITUDE);
    }

    public PulseTimer(Project p, double TIME, double DURATION, double AMPLITUDE, String AMPUNITS, int numPulses, double pulseInterval) {
        super(p);
        addTrain(TIME, DURATION, AMPLITUDE, numPulses, pulseInterval);
        ampUnits = AMPUNITS;
    }

    public final int add(double TIME) {
        return addPulse(new Pulse(project, TIME, 0, 1));
    }

    public final int add(double TIME, double DURATION) {
        return addPulse(new Pulse(project, TIME, DURATION, 1));
    }

    public final int add(double TIME, double DURATION, double AMPLITUDE) {
        return addPulse(new Pulse(project, TIME, DURATION, AMPLITUDE));
    }

    public final int addTrain(double TIME, double DURATION, double AMPLITUDE, int numPulses, double pulseInterval) {

        int j = -1;

        for (int i = 0; i < numPulses; i++) {
            j = addPulse(new Pulse(project, TIME + i * pulseInterval, DURATION, AMPLITUDE));
        }

        return j;

    }

    private int addPulse(Pulse newPulse) {

        int i = 0;

        if (pulses != null) {
            i = pulses.length;
        }

        Pulse[] newArray = new Pulse[i + 1];

        if (i > 0) {
            System.arraycopy(pulses, 0, newArray, 0, i);
        }

        newArray[i] = newPulse;

        pulses = newArray; // replace old array with new one

        numPulses = pulses.length;
        
        for (int j = 0 ; j < pulses.length; j++) {
            pulses[j].name = "pulse" + j;
        }

        //Master.log("Added pulse #" + i);
        return i;

    }

    public final void initTimer(boolean normalizeAmplitudePerDT, boolean updatePulses) {

        int itmax, itbgn, itend;
        
        timer = null;

        if ((pulses == null) || (pulses.length == 0)) {
            return;
        }

        if (updatePulses) {
            for (Pulse p : pulses) {
                p.update_it_params(normalizeAmplitudePerDT); // update it_bgn, it_end, it_amp
            }
        }

        if (!useTimerArray) {
            return;
        }
        
        itmax = project.simPoints();

        if (itmax <= 0) {
            return;
        }

        timer = new double[itmax];

        for (Pulse p : pulses) {

            itbgn = Math.max(p.it_bgn, 0);
            itend = Math.min(p.it_end, timer.length - 1);

            for (int it = itbgn; it <= itend; it++) {
                timer[it] += p.it_amp;
            }

        }

    }

    public double high(int itime) {
        double sumAmp = 0;

        if ((timer != null) && (itime < timer.length)) {
            return timer[itime];
        }

        if ((pulses == null) || (pulses.length == 0)) {
            return 0;
        }

        for (Pulse p : pulses) {
            if ((itime >= p.it_bgn) && (itime <= p.it_end)) {
                sumAmp += p.it_amp;
            }
        }

        return sumAmp;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (pulses != null) {
            for (Pulse p : pulses) {
                if (p != null) {
                    p.addUser(pv);
                }
            }
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (pulses != null) {
            for (Pulse p : pulses) {
                if (p != null) {
                    addBlankParam();
                    p.createVector(true);
                    addVector(p.getVector());
                    p.addUser(this);
                }
            }
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (pulses != null) {
            for (Pulse p : pulses) {
                if (p != null) {
                    p.updateVector(v);
                }
            }
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("useTimerArray")) {
            useTimerArray = (v == 1);
            return true;
        }
        //if (n.equalsIgnoreCase("normalizeAmplitudesPerDT")) {
            //normalizeAmplitudesPerDT = (v == 1);
            //return true;
        //}

        return super.setMyParams(o, v);

    }
}
