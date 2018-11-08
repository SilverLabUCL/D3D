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
public class PulseTimer extends ParamVector {

    public int numPulses = 0;

    public boolean useTimerArray = true; // create timer array to be used during simulations
    // this is faster than looping thru each pulse on each time step
    // but creates an array of data equal to number of simulation points

    public Pulse[] pulses = null;

    public double[] timer = null; // pulse array used during simulations

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

    public PulseTimer(Project p, double TIME, double DURATION, double AMPLITUDE) {
        super(p);
        add(TIME, DURATION, AMPLITUDE);
    }

    public PulseTimer(Project p, double TIME, double DURATION, double AMPLITUDE, int numPulses, double pulseInterval) {
        super(p);
        addTrain(TIME, DURATION, AMPLITUDE, numPulses, pulseInterval);
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

        //Master.log("Added pulse #" + i);
        return i;

    }

    public double sumAmplitude() {
        double sum = 0;

        if (pulses == null) {
            return Double.NaN;
        }

        for (Pulse p : pulses) {
            sum += p.amplitude;
        }

        return sum;

    }

    public double sumDuration() {
        double sum = 0;

        if (pulses == null) {
            return Double.NaN;
        }

        for (Pulse p : pulses) {
            sum += p.duration;
        }

        return sum;

    }

    public final void initTimer() {

        int itmax = project.simPoints();
        timer = null;

        if (pulses == null) {
            return;
        }

        for (Pulse p : pulses) {
            p.update();
        }

        if (!useTimerArray) {
            return;
        }

        if (itmax <= 0) {
            return;
        }

        timer = new double[itmax];

        for (Pulse p : pulses) {
            for (int it = p.itbgn; it <= p.itend; it++) {
                timer[it] += p.amp;
            }
        }

    }

    public double high(int itime) {
        double sumAmp = 0;

        if ((timer != null) && (itime < timer.length)) {
            return timer[itime];
        }

        if (pulses == null) {
            return 0;
        }

        for (Pulse p : pulses) {
            if ((itime >= p.itbgn) && (itime <= p.itend)) {
                sumAmp += p.amp;
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

        return super.setMyParams(o, v);

    }
}
