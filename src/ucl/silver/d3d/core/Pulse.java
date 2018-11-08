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
public class Pulse extends ParamVector {

    public double time = 0;
    public double duration = 0;
    public double amplitude = 1;

    public int itbgn, itend; // used during simulation
    public double amp; // used during simulation

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("itbgn")) {
            return false;
        }
        if (name.equalsIgnoreCase("itend")) {
            return false;
        }
        if (name.equalsIgnoreCase("amp")) {
            return false;
        }
        return super.canEdit(name);
    }

    public Pulse(Project p, double TIME, double DURATION, double AMPLITUDE) {
        super(p);
        time = TIME;
        duration = DURATION;
        amplitude = AMPLITUDE;
    }

    public void update() {
        double pnts;

        int itmax = project.simPoints();

        itbgn = (int) (time / project.dt);

        if (duration < 0) {
            itend = itbgn - 1; // off
        } else if (duration == 0) {
            itend = itbgn; // impulse
        } else {
            itend = itbgn + (int) (duration / project.dt);
        }

        itend = Math.min(itend, itmax - 1);

        pnts = itend - itbgn + 1;
        amp = amplitude / pnts; // spread evenly across pulse

    }

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("time")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("duration")) {
            return project.timeUnits;
        }
        return super.units(name);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("time")) {
            time = v;
            return true;
        }
        if (n.equalsIgnoreCase("duration")) {
            duration = v;
            return true;
        }
        if (n.equalsIgnoreCase("amplitude")) {
            amplitude = v;
            return true;
        }
        return super.setMyParams(o, v);
    }
}
