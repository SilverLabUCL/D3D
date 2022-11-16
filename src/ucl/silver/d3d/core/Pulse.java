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
public class Pulse extends ParamVector {

    public double time = 0;
    public double duration = 0;
    public double amplitude = 0;

    public int it_bgn = -1, it_end = -1; // used during simulation
    public double it_amp = Double.NaN; // used during simulation

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("it_bgn")) {
            return false;
        }
        if (name.equalsIgnoreCase("it_end")) {
            return false;
        }
        if (name.equalsIgnoreCase("it_amp")) {
            return false;
        }
        return super.canEdit(name);
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

    public Pulse(Project p, double TIME, double DURATION, double AMPLITUDE) {
        super(p);
        time = TIME;
        duration = DURATION;
        amplitude = AMPLITUDE;
    }

    public void update_it_params(boolean normalizeAmplitudePerDT) {
        double pnts;

        int itmax = project.simPoints();

        it_bgn = (int) (time / project.dt);

        if (duration < 0) {
            it_end = it_bgn - 1; // off
        } else if (duration == 0) {
            it_end = it_bgn; // impulse
        } else {
            it_end = it_bgn + (int) (duration / project.dt);
        }

        it_end = Math.min(it_end, itmax - 1);
        
        if (normalizeAmplitudePerDT) {
            pnts = it_end - it_bgn + 1;
            it_amp = amplitude / pnts; // spread across all points
        } else {
            it_amp = amplitude;
        }

    }
    
    public static int pulseNum(String pulseName) { // e.g. "Pulse0"
        
        if (!pulseName.startsWith("Pulse")) {
            return -1; // name must begin with "Pulse"
        }
        
        String istr = pulseName.substring(5, pulseName.length());

        if (istr.length() >= 1) {
            try {
                return Integer.parseInt(istr);
            } catch (NumberFormatException e) {
                return -1;
            }
        }
        
        return -1;
        
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
