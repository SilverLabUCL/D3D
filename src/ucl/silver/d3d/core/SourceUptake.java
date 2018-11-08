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
public class SourceUptake extends Source {

    private double dummyvar; // need this here so krate is displayed on Parameter tab
    public double krate = 0; // rate of uptake (1/ms)
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("krate")) {
            return "1/" + project.timeUnits;
        }
        return super.units(name);
    }
    
    public SourceUptake(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt) {
        super(p, NAME, DiffusantNum, c, pt);
        createVector(true); // sets ParamVector
    }

    public SourceUptake(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double Krate) {
        super(p, NAME, DiffusantNum, c, null);
        krate = Krate;
        createVector(true); // sets ParamVector
    }

    @Override
    public double release(RunFiniteDifference fd, Geometry geometry) {

        double k, c0, c1, c2, newconc, cdelta = 0;
        boolean timerHigh;

        if (spaceVoxels == 0) {
            return 0;
        }

        if ((fd.diffus == null) || (fd.it >= fd.itmax)) {
            return 0;
        }

        if ((pulseTimer == null) || (pulseTimer.timer == null)) {
            timerHigh = true; // on for all time
            k = Math.abs(krate);
        } else {
            k = Math.abs(pulseTimer.timer[fd.it]);
            timerHigh = k > 0;
        }
        
        if (timerHigh) {
            
            c0 = project.diffusants[diffusantNum].C0;

            if (fd.diffus[0].length == 1) { // single compartment

                c1 = fd.diffus[diffusantNum][0];
                c2 = (c1 - c0) * k * project.dt;
                newconc = c1 - c2;
                
                if (newconc >= c0) {
                    fd.diffus[diffusantNum][0] = newconc;
                    cdelta = c2;
                }

            } else {

                //if (indexRFD == null) {
                    //Master.exit("Source error: coordinates index array has not be initialized.");
                    //return 0;
                //}

                for (int i : indexRFD) {

                    c1 = fd.diffus[diffusantNum][i];
                    c2 = (c1 - c0) * k * project.dt;
                    newconc = c1 - c2;

                    if (newconc >= c0) {
                        fd.diffus[diffusantNum][i] = newconc;
                        cdelta += c2;
                    }

                }

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

        if (!(o.paramVector instanceof SourceUptake)) {
            return false;
        }

        if (n.equalsIgnoreCase("krate")) {
            if (v > 0) {
                krate = (int) v;
                return true;
            }
        }
        return super.setMyParams(o, v);
    }

}
