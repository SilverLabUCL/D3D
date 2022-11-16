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
public class SourceUptake extends Source {
    
    public double Kuptake = Double.NaN; // 1/ms (e.g. uptake)
    
    public Q10 Q10_Kuptake = null; // Q10 temperature scaling for Kuptake
    
    public SourceUptake(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PulseTimer pt, double KUPTAKE) {
        super(p, NAME, DiffusantNum, c, pt, "K");
        Kuptake = KUPTAKE;
        createVector(true); // sets ParamVector
    }

    public SourceUptake(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double KUPTAKE) {
        super(p, NAME, DiffusantNum, c, null, "K");
        Kuptake = KUPTAKE;
        createVector(true); // sets ParamVector
    }
    
    @Override
    public boolean canEdit(String name) {
        boolean pulseTimerExists = (pulseTimer != null);
        if (name.equalsIgnoreCase("Kuptake")) {
            return type.equalsIgnoreCase("K") && !pulseTimerExists;
        }
        if (name.equalsIgnoreCase("clamp")) {
            return false;
        }
        return super.canEdit(name);
    }
    
    @Override
    public void scaleQ10() {
        double k;
        super.scaleQ10();
        if ((Q10_Kuptake != null) && Q10_Kuptake.on) {
            k = Q10_Kuptake.getScaledValue();
            if (Double.isFinite(k)) {
                Kuptake = k;
            }
        }
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
        
        if (pulseTimer != null) {
            k = pulseTimer.high(fd.it);
            k = Math.abs(k);
            timerHigh = k > 0;
        } else {
            timerHigh = true; // on for all time
            k = Math.abs(Kuptake);
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
    
    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);
        
        if (Q10_Kuptake != null) {
            Q10_Kuptake.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }
        
        if (Q10_Kuptake != null) {
            addBlankParam();
            Q10_Kuptake.createVector(true);
            addVector(Q10_Kuptake.getVector());
            Q10_Kuptake.addUser(this);
        }

        if (close){
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);
        
        if (Q10_Kuptake != null) {
            Q10_Kuptake.updateVector(v);
        }

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

        if (n.equalsIgnoreCase("Kuptake")) {
            if (v >= 0) {
                Kuptake = v;
                if (Q10_Kuptake != null) {
                    Q10_Kuptake.value = v;
                }
                return true;
            }
        }
        
        return super.setMyParams(o, v);
        
    }

}
