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
public class Q10 extends ParamVector {
    
    public double value = Double.NaN; // parameter value to scale
    public double temp = Double.NaN; // temperature parameter value was measured
    public double Q10 = Double.NaN; // Q10 factor for parameter
    
    public boolean on = true;
    
    public static double Q10_D = 1.3;
    public static double Q10_K = 2.0;
    
    public Q10(Project p, String paramName, double paramValue, double TEMP, double q10) {
        super(p);
        name = paramName;
        value = paramValue;
        temp = TEMP;
        Q10 = q10;
        createVector(true);
    }
    
    public double getScaledValue() {
        if (!on) {
            return 1;
        }
        if (Double.isNaN(Q10) || (Q10 <= 0)) {
            error("getScaledValue", "Q10", "bad value");
            return Double.NaN;
        }
        if (Double.isNaN(value)) {
            error("getScaledValue", "value", "NaN");
            return Double.NaN;
        }
        if (Double.isNaN(temp)) {
            error("getScaledValue", "temp", "NaN");
            return Double.NaN; // error
        }
        if (Double.isNaN(project.simTemp)) {
            error("getScaledValue", "project.simTemp", "NaN");
            return Double.NaN; // error
        }
        if (temp == project.simTemp) {
            return value;
        }
        return value * Math.pow(Q10, (project.simTemp - temp) / 10.0);
    }
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("temp")) {
            return project.tempUnits;
        }
        return super.units(name);
    }
    
    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("value")) {
            return false;
        }
        return super.canEdit(name);
    }
    
    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof Q10)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("Q10")) {
            if (v <= 0) {
                return false;
            }
            Q10 = v;
            return true;
        }
        if (n.equalsIgnoreCase("temp")) {
            if (v < 0) {
                return false;
            }
            temp = v;
            return true;
        }
        if (n.equalsIgnoreCase("on")) {
            on = (v == 1);
            return true;
        }
        return super.setMyParams(o, v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof Q10)) {
            return false;
        }

        String n = o.getName();

        //if (n.equalsIgnoreCase("name")) {
            //name = s;
            //return true;
        //}
        
        return super.setMyParams(o, s);
    }
    
}
