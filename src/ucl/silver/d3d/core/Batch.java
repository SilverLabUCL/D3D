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
public class Batch extends ParamVector {

    public String folder = ""; // subfolder where outputs are saved

    public int batchNum;

    public String parameter = ""; // long variable name (i.e. "Project.simTime" or "Ca.D")
    public String strValue = ""; // parameter value - if it's a string

    public double value; // parameter value - if it's a number

    public boolean finished = false;
    public boolean execute = true;

    public Batch(int BATCHNUM, String longName, double VALUE, boolean EXECUTE) {
        super(null);
        batchNum = BATCHNUM;
        parameter = longName;
        value = VALUE;
        strValue = "";
        execute = EXECUTE;
        createVector(true);
    }

    public Batch(int BATCHNUM, String longName, String STRVALUE, boolean EXECUTE) {
        super(null);
        batchNum = BATCHNUM;
        parameter = longName;
        strValue = STRVALUE;
        value = Double.NaN;
        execute = EXECUTE;
        createVector(true);
    }

    public boolean isString() {
        return (strValue.length() > 0);
    }

    public void init() {
        // nothing to do
    }

    public void setFolderName() {

        String fname = "Batch" + Integer.toString(batchNum) + "_" + parameter;

        if (isString()) {
            fname += "_" + strValue;
        } else {
            fname += "_" + Float.toString((float) value);
        }

        fname = fname.replace('.', '_');
        fname = fname.replace(';', '_');
        fname = fname.replace(',', '_');

        if (fname.length() > 0) {
            folder = fname;
            setParamObject("folder", folder);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("batchNum")) {
            if (v < 0) {
                return false;
            }
            batchNum = (int) v;
            setFolderName();
            return true;
        }
        if (n.equalsIgnoreCase("value")) {
            value = v;
            setFolderName();
            return true;
        }
        if (n.equalsIgnoreCase("finished")) {
            finished = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("execute")) {
            execute = (v == 1);
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

        String n = o.getName();

        if (n.equalsIgnoreCase("parameter")) {
            parameter = s;
            setFolderName();
            return true;
        }
        if (n.equalsIgnoreCase("strValue")) {
            strValue = s;
            return true;
        }
        if (n.equalsIgnoreCase("folder")) {
            folder = s;
            return true;
        }
        return super.setMyParams(o, s);
    }

}
