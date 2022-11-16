package ucl.silver.d3d.core;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.io.Serializable;

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
 *
 * Description: a class that holds parameters (ParamObjects) to be displayed in
 * the front table of D3D (PanelParams). This class is to be extended by other
 * classes (Diffusant, Detector, etc) thereby inhereting the functions defined
 * below.
 *
 * The ParamVector will include declared fields that are not private (integers,
 * floats, doubles, booleans and Strings).
 *
 * In order to change the value of a parameter, one must overwrite the
 * setMyParams function in subclass, which accepts string name and double value
 * (or string name and string value), and then change the value of the
 * appropriate parameter, returning true if the value has been set. Otherwise,
 * false is returned, and the parameter is un-editable.
 *
 * Finally, use the "set" functions to set parameter values (e.g.
 * DiffusantSimple.set("d", 0.33))
 */
public class ParamVector implements Serializable {

    public String name = null;

    private String D3Dversion = null;

    public Project project;
    
    private transient ParamObject[] vector = null;
    private transient ParamVector[] userRegistry = null;

    public ParamVector(Project p) {
        project = p;
        D3Dversion = Master.D3D_VERSION;
    }

    public void error(String message) {

        if (message == null) {
            return;
        }

        System.err.println(message + " (" + toString() + ")");

        if (project != null) {
            Master.writeToLogFile(message + " (" + toString() + ")", false);
            project.addError(message);
        }

    }

    public void error(String function, String parameter, String message) {

        if ((function == null) || (message == null)) {
            return;
        }

        setParamError(parameter, message);

        D3Derror e = new D3Derror(toString(), function, parameter, message);

        Master.log(function + ", \"" + parameter + "\" " + message + " (" + toString() + ")");

        if (project != null) {
            project.addError(e);
        } else {
            System.err.println(message + " (" + toString() + ")");
        }

    }

    public String getD3Dversion() {
        return D3Dversion;
    }

    public ParamObject[] getVector() {
        return vector;
    }

    public void addVector(ParamObject[] v) {

        if (v == null) {
            return;
        }

        for (ParamObject o : v) {
            addElement(o);
        }

    }

    private boolean addElement(ParamObject o) {

        int ipnts;

        ParamObject[] temp;

        if (o == null) {
            return false;
        }

        if (vector == null) {
            vector = new ParamObject[1];
            vector[0] = o;
            return true;
        }

        ipnts = vector.length;

        temp = new ParamObject[ipnts + 1];

        System.arraycopy(vector, 0, temp, 0, ipnts);

        temp[ipnts] = o;

        vector = temp;

        return false;

    }

    public boolean createVector(boolean close) {
        
        int j = 0;
        int maxNumSuperClasses = 20; // can increase if necessary
        
        Class[] c = new Class[maxNumSuperClasses];

        vector = null;

        String NAME = getClassName(getClass());

        addClassParam(NAME);
        
        c[0] = getClass();
        
        if (c[0] == null) {
            return false;
        }
        
        for (int i = 1; i < c.length; i++) {
            if (c[i-1] == null) {
                break;
            } else {
                c[i] = c[i-1].getSuperclass();
                j = i;
            }
        }
        
        for (int i = j; i >= 0; i--) {
            addParams(c[i]);
        }

        if (close) {
            addEndClassParam(NAME);
        }

        return true;

    }

    public void clearVector() {
        vector = null;
    }

    public void updateVectors() {
        updateVector(vector);
        updateUsers();
    }

    public void updateVector() {
        updateVector(vector);
    }

    public void updateVector(ParamObject[] v) {
        updateParams(v, this);
    }

    public void closeVector() {
        addEndClassParam(getClassName(getClass()));
    }

    private int addParams(Class c) {

        int modifier, count = 0;
        String type, NAME, units, help;

        if (c == null) {
            return 0;
        }

        Field[] fields = c.getDeclaredFields(); // create an array of all parameters in the object

        if (fields.length == 0) {
            return count;
        }

        for (Field f : fields) {

            if (f == null) {
                continue;
            }

            try {
                if (f.get(this) == null) {
                    continue;
                }
            } catch (IllegalAccessException e) {
            }

            type = f.getType().toString();
            NAME = f.getName();
            units = units(NAME);
            help = help(NAME);
            modifier = f.getModifiers();

            if (Modifier.isPrivate(modifier)) {
                continue;
            }

            if (hide(NAME)) {
                continue;
            }

            try {

                if (type.equalsIgnoreCase("int")) {
                    addParam(this, NAME, f.getInt(this), units, help);
                    count++;
                } else if (type.equalsIgnoreCase("double")) {
                    addParam(this, NAME, f.getDouble(this), units, help);
                    count++;
                } else if (type.equalsIgnoreCase("boolean")) {
                    addParam(this, NAME, f.getBoolean(this), help);
                    count++;
                } else if (type.equals("class java.lang.String")) {
                    addParam(this, NAME, f.get(this).toString(), units, help);
                    count++;
                } else {
                    // skip everything else
                    //Master.log(f);
                }

            } catch (IllegalAccessException e) {
            }

        }

        return count;

    }

    private void updateParams(ParamObject[] v, ParamVector pv) {

        if (pv == null) {
            return;
        }

        updateParams(v, pv, pv.getClass().getSuperclass().getSuperclass());
        updateParams(v, pv, pv.getClass().getSuperclass());
        updateParams(v, pv, pv.getClass());

    }

    private void updateParams(ParamObject[] v, ParamVector pv, Class c) {

        String type, NAME;

        ParamObject o;

        if ((pv == null) || (c == null)) {
            return;
        }

        Field[] fields = c.getDeclaredFields(); // create an array of all parameters in the object

        if (fields.length == 1) {
            return;
        }

        for (Field f : fields) {

            if (f == null) {
                continue;
            }

            try {
                if (f.get(pv) == null) {
                    continue;
                }
            } catch (IllegalAccessException e) {
                continue; //ignore
            }

            type = f.getType().toString();
            NAME = f.getName();

            o = getParamObject(v, pv, NAME);

            if (o == null) {
                continue;
            }

            try {

                if (type.equalsIgnoreCase("int")) {
                    o.setValue(f.getInt(pv));
                } else if (type.equalsIgnoreCase("double")) {
                    o.setValue(f.getDouble(pv));
                } else if (type.equalsIgnoreCase("boolean")) {
                    if (f.getBoolean(pv)) {
                        o.setValue(1);
                    } else {
                        o.setValue(0);
                    }
                } else if (type.equals("class java.lang.String")) {
                    o.setValue(f.get(pv).toString());
                } else {
                    // skip everything else
                    //Master.log(f);
                }

            } catch (IllegalAccessException e) {
            }

        }

    }

    public String getClassName(Class c) {

        String pack, cstr;

        if (c == null) {
            return null;
        }

        pack = c.getPackage().toString(); // package name
        pack = pack.replaceAll("package ", "");

        cstr = c.getName();
        cstr = cstr.replaceAll(pack + ".", ""); // remove package name
        //cstr = cstr.toUpperCase(); // this is the class name

        return cstr;

    }

    public boolean canEdit(String name) {
        return true;
    }

    public boolean hide(String name) {
        return false;
    }

    public String units(String name) {
        return "";
    }

    public String help(String name) {
        return "";
    }

    public void addBlankParam() {
        addParam(null, "", "", "", "");
    }

    public void addClassParam(String name) {
        addParam(null, "Class", name, "", "");
    }

    public void addEndClassParam(String name) {
        addParam(null, "EndClass", name, "", "");
    }

    public void addParam(ParamVector pv, String name, double value, String units, String help) {
        addElement(new ParamObject(pv, name, value, units, help));
    }

    public void addParam(ParamVector pv, String name, int value, String units, String help) {
        addElement(new ParamObject(pv, name, value, units, help));
    }

    public void addParam(ParamVector pv, String name, boolean value, String help) {
        addElement(new ParamObject(pv, name, value, help));
    }

    public void addParam(ParamVector pv, String name, String value, String units, String help) {
        addElement(new ParamObject(pv, name, value, units, help));
    }

    // get the ParamObject class for particular variable
    public ParamObject getParamObject(ParamObject[] v, ParamVector pv, String name) {

        if (name == null) {
            return null;
        }

        if (vector == null) {
            return null;
        }

        for (ParamObject po : vector) {
            if ((po.paramVector == pv) && name.equalsIgnoreCase(po.getName())) {
                return po;
            }
        }

        //Master.log("Failed to find parameter: " + name + " (" + toString() + ")");
        return null;

    }

    public double getVar(String name) {
        return getVar(getParamObject(vector, this, name));
    }

    public double getVar(ParamObject o) {

        double dvalue = Double.NaN;

        if (o == null) {
            return dvalue;
        }

        try {
            dvalue = Double.parseDouble(o.getValue());
        } catch (NumberFormatException e) {
        }

        return dvalue;

    }

    public String getStr(String name) {
        return getStr(getParamObject(vector, this, name));
    }

    public String getStr(ParamObject o) {

        if (o == null) {
            return null;
        }

        return o.getValue();

    }

    public void setParamError(String name, String errorStr) {

        ParamObject o = getParamObject(vector, this, name);

        if (o == null) {
            return;
        }

        o.ERROR = errorStr;

    }

    public boolean set(String name, double d) {

        ParamObject o = getParamObject(vector, this, name);

        return set(o, d);

    }

    public boolean set(String name, boolean boo) {

        ParamObject o = getParamObject(vector, this, name);

        if (boo) {
            return set(o, 1);
        } else {
            return set(o, 0);
        }

    }

    public boolean set(String name, String str) {

        ParamObject o = getParamObject(vector, this, name);

        return set(o, str);

    }

    public boolean set(ParamObject o, double d) {

        if (o == null) {
            return false;
        }

        if (setMyParams(o, d)) {
            return o.setValue(d);
        }

        return false;

    }

    public boolean set(ParamObject o, String s) {

        Double d;

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!o.testValue(s)) {
            return false;
        }

        if (o.getType().equalsIgnoreCase("int") || o.getType().equalsIgnoreCase("double")) {

            d = Double.valueOf(s);

            if (setMyParams(o, d)) {
                return o.setValue(d);
            }

        } else if (o.getType().equalsIgnoreCase("boolean")) {

            if (s.equalsIgnoreCase("true") || s.equals("1")) {
                if (setMyParams(o, 1)) {
                    return o.setValue(1);
                }
            } else if (s.equalsIgnoreCase("false") || s.equals("0")) {
                if (setMyParams(o, 0)) {
                    return o.setValue(0);
                }
            }

        } else if (o.getType().equalsIgnoreCase("String")) {

            if (setMyParams(o, s)) {
                return o.setValue(s);
            }

        }

        return false;

    }

    public boolean setParamObject(String name, double v) {

        if (name == null) {
            return false;
        }

        ParamObject o = getParamObject(vector, this, name);

        if (o == null) {
            return false;
        }

        return o.setValue(v);

    }

    public boolean setParamObject(String name, String str) {

        if (name == null) {
            return false;
        }

        ParamObject o = getParamObject(vector, this, name);

        if (o == null) {
            return false;
        }

        return o.setValue(str);

    }

    public boolean setParamObject(String name, boolean b) {

        if (name == null) {
            return false;
        }

        ParamObject o = getParamObject(vector, this, name);

        if (o == null) {
            return false;
        }

        if (b) {
            return o.setValue(1);
        } else {
            return o.setValue(0);
        }

    }

    // this function needs to be overwritten in extended class (see Project class)
    public boolean setMyParams(ParamObject o, double value) {
        return false;
    }

    // this function needs to be overwritten in extended class (see Project class)
    public boolean setMyParams(ParamObject o, String s) {

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        return false;
    }

    public boolean addUser(ParamVector pv) {

        int ipnts;

        ParamVector[] temp;

        if (pv == null) {
            return false;
        }

        if (userRegistry == null) {
            userRegistry = new ParamVector[1];
            userRegistry[0] = pv;
            return true;
        }

        ipnts = userRegistry.length;

        temp = new ParamVector[ipnts + 1];

        System.arraycopy(userRegistry, 0, temp, 0, ipnts);

        temp[ipnts] = pv;

        userRegistry = temp;

        return false;

    }

    public void updateUsers() {

        if (userRegistry == null) {
            return;
        }

        for (ParamVector pv : userRegistry) {
            if (pv.getVector() != null) {
                updateVector(pv.getVector());
                //System.out.println("updated user " + userRegistry[i]);
            }
        }

    }

}
