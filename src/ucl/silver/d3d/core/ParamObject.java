package ucl.silver.d3d.core;

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
 */
public class ParamObject implements Serializable {

    private String name = ""; // parameter name
    private String type = ""; // "double" or "int" or "boolean" or "String" or "Class" or "EndClass"
    private String strValue = ""; // parameter value, saved as a string
    private String units = ""; // parameter units
    private String help = ""; // parameter help/description

    public String ERROR = null;

    public ParamVector paramVector = null; // pointer to the vector where parameter object resides

    public ParamObject(ParamVector PV, String NAME, double value, String UNITS, String HELP) { // create a double-type object
        if (NAME.length() > 0) {
            paramVector = PV;
            name = NAME;
            strValue = Double.toString(value);
            type = "double";
            units = UNITS;
            help = HELP;
        }
    }

    public ParamObject(ParamVector PV, String NAME, int ivalue, String UNITS, String HELP) { // create an integer-type object
        if (NAME.length() > 0) {
            paramVector = PV;
            name = NAME;
            strValue = Integer.toString(ivalue);
            type = "int";
            units = UNITS;
            help = HELP;
        }
    }

    public ParamObject(ParamVector PV, String NAME, boolean bvalue, String HELP) { // create a boolean-type object
        if (NAME.length() > 0) {
            paramVector = PV;
            name = NAME;
            strValue = Boolean.toString(bvalue);
            type = "boolean";
            help = HELP;
        }
    }

    public ParamObject(ParamVector PV, String NAME, String StrValue, String UNITS, String HELP) { // create a string-type object

        if (NAME.length() > 0) {

            strValue = StrValue;

            if (NAME.equalsIgnoreCase("Class")) {
                type = "Class";
            } else if (NAME.equalsIgnoreCase("EndClass")) {
                type = "EndClass";
            } else {
                paramVector = PV;
                type = "String";
                name = NAME;
                units = UNITS;
                help = HELP;
            }

        }

    }

    public boolean setValue(String sValue) {

        if (testValue(sValue)) {
            strValue = sValue;
            return true;
        }

        return false;

    }

    public boolean setValue(double value) {

        if (type.equalsIgnoreCase("int")) {
            strValue = Integer.toString((int) value);
        } else if (type.equalsIgnoreCase("double")) {
            strValue = Double.toString(value);
        } else if (type.equalsIgnoreCase("boolean")) {
            if (value == 0) {
                strValue = "false";
            } else {
                strValue = "true";
            }
        }

        return true;

    }

    public boolean testValue(String s) {

        if (s == null) {
            return false;
        }

        try {
            if (type.equalsIgnoreCase("int")) {
                Integer.valueOf(s); // may throw NumberFormatException
                return true;
            } else if (type.equalsIgnoreCase("double")) {
                Double.valueOf(s); // may throw NumberFormatException
                return true;
            } else if (type.equalsIgnoreCase("boolean")) {
                if (s.equalsIgnoreCase("true") || s.equalsIgnoreCase("false")
                        || s.equals("1") || s.equals("0")) {
                    return true;
                }
            } else if (type.equalsIgnoreCase("String")) {
                return true;
            }
        } catch (NumberFormatException e) {
            return false;
        }

        return false;

    }

    public String getName() {
        return name;
    }

    public String getType() {
        return type;
    }

    public String getValue() {
        return strValue;
    }

    public String getUnits() {
        return units;
    }

    public String getHelp() {
        return help;
    }

    public boolean error() {
        return (ERROR != null);
    }

}
