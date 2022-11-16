package ucl.silver.d3d.core;

import java.awt.*;
import javax.swing.*;
import ucl.silver.d3d.utils.*;

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
public class ColorD3D extends ParamVector {

    public Color color = null; // main RGB color
    public Color gradientColor = null; // second color for 2-color gradient
    
    public String colorRGB = null; // string representation of color, or color scale name
    public String gradientRGB = ""; // string representation of second color

    public boolean invert = false; // invert color scale or gradient

    public static String[] colorScaleList = {"BTC", "BTY", "Gray", "LinGray", "Heat", "Optimal", "LinOptimal", "Magenta", "Rainbow"};

    public static ColorScaleBTC colorScaleBTC = null;
    public static ColorScaleBTY colorScaleBTY = null;
    public static ColorScaleGray colorScaleGray = null;
    public static ColorScaleLinGray colorScaleLinGray = null;
    public static ColorScaleHeat colorScaleHeat = null;
    public static ColorScaleOptimal colorScaleOptimal = null;
    public static ColorScaleLinOptimal colorScaleLinOptimal = null;
    public static ColorScaleMagenta colorScaleMagenta = null;
    public static ColorScaleRainbow colorScaleRainbow = null;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("colorRGB")) {
            return "RGB";
        }
        if (name.equalsIgnoreCase("gradientRGB")) {
            if (isColorScale(colorRGB)) {
                return "";
            } else {
                return "RGB";
            }
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("colorRGB")) {
            return true;
        }
        if (name.equalsIgnoreCase("gradientRGB")) {
            if (isColorScale(colorRGB)) {
                return false;
            } else {
                return true;
            }
        }
        if (name.equalsIgnoreCase("invert")) {
            return true;
        }
        return false;
    }

    public ColorD3D(String NAME, Color COLOR) {
        super(null);
        name = NAME;
        color = COLOR;
        colorRGB = COLOR.toString();
        createVector(true);
    }

    public ColorD3D(String NAME, String colorScale) {
        super(null);
        name = NAME;
        colorRGB = colorScale;
        createVector(true);
    }

    public ColorD3D(String NAME, Color COLOR, Color GRADIENT) {
        super(null);
        name = NAME;
        color = COLOR;
        colorRGB = COLOR.toString();
        gradientColor = GRADIENT;
        gradientRGB = GRADIENT.toString();
        createVector(true);

    }

    public boolean isColorScale() {
        return isColorScale(colorRGB);
    }

    public boolean isColorGradient() {

        if (isColorScale(colorRGB)) {
            return false;
        }
        
        return ((color != null) && (gradientColor != null));

    }

    public Color getColor(double value, double minValue, double maxValue) {

        if (isColorGradient()) {
            return getColorGradient(value, minValue, maxValue);
        }

        if (isColorScale(colorRGB)) {
            return getColorScale(colorRGB, value, minValue, maxValue, invert);
        }

        if (color == null) {
            return Color.white;
        } else {
            return color;
        }

    }

    public Color getColorGradient(double value, double minValue, double maxValue) {

        int red, green, blue;
        double normRatio;

        if (gradientColor == null) {
            if (color == null) {
                return Color.white;
            } else {
                return color;
            }
        }
        
        normRatio = (value - minValue) / (maxValue - minValue);

        if (normRatio < 0) {
            normRatio = 0;
        } else if (normRatio > 1.0) {
            normRatio = 1.0;
        }

        if (invert) {
            red = (int) (gradientColor.getRed() * normRatio + color.getRed() * (1.0 - normRatio));
            green = (int) (gradientColor.getGreen() * normRatio + color.getGreen() * (1.0 - normRatio));
            blue = (int) (gradientColor.getBlue() * normRatio + color.getBlue() * (1.0 - normRatio));
        } else {
            red = (int) (color.getRed() * normRatio + gradientColor.getRed() * (1.0 - normRatio));
            green = (int) (color.getGreen() * normRatio + gradientColor.getGreen() * (1.0 - normRatio));
            blue = (int) (color.getBlue() * normRatio + gradientColor.getBlue() * (1.0 - normRatio));
        }
        
        return new Color(red, green, blue);

    }
    
    public boolean setColor(String newColorStr) {

        Color newColor;

        if (isColorScale(newColorStr)) {

            colorRGB = newColorStr;
            color = null;

        } else {

            newColor = string2color(newColorStr);

            if (newColor == null) {
                return false;
            }

            colorRGB = newColorStr;
            color = newColor;

        }

        return true;

    }

    public boolean setGradientColor(String newColorStr) {

        Color newColor;

        if (isColorScale(newColorStr)) {

            return false;

        } else {

            newColor = string2color(newColorStr);

            if (newColor == null) {
                return false;
            }

            gradientRGB = newColorStr;
            gradientColor = newColor;

            

        }

        return true;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     General Color Scale Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static boolean isColorScale(String scaleName) {
        int i = Utility.whichListItem(colorScaleList, scaleName, true);
        return (i >= 0);
    }

    public static Color getColorScale(String scaleName, double value, double minValue, double maxValue, boolean inverted) {

        if (scaleName == null) {
            return null;
        }

        if (scaleName.equalsIgnoreCase("BTC")) {

            if (colorScaleBTC == null) {
                colorScaleBTC = new ColorScaleBTC();
            }

            return colorScaleBTC.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("BTY")) {

            if (colorScaleBTY == null) {
                colorScaleBTY = new ColorScaleBTY();
            }

            return colorScaleBTY.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("Gray")) {

            if (colorScaleGray == null) {
                colorScaleGray = new ColorScaleGray();
            }

            return colorScaleGray.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("LinGray")) {

            if (colorScaleLinGray == null) {
                colorScaleLinGray = new ColorScaleLinGray();
            }

            return colorScaleLinGray.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("Heat")) {

            if (colorScaleHeat == null) {
                colorScaleHeat = new ColorScaleHeat();
            }

            return colorScaleHeat.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("Optimal")) {

            if (colorScaleOptimal == null) {
                colorScaleOptimal = new ColorScaleOptimal();
            }

            return colorScaleOptimal.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("LinOptimal")) {

            if (colorScaleLinOptimal == null) {
                colorScaleLinOptimal = new ColorScaleLinOptimal();
            }

            return colorScaleLinOptimal.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("Magenta")) {

            if (colorScaleMagenta == null) {
                colorScaleMagenta = new ColorScaleMagenta();
            }

            return colorScaleMagenta.getColor(value, minValue, maxValue, inverted);

        } else if (scaleName.equalsIgnoreCase("Rainbow")) {

            if (colorScaleRainbow == null) {
                colorScaleRainbow = new ColorScaleRainbow();
            }

            return colorScaleRainbow.getColor(value, minValue, maxValue, inverted);

        }

        return null;

    }

    public static String promptColorScale(String defaultColor) {

        ImageIcon icon = null;

        String s = (String) JOptionPane.showInputDialog(
                Master.mainframe,
                "Choose Color Scale:",
                "Color Scale",
                JOptionPane.PLAIN_MESSAGE,
                icon,
                colorScaleList,
                defaultColor);

        if ((s == null) || (s.length() == 0)) {
            return null; // error
        }

        return s;

    }

    ///////////////////////////////////////////////////////////////////
    //
    //     General Color Functions
    //
    ///////////////////////////////////////////////////////////////////

    public static String promptColorOrScale(String currentColor) {

        Color newColor;

        String defaultSelect = "One Color";

        if (isColorScale(currentColor)) {
            defaultSelect = "Color Scale";
        }

        String[] possibilities = {"One Color", "Color Scale"};

        ImageIcon icon = null;

        String s = (String) JOptionPane.showInputDialog(
                Master.mainframe,
                "choose type of color:",
                "Change Color",
                JOptionPane.PLAIN_MESSAGE,
                icon,
                possibilities,
                defaultSelect);

        if ((s == null) || (s.length() == 0)) {
            return null; // error
        }

        if (s.equalsIgnoreCase("Color Scale")) {

            return promptColorScale(currentColor);

        } else {

            newColor = promptColor(currentColor);

            if (newColor == null) {
                return "";
            } else {
                return newColor.toString();
            }
        }

    }

    public static Color promptColor(String currentColor) {

        Color color = string2color(currentColor);

        if (color == null) {
            color = Color.BLACK;
        }

        return promptColor(color);

    }

    public static Color promptColor(Color currentColor) {
        return JColorChooser.showDialog( null, "Choose Color", currentColor);
    }

    public static Color string2color(String colorStr) {

        if (colorStr.length() == 0) {
            return Color.white;
        }

        if (isColorScale(colorStr)) {
            return getColorScale(colorStr, 0.9, 0, 1, false);
        }

        int ired = getColorRGB(colorStr, "r=" );
        int igreen = getColorRGB(colorStr, "g=" );
        int iblue = getColorRGB(colorStr, "b=" );

        if ((ired < 0) || (ired > 255)) {
            return Color.white;
        }

        if ((igreen < 0) || (igreen > 255)) {
            return Color.white;
        }

        if ((iblue < 0) || (iblue > 255)) {
            return Color.white;
        }

        return new Color(ired, igreen, iblue);

    }

    public static Color getColorForeground(String colorStr) {

        int ired, igreen, iblue;

        if (colorStr.length() == 0) {
            return Color.black;
        }

        if (isColorScale(colorStr)) {
            return getColorScale(colorStr, 0.1, 0, 1, false);
        }

        ired = getColorRGB(colorStr, "r=" );
        igreen = getColorRGB(colorStr, "g=" );
        iblue = getColorRGB(colorStr, "b=" );

        double avg = (ired + igreen + iblue) / 3.0;

        if (avg < 255 / 2) {
            return Color.WHITE;
        } else {
            return Color.BLACK;
        }

    }

    public static int getColorRGB(String colorStr, String findRGB) {

        int ibgn, iend, itest;
        String subStr;

        if (colorStr == null) {
            return -1;
        }

        ibgn = colorStr.indexOf(findRGB);

        if (ibgn < 0) {
            return -1;
        }

        subStr = colorStr.substring(ibgn+findRGB.length());

        iend = -1;

        for (int i=0;i<subStr.length();i++) {
            try {
                Integer.parseInt(subStr.substring(i, i+1));
            } catch (NumberFormatException e) {
                iend = i;
                break;
            }
        }

        if (iend <= 0 ) {
            return -1;
        }

        try {
            itest = Integer.parseInt(subStr.substring(0, iend));
        } catch (NumberFormatException e) {
            return -1;
        }

        return itest;

    }

     // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof ColorD3D)) {
            return false;
        }

        if (n.equalsIgnoreCase("invert")) {
            invert = (v == 1);
            return true;
        }

        return false;
        
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof ColorD3D)) {
            return false;
        }

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        if (n.equalsIgnoreCase("colorRGB")) {
            return setColor(s);
        }
        if (n.equalsIgnoreCase("gradientRGB")) {
            return setGradientColor(s);
        }
        return false;
    }

}
