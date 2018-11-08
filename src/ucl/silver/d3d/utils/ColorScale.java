package ucl.silver.d3d.utils;

import java.awt.*;

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
public class ColorScale {

    static int size = 256; // 256 defined RGB Colors
    static int lastPnt = size - 1;
    public Color[] rgb;
    //Test 4

    public ColorScale() {
        rgb = new Color[256];
    }

    public Color getColor(double value, boolean inverted) { // min = 0, max = 1
        if (inverted) {
            return getColorInverted(value);
        } else {
            return getColor(value);
        }
    }

    public Color getColor(double value) { // min = 0, max = 1
        int i = (int) (value * lastPnt);
        i = Math.min(Math.max(i, 0), lastPnt);
        return rgb[i];
    }

    public Color getColorInverted(double value) { // min = 0, max = 1
        int i = (int) (value * lastPnt);
        i = Math.min(Math.max(i, 0), lastPnt);
        i = Math.abs(i - lastPnt);
        return rgb[i];
    }

    public Color getColor(double value, double min, double max, boolean inverted) {
        if (inverted) {
            return getColorInverted(value, min, max);
        } else {
            return getColor(value, min, max);
        }
    }

    public Color getColor(double value, double min, double max) {
        value = (value - min) * lastPnt / max;
        int i = Math.min(Math.max((int) value, 0), lastPnt);
        return rgb[i];
    }

    public Color getColorInverted(double value, double min, double max) {
        value = (value - min) * lastPnt / max;
        int i = Math.min(Math.max((int) value, 0), lastPnt);
        i = Math.abs(i - lastPnt);
        return rgb[i];
    }
}
