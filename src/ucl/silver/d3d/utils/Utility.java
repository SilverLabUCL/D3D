package ucl.silver.d3d.utils;

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
public final class Utility {
   
    public static final double AVOGADRO = 6.022140857e23; // molecules per mole
    public static final double FARADAY = 96485.3329; // coulombs per mole
    public static final double ELEMENTARYCHARGE = 1.6021766208e-19; // coulombs per molecule // Faraday / Avogadro
    private final static double SQRT2 = Math.sqrt(2);
    private final static double SQRT2PI = Math.sqrt(2 * Math.PI);
    
    public static final double D_GLUTAMATE = 0.33; // um^2/ms
    public static final double D_CALCIUM = 0.22; // um^2/ms

    private Utility() {
        // cannot be instantiated
    }
    
    public static double getCharge(String name) {
        
        if (name.equalsIgnoreCase("glutamate")) {
            return -1;
        }
        
        if (name.equalsIgnoreCase("calcium") || name.equalsIgnoreCase("Ca")) {
            return 2;
        }
        
        if (name.equalsIgnoreCase("potassium") || name.equalsIgnoreCase("K")) {
            return 1;
        } 
        
        if (name.equalsIgnoreCase("sodium") || name.equalsIgnoreCase("Na")) {
            return 1;
        }
        
        if (name.equalsIgnoreCase("chloride") || name.equalsIgnoreCase("Cl")) {
            return 1;
        }
        
        if (name.equalsIgnoreCase("magnesium") || name.equalsIgnoreCase("Mg")) {
            return 2;
        }
        
        if (name.equalsIgnoreCase("ATP")) {
            return -3.5;
        }
        
        return 0;
        
    }
    
    public static double litersPerVoxel(double cubeWidth_um) { // cube width (um)
        return Math.pow(cubeWidth_um, 3) * 1e-18 * 1000; // 1L = 1e-3 m^3
    }
    
    public static double litersPerUM3(double volume_um3) { // volume (um^3)
        return volume_um3 * 1e-18 * 1000; // 1L = 1e-3 m^3
    }
    
    public static double C2N(double c_mM, double volume_um3) { // convert mM to molecules
        // volume is um^3
        
        double liters = Utility.litersPerUM3(volume_um3);
        double moles = (c_mM * 1e-3) * liters;
        
        return moles * Utility.AVOGADRO;

    }
    
    public static double C2N(double c_mM, double dx_um, int numVoxels) { // convert mM to molecules
        
        double liters = numVoxels * Utility.litersPerVoxel(dx_um);
        double moles = (c_mM * 1e-3) * liters;
        
        return moles * Utility.AVOGADRO;

    }
    
    public static double N2C(double n, double volume_um3) { // convert molecules to mM
        
        double liters = Utility.litersPerUM3(volume_um3);
        double moles = n / Utility.AVOGADRO;
        
        return moles * 1e3 / liters; // mM

    }
    
    public static double N2C(double n, double dx_um, int numVoxels) { // convert molecules to mM
        
        double liters = numVoxels * Utility.litersPerVoxel(dx_um);
        double moles = n / Utility.AVOGADRO;
        
        return moles * 1e3 / liters; // mM

    }
    
    public static double gaussFWHM2STDV(double FWHM) {
        double cfactor = 2 * Math.sqrt(2 * Math.log(2)); // 2.35482
        if (FWHM < 0) {
            return -1;
        } else {
            return FWHM / cfactor;
            // http://hyperphysics.phy-astr.gsu.edu/hbase/math/gaufcn2.html
        }
    }

    public static double gaussSTDV2FWHM(double STDV) {
        double cfactor = 2 * Math.sqrt(2 * Math.log(2)); // 2.35482
        if (STDV < 0) {
            return -1;
        } else {
            return STDV * cfactor;
        }
    }
    
    public static double gaussNormArea(double t, double tpeak, double stdv) { // area = 1.0
        return (Math.exp(-Math.pow((t - tpeak) / (SQRT2 * stdv), 2)) / (stdv * SQRT2PI));
    }
    
    public static double gaussNormPeak(double t, double tpeak, double stdv) { // peak = 1.0
        return Math.exp(-Math.pow((t - tpeak) / (SQRT2 * stdv), 2));
    }
    
    public static double gamma(double t, double tonset, int alpha, double beta) { // area = 1.0
        double gamma = 1;
        int j = alpha - 1;
        
        if ((alpha <= 0) || (beta <= 0)) {
            return Double.NaN;
        }
        
        switch (alpha) {
            case 1:
            case 2:
                break;
            case 3:
                gamma = 2;
                break;
            case 4:
                gamma = 6;
                break;
            default:
                for (int i = j; i <= 1; i--) {
                    gamma *= i; // factorial
                }
        }
        
        if (t > tonset) {
            return Math.pow(beta, alpha) * Math.pow(t - tonset, alpha - 1) * Math.exp(-beta * (t - tonset)) / gamma;
        } else {
            return 0;
        }
        
    }

    public static int makeOdd(int i) {
        if (Math.IEEEremainder(i, 2) == 0) {
            return i + 1;
        }
        return i;
    }

    public static int makeEven(int i) {
        if (Math.IEEEremainder(i, 2) == 0) {
            return i;
        }
        return i + 1;
    }

    public static String parentClass(String longNameWithPeriod) { // e.g. "Detector0.PSF.numericalAperture"

        if (longNameWithPeriod == null) {
            return null;
        }

        String[] splitStr = longNameWithPeriod.split("\\.");

        if (splitStr.length >= 2) {
            return splitStr[0]; // e.g. "Detector0"
        }

        return null;

    }

    public static String childClass(String longNameWithPeriod) { // e.g. "Detector0.PSF.numericalAperture"

        if (longNameWithPeriod == null) {
            return null;
        }

        String[] splitStr = longNameWithPeriod.split("\\.");

        if (splitStr.length == 3) {
            return splitStr[1]; // e.g. "PSF"
        }

        return null;

    }

    public static String parameterName(String longName) { // e.g. "Detector0.PSF.numericalAperture"

        if (longName == null) {
            return null;
        }

        String[] splitStr = longName.split("\\.");

        if (splitStr.length > 1) {
            return splitStr[splitStr.length - 1]; // e.g. "numericalAperture"
        }

        return null;

    }

    public static int itemsInList(String strList, char stop) {

        int count;

        if ((strList == null) || (strList.length() == 0)) {
            return 0;
        }

        count = 1;

        for (int i = 0; i < strList.length() - 1; i += 1) {
            if (strList.charAt(i) == stop) {
                count += 1;
            }
        }

        return count;

    }

    public static String itemFromList(String strList, int item) {

        int count = 0, start, end = -1;

        if ((strList == null) || (strList.length() == 0)) {
            return "";
        }

        if (strList.indexOf(",") < 0) {
            return strList;
        }

        for (int i = 0; i < strList.length(); i += 1) {
            if (strList.charAt(i) == ',') {
                count += 1;
                start = end + 1;
                end = i;
                if (item == count - 1) {
                    return strList.substring(start, end);
                }
            }
        }

        start = end + 1;
        end = strList.length();

        return strList.substring(start, end);

    }

    public static int whichListItem(String[] strList, String itemStr, boolean ignoreCase) {

        if ((strList == null) || (strList.length == 0)) {
            return -1;
        }

        if (itemStr == null) {
            return -1;
        }

        for (int i = 0; i < strList.length; i++) {

            if (strList[i] == null) {
                continue;
            }

            if (ignoreCase) {
                if (strList[i].equalsIgnoreCase(itemStr)) {
                    return i;
                }
            } else {
                if (strList[i].equals(itemStr)) {
                    return i;
                }
            }

        }

        return -1;

    }

    // limit integer value
    public static int limit(int value, int imin, int imax) {
        return Math.min(imax, Math.max(value, imin));
    }

    // find min value of double array
    public static double getMin(double array[]) {

        double min = Double.POSITIVE_INFINITY;
        
        for (double a : array) {
            min = Math.min(min, a);
        }

        return min;

    }

    // find min value of double array
    public static double getMin(double array[][]) {

        double min = Double.POSITIVE_INFINITY;
        
        for (double[] i : array) {
            for (double a : i ) {
                min = Math.min(min, a);
            }
        }

        return min;

    }

    // find max value of double array
    public static double getMax(double array[]) {

        double max = Double.NEGATIVE_INFINITY;
        
        for (double a : array) {
            max = Math.max(max, a);
        }

        return max;

    }

    // find max value of double array
    public static double getMax(double array[][]) {

        double max = Double.NEGATIVE_INFINITY;
        
        for (double[] i : array) {
            for (double a : i ) {
                max = Math.max(max, a);
            }
        }

        return max;

    }

    public static int[] resizeArray(int[] array, int iNew) {

        int iNum = 0;

        if (array != null) {
            iNum = array.length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (iNew <= 0) {
            return null;
        }

        int[] temp = new int[iNew];

        for (int i = 0; i < iNew; i++) {
            if (i >= iNum) {
                break;
            }
            temp[i] = array[i];
        }

        return temp;

    }

    public static byte[] resizeArray(byte[] array, int iNew) {

        int iNum = 0;

        if (array != null) {
            iNum = array.length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (iNew <= 0) {
            return null;
        }

        byte[] temp = new byte[iNew];

        for (int i = 0; i < iNew; i++) {
            if (i >= iNum) {
                break;
            }
            temp[i] = array[i];
        }

        return temp;

    }

    public static double[] resizeArray(double[] array, int iNew) {

        int iNum = 0;

        if (array != null) {
            iNum = array.length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (iNew <= 0) {
            return null;
        }

        double[] temp = new double[iNew];

        for (int i = 0; i < iNew; i++) {
            if (i >= iNum) {
                break;
            }
            temp[i] = array[i];
        }

        return temp;

    }

    public static double[][][] resizeArray(double[][][] array, int iNew, int jNew, int kNew) {

        int i, j, k;

        int iNum = 0, jNum = 0, kNum = 0;

        if (array != null) {
            iNum = array.length;
            jNum = array[0].length;
            kNum = array[0][0].length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (jNew < 0) {
            jNew = jNum;
        }

        if (kNew < 0) {
            kNew = kNum;
        }

        if ((iNew <= 0) || (jNew <= 0) || (kNew <= 0)) {
            return null;
        }

        double[][][] temp = new double[iNew][jNew][kNew];

        for (k = 0; k < kNew; k++) {
            if (k >= kNum) {
                break;
            }
            for (j = 0; j < jNew; j++) {
                if (j >= jNum) {
                    break;
                }
                for (i = 0; i < iNew; i++) {
                    if (i >= iNum) {
                        break;
                    }
                    temp[i][j][k] = array[i][j][k];
                }
            }
        }

        return temp;

    }

    // resize byte array matrix
    public static int[][][] resizeArray(int[][][] array, int iNew, int jNew, int kNew) {

        int i, j, k;

        int iNum = 0, jNum = 0, kNum = 0;

        if (array != null) {
            iNum = array.length;
            jNum = array[0].length;
            kNum = array[0][0].length;
        }

        if (iNew < 0) {
            iNew = iNum;
        }

        if (jNew < 0) {
            jNew = jNum;
        }

        if (kNew < 0) {
            kNew = kNum;
        }

        if ((iNew <= 0) || (jNew <= 0) || (kNew <= 0)) {
            return null;
        }

        int[][][] temp = new int[iNew][jNew][kNew];

        for (k = 0; k < kNew; k++) {
            if (k >= kNum) {
                break;
            }
            for (j = 0; j < jNew; j++) {
                if (j >= jNum) {
                    break;
                }
                for (i = 0; i < iNew; i++) {
                    if (i >= iNum) {
                        break;
                    }
                    temp[i][j][k] = array[i][j][k];
                }
            }
        }

        return temp;

    }
}
