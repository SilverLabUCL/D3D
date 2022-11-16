package ucl.silver.d3d.utils;
import ucl.silver.d3d.utils.MoreMath;

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
public final class Utility {
   
    public static final double AVOGADRO = 6.022140857e23; // molecules per mole
    public static final double FARADAY = 96485.3321233100184; // coulombs per mole
    public static final double R = 8.31446; // gas constant // J / ( mol * K )
    public static final double ELEMENTARYCHARGE = 1.6021766208e-19; // coulombs per molecule // Faraday / Avogadro
    private final static double SQRT2 = Math.sqrt(2);
    private final static double SQRT2PI = Math.sqrt(2 * Math.PI);
    
    public static final double D_GLUTAMATE_37C = 0.33; // um^2/ms
    // Neuron. 2004 Jun 10;42(5):757-71.
    // Modulation of glutamate mobility reveals the mechanism underlying slow-rising AMPAR EPSCs and the diffusion coefficient in the synaptic cleft.
    // Nielsen TA, DiGregorio DA, Silver RA.
    
    public static final double D_CALCIUM_20C = 0.223; // um^2/ms
    // Science. 1992 Dec 11;258(5089):1812-5.
    // Range of messenger action of calcium ion and inositol 1,4,5-trisphosphate.
    // Allbritton NL, Meyer T, Stryer L.
    // D = 0.223, Temp = 20C
    
    public static final double D_CALCIUM_25C = 0.247; // um^2/ms
    // J Exp Biol. 1987 May;129:191-203.
    // Temperature affects the diffusion of small molecules through cytosol of fish muscle.
    // Sidell BD, Hazel JR.
    // D = 0.247, Temp = 25C, Q10 = 2.04 (Temp1 = 5C, Temp2 = 25C)
    // Q10 = 1.5 computed using mean values

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
    
    public static double Q2I(double Q, double charge, double liters) { // convert mM/ms (same as M/s) to pA
        double moles_per_second = Q * liters;
        double I_A = moles_per_second * FARADAY * charge;
        return I_A * 1e12; // pA
    }
    
    public static double I2Q(double I_pA, double charge, double liters) { // convert pA to mM/ms
        double I_A = I_pA * 1e-12; // coulombs/s
        double moles_per_second = I_A / ( FARADAY * charge); // m/s
        return moles_per_second / liters; // M/s (same as mM/ms)
    }
    
    public static double C2N(double c_mM, double volume_um3) { // convert mM to molecules
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
    
    public static double expNormArea(double t, double tonset, double tau) { // area = 1.0
        if (t < tonset) {
            return 0;
        }
        return (1 / tau) * Math.exp(-(t - tonset) / tau);
    }
    
    public static double expCtotal2Qpeak(double Ctotal_mM, double tau_ms) {
        return Ctotal_mM / tau_ms; // mM/ms
    }
    
    public static double expQpeak2Ctotal(double Qpeak, double tau_ms) {
        return Qpeak * tau_ms; // mM
    }
    
    public static double expCtotal2Ipeak(double Ctotal_mM, double tau_ms, double charge, double volume_um3) {
        double liters = Utility.litersPerUM3(volume_um3);
        double Qpeak = Ctotal_mM / tau_ms; // mM/ms
        return Utility.Q2I(Qpeak, charge, liters); // pA
    }

    public static double expIpeak2Ctotal(double Ipeak_pA, double tau_ms, double charge, double volume_um3) {
        double liters = Utility.litersPerUM3(volume_um3);
        double Qpeak = Utility.I2Q(Ipeak_pA, charge, liters); // mM/ms
        return Qpeak * tau_ms; // mM
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
    
    public static double gaussPDF(double t, double tpeak, double stdv) { // area = 1.0
        return (Math.exp(-Math.pow((t - tpeak) / (SQRT2 * stdv), 2)) / (stdv * SQRT2PI));
    }

    public static double gaussCDF(double t, double tpeak, double stdv) {
        double x = (t - tpeak) / (SQRT2 * stdv);
        return 0.5 * (1 + MoreMath.erf(x)); // mM
    }
    
    public static double gaussNormPeak(double t, double tpeak, double stdv) { // peak = 1.0
        return Math.exp(-Math.pow((t - tpeak) / (SQRT2 * stdv), 2));
    }
    
    public static double gaussCtotal2Qpeak(double Ctotal_mM, double stdv_ms) {
        return Ctotal_mM / (stdv_ms * SQRT2PI); // mM/ms
    }
    
    public static double gaussQpeak2Ctotal(double Qpeak, double stdv_ms) {
        return Qpeak * stdv_ms * SQRT2PI; // mM
    }
    
    public static double gaussCtotal2Ipeak(double Ctotal_mM, double stdv_ms, double charge, double volume_um3) {
        double liters = Utility.litersPerUM3(volume_um3);
        double Qpeak = Ctotal_mM / (stdv_ms * SQRT2PI); // mM/ms
        return Utility.Q2I(Qpeak, charge, liters); // pA
    }

    public static double gaussIpeak2Ctotal(double Ipeak_pA, double stdv_ms, double charge, double volume_um3) {
        double liters = Utility.litersPerUM3(volume_um3);
        double Qpeak = Utility.I2Q(Ipeak_pA, charge, liters); // mM/ms
        return Qpeak * stdv_ms * SQRT2PI; // mM
    }
    
    // Gamma distribution (PDF) with parameters alpha and beta
    // https://en.wikipedia.org/wiki/Gamma_distribution
    public static double gammaPDF(double t, double tonset, int alpha, double beta) { // area = 1.0, beta = 1/tau
        
        if ((alpha <= 0) || (beta <= 0)) {
            return Double.NaN;
        }
        
        if (t >= tonset) {
            return Math.pow(beta, alpha) * Math.pow(t - tonset, alpha - 1) * Math.exp(-beta * (t - tonset)) / MoreMath.gamma(alpha);
        } else {
            return 0;
        }
        
    }
    
    // Gamma distribution (CDF) with parameters alpha and beta
    // https://en.wikipedia.org/wiki/Gamma_distribution
    public static double gammaCDF(double t, double tonset, int alpha, double beta) { // area = 1.0, beta = 1/tau

        if ((alpha <= 0) || (beta <= 0)) {
            return Double.NaN;
        }

        if (t >= tonset) {
            return MoreMath.igam(alpha, beta * (t - tonset)) / MoreMath.gamma(alpha);
        } else {
            return 0;
        }

    }
    
    public static double gammaTimeOfPeak(int alpha, double beta) {
        return (alpha - 1) / beta;
    }
    
    public static double gammaCtotal2Qpeak(double Ctotal_mM, int alpha, double beta) {
        double timeOfPeak = Utility.gammaTimeOfPeak(alpha, beta);
        return Ctotal_mM * Utility.gammaPDF(timeOfPeak, 0, alpha, beta); // mM/ms
    }
    
    public static double gammaQpeak2Ctotal(double Qpeak, int alpha, double beta) {
        double timeOfPeak = Utility.gammaTimeOfPeak(alpha, beta);
        return Qpeak / Utility.gammaPDF(timeOfPeak, 0, alpha, beta); // mM
    }
    
    public static double gammaCtotal2Ipeak(double Ctotal_mM, int alpha, double beta, double charge, double volume_um3) {
        double timeOfPeak = Utility.gammaTimeOfPeak(alpha, beta);
        double liters = Utility.litersPerUM3(volume_um3);
        double Qpeak = Ctotal_mM * Utility.gammaPDF(timeOfPeak, 0, alpha, beta); // mM/ms
        return Utility.Q2I(Qpeak, charge, liters); // pA
    }

    public static double gammaIpeak2Ctotal(double Ipeak_pA, int alpha, double beta, double charge, double volume_um3) {
        double timeOfPeak = Utility.gammaTimeOfPeak(alpha, beta);
        double liters = Utility.litersPerUM3(volume_um3);
        double Qpeak = Utility.I2Q(Ipeak_pA, charge, liters); // mM/ms
        return Qpeak / Utility.gammaPDF(timeOfPeak, 0, alpha, beta); // mM
    }
    
    public static double chiMean(double df, double beta) {
        double gratio = MoreMath.gamma(0.5 * (df + 1.0)) / MoreMath.gamma(0.5 * df);
        return Math.sqrt(2.0 * beta) * gratio;
    }

    public static double chiStdv(double df, double beta) {
        double gratio = MoreMath.gamma(0.5 * (df + 1.0)) / MoreMath.gamma(0.5 * df);
        return Math.sqrt(beta * (df - 2.0 * gratio * gratio));
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
    
    public static double[][] histogram(double array[], double min_value, double max_value, double bin_width, int norm_flag) {
        // norm_flag = 0: count
        // norm_flag = 1: probability
        // norm_flag = 2: frequency
        // returns histo: i-bins, j=2
        //      j[0] = x-axis
        //      j[1] = y-axis
        int ibin, nbins, skipped = 0;
        double bin_half_width, count = 0;
        
        if (Double.isNaN(min_value)) {
            min_value = getMin(array);
        }
        
        if (Double.isNaN(max_value)) {
            max_value = getMin(array);
        }
        
        if (Double.isNaN(bin_width) || (bin_width <= 0)) {
            bin_width = (max_value - min_value) / 100;
        }
        
        nbins = 1 + (int)((max_value - min_value) / bin_width);
        
        if (nbins <= 0) {
            return null;
        }
        
        double[][] histo = new double[nbins][2];
        
        ibin = 0;
        bin_half_width = 0.5 * bin_width; // center x-scale on middle of bins
        
        for (double[] h : histo) {
            h[0] = min_value + bin_half_width + ibin * bin_width;
            ibin++;
        }
        
        for (double a : array) {
            
            if (!Double.isFinite(a) ) {
                continue;
            }
            
            if ((a < min_value) || (a > max_value)) {
                skipped++;
                continue;
            }
            
            ibin = (int) Math.floor((a - min_value) / bin_width);
            
            if ((ibin >= 0) && (ibin < nbins)) {
                histo[ibin][1]++;
                count++;
            } else {
                skipped++;
            }
            
        }
        
        if (skipped > 0) {
            System.out.println("warning: histogram: skipped data points N = " + skipped);
        }
        
        if (count == 0) {
            return histo;
        }
        
        if (norm_flag == 1) { // Probability
            for (double[] h : histo) {
                h[1] /= bin_width * count * 1.0;
            }
        } else if (norm_flag == 2) { // Frequency
            for (double[] h : histo) {
                h[1] /= count * 1.0;
            }
        }
        
        return histo;
        
    }
    
    public static double[][] histogramIJ(double array[][], double min_value, double max_value, double bin_width) {
        
        int ibin, nbins;
        int skipped = 0;
        double bin_half_width;
        
        if (Double.isNaN(min_value)) {
            min_value = getMinI(array, 0);
        }
        
        if (Double.isNaN(max_value)) {
            max_value = getMinI(array, 0);
        }
        
        if (Double.isNaN(bin_width) || (bin_width <= 0)) {
            bin_width = (max_value - min_value) / 100;
        }
        
        nbins = 1 + (int)((max_value - min_value) / bin_width);
        
        if (nbins <= 0) {
            return null;
        }
        
        double[][] histo = new double[nbins][4];
        
        ibin = 0;
        
        bin_half_width = 0.5 * bin_width; // center x-scale on middle of bins
        
        for (double[] h : histo) {
            h[0] = min_value + bin_half_width + ibin * bin_width;
            ibin++;
        }
        
        for (double[] a : array) {
            
            if (!Double.isFinite(a[0]) || !Double.isFinite(a[1])) {
                continue;
            }
            
            if ((a[0] < min_value) || (a[0] > max_value)) {
                skipped++;
                continue;
            }
            
            ibin = (int) Math.floor((a[0] - min_value) / bin_width);

            if ((ibin >= 0) && (ibin < nbins)) {
                histo[ibin][1] += a[1];
                histo[ibin][2] += a[1]* a[1];
                histo[ibin][3]++;
            } else {
                skipped++;
            }
            
        }
        
        double[][] histo2 = new double[nbins][4];
        
        for (int i = 0; i < histo.length; i++) {
            histo2[i][0] = histo[i][0];
            if (histo[i][3] > 0) {
                histo2[i][1] = histo[i][1] / histo[i][3]; // mean
                histo2[i][2] = Math.sqrt((histo[i][2] - histo[i][1] * histo[i][1] / histo[i][3]) / (histo[i][3] - 1)); // stdv
                histo2[i][3] = histo[i][3];
            }
        }
        
        if (skipped > 0) {
            System.out.println("warning: histogramIJ: out-of-range data points N = " + skipped);
        }
        
        return histo2;
        
    }
    
    public static int countFinite(double array[]) {

        int count = 0;

        for (double a : array) {
            if (Double.isFinite(a)) {
                count++;
            }
        }

        return count;

    }
    
    public static double average(double array[]) {

        double sum = 0, count = 0;

        for (double a : array) {
            if (Double.isFinite(a)) {
                sum += a;
                count++;
            }
        }

        if (count == 0) {
            return Double.NaN;
        }

        return (sum / count);

    }
    
    public static double stdev(double array[]) {

        double sum = 0, count = 0;
        
        double avg = Utility.average(array);
        
        if (Double.isNaN(avg)) {
            return Double.NaN;
        }

        for (double a : array) {
            if (Double.isFinite(a)) {
                sum += Math.pow(a - avg, 2);
                count++;
            }
        }

        if (count == 0) {
            return Double.NaN;
        }

        return Math.sqrt((sum / (count-1)));

    }

    // find min value of double array
    public static double getMin(double array[]) {

        double min = Double.POSITIVE_INFINITY;
        
        for (double a : array) {
            if (Double.isFinite(a)) {
                min = Math.min(min, a);
            }
        }

        return min;

    }

    // find min value of double array
    public static double getMin(double array[][]) {

        double min = Double.POSITIVE_INFINITY;
        
        for (double[] i : array) {
            for (double a : i) {
                if (Double.isFinite(a)) {
                    min = Math.min(min, a);
                }
            }
        }

        return min;

    }
    
    // find min value of double array
    public static double getMinI(double array[][], int j) {

        double min = Double.POSITIVE_INFINITY;
        
        if ((j < 0) || (j >= array[0].length)) {
            return Double.NaN;
        }
        
        for (int i = 0; i < array.length; i++) {
            if (Double.isFinite(array[i][j])) {
                min = Math.min(min, array[i][j]);
            }
        }

        return min;

    }
    
    // find min value of double array
    public static double getMinJ(double array[][], int i) {

        double min = Double.POSITIVE_INFINITY;
        
        if ((i < 0) || (i >= array.length)) {
            return Double.NaN;
        }
        
        for (int j = 0; j < array[i].length; j++) {
            if (Double.isFinite(array[i][j])) {
                min = Math.min(min, array[i][j]);
            }
        }

        return min;

    }

    // find max value of double array
    public static double getMax(double array[]) {

        double max = Double.NEGATIVE_INFINITY;
        
        for (double a : array) {
            if (Double.isFinite(a)) {
                max = Math.max(max, a);
            }
        }

        return max;

    }

    // find max value of double array
    public static double getMax(double array[][]) {

        double max = Double.NEGATIVE_INFINITY;
        
        for (double[] i : array) {
            for (double a : i) {
                if (Double.isFinite(a)) {
                    max = Math.max(max, a);
                }
            }
        }

        return max;

    }
    
    // find max value of double array
    public static double getMaxI(double array[][], int j) {

        double max = Double.NEGATIVE_INFINITY;
        
        for (int i = 0; i < array.length; i++) {
            if (Double.isFinite(array[i][j])) {
                max = Math.max(max, array[i][j]);
            }
        }

        return max;

    }
    
    // find max value of double array
    public static double getMaxJ(double array[][], int i) {

        double max = Double.NEGATIVE_INFINITY;
        
        for (int j = 0; j < array[i].length; j++) {
            if (Double.isFinite(array[i][j])) {
                max = Math.max(max, array[i][j]);
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
