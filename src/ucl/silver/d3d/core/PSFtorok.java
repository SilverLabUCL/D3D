package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;
import java.io.*;

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
public class PSFtorok extends PSF {

    public int integrationSteps = 1000; // numerical integration steps
    public double numericalAperture = 1.0; // numerical aperture
    public double waveLength = 0.351; // light wavelength
    public double refractiveIndex = 1.338; // refractive index of medium
    public double ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8, ab9, ab10, abZoff, zFxnOffset; // aberration coefficients

    public double rscale = 1.0;
    // variables for computing Torok fxn

    private double alpha, dtheta, tk0, tk1;
    private Complex c_dtheta, c_2, i0, i1, i2, I0, I1, I2;

    // variables for reading/writing random access (ra) files

    boolean randomAccess = false;
    private int BYTES_HEADER = 107; // bytes in PSF file header
    int ra_xVoxels, ra_yVoxels, ra_zVoxels; // size of random access file matrix
    double ra_x0, ra_y0, ra_z0; // psf center point
    double ra_ioff, ra_joff, ra_koff, ra_scale = 1; // offset reading values
    
    private RandomAccessFile raFile = null;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("waveLength")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    public PSFtorok(Project p, CoordinatesVoxels c) {
        super(p, c);
        //xySymmetric = true; // should manually set to true
        name = "Torok PSF";
        createVector(true);
    }

    public void setFitDavid() { // This fit turns out to be wrong, DO NOT USE
        waveLength = 0.351;
        numericalAperture = 0.915;
        ab1 = -4.7709;
        ab2 = 31.997;
        ab3 = -234.06;
        ab4 = 959.61;
        ab5 = -2459.8;
        ab6 = 4008.5;
        ab7 = -3847.9;
        ab8 = 1521.4;
        zFxnOffset = 0.065; // fixes small z-offset of David's PSF data
        abZoff = 4.06 + zFxnOffset; // adjusted to center z axis
    }

    public void setFitJason() { // recomputed fit with 1000 integration Steps
        waveLength = 0.351;
        numericalAperture = 0.922393;
        ab1 = -2.106120;
        ab2 = 51.335700;
        ab3 = -387.978;
        ab4 = 1530.59;
        ab5 = -3020.46;
        ab6 = 1944.49;
        ab7 = 2155.03;
        ab8 = -3005.869809;
        zFxnOffset = 0.065; // fixes small z-offset of David's PSF data
        abZoff = -1.278070 + zFxnOffset; // adjusted to center z axis
    }

    public void setFitZoltan() {
        waveLength = 0.488;
        numericalAperture = 0.86782;
        ab1 = 0.30314;
        ab2 = -0.19652;
        ab3 = -11.556;
        ab4 = 123.73;
        ab5 = -553.8;
        ab6 = 738.29;
        ab7 = 1274.7;
        ab8 = -3099;
        abZoff = 0;
    }

    public void setFitZoltan2() {
        waveLength = 0.488;
        numericalAperture = 0.91182;
        ab1 = 3.62;
        ab2 = -52.092;
        ab3 = 372.71;
        ab4 = -1234.4;
        ab5 = 939.45;
        ab6 = 4736.9;
        ab7 = -11943;
        ab8 = 8010.6;
        abZoff = 0;
    }

    public void setFitZoltan3() { // this is best fit to Zoltan's FRAP PSF
        waveLength = 0.488;
        numericalAperture = 0.90643;
        ab1 = -0.15507;
        ab2 = 13.488;
        ab3 = -154.08;
        ab4 = 839.65;
        ab5 = -2324.9;
        ab6 = 2568.2;
        ab7 = 1049.8;
        ab8 = -3276.1;
        abZoff = 0;
    }

    public void setFitZoltan4() {
        waveLength = 0.488;
        numericalAperture = 0.9523;
        ab1 = 3.8998;
        ab2 = -50.438;
        ab3 = 319.63;
        ab4 = -801.59;
        ab5 = -749.57;
        ab6 = 8225.6;
        ab7 = -15721;
        ab8 = 9855.6;
        abZoff = 0;
    }

    public void setTest() {
        waveLength = 0.488;
        numericalAperture = 0.91182;
    }

    public void check(boolean recompute) {
        if (recompute || (array == null)) {
            if ((fileName.length() > 0) && randomAccess) {
                array = readRandomAccess(fileName, coordinates().xVoxels, coordinates().yVoxels, coordinates().zVoxels);
                return;
            }
            compute();
        }
    }

    @Override
    public double getValue(int i, int j, int k) {
        if (exists()) {
            return getArrayValue(i, j, k); // already computed, get from 3D array
        } else if (randomAccess && (fileName.length() > 0)) {
            if (raFile == null) {
                initRandomAccess(fileName);
            }
            return readRandomAccessVoxelOffset(i, j, k); // get from random access file
        } else {
            return computeValue(i, j, k, normalize); // else compute the value
        }
    }

    @Override
    public void init() {

        alpha = Math.asin(numericalAperture / refractiveIndex);
        dtheta = alpha / integrationSteps;
        tk0 = 2 * Math.PI / waveLength;
        tk1 = 2 * Math.PI * refractiveIndex / waveLength;
        c_dtheta = new Complex(dtheta, 0);
        c_2 = new Complex(2.0, 0);
        //set("xySymmetric", true);

        super.init();

    }

    @Override
    public double computeVoxel(double i, double j, double k) {

        double ii, jj, kk;
        double r, z, phi, rxz, axz, ryz, ayz, zoff = abZoff;
        double rotatexz = 0;
        double rotateyz = 0;

        dx = project.dx;

        //rotatexz = 2 * Math.PI / 8; // 45 degrees
        //rotateyz = 2 * Math.PI / 8; // 45 degrees

        ii = i - xVoxelCenter();
        jj = j - yVoxelCenter();
        kk = k - zVoxelCenter();

        if (rotatexz > 0) {

            // convert to cylindrical coordinates
            rxz = Math.sqrt(Math.pow(ii, 2) + Math.pow(kk, 2));
            axz = Math.atan2(kk, ii);

            axz += rotatexz; // rotate 45 degrees in xz plane

            // now convert back
            ii = rxz * Math.cos(axz);
            kk = rxz * Math.sin(axz);

        }

        if (rotateyz > 0) {

            // convert to cylindrical coordinates
            ryz = Math.sqrt(Math.pow(jj, 2) + Math.pow(kk, 2));
            ayz = Math.atan2(kk, jj);

            ayz += rotateyz; // rotate 45 degrees in yz plane

            // now convert back
            jj = ryz * Math.cos(ayz);
            kk = ryz * Math.sin(ayz);

        }

        r = Math.sqrt(Math.pow(ii, 2) + Math.pow(jj, 2)) * dx;

        z = (kk * dx) - zoff;

        phi = getPhi(ii, jj);

        integrate_Torok_New(r * rscale, z); // trapezoid method

        return densityE(I0, I1, I2, phi);
        //return densityEB(I0, I1, I2);
        //return PoyntingVector(I0, I1, I2);

    }

    double getPhi(double i, double j) {
        double phi = 0;

        if (i == 0) {
            if (j == 0) {
                return 0;
            } else if (j > 0) {
                return (Math.PI / 2.0);
            } else if (j < 0) {
                return (3 * Math.PI / 2.0);
            }
        } else if (i > 0) {
            return Math.atan(j / i);
        } else if (i < 0) {
            return (Math.PI + Math.atan(j / i));
        }

        return phi;

    }

    private void integrate_Torok_Old(double r, double z) {

        I0 = new Complex(0, 0);
        I1 = new Complex(0, 0);
        I2 = new Complex(0, 0);

        for (int itheta = 0; itheta <= integrationSteps; itheta += 1) {
            Torok(r, z, itheta);
            I0 = I0.plus(i0);
            I1 = I1.plus(i1);
            I2 = I2.plus(i2);
        }

    }

    private void integrate_Torok_New(double r, double z) { // trapezoid integration

        Complex c0, c1, c2;

        I0 = new Complex(0, 0);
        I1 = new Complex(0, 0);
        I2 = new Complex(0, 0);

        Torok(r, z, 0);

        c0 = i0;
        c1 = i1;
        c2 = i2;

        for (int itheta = 0; itheta < integrationSteps; itheta += 1) {

            Torok(r, z, itheta + 1);

            c0 = c0.plus(i0);
            c0 = c0.div(c_2);

            c1 = c1.plus(i1);
            c1 = c1.div(c_2);

            c2 = c2.plus(i2);
            c2 = c2.div(c_2);

            I0 = I0.plus(c0);
            I1 = I1.plus(c1);
            I2 = I2.plus(c2);

            c0 = i0;
            c1 = i1;
            c2 = i2;

        }

    }

    // basic computation of Torok PSF
    private void Torok(double r, double z, int itheta) {
        double theta, sinT, cosT, psi;

        Complex arg, arg1, arg2;

        theta = itheta * dtheta;

        sinT = Math.sin(theta);
        cosT = Math.cos(theta);

        //taus = 2 * sinT * cosT / (Math.sin(2 * theta));
        //taup = 2 * sinT * cosT / (Math.sin(2 * theta) * Math.cos(0));

        psi = ab1 * Math.pow(sinT, 2) +
                ab2 * Math.pow(sinT, 4) +
                ab3 * Math.pow(sinT, 6) +
                ab4 * Math.pow(sinT, 8) +
                ab5 * Math.pow(sinT, 10) +
                ab6 * Math.pow(sinT, 12) +
                ab7 * Math.pow(sinT, 14) +
                ab8 * Math.pow(sinT, 16) +
                ab9 * Math.pow(sinT, 18) +
                ab10 * Math.pow(sinT, 20);

        // from Sheppard & Torok, 1996

        arg = new Complex(Math.sqrt(cosT) * sinT, 0);
        arg1 = new Complex(0, tk0 * psi);
        arg1 = arg1.exp();
        arg2 = new Complex(0, tk1 * z * cosT);
        arg2 = arg2.exp();
        arg = c_dtheta.times(arg.times(arg1.times(arg2)));

        arg1 = new Complex((1 + 1 * cosT) * MoreMath.jn(0, tk1 * r * sinT), 0);
        i0 = arg.times(arg1);

        arg1 = new Complex((1 * sinT) * MoreMath.jn(1, tk1 * r * sinT), 0);
        i1 = arg.times(arg1);

        arg1 = new Complex((1 - cosT) * MoreMath.jn(2, tk1 * r * sinT), 0);
        i2 = arg.times(arg1);

    }

    double densityEB(Complex I0, Complex I1, Complex I2) {

        return (Math.pow(I0.mod(), 2) + 2 * Math.pow(I1.mod(), 2) +
                Math.pow(I2.mod(), 2));

    }

    double densityE(Complex I0, Complex I1, Complex I2, double phi) {

        double d;

        Complex c = I2.conj();

        c = c.times(I0);

        d = Math.pow(I0.mod(), 2) + Math.pow(I2.mod(), 2);

        d += 2 * (1 + Math.cos(2 * phi)) * Math.pow(I1.mod(), 2);

        d += 2 * Math.cos(2 * phi) * c.real();

        return d;

    }

    double PoyntingVector(Complex I0, Complex I1, Complex I2) {

        double pp, q;

        pp = (Math.pow(I0.mod(), 2) - Math.pow(I2.mod(), 2)) / 2;

        Complex cc = I2.conj();
        cc = cc.minus(I0.conj());
        cc = I1.times(cc);

        q = cc.imag();

        return (Math.sqrt(Math.pow(pp, 2) + Math.pow(q, 2)));

    }

    // write PSF to a file
    public void writeRandomAccess(String fName, int iNum, int jNum, int kNum) {

        double v;
        double xvCenter, yvCenter, zvCenter;

        fileName = fName;

        StopWatch timer = new StopWatch();

        if (fileName.length() == 0) {
            return;
        }

        xvCenter = xVoxelCenter();
        yvCenter = yVoxelCenter();
        zvCenter = zVoxelCenter();

        try {

            raFile = new RandomAccessFile(new File(fileName), "rw");

            Master.log("Writing to file: " + fileName);
            //Master.log("D3D version: " + m.D3Dversion[0] +
            //                   m.D3Dversion[1] +
            //                   m.D3Dversion[2]);

            // header variables

            raFile.writeChar('D');
            raFile.writeChar('3');
            raFile.writeChar('D');
            raFile.writeChar('v');
            raFile.writeChar('1');
            raFile.writeChar('.');
            raFile.writeChar('2');

            raFile.writeInt(iNum);
            raFile.writeInt(jNum);
            raFile.writeInt(kNum);
            raFile.writeFloat((float) project.dx);

            raFile.writeFloat((float) numericalAperture);
            raFile.writeFloat((float) refractiveIndex);
            raFile.writeFloat((float) waveLength);

            raFile.writeFloat((float) xvCenter);
            raFile.writeFloat((float) yvCenter);
            raFile.writeFloat((float) zvCenter);

            raFile.writeFloat((float) ab1);
            raFile.writeFloat((float) ab2);
            raFile.writeFloat((float) ab3);
            raFile.writeFloat((float) ab4);
            raFile.writeFloat((float) ab5);
            raFile.writeFloat((float) ab6);
            raFile.writeFloat((float) ab7);
            raFile.writeFloat((float) ab8);
            raFile.writeFloat((float) ab9);
            raFile.writeFloat((float) ab10);
            raFile.writeFloat((float) abZoff);

            raFile.writeInt(integrationSteps);

            raFile.writeBoolean(xySymmetric());

            //Master.log("Random access header bytes: " + raFile.getFilePointer());

            // compute and write PSF

            max = computeValue(xvCenter, yvCenter, zvCenter, false);

            Master.log("PSF peak intensity: " + max);

            timer.start();

            for (int k = 0; k < kNum; k++) {

                if ((k > 0) && Math.IEEEremainder(k, zCount) == 0) {
                    timer.stop();
                    Master.log("writing k: " + k + ", t:" + timer.toString());
                    timer.start();
                }

                for (int j = 0; j < jNum; j++) {
                    for (int i = 0; i < iNum; i++) {
                        if (array != null) {
                            raFile.writeFloat((float) array[i][j][k]);
                        } else {
                            v = computeValue(i, j, k, normalize);
                            raFile.writeFloat((float) v);
                        }
                    }
                }
            }

            timer.stop();

            raFile.close();

            Master.log("Finished writing PSF");

        } catch (IOException e) {
            System.err.println("unable to write file: " + fileName);
        }

    }

    public boolean initRandomAccess(String file) {

        if (!readRandomAccessHeader(file)) {
            return false;
        }

        try {
            raFile = new RandomAccessFile(new File(file), "r"); // open file
        } catch (IOException e) {
            System.err.println("unable to read file: " + file);
            return false;
        }

        ra_ioff = (int) (ra_x0 - xVoxelCenter()); // offsets
        ra_joff = (int) (ra_y0 - yVoxelCenter());
        ra_koff = (int) (ra_z0 - zVoxelCenter());

        Master.log("initializaed random access file: " + file);

        return true;

    }

    public double[][][] readRandomAccess(String file, int imax, int yVoxels,
            int zVoxels) {
        int ii, jj;
        double pmax = 0, im = 0, jm = 0, km = 0;
        boolean error = false;

        if (file.length() == 0) {
            return null;
        }

        if (!initRandomAccess(file)) {
            return null; // error
        }

        Master.log("Reading Random Access...");

        double[][][] pp = new double[imax][yVoxels][zVoxels];

        for (int k = 0; k < zVoxels; k++) {

            if (error) {
                break;
            }

            if ((k > 0) && Math.IEEEremainder(k, zCount) == 0) {
                Master.log("reading k: " + k);
            }

            for (int j = 0; j < yVoxels; j++) {

                if (error) {
                    break;
                }

                for (int i = 0; i < imax; i++) {

                    if (xySymmetric()) {

                        ii = (int) symI(i);
                        jj = (int) symJ(j);

                        if ((ii >= 0) && (jj >= 0) && (k >= 0) && (pp[ii][jj][k] > 0)) {
                            pp[i][j][k] = pp[ii][jj][k]; // x-y symmetrical, already computed
                            //p[i][j][k] = 0; // testing only
                            continue;
                        }

                    }

                    pp[i][j][k] = readRandomAccessVoxelOffset(i, j, k); // otherwise we need to compute

                    if (pp[i][j][k] < 0) {
                        error = true;
                    }

                    if (pp[i][j][k] > pmax) {
                        pmax = pp[i][j][k];
                        im = i;
                        jm = j;
                        km = k;
                    }

                }
            }
        }

        Master.log("pmax: " + pmax + " at: " + im + ", " + jm + ", " + km);

        closeRandomAccess();

        return pp;

    }

    public double readRandomAccessVoxelOffset(int i, int j, int k) {
        int ii, jj;

        if ((i < 0) || (j < 0) || (k < 0)) {
            closeRandomAccess();
            return -1;
        }

        if (raFile == null) {
            return -1;
        }

        ii = i;
        jj = j;

        if (xySymmetric()) {
            ii = (int) symI(i);
            jj = (int) symJ(j);
        }

        ii += ra_ioff;
        jj += ra_joff;
        k += ra_koff;

        if ((ii < 0) || (ii >= ra_xVoxels)) {
            Master.log("Attempted to read Out Of Bounds i:" + ii);
            return -1;
        }

        if ((jj < 0) || (jj >= ra_yVoxels)) {
            Master.log("Attempted to read Out Of Bounds j:" + jj);
            return -1;
        }

        if ((k < 0) || (k >= ra_zVoxels)) {
            Master.log("Attempted to read Out Of Bounds k:" + k);
            return -1;
        }

        try {
            raFile.seek(getRandomAccessPointer(ii, jj, k));
            return (raFile.readFloat() * ra_scale);
        } catch (IOException e) {
            return -1;
        }

    }

    private long getRandomAccessPointer(int i, int j, int k) {
        return (BYTES_HEADER +
                4 * (i + j * ra_xVoxels + k * ra_xVoxels * ra_yVoxels));
    }

    public void closeRandomAccess() {
        try {
            raFile.close();
            raFile = null;
            Master.log("Closed PSF file: " + fileName);
        } catch (IOException e) {
        }
    }

    public boolean readRandomAccessHeader(String file) {

        double dtemp;
        double xvCenter, yvCenter, zvCenter;
        char c;
        char[] v = {'0', '0', '0'};
        boolean newVersion = true;

        String vname;

        if (file.length() == 0) {
            return false;
        }

        xvCenter = xVoxelCenter();
        yvCenter = yVoxelCenter();
        zvCenter = zVoxelCenter();

        try {

            raFile = new RandomAccessFile(new File(file), "r"); // open file
            Master.log("Reading PSF Header: " + file);

            c = raFile.readChar();

            if (c == 'D') { // new version
                c = raFile.readChar(); // 3
                c = raFile.readChar(); // D
                c = raFile.readChar(); // v
                v[0] = raFile.readChar();
                v[1] = raFile.readChar();
                v[2] = raFile.readChar();
                BYTES_HEADER = 107;
                Master.log("Encountered new random access file format.");
            } else { // older version
                newVersion = false;
                raFile.seek(0); // start at beginning
                BYTES_HEADER = 92;
                Master.log("Encountered old random access file format.");
            }

            Master.log("PSF D3Dv" + v[0] + v[1] + v[2]);

            ra_xVoxels = raFile.readInt();
            ra_yVoxels = raFile.readInt();
            ra_zVoxels = raFile.readInt();

            Master.log("PSF ra_xVoxels: " + ra_xVoxels);
            Master.log("PSF ra_yVoxels: " + ra_yVoxels);
            Master.log("PSF ra_zVoxels: " + ra_zVoxels);

            vname = "dx";
            dtemp = raFile.readFloat();
            dx = dtemp;
            Master.log("PSF " + vname + ": " + dtemp);

            vname = "numericalAperture";
            dtemp = raFile.readFloat();
            numericalAperture = dtemp;
            //Master.log(vname + ": " + dtemp);

            vname = "refractiveIndex";
            dtemp = raFile.readFloat();
            refractiveIndex = dtemp;
            //Master.log(vname + ": " + dtemp);

            vname = "waveLength";
            dtemp = raFile.readFloat();
            waveLength = dtemp;
            //Master.log(vname + ": " + dtemp);

            ra_x0 = raFile.readFloat();
            ra_y0 = raFile.readFloat();
            ra_z0 = raFile.readFloat();

            //rotZ = (int) raFile.readFloat();
            raFile.readFloat();

            Master.log("PSF X0: " + ra_x0);
            Master.log("PSF Y0: " + ra_y0);
            Master.log("PSF Z0: " + ra_z0);
            //Master.log("PSF RotZ: " + rotZ);

            for (int i = 1; i <= 10; i++) {
                vname = "ab" + Integer.toString(i);
                dtemp = raFile.readFloat();
                set(vname, dtemp);
                //Master.log(vname + ": " + dtemp);
            }

            vname = "abZoff";
            dtemp = raFile.readFloat();
            abZoff = dtemp;
            Master.log(vname + ": " + dtemp);

            dtemp = raFile.readInt();
            integrationSteps = (int) dtemp;
            //Master.log("integrationSteps: " + raFile.readInt());

            if (newVersion) {
                set("xySymmetric", raFile.readBoolean());
                //Master.log("xy symmetry : " + xySymmetric);
            }

            ra_scale = 1;

            if ((abZoff > 4.04) && (abZoff < 4.08)) {

                // this is David's PSF, abZoff = 4.06
                // need to shift slightly to center at peak

                ra_scale = computeValue(xvCenter, yvCenter, zvCenter, normalize); // old peak value

                abZoff += zFxnOffset; // adjust z-offset

                Master.log("Random access adjusted z-offset: " + zFxnOffset);

                ra_scale /= computeValue(xvCenter, yvCenter, zvCenter, normalize); // new peak value

                ra_scale = 1.009; // THIS IS OLD SCALE WHICH IS WRONG (Keep for now)

                ra_z0 -= (int) (zFxnOffset / dx); // adjust z-offset

            }

            //Master.log("Random access header bytes: " + raFile.getFilePointer());

            raFile = null;

            return true;

        } catch (IOException e) {
            System.err.println("unable to read file: " + file);
            return false;
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }
        
        if (!(o.paramVector instanceof PSFtorok)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("integrationSteps")) {
            if (v < 5) {
                return false;
            }
            integrationSteps = (int) v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("numericalAperture")) {
            if (v < 0) {
                return false;
            }
            numericalAperture = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("waveLength")) {
            if (v < 0) {
                return false;
            }
            waveLength = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("refractiveIndex")) {
            if (v < 0) {
                return false;
            }
            refractiveIndex = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab1")) {
            ab1 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab2")) {
            ab2 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab3")) {
            ab3 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab4")) {
            ab4 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab5")) {
            ab5 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab6")) {
            ab6 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab7")) {
            ab7 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab8")) {
            ab8 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab9")) {
            ab9 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("ab10")) {
            ab10 = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("abZoff")) {
            abZoff = v;
            array = null;
            return true;
        }
        return super.setMyParams(o, v);
    }
}
