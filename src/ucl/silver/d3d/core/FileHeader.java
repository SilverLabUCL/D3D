package ucl.silver.d3d.core;

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
public class FileHeader {

    public double D3Dversion;
    public String name;
    public int x1, y1, z1, x2, y2, z2;
    public int pointsPerSample;
    public int samples2save;
    public double dx, dt;
    public String xdim;
    public String ydim;
    public String format;
    public int bytes;

    public static boolean readTextHeader(String file, boolean print, FileHeader header) {

        RandomAccessFile raFile = null;

        int counter = 0;
        int maxLines = 1000;
        String strValue;
        boolean readFile = true;

        try {

            raFile = new RandomAccessFile(new File(file), "r"); // open file

            Master.log("reading File: " + file);

            do {

                strValue = raFile.readLine();

                if (strValue == null) {
                    break;
                }

                if (print) {
                    Master.log(strValue);
                }

                counter++;

                if (strValue.startsWith("D3D=")) {
                    strValue = strValue.replaceAll("D3D=", "");
                    header.D3Dversion = Double.parseDouble(strValue);
                }

                if (strValue.startsWith("name=")) {
                    strValue = strValue.replaceAll("name=", "");
                    header.name = strValue;
                }

                if (strValue.startsWith("x1=")) {
                    strValue = strValue.replaceAll("x1=", "");
                    header.x1 = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("y1=")) {
                    strValue = strValue.replaceAll("y1=", "");
                    header.y1 = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("z1=")) {
                    strValue = strValue.replaceAll("z1=", "");
                    header.z1 = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("x2=")) {
                    strValue = strValue.replaceAll("x2=", "");
                    header.x2 = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("y2=")) {
                    strValue = strValue.replaceAll("y2=", "");
                    header.y2 = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("z2=")) {
                    strValue = strValue.replaceAll("z2=", "");
                    header.z2 = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("pointsPerSample=")) {
                    strValue = strValue.replaceAll("pointsPerSample=", "");
                    header.pointsPerSample = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("samples2save=")) {
                    strValue = strValue.replaceAll("samples2save=", "");
                    header.samples2save = Integer.parseInt(strValue);
                }

                if (strValue.startsWith("dx=")) {
                    strValue = strValue.replaceAll("dx=", "");
                    header.dx = Double.parseDouble(strValue);
                }

                if (strValue.startsWith("dt=")) {
                    strValue = strValue.replaceAll("dt=", "");
                    header.dt = Double.parseDouble(strValue);
                }

                if (strValue.startsWith("format=")) {
                    strValue = strValue.replaceAll("format=", "");
                    header.format = strValue;
                    readFile = false;
                }

                if (counter >= maxLines) {
                    readFile = false;
                }

            } while (readFile);

        } catch (IOException e) {
            System.err.println("unable to read file: " + file);
            return false;
        }

        try {
            header.bytes = (int) raFile.getFilePointer();
            raFile.close();
            raFile = null;
            Master.log("closed file: " + file);
        } catch (IOException e) {
        }

        return true;

    }

    public void print() {
        Master.log("D3D=" + D3Dversion);
        Master.log("name=" + name);
        Master.log("x1=" + x1);
        Master.log("y1=" + y1);
        Master.log("z1=" + z1);
        Master.log("x2=" + x2);
        Master.log("y2=" + y2);
        Master.log("z2=" + z2);
        Master.log("pointsPerSample=" + pointsPerSample);
        Master.log("samples2save=" + samples2save);
        Master.log("dx=" + dx);
        Master.log("dt=" + dt);
        Master.log("xdim=" + xdim);
        Master.log("ydim=" + ydim);
        Master.log("format=" + format);
        Master.log("dataOffsetBytes=" + bytes);
    }

}
