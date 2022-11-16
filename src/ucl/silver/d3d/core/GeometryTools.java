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
public final class GeometryTools {

    public static CoordinatesVoxels coordinates[];
    public static boolean useCoordinates = false;

    // cannot instantiate
    private GeometryTools() {
    }

    public static int spaceVoxels(int[][][] array) {

        int spaceVoxels = 0;
        
        for (int[][] i : array) {
            for (int[] j :  i) {
                for (int a : j) {
                    if (a < 0) {
                        //
                    } else {
                        spaceVoxels += 1;
                    }
                    
                }
            }
        }

        return spaceVoxels;

    }
    
    public static double surfaceArea(int[][][] a, double dx) {

        int faces = 0;
        double surfaceArea;

        for (int k = 0; k < a[0][0].length; k++) {
            for (int j = 0; j < a[0].length; j++) {
                for (int i = 0; i < a.length; i++) {

                    if (a[i][j][k] < 0) {
                        // non-space, do nothing
                    } else { // space

                        if ((i == 0) || (i == a.length - 1)) {
                            faces++;
                        } else {
                            if (a[i - 1][j][k] < 0) {
                                faces++;
                            }
                            if (a[i + 1][j][k] < 0) {
                                faces++;
                            }
                        }

                        if ((j == 0) || (j == a[0].length - 1)) {
                            faces++;
                        } else {
                            if (a[i][j - 1][k] < 0) {
                                faces++;
                            }
                            if (a[i][j + 1][k] < 0) {
                                faces++;
                            }
                        }

                        if ((k == 0) || (k == a[0][0].length - 1)) {
                            faces++;
                        } else {
                            if (a[i][j][k - 1] < 0) {
                                faces++;
                            }
                            if (a[i][j][k - 1] < 0) {
                                faces++;
                            }
                        }

                    }

                }
            }
        }
        
        surfaceArea = faces * dx * dx;

        return surfaceArea;

    }

    // Scale values are for circular variation (0) no (1) yes
    public static int ellipsoid(double[][][] a, CoordinatesVoxels c, double value) { // create elliptical space (center a0, b0, c0)

        int count = 0;
        double xd, yd, zd;

        if ((a == null) || (c == null)) {
            return 0;
        }

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

                    if ((i < 0) || (i >= a.length)) {
                        continue;
                    }

                    if ((j < 0) || (j >= a[0].length)) {
                        continue;
                    }

                    if ((k < 0) || (k >= a[0][0].length)) {
                        continue;
                    }

                    xd = (i - c.xVoxelCenter) / (c.xVoxels / 2.0);
                    yd = (j - c.yVoxelCenter) / (c.yVoxels / 2.0);
                    zd = (k - c.zVoxelCenter) / (c.zVoxels / 2.0);

                    if ((xd * xd * c.xScale + yd * yd * c.yScale + zd * zd * c.zScale) <= 1) {
                        a[i][j][k] = value;
                        count++;
                    }

                }
            }
        }
        
        return count;

    }

    public static int ellipsoid(int[][][] a, CoordinatesVoxels c, int value, boolean inside) { // create elliptical space (center a0, b0, c0)

        int count = 0;
        double xd, yd, zd;

        if ((a == null) || (c == null)) {
            return 0;
        }

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

                    if ((i < 0) || (i >= a.length)) {
                        continue;
                    }

                    if ((j < 0) || (j >= a[0].length)) {
                        continue;
                    }

                    if ((k < 0) || (k >= a[0][0].length)) {
                        continue;
                    }

                    xd = (i - c.xVoxelCenter) / (c.xVoxels / 2.0);
                    yd = (j - c.yVoxelCenter) / (c.yVoxels / 2.0);
                    zd = (k - c.zVoxelCenter) / (c.zVoxels / 2.0);

                    if (inside) {
                        if ((xd * xd * c.xScale + yd * yd * c.yScale + zd * zd * c.zScale) <= 1) {
                            a[i][j][k] = value;
                            count++;
                        }
                    } else if ((xd * xd * c.xScale + yd * yd * c.yScale + zd * zd * c.zScale) > 1) {
                        a[i][j][k] = value;
                        count++;
                    }
                    

                }
            }
        }
        
        return count;

    }

    public static boolean addCuboids(Geometry geometry, CoordinatesVoxels c, int N, double xdim, double ydim, double zdim) {

        int i, j, k, success = 0, failure = 0, maxTrials = 200000;
        double rani, ranj, rank;

        double dx = geometry.project.dx;

        CoordinatesVoxels cuboid = new CoordinatesVoxels(geometry.project);

        int xv = (int) (xdim / dx);
        int yv = (int) (ydim / dx);
        int zv = (int) (zdim / dx);

        if (c == null) {
            return true;
        }

        while (success < N) {

            // generate random numbers and map to available range of points in space
            rani = c.xVoxel1 + Master.mt.nextDouble() * (c.xVoxel2 - c.xVoxel1);
            ranj = c.yVoxel1 + Master.mt.nextDouble() * (c.yVoxel2 - c.yVoxel1);
            rank = c.zVoxel1 + Master.mt.nextDouble() * (c.zVoxel2 - c.zVoxel1);

            // change to integers
            i = (int) rani;
            j = (int) ranj;
            k = (int) rank;

            cuboid.setVoxels(i, j, k, i + xv - 1, j + yv - 1, k + zv - 1);

            if ((success < N) && geometry.isSpace(cuboid)) {
                geometry.setSpace(cuboid, -1);
                success += 1;
            } else {
                failure += 1;
            }

            if ((success == N) || (failure > maxTrials)) {
                System.out.println("Finished placing cuboids. N = " + success);
                return false;
            }

        }

        Master.log("finished placing cuboids. N = " + success);
        
        return (success < N);

    }

    public static double addEllipsoids(Geometry geometry, CoordinatesVoxels c, double radius, double axial_ratio, double volumeFraction, boolean exact, int ijkSelect, boolean cylinder, int linkNum, int extraVoxelsXYZ) {

        int i, j, k, iradius, ihalfLength, idiameter, numEllipsoids = 0;
        int oldVoxels, newVoxels, ellipsoidVoxels, spaceVoxels, iLast;
        int maxEllipsoids = 1000000, maxTrials = 1000, trialCount;
        double d, rani, ranj, rank;
        boolean overlap;
        boolean randomIJK = false;
        boolean even;
        boolean link = false;
        int linkCounter = 0;

        double dx = geometry.project.dx;
        
        int[][][] spaceCopy = geometry.getSpaceCopy();

        int xv1 = c.xVoxel1 - extraVoxelsXYZ;
        int yv1 = c.yVoxel1 - extraVoxelsXYZ;
        int zv1 = c.zVoxel1 - extraVoxelsXYZ;
        int xv2 = c.xVoxel2 + extraVoxelsXYZ;
        int yv2 = c.xVoxel2 + extraVoxelsXYZ;
        int zv2 = c.xVoxel2 + extraVoxelsXYZ;

        int cx1, cy1, cz1, cx2, cy2, cz2;

        CoordinatesVoxels ctemp = new CoordinatesVoxels(geometry.project);

        if (useCoordinates) {
            maxEllipsoids = coordinates.length;
        } else {
            coordinates = new CoordinatesVoxels[maxEllipsoids];
        }

        iradius = (int) Math.round(radius / dx);
        idiameter = (int) Math.round(2 * radius / dx);
        ihalfLength = (int) Math.round(0.5 * axial_ratio * idiameter);
        
        //Master.log("" + iradius + "," + ihalfLength);

        if (Math.IEEEremainder(idiameter, 2) == 0) {
            even = true;
        } else {
            even = false;
            iradius -= 1;
            ihalfLength -= 1;
        }

        spaceVoxels = geometry.spaceVoxels;
        oldVoxels = geometry.spaceVoxels;
        newVoxels = (int) (oldVoxels * (1.0 - volumeFraction));

        if (ijkSelect < 0) {
            randomIJK = true;
        }

        if (linkNum > 1) {
            link = true;
        }

        for (int ii = 0; ii < maxEllipsoids; ii++) {

            overlap = false;

            if (useCoordinates) {

                if (coordinates[ii] == null) {
                    continue;
                }
                
            } else {

                coordinates[ii] = new CoordinatesVoxels(geometry.project);

                if (randomIJK) {

                    d = Master.mt.nextDouble();

                    if (d < 0.3333333333) {
                        ijkSelect = 0; // ellipsoid-x
                    } else if (d < 0.6666666666) {
                        ijkSelect = 1; // ellipsoid-y
                    } else {
                        ijkSelect = 2; // ellipsoid-z
                    }

                }

                trialCount = 0;

                for (int jj = 0; jj < maxTrials; jj++) {

                    if (link && (linkCounter > 0)) {

                        switch (ijkSelect) {
                            case 0: // ellipsoid-x
                                cx1 = coordinates[ii - linkCounter].xVoxel1;
                                cy1 = coordinates[ii - linkCounter].yVoxel1 - iradius;
                                cz1 = coordinates[ii - linkCounter].zVoxel1 - iradius;
                                cx2 = coordinates[ii - linkCounter].xVoxel2;
                                cy2 = coordinates[ii - linkCounter].yVoxel2 + iradius;
                                cz2 = coordinates[ii - linkCounter].zVoxel2 + iradius;
                                cy1 = Math.max(cy1, yv1 + iradius);
                                cz1 = Math.max(cz1, zv1 + iradius);
                                cy2 = Math.min(cy2, yv2 - iradius);
                                cz2 = Math.min(cz2, zv2 - iradius);
                                break;
                            case 1: // ellipsoid-y
                                cx1 = coordinates[ii - linkCounter].xVoxel1 - iradius;
                                cy1 = coordinates[ii - linkCounter].yVoxel1;
                                cz1 = coordinates[ii - linkCounter].zVoxel1 - iradius;
                                cx2 = coordinates[ii - linkCounter].xVoxel2 + iradius;
                                cy2 = coordinates[ii - linkCounter].yVoxel2;
                                cz2 = coordinates[ii - linkCounter].zVoxel2 + iradius;
                                cx1 = Math.max(cx1, xv1 + iradius);
                                cz1 = Math.max(cz1, zv1 + iradius);
                                cx2 = Math.min(cx2, xv2 - iradius);
                                cz2 = Math.min(cz2, zv2 - iradius);
                                break;
                            case 2: // ellipsoid-z
                                cx1 = coordinates[ii - linkCounter].xVoxel1 - iradius;
                                cy1 = coordinates[ii - linkCounter].yVoxel1 - iradius;
                                cz1 = coordinates[ii - linkCounter].zVoxel1;
                                cx2 = coordinates[ii - linkCounter].xVoxel2 + iradius;
                                cy2 = coordinates[ii - linkCounter].yVoxel2 + iradius;
                                cz2 = coordinates[ii - linkCounter].zVoxel2;
                                cx1 = Math.max(cx1, xv1 + iradius);
                                cy1 = Math.max(cy1, yv1 + iradius);
                                cx2 = Math.min(cx2, xv2 - iradius);
                                cy2 = Math.min(cy2, yv2 - iradius);
                                break;
                            default:
                                return 0;
                        }

                        rani = cx1 + Master.mt.nextDouble() * (cx2 - cx1);
                        ranj = cy1 + Master.mt.nextDouble() * (cy2 - cy1);
                        rank = cz1 + Master.mt.nextDouble() * (cz2 - cz1);

                    } else {
                        
                        switch (ijkSelect) {
                            case 0: // ellipsoid-x
                                cx1 = xv1;
                                cy1 = yv1 + iradius;
                                cz1 = zv1 + iradius;
                                cx2 = xv2;
                                cy2 = yv2 - iradius;
                                cz2 = zv2 - iradius;
                                break;
                            case 1: // ellipsoid-y
                                cx1 = xv1 + iradius;
                                cy1 = yv1;
                                cz1 = zv1 + iradius;
                                cx2 = xv2 - iradius;
                                cy2 = yv2;
                                cz2 = zv2 - iradius;
                                break;
                            case 2: // ellipsoid-z
                                cx1 = xv1 + iradius;
                                cy1 = yv1 + iradius;
                                cz1 = zv1;
                                cx2 = xv2 - iradius;
                                cy2 = yv2 - iradius;
                                cz2 = zv2;
                                break;
                            default:
                                return 0;
                        }

                        rani = cx1 + Master.mt.nextDouble() * (cx2 - cx1);
                        ranj = cy1 + Master.mt.nextDouble() * (cy2 - cy1);
                        rank = cz1 + Master.mt.nextDouble() * (cz2 - cz1);

                    }

                    i = (int) Math.round(rani);
                    j = (int) Math.round(ranj);
                    k = (int) Math.round(rank);

                    coordinates[ii].xScale = 1;
                    coordinates[ii].yScale = 1;
                    coordinates[ii].zScale = 1;

                    switch (ijkSelect) {

                        case 0: // ellipsoid-x

                            if (even) {
                                coordinates[ii].setVoxels(i - ihalfLength + 1, j - iradius + 1, k - iradius + 1, i + ihalfLength, j + iradius, k + iradius);
                            } else {
                                coordinates[ii].setVoxels(i - ihalfLength, j - iradius, k - iradius, i + ihalfLength, j + iradius, k + iradius);
                            }

                            if (cylinder) {
                                coordinates[ii].xScale = 0;
                            }

                            break;

                        case 1: // ellipsoid-y

                            if (even) {
                                coordinates[ii].setVoxels(i - iradius + 1, j - ihalfLength + 1, k - iradius + 1, i + iradius, j + ihalfLength, k + iradius);
                            } else {
                                coordinates[ii].setVoxels(i - iradius, j - ihalfLength, k - iradius, i + iradius, j + ihalfLength, k + iradius);
                            }

                            if (cylinder) {
                                coordinates[ii].yScale = 0;
                            }

                            break;

                        case 2: // ellipsoid-z

                            if (even) {
                                coordinates[ii].setVoxels(i - iradius + 1, j - iradius + 1, k - ihalfLength + 1, i + iradius, j + iradius, k + ihalfLength);
                            } else {
                                coordinates[ii].setVoxels(i - iradius, j - iradius, k - ihalfLength, i + iradius, j + iradius, k + ihalfLength);
                            }

                            if (cylinder) {
                                coordinates[ii].zScale = 0;
                            }

                            break;

                    }

                    overlap = false;

                    for (int kk = 0; kk < ii; kk++) {
                        if (coordinates[ii].isInside(coordinates[kk])) {
                            overlap = true;
                        }
                    }


                    trialCount++;

                    if (!overlap) {
                        break;
                    }

                }

                if (overlap) {
                    continue;
                }

            }

            if (overlap) {
                
                coordinates[ii] = null;
                
            } else {

                GeometryTools.ellipsoid(spaceCopy, coordinates[ii], -1, true);
                numEllipsoids++;
                
                //Master.log("added ellipsoid = " + numEllipsoids);

                if (link) {

                    linkCounter++;

                    if (linkCounter >= linkNum) {
                        linkCounter = 0;
                    }

                }
                
            }
            
            spaceVoxels = spaceVoxels(spaceCopy);

            if (spaceVoxels <= newVoxels) {
                break; // finished
            }

        }
        
        if ((exact && (spaceVoxels < newVoxels))) {

            iLast = -1;

            for (int ii = 0; ii < numEllipsoids; ii++) {

                if (coordinates[ii] == null) {
                    continue;
                }

                iLast = ii;

            }

            for (int jj = 0; jj < geometry.xVoxels; jj++) {

                if (iLast < 0) {
                    break;
                }

                spaceCopy = geometry.getSpaceCopy();

                ctemp.matchVoxels(coordinates[iLast]);

                if (ctemp.zScale == 0) {
                    ctemp.zVoxel1 = geometry.zVoxels - jj - 1;
                    ctemp.update();
                } else if (ctemp.xScale == 0) {
                    ctemp.xVoxel1 = geometry.xVoxels - jj - 1;
                    ctemp.update();
                } else if (ctemp.yScale == 0) {
                    ctemp.yVoxel1 = geometry.yVoxels - jj - 1;
                    ctemp.update();
                }

                GeometryTools.ellipsoid(spaceCopy, ctemp, -1, true);

                spaceVoxels = spaceVoxels(spaceCopy);

                if (spaceVoxels <= newVoxels) {
                    break; // finished
                }

            }

            geometry.setSpace(spaceCopy);
            
        } else {

            geometry.setSpace(spaceCopy);

        }

        ellipsoidVoxels = oldVoxels - geometry.spaceVoxels;

        d = ( 1.0 * ellipsoidVoxels ) / ( 1.0  * oldVoxels); // volume fraction

        Master.log("added ellipsoids = " + numEllipsoids + ", space voxels = " + spaceVoxels + ", density = " + d);

        return d;

    }

    // function to write shape matrix in char A's and B's
    // file begins with 3 integers for i, j, k dimensions
    public static boolean exportGeometry(String fileName, int[][][] array, int iMax,
            int jMax, int kMax, boolean writeDim) {

        DataOutputStream dout;

        try {

            dout = new DataOutputStream(new FileOutputStream(fileName));

            if (writeDim) {
                dout.writeBytes(Integer.toString(iMax));
                dout.writeByte(' ');
                dout.writeBytes(Integer.toString(jMax));
                dout.writeByte(' ');
                dout.writeBytes(Integer.toString(kMax));
                dout.writeByte(' ');
            }

            for (int k = 0; k < kMax; k++) {
                for (int j = 0; j < jMax; j++) {
                    for (int i = 0; i < iMax; i++) {
                        if (array[i][j][k] >= 0) {
                            dout.writeByte(66); // B, space
                        } else {
                            dout.writeByte(65); // A, non-space
                        }
                    }
                }
            }

            dout.close();

            //System.err.println("Exported shape to file: " + fileName);

        } catch (IOException e) {
            System.err.println("unable to write file: " + fileName);
            return false;
        }

        return true;

    }

    // function to read shape matrix of char A's and B's
    // pass -1 to read all i, j, k voxels
    public static boolean importGeometry(Geometry g, String fileName, boolean fileContainsNonSpaceBorder) {

        int i, j, k, iMax, jMax, kMax, space, count = 0;
        int ii, jj, kk;
        char c;

        int dim[] = determineGeometrySize(fileName);

        DataInputStream din;

        if ((dim == null) || (dim.length != 3)) {
            System.err.println("encountered bad matrix dimensions for " + fileName);
            return false;
        }

        iMax = dim[0];
        jMax = dim[1];
        kMax = dim[2];

        if ((iMax < 0) || (jMax < 0) || (kMax < 0)) {
            System.err.println("encountered bad matrix dimensions for " + fileName);
            return false;
        }

        if (fileContainsNonSpaceBorder) {
            g.resizeWithSpace(iMax - 2, jMax - 2, kMax - 2);
        } else {
            g.resizeWithSpace(iMax, jMax, kMax);
        }

        g.clear();

        try {

            din = new DataInputStream(new FileInputStream(fileName));

            do { // read header until last space character
                c = (char) din.readByte();
                if (c == ' ') {
                    count++;
                }
                if ((c == 'A') || (c == 'B')) {
                    return false; // error, read too far
                }
            } while (count < 3);

            for (k = 0; k < kMax; k++) {
                for (j = 0; j < jMax; j++) {
                    for (i = 0; i < iMax; i++) {

                        c = (char) din.readByte();

                        space = -1; // non-space

                        if (c == 'B') {
                            space = 0;
                        }

                        ii = i;
                        jj = j;
                        kk = k;

                        if (fileContainsNonSpaceBorder) {

                            if ((i == 0) || (j == 0) || (k == 0)) {
                                continue; // skip border of shape file
                            }

                            if ((i == iMax - 1) || (j == jMax - 1) || (k == kMax - 1)) {
                                continue; // skip border of shape file
                            }

                            ii = i - 1;
                            jj = j - 1;
                            kk = k - 1;

                        }

                        //if ((i == iMax - 1) || (j == jMax - 1) || (k == kMax - 1)) {
                        //    continue; // skip border of shape file
                        //}

                        //if ((i == 0) || (j == 0) || (k == 0)) {
                        //    space = 0; // border of space matrix
                        //}

                        //if ((i == g.xVoxels - 1) || (j == g.yVoxels - 1) ||
                        //        (k == g.zVoxels - 1)) {
                        //    space = 0; // border of space matrix
                        //}

                        g.setSpace(ii, jj, kk, space);

                    }

                }

            }

            din.close();

        } catch (EOFException e) {
            System.err.println("Unexpectedly encountered End Of File when reading "
                    + fileName);
            return false;
        } catch (IOException e) {
            System.err.println("unable to read file: " + fileName);
            return false;
        }

        return true;

    }

    public static int[] determineGeometrySize(String fileName) {
        int count = 0;
        char c;
        boolean read = true;
        String snum = "";

        int dim[] = {-1, -1, -1};

        DataInputStream din;

        try {

            din = new DataInputStream(new FileInputStream(fileName));

            do {
                c = (char) din.readByte();
                switch (c) {
                    case 'A':
                    case 'B':
                        read = false;
                        break;
                    case ' ':
                        dim[count] = Integer.parseInt(snum);
                        count++;
                        snum = "";
                        break;
                    default:
                        snum = snum + c;
                }
            } while (read);

            din.close();

        } catch (EOFException e) {
            System.err.println("unexpectedly encountered End Of File when reading " +
                    fileName);
        } catch (IOException e) {
            System.err.println("unable to read file: " + fileName);
        }

        if ((count == 3) && (dim[0] >= 0) && (dim[1] >= 0) && (dim[2] >= 0)) {
            return dim;
        } else {
            return null;
        }

    }

    // function to write shape matrix as packed 1's and 0's
    // file begins with length of integer array
    public static void exportPacked(String fileName, byte[][][] array) {
        int iNum;

        int[] packed = pack3D(array); // pack as integer array

        iNum = packed.length;

        DataOutputStream dout;

        try {

            dout = new DataOutputStream(new FileOutputStream(fileName));

            dout.writeInt(iNum); // first integer is length of array

            for (int i = 0; i < iNum; i++) {
                dout.writeInt(packed[i]);
            }

            dout.close();

        } catch (IOException e) {
            System.err.println("unable to write file: " + fileName);
        }

    }

    // function to read shape matrix of packed 1's and 0's
    public static byte[][][] importPacked(String fileName) {

        int iNum;
        int[] packed;

        DataInputStream din;

        try {

            din = new DataInputStream(new FileInputStream(fileName));
            iNum = din.readInt();

            packed = new int[iNum];

            for (int i = 0; i < iNum; i++) {
                try {
                    packed[i] = din.readInt();
                } catch (EOFException e) {
                    System.err.println("Unexpectedly encountered End Of File when reading " + fileName);
                    break;
                }

            }

            din.close();

            return unpack3D(packed);

        } catch (EOFException e) {
            System.err.println("Unexpectedly encountered End Of File when reading " + fileName);
        } catch (IOException e) {
            System.err.println("unable to read file: " + fileName);
        }

        return null;

    }

    // function to pack 3D shape matrix into
    // 1D integer array. reduces bits by factor of 8.
    public static int[] pack3D(byte[][][] array3D) {

        int i, j, k, p, s, pack;

        int iNum = array3D.length;
        int jNum = array3D[0].length;
        int kNum = array3D[0][0].length;
        int pNum = 4 + (iNum * jNum * kNum / 32);

        int packed[] = new int[pNum];

        packed[0] = iNum;
        packed[1] = jNum;
        packed[2] = kNum;

        p = 3;
        s = 0;
        pack = 0;

        for (k = 0; k < kNum; k++) {
            for (j = 0; j < jNum; j++) {
                for (i = 0; i < iNum; i++) {

                    pack = array3D[i][j][k] | pack;
                    s += 1;

                    if (s < 32) {

                        pack = pack << 1;

                    } else {

                        packed[p] = pack;
                        p += 1;
                        s = 0;
                        pack = 0;

                        if (p == pNum) {
                            return packed;
                        }

                    }
                }
            }
        }

        packed[p] = pack; // make sure to save last array value

        return packed;

    }

    // function to unpack 1D integer array into 3D shape matrix
    public static byte[][][] unpack3D(int[] packed) {

        int i, j, k, p, s, ipack, unpack, pNum = packed.length;

        int iNum = packed[0];
        int jNum = packed[1];
        int kNum = packed[2];

        byte array3D[][][] = new byte[iNum][jNum][kNum];

        p = 3;
        s = 32 - 1;
        ipack = packed[p];

        for (k = 0; k < kNum; k++) {
            for (j = 0; j < jNum; j++) {
                for (i = 0; i < iNum; i++) {

                    unpack = (ipack >> s) & 1;
                    array3D[i][j][k] = (byte) unpack;
                    s -= 1;

                    if (s < 0) {

                        p += 1;
                        s = 32 - 1;

                        if (p == pNum) {
                            Master.log("alert: incomplete packed array.");
                            return array3D; // shouldnt happen
                        }

                        ipack = packed[p];

                    }
                }
            }
        }

        return array3D;

    }
}
