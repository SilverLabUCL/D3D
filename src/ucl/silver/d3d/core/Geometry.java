package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;
import java.awt.Color;

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
public class Geometry extends CoordinatesVoxels {

    public String file = null;

    public boolean simpleCuboid = true; // shape is a simple cuboid of all space voxels
    public boolean forceEvenVoxels = false;
    public boolean forceOddVoxels = false;

    public double cubeWidth = 0; // um

    public ColorD3D colorSpace = new ColorD3D("space", Color.lightGray);
    //public ColorD3D colorNonSpace = new ColorD3D("non-space", Color.WHITE);
    public ColorD3D colorNonSpace = new ColorD3D("non-space", Color.gray);

    public int[][][] space = null; // the space array (-1) nonspace (>= 0) space
    
    static public int SPACEVALUE = 1;
    static public int NONSPACEVALUE = -1;
    public int NEWSPACEVALUE = 1; // default space value when resizing
    
    public double surfaceArea; // um^2

    public Voxel[][][] voxelSpace = null; // space array as objects (used in Monte Carlo simulations)
    public VoxelPBC[][][] voxelSpacePBC = null; // space array as objects (used in Monte Carlo simulations)

    public PSF[] PSFs = null; // array of PSFs
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("cubeWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("surfaceArea")) {
            return project.spaceUnits + "^2";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {

        if (name.equalsIgnoreCase("xVoxels")) {
            return true;
        }
        if (name.equalsIgnoreCase("yVoxels")) {
            return true;
        }
        if (name.equalsIgnoreCase("zVoxels")) {
            return true;
        }
        if (name.equalsIgnoreCase("xWidth")) {
            return true;
        }
        if (name.equalsIgnoreCase("yWidth")) {
            return true;
        }
        if (name.equalsIgnoreCase("zWidth")) {
            return true;
        }
        if (name.equalsIgnoreCase("cubeWidth")) {
            return true;
        }
        if (name.equalsIgnoreCase("displayPSF")) {
            return true;
        }
        if (name.equalsIgnoreCase("forceEvenVoxels")) {
            return true;
        }
        if (name.equalsIgnoreCase("forceOddVoxels")) {
            return true;
        }
        return false;
        
    }

    public Geometry(Project p, double Xwidth, double Ywidth, double Zwidth) {

        super(p);

        name = "Cuboid";

        xVoxels = (int) (Xwidth / project.dx);
        yVoxels = (int) (Ywidth / project.dx);
        zVoxels = (int) (Zwidth / project.dx);

        xVoxels = Math.max(1, xVoxels);
        yVoxels = Math.max(1, yVoxels);
        zVoxels = Math.max(1, zVoxels);

        setVoxels(0, 0, 0, xVoxels-1, yVoxels-1, zVoxels-1);

        space = new int[xVoxels][yVoxels][zVoxels];

        createVector(true);

        setSpace(1); // empty space

    }

    public Geometry(Project p, int Xvoxels, int Yvoxels, int Zvoxels) {

        super(p);

        name = "Cuboid";

        xVoxels = Math.max(1, Xvoxels);
        yVoxels = Math.max(1, Yvoxels);
        zVoxels = Math.max(1, Zvoxels);

        setVoxels(0, 0, 0, xVoxels-1, yVoxels-1, zVoxels-1);

        space = new int[xVoxels][yVoxels][zVoxels];

        createVector(true);

        setSpace(1); // empty space

    }

    public void init() {

        if (forceEvenVoxels) {
            forceOddVoxels = false;
        }

        checkSpace();

        if (PSFs != null) {
            for (PSF p : PSFs) {
                p.init();
                p.useGeometryCoordinates = true; // force this option
            }
        }

    }

    @Override
    public void update() {

        if (project == null) {
            return;
        }

        double dx = project.dx;

        xVoxel1 = 0;
        yVoxel1 = 0;
        zVoxel1 = 0;
        xVoxel2 = xVoxels - 1;
        yVoxel2 = yVoxels - 1;
        zVoxel2 = zVoxels - 1;

        xVoxelCenter = xVoxel1 + (xVoxels - 1) / 2.0;
        yVoxelCenter = yVoxel1 + (yVoxels - 1) / 2.0;
        zVoxelCenter = zVoxel1 + (zVoxels - 1) / 2.0;

        voxels = xVoxels * yVoxels * zVoxels;
        volume = voxels * dx * dx * dx;

        spaceVoxels = countSpaceVoxels(this);
        spaceVolume = spaceVoxels * dx * dx * dx;

        x1 = computeXedge(xVoxel1);
        y1 = computeYedge(yVoxel1);
        z1 = computeZedge(zVoxel1);

        x2 = computeXedge(xVoxel2);
        y2 = computeYedge(yVoxel2);
        z2 = computeZedge(zVoxel2);

        xWidth = x2 - x1;
        yWidth = y2 - y1;
        zWidth = z2 - z1;

        xCenter = x1 + xWidth / 2.0;
        yCenter = y1 + yWidth / 2.0;
        zCenter = z1 + zWidth / 2.0;

        if ((xWidth == yWidth) && (xWidth == zWidth)) {
            cubeWidth = xWidth;
        } else {
            cubeWidth = 0;
        }
        
        surfaceArea = GeometryTools.surfaceArea(space, project.dx);

        updateVectors();

    }

    public final double computeVoxelX(double x) {
        return (x / project.dx) + (xVoxels / 2.0);
    }

    public final double computeVoxelY(double y) {
        return (y / project.dx) + (yVoxels / 2.0);
    }

    public final double computeVoxelZ(double z) {
        return (z / project.dx) + (zVoxels / 2.0);
    }

    public final double computeX(double xVoxel) {
        return (xVoxel - xVoxelCenter) * project.dx; // um (voxel center)
    }

    public final double computeY(double yVoxel) {
        return (yVoxel - yVoxelCenter) * project.dx; // um (voxel center)
    }

    public final double computeZ(double zVoxel) {
        return (zVoxel - zVoxelCenter) * project.dx; // um (voxel center)
    }

    public final double computeXedge(int xVoxel) {

        double dx = project.dx;
        double x = (xVoxel - xVoxelCenter) * dx;

        if (x < 0) {
            x -= dx / 2.0;
        } else {
            x += dx / 2.0;
        }

        return x; // um (voxel center)

    }

    public final double computeYedge(int yVoxel) {

        double dx = project.dx;
        double y = (yVoxel - yVoxelCenter) * dx;

        if (y < 0) {
            y -= dx / 2.0;
        } else {
            y += dx / 2.0;
        }

        return y; // um (voxel center)

    }

    public final double computeZedge(int zVoxel) {

        double dx = project.dx;
        double z = (zVoxel - zVoxelCenter) * dx;

        if (z < 0) {
            z -= dx / 2.0;
        } else {
            z += dx / 2.0;
        }

        return z; // um (voxel center)

    }

    public double displayValue(int xVoxel, int yVoxel, int zVoxel) {

        if (!checkBounds(xVoxel, yVoxel, zVoxel)) {
            return 0; // non-space
        }

        if (space[xVoxel][yVoxel][zVoxel] < 0) {
            return 0; // non-space
        }

        //return getSpacePSFweight(xVoxel, yVoxel, zVoxel);

        //return space[xVoxel][yVoxel][zVoxel]; // space

        return 1;

    }

    public final int getSpace(int xVoxel, int yVoxel, int zVoxel) {
        if (checkBounds(xVoxel, yVoxel, zVoxel)) {
            return space[xVoxel][yVoxel][zVoxel];
        } else {
            return -1; // non-space
        }
    }

    public final double getSpacePSFweight(int xVoxel, int yVoxel, int zVoxel) {

        if (voxelSpace == null) {
            return -1;
        }

        if (checkBounds(xVoxel, yVoxel, zVoxel)) {
            return voxelSpace[xVoxel][yVoxel][zVoxel].PSFweight;
        } else {
            return -1; // non-space
        }

    }

    public final int[][][] getSpaceCopy() {

        int[][][] spaceCopy = new int[xVoxels][yVoxels][zVoxels];

        for (int k = zVoxel1; k <= zVoxel2; k++) {
            for (int j = yVoxel1; j <= yVoxel2; j++) {
                for (int i = xVoxel1; i <= xVoxel2; i++) {
                    spaceCopy[i][j][k] = space[i][j][k];
                }
            }
        }

        return spaceCopy;

    }

    public final boolean isSpace(CoordinatesVoxels c) {
        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                    if (!isSpace(i, j, k)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    public final boolean isSpace(double xVoxel, double yVoxel, double zVoxel) {
        return isSpace((int) xVoxel, (int) yVoxel, (int) zVoxel);
    }

    public final boolean isSpace(int xVoxel, int yVoxel, int zVoxel) {
        return ((checkBounds(xVoxel, yVoxel, zVoxel)) && (space[xVoxel][yVoxel][zVoxel] >= 0));
    }

    public final boolean checkBounds(double xVoxel, double yVoxel, double zVoxel) {
        return checkBounds((int) xVoxel, (int) yVoxel, (int) zVoxel);
    }

    public final boolean checkBounds(int xVoxel, int yVoxel, int zVoxel) {
        if (space == null) {
            return false;
        }
        return ((xVoxel >= 0) && (yVoxel >= 0) && (zVoxel >= 0) && (xVoxel < space.length) && (yVoxel < space[0].length) && (zVoxel < space[0][0].length));
    }

    public final int checkSpace() {

        spaceVoxels = 0;
        simpleCuboid = true;
        
        for (int[][] i : space) {
            for (int[] j : i ) {
                for (int s : j) {
                    if (s < 0) {
                        simpleCuboid = false;
                        s = NONSPACEVALUE;
                    } else {
                        spaceVoxels += 1;
                        s = SPACEVALUE;
                    }
                }
            }
        }

        update();
        
        //Master.log("spaceVoxels = " + spaceVoxels);

        return spaceVoxels;

    }

    public final int countSpaceVoxels(CoordinatesVoxels c) {

        int numSpaceVoxels = 0;

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                    if ((checkBounds(i, j, k)) && (space[i][j][k] >= 0)) {
                        numSpaceVoxels += 1;
                    }
                }
            }
        }

        return numSpaceVoxels;

    }

    public final boolean setSpace(int[][][] a) {

        if (a == null ) {
            return false;
        }

        for (int k = 0; k < space[0][0].length; k++) {
            for (int j = 0; j < space[0].length; j++) {
                for (int i = 0; i < space.length; i++) {

                    if (i >= a.length) {
                        continue;
                    }

                    if (j >= a[0].length) {
                        continue;
                    }

                    if (k >= a[0][0].length) {
                        continue;
                    }

                    space[i][j][k] = a[i][j][k];

                }
            }
        }

        checkSpace();

        return true;

    }

    public final boolean setSpace(double xVoxel, double yVoxel, double zVoxel, int v) {
        return setSpace((int) xVoxel, (int) yVoxel, (int) zVoxel, v);
    }

    public final boolean setSpace(int xVoxel, int yVoxel, int zVoxel, int v) {

        if (checkBounds(xVoxel, yVoxel, zVoxel)) {
            space[xVoxel][yVoxel][zVoxel] = v;
            return true;
        }
        return false;
    }

    public final void clear() {
        setSpace(SPACEVALUE);
    }

    public final void fill() {
        setSpace(NONSPACEVALUE);
    }

    public final void setSpace(int what) {
        CoordinatesVoxels c = new CoordinatesVoxels(project, 0, 0, 0, xVoxels - 1, yVoxels - 1, zVoxels - 1);
        setSpace(c, what);
    }

    public final int setSpace(CoordinatesVoxels c, int what) {
        
        int count = 0;

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {
                    if (checkBounds(i, j, k) && c.isInside(i, j, k)) {
                        space[i][j][k] = what;
                        count++;
                    }
                }
            }
        }

        checkSpace();
        
        return count;

    }

    public final void ellipsoid() {
        ellipsoid(1, 1, 1);
        name = "Ellipsoid";
        updateVectors();
    }

    public final void cylinder() {
        //ellipsoid(this, 0, 1, 1);
        //ellipsoid(this, 1, 0, 1);
        ellipsoid(1, 1, 0);
        name = "Cylinder";
        updateVectors();
    }

    public final void ellipsoid(int Xscale, int Yscale, int Zscale) { // create elliptical space
        fill();
        xScale = Xscale;
        yScale = Yscale;
        zScale = Zscale;
        GeometryTools.ellipsoid(space, this, 1, true);
        checkSpace();
    }

    public final void ellipsoid(CoordinatesVoxels c) { // create elliptical space
        fill();
        GeometryTools.ellipsoid(space, c, 1, true);
        checkSpace();
    }
    
    public final boolean resizeWithSpace(double Xwidth, double Ywidth, double Zwidth) {
        return resize(Xwidth, Ywidth, Zwidth, SPACEVALUE);
    }
    
    public final boolean resizeWithNonSpace(double Xwidth, double Ywidth, double Zwidth) {
        return resize(Xwidth, Ywidth, Zwidth, NONSPACEVALUE);
    }

    private boolean resize(double Xwidth, double Ywidth, double Zwidth, int newSpaceValue) {

        int i = -1, j = -1, k = -1;

        double dx = project.dx;

        if (Xwidth >= 0) {
            i = (int) (Xwidth / dx);
        }

        if (Ywidth >= 0) {
            j = (int) (Ywidth / dx);
        }

        if (Zwidth >= 0) {
            k = (int) (Zwidth / dx);
        }

        return resize(i, j, k, newSpaceValue);

    }
    
    public final void resizeAndFill(int Xvoxels, int Yvoxels, int Zvoxels, int fillValue) {
        resize(Xvoxels, Yvoxels, Zvoxels, NONSPACEVALUE);
        setSpace(fillValue);
    }
    
    public final boolean resizeWithSpace(int Xvoxels, int Yvoxels, int Zvoxels) {
        return resize(Xvoxels, Yvoxels, Zvoxels, SPACEVALUE);
    }
    
    public final boolean resizeWithNonSpace(int Xvoxels, int Yvoxels, int Zvoxels) {
        return resize(Xvoxels, Yvoxels, Zvoxels, NONSPACEVALUE);
    }

    private boolean resize(int Xvoxels, int Yvoxels, int Zvoxels, int newSpaceValue) {

        if (Xvoxels < 1) {
            Xvoxels = xVoxels;
        }

        if (Yvoxels < 1) {
            Yvoxels = yVoxels;
        }

        if (Zvoxels < 1) {
            Zvoxels = zVoxels;
        }

        if (forceEvenVoxels) {
            Xvoxels = Utility.makeEven(Xvoxels);
            Yvoxels = Utility.makeEven(Yvoxels);
            Zvoxels = Utility.makeEven(Zvoxels);
        } else if (forceOddVoxels) {
            Xvoxels = Utility.makeOdd(Xvoxels);
            Yvoxels = Utility.makeOdd(Yvoxels);
            Zvoxels = Utility.makeOdd(Zvoxels);
        }

        int[][][] temp = new int[Xvoxels][Yvoxels][Zvoxels];

        for (int k = 0; k < Zvoxels; k++) {
            for (int j = 0; j < Yvoxels; j++) {
                for (int i = 0; i < Xvoxels; i++) {

                    if ((space != null) && (i < space.length) &&
                            (j < space[0].length) && (k < space[0][0].length)) {
                        temp[i][j][k] = space[i][j][k]; // old space
                    } else {
                        temp[i][j][k] = newSpaceValue; // new space
                    }

                }
            }
        }

        setVoxels(0, 0, 0, Xvoxels-1, Yvoxels-1, Zvoxels-1);

        space = temp;

        checkSpace();

        if (voxelSpace != null) {
            voxelSpace = null; // force re-creation
        }

        resizePSFs();

        return true;

    }

    //
    // PSF functions
    //
    public void checkPSFs() {

        if (PSFs == null) {
            return;
        }

        for (PSF p : PSFs) {
            p.checkExists();
        }

    }

    public void resizePSFs() {

        if (PSFs == null) {
            return;
        }

        for (PSF p : PSFs) {

            p.useGeometryCoordinates = true;

            if (p.array != null) {
                p.array = null;
                p.checkExists();
            }

        }

    }

    public int addPSF(PSF p) {

        int i, j;

        if (PSFs == null) {
            i = 1;
        } else {
            i = PSFs.length + 1;
        }

        PSF[] pnew = new PSF[i]; // new array

        if (i > 1) {
            for (j = 0; j < i - 1; j++) {
                pnew[j] = PSFs[j]; // copy existing array to new array
            }
        }

        pnew[i - 1] = p;

        PSFs = pnew; // replace old array with new one

        i = pnew.length - 1;

        return i;

    }

    public boolean checkPSFArrayNum(int arrayNum) {
        return ((PSFs != null) && (arrayNum >= 0) && (arrayNum < PSFs.length));
    }

    public boolean killPSF(int arrayNum) {

        int k = 0;

        if (!checkPSFArrayNum(arrayNum)) {
            return false;
        }

        if (PSFs.length == 1) {
            PSFs = null;
            return true;
        }

        PSF[] p = new PSF[PSFs.length - 1]; // new array

        for (int j = 0; j < PSFs.length; j++) {
            if (j == arrayNum) {
                continue;
            }
            p[k] = PSFs[j];
            k++;
        }

        PSFs = p; // replace old array with new one

        return true;

    }

    //
    // Voxel functions
    //
    public Voxel[][][] getVoxelSpace() {
        return voxelSpace;
    }

    public boolean checkVoxelBounds(boolean PBC, double xVoxel, double yVoxel, double zVoxel) {
        return checkVoxelBounds(PBC, (int) xVoxel, (int) yVoxel, (int) zVoxel);
    }

    public boolean checkVoxelBounds(boolean PBC, int xVoxel, int yVoxel, int zVoxel) {

        if (PBC) {

            if (voxelSpacePBC == null) {
                return false;
            }
            
            return ((xVoxel >= 0) && (yVoxel >= 0) && (zVoxel >= 0) && (xVoxel < voxelSpacePBC.length) && (yVoxel < voxelSpacePBC[0].length) &&
                    (zVoxel < voxelSpacePBC[0][0].length));

        } else {

            if (voxelSpace == null) {
                return false;
            }
            
            return ((xVoxel >= 0) && (yVoxel >= 0) && (zVoxel >= 0) && (xVoxel < voxelSpace.length) && (yVoxel < voxelSpace[0].length) &&
                    (zVoxel < voxelSpace[0][0].length));


        }

    }

    public boolean initVoxels(int iPSF, int dPSF) {

        int index;
        int iscale, jscale, kscale;
        int itemp, jtemp, ktemp;
        boolean foundNonSpace;

        boolean iPSFOK = checkPSFArrayNum(iPSF);
        boolean dPSFOK = checkPSFArrayNum(dPSF);

        Master.log("initializing voxel arrays ( " + xVoxels + " x " + yVoxels +" x " + zVoxels + " )...");

        voxelSpace = new Voxel[xVoxels][yVoxels][zVoxels];

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    voxelSpace[i][j][k] = new Voxel();
                    voxelSpace[i][j][k].isSpace = isSpace(i, j, k);
                    voxelSpace[i][j][k].i = i;
                    voxelSpace[i][j][k].j = j;
                    voxelSpace[i][j][k].k = k;
                    voxelSpace[i][j][k].x = computeX(i);
                    voxelSpace[i][j][k].y = computeY(j);
                    voxelSpace[i][j][k].z = computeZ(k);

                    if (iPSFOK && (PSFs[iPSF] != null)) {
                        voxelSpace[i][j][k].PSFi = PSFs[iPSF].getValue(i, j, k);
                    }

                    if (dPSFOK && (PSFs[dPSF] != null)) {
                        voxelSpace[i][j][k].PSFd = PSFs[dPSF].getValue(i, j, k);
                    }

                }
            }
        }

        // compute array of neighbors

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    if (!voxelSpace[i][j][k].isSpace) {
                        continue;
                    }

                    voxelSpace[i][j][k].neighbors[0] = voxelSpace[i][j][k]; // itself
                    voxelSpace[i][j][k].numNeighbors = 1;

                    voxelSpace[i][j][k].numNonSpaceNeighbors = 0;

                    for (int ii = i - 1; ii <= i + 1; ii++) {
                        for (int jj = j - 1; jj <= j + 1; jj++) {
                            for (int kk = k - 1; kk <= k + 1; kk++) {

                                if ((ii < 0) || (ii >= xVoxels)) {
                                    continue; // out of bounds
                                }
                                if ((jj < 0) || (jj >= yVoxels)) {
                                    continue; // out of bounds
                                }
                                if ((kk < 0) || (kk >= zVoxels)) {
                                    continue; // out of bounds
                                }

                                if ((ii == i) && (jj == j) && (kk == k)) {
                                    continue; // current voxel cannot be a neighbor
                                }

                                if (!voxelSpace[ii][jj][kk].isSpace) {
                                    index = voxelSpace[i][j][k].numNonSpaceNeighbors;
                                    voxelSpace[i][j][k].nonSpaceNeighbors[index] = voxelSpace[ii][jj][kk];
                                    voxelSpace[i][j][k].numNonSpaceNeighbors++;
                                    continue;
                                }

                                index = voxelSpace[i][j][k].numNeighbors;

                                voxelSpace[i][j][k].neighbors[index] = voxelSpace[ii][jj][kk];
                                voxelSpace[i][j][k].numNeighbors++; // increment counter

                            }
                        }
                    }

                }
            }
        }

        if (!iPSFOK && !dPSFOK) {

            Master.log("finished initializing voxel arrays");

            return false;
            
        }

        // compute PSF weight

        //double minWeight = Double.POSITIVE_INFINITY, maxWeight = Double.NEGATIVE_INFINITY;

        //int counter;

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    if (!voxelSpace[i][j][k].isSpace) {
                        continue;
                    }

                    voxelSpace[i][j][k].PSFweight = 0;

                    for (int ii = 0; ii < 2; ii++) {
                        for (int jj = 0; jj < 2; jj++) {
                            for (int kk = 0; kk < 2; kk++) {

                                if (ii == 0) {
                                    iscale = -1;
                                } else {
                                    iscale = 1;
                                }

                                if (jj == 0) {
                                    jscale = -1;
                                } else {
                                    jscale = 1;
                                }

                                if (kk == 0) {
                                    kscale = -1;
                                } else {
                                    kscale = 1;
                                }

                                foundNonSpace = false;

                                //counter = 0;

                                //if ((i==10)&&(j==10)&&(k==10)) {
                                //Master.log("i=" + iscale + ", j=" + jscale + ", k=" + kscale);
                                //}
                                
                                for (int iii = 0; iii < 2; iii++) {
                                    for (int jjj = 0; jjj < 2; jjj++) {
                                        for (int kkk = 0; kkk < 2; kkk++) {

                                            if (iii + jjj + kkk == 0) {
                                                continue;
                                            }

                                            itemp = i + iii * iscale;
                                            jtemp = j + jjj * jscale;
                                            ktemp = k + kkk * kscale;

                                            if ((itemp < 0) || (itemp >= xVoxels)) {
                                                foundNonSpace = true;
                                                continue;
                                            }

                                            if ((jtemp < 0) || (jtemp >= yVoxels)) {
                                                foundNonSpace = true;
                                                continue;
                                            }

                                            if ((ktemp < 0) || (ktemp >= zVoxels)) {
                                                foundNonSpace = true;
                                                continue;
                                            }

                                            //if ((i==10)&&(j==10)&&(k==10)) {
                                            //if ((ii==1)&&(jj==1)&&(kk==1)) {
                                            //    Master.log("i=" + itemp + ", j=" + jtemp + ", k=" + ktemp);
                                            //}
                                            //}

                                            if (!voxelSpace[itemp][jtemp][ktemp].isSpace) {
                                                foundNonSpace = true;
                                            }

                                        }
                                    }
                                }

                                //Master.log("counter = " + counter);

                                if (!foundNonSpace) {
                                    voxelSpace[i][j][k].PSFweight += 0.125;
                                }

                            }
                        }
                    }

                    //minWeight = Math.min(minWeight, voxelSpace[i][j][k].PSFweight);
                    //maxWeight = Math.max(maxWeight, voxelSpace[i][j][k].PSFweight);

                    //Master.log("weight = " + voxelSpace[i][j][k].PSFweight);

                }

            }

        }

        //Master.log("minWeight = " + minWeight);
        //Master.log("maxWeight = " + maxWeight);

        //Master.log("finished initializing voxel arrays");

        return false;

    }

    public boolean initVoxelsPBC(int iPSF, int dPSF) {

        // PBC = Periodic Boundary Conditions for Monte Carlo simulations

        int index;
        int itemp, jtemp, ktemp;
        boolean foundPBC;

        boolean iPSFOK = checkPSFArrayNum(iPSF);
        boolean dPSFOK = checkPSFArrayNum(dPSF);

        voxelSpacePBC = new VoxelPBC[xVoxels][yVoxels][zVoxels];

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    voxelSpacePBC[i][j][k] = new VoxelPBC();
                    voxelSpacePBC[i][j][k].isSpace = isSpace(i, j, k);
                    voxelSpacePBC[i][j][k].i = i;
                    voxelSpacePBC[i][j][k].j = j;
                    voxelSpacePBC[i][j][k].k = k;
                    voxelSpacePBC[i][j][k].x = computeX(i);
                    voxelSpacePBC[i][j][k].y = computeY(j);
                    voxelSpacePBC[i][j][k].z = computeZ(k);

                    if (iPSFOK && (PSFs[iPSF] != null)) {
                        voxelSpacePBC[i][j][k].PSFi = PSFs[iPSF].getValue(i, j, k);
                    }

                    if (dPSFOK && (PSFs[dPSF] != null)) {
                        voxelSpacePBC[i][j][k].PSFd = PSFs[dPSF].getValue(i, j, k);
                    }

                }
            }
        }

        // compute array of neighbors

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    if (!voxelSpacePBC[i][j][k].isSpace) {
                        continue;
                    }

                    voxelSpacePBC[i][j][k].neighbors[0] = voxelSpacePBC[i][j][k]; // itself
                    voxelSpacePBC[i][j][k].numNeighbors = 1;

                    voxelSpacePBC[i][j][k].numNonSpaceNeighbors = 0;

                    voxelSpacePBC[i][j][k].numPBCneighbors = 0;

                    for (int ii = i - 1; ii <= i + 1; ii++) {
                        for (int jj = j - 1; jj <= j + 1; jj++) {
                            for (int kk = k - 1; kk <= k + 1; kk++) {

                                foundPBC = false;

                                if ((ii < 0) || (ii >= xVoxels)) {
                                    foundPBC = true;
                                }
                                if ((jj < 0) || (jj >= yVoxels)) {
                                    foundPBC = true;
                                }
                                if ((kk < 0) || (kk >= zVoxels)) {
                                    foundPBC = true;
                                }

                                itemp = voxelPBC(ii, xVoxels);
                                jtemp = voxelPBC(jj, yVoxels);
                                ktemp = voxelPBC(kk, zVoxels);

                                if ((itemp == i) && (jtemp == j) && (ktemp == k)) {
                                    continue; // current voxel cannot be a neighbor
                                }

                                if (!voxelSpacePBC[itemp][jtemp][ktemp].isSpace) {
                                    index = voxelSpacePBC[i][j][k].numNonSpaceNeighbors;
                                    voxelSpacePBC[i][j][k].nonSpaceNeighbors[index] = voxelSpacePBC[itemp][jtemp][ktemp];
                                    voxelSpacePBC[i][j][k].numNonSpaceNeighbors++;
                                    continue;
                                }

                                if (foundPBC) {
                                    index = voxelSpacePBC[i][j][k].numPBCneighbors;
                                    voxelSpacePBC[i][j][k].PBCneighbors[index] = voxelSpacePBC[itemp][jtemp][ktemp];
                                    voxelSpacePBC[i][j][k].PBCi[index] = indexPBC(ii, xVoxels); // save ijk positions
                                    voxelSpacePBC[i][j][k].PBCj[index] = indexPBC(jj, yVoxels); // save ijk positions
                                    voxelSpacePBC[i][j][k].PBCk[index] = indexPBC(kk, zVoxels); // save ijk positions
                                    voxelSpacePBC[i][j][k].numPBCneighbors++; // increment counter
                                } else {
                                    index = voxelSpacePBC[i][j][k].numNeighbors;
                                    voxelSpacePBC[i][j][k].neighbors[index] = voxelSpacePBC[itemp][jtemp][ktemp];
                                    voxelSpacePBC[i][j][k].numNeighbors++; // increment counter
                                }

                            }
                        }
                    }

                }
            }
        }

        return false;

    }

    private int voxelPBC(int voxelNum, int limit) {

        
        if (voxelNum < 0) {
            return limit - 1;
        }

        if (voxelNum >= limit) {
            return 0;
        }

        return voxelNum;

    }

    private int indexPBC(int voxelNum, int limit) {

        
        if (voxelNum < 0) {
            return -1;
        }

        if (voxelNum >= limit) {
            return 1;
        }

        return 0;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (colorSpace != null) {
            colorSpace.addUser(pv);
        }

        if (colorNonSpace != null) {
            colorNonSpace.addUser(pv);
        }

        if (PSFs != null) {
            for (PSF p : PSFs) {
                if (p != null) {
                    p.addUser(pv);
                }
            }
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (colorSpace != null) {
            addBlankParam();
            addVector(colorSpace.getVector());
        }

        if (colorNonSpace != null) {
            addBlankParam();
            addVector(colorNonSpace.getVector());
        }

        if (PSFs != null) {
            for (PSF p : PSFs) {
                if (p != null) {
                    addBlankParam();
                    p.createVector(true);
                    addVector(p.getVector());
                    p.addUser(this);
                }
            }
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (colorSpace != null) {
            colorSpace.updateVector(v);
        }

        if (colorNonSpace != null) {
            colorNonSpace.updateVector(v);
        }

        if (PSFs != null) {
            for (PSF p : PSFs) {
                if (p != null) {
                    p.updateVector(v);
                }
            }
        }

    }

    public boolean exportAB(String fileName) {
        return GeometryTools.exportGeometry(fileName, space, xVoxels, yVoxels, zVoxels, true);
    }

    public boolean importAB(String fileName, boolean fileContainsNonSpaceBorder) {
        return GeometryTools.importGeometry(this, fileName, fileContainsNonSpaceBorder);
    }

    // do not use this function directly, use ParamVector.set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof Geometry)) {
            return false;
        }

        if (n.equalsIgnoreCase("xVoxels")) {
            if (v < 0) {
                return false;
            }
            if ((v > 0) && (v == xVoxels)) {
                return true;
            }
            resize((int) v, -1, -1, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("yVoxels")) {
            if (v < 0) {
                return false;
            }
            if ((v > 0) && (v == yVoxels)) {
                return true;
            }
            resize(-1, (int) v, -1, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("zVoxels")) {
            if (v < 0) {
                return false;
            }
            if ((v > 0) && (v == zVoxels)) {
                return true;
            }
            resize(-1, -1, (int) v, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("xWidth")) {
            if (v < 0) {
                return false;
            }
            if ((v > 0) && (v == xWidth)) {
                return true;
            }
            resize(v, -1, -1, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("yWidth")) {
            if (v < 0) {
                return false;
            }
            if ((v > 0) && (v == yWidth)) {
                return true;
            }
            resize(-1, v, -1, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("zWidth")) {
            if (v < 0) {
                return false;
            }
            if ((v > 0) && (v == zWidth)) {
                return true;
            }
            resize(-1, -1, v, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("cubeWidth")) {
            if (v < 0) {
                return false;
            }
            if ((v > 0) && (v == cubeWidth)) {
                return true;
            }
            resize(v, v, v, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("forceEvenVoxels")) {
            if (v == 1) {
                forceEvenVoxels = true;
                forceOddVoxels = false;
            } else {
                forceEvenVoxels = false;
            }
            resize(-1, -1, -1, NEWSPACEVALUE);
            return true;
        }
        if (n.equalsIgnoreCase("forceOddVoxels")) {
            if (v == 1) {
                forceOddVoxels = true;
                forceEvenVoxels = false;
            } else {
                forceOddVoxels = false;
            }
            resize(-1, -1, -1, NEWSPACEVALUE);
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

        if (!(o.paramVector instanceof Geometry)) {
            return false;
        }

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        return super.setMyParams(o, s);
    }
}
