package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;
import java.util.Random;

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
public final class InitMC_GeometryWithMito_Jason extends ParamVector {
    
    public double xdim = 8;
    public double ydim = 8;
    public double zdim = 10;
    public boolean ellipsoid = false;
    public boolean eighth_geometry = false;
    public boolean nob = false;
    
    public double xvcenter, yvcenter, zvcenter;
    
    public boolean mitoNonSpaceVoxels = false;
    public double mitoRadius = (0.25 / 0.89) / 2.0; // um Palay
    public double mitoAxialRatio = 2.0 / 0.25; // Palay
    public double mitoVolumeFraction = 0.28; // Zoltan average
    public double mitoVolumeFractionTolerance = 0.02;
    public double mitoVolumeFractionActual = 0;
    public boolean mitoVolumeFractionExact = false;
    public boolean mitoPreset = false;
    public boolean mitoCylinder = true;
    public int mito_ijkselect = -1;
    
    public long seed = 8682522807148012L + System.nanoTime();
    
    Random ran = new Random(seed);

    public InitMC_GeometryWithMito_Jason(Project p) {
        super(p);
        createVector(true);
    }

    @Override
    public boolean canEdit(String name) {
        return false;
    }

    public boolean initGeometryMC(int iPSFselect) {

        boolean error = false;

        double xdim1 = xdim;
        double ydim1 = ydim;
        double zdim1 = zdim;

        Geometry geometry = project.geometry;

        geometry.forceEvenVoxels = true;

        if (eighth_geometry) {

            xdim1 = xdim / 2;
            ydim1 = ydim / 2;

            if (iPSFselect < 2) {
                zdim1 = zdim / 2;
            }

        }

        geometry.resizeWithSpace(xdim1, ydim1, zdim1);

        initGeometryMC2();

        if (mitoNonSpaceVoxels) {
            if (mitoPreset) {
                error = initGeometryMitochondriaPreset();
            } else {
                error = initGeometryMitochondria();
            }
        }

        xvcenter = geometry.xVoxelCenter;
        yvcenter = geometry.yVoxelCenter;
        zvcenter = geometry.zVoxelCenter;

        if (eighth_geometry) {

            xvcenter = geometry.xVoxels - 0.5;
            yvcenter = geometry.yVoxels - 0.5;

            if (iPSFselect < 2) {
                zvcenter = geometry.zVoxels - 0.5;
            }

        }

        geometry.update();

        return error;

    }

    public void initGeometryMC2() {

        Geometry geometry = project.geometry;

        geometry.clear();

        if (ellipsoid) {
            geometry.ellipsoid();
        }

        if (nob) {

            double zwidth = 0.4;
            double xywidth = 0.5;

            int iz = (int) (zwidth / project.dx) - 1;
            int xyv = (int) (xywidth / project.dx);
            int ixy = (int) ((geometry.xVoxels - xyv) / 2.0);

            CoordinatesVoxels c = new CoordinatesVoxels(project, geometry);

            c.setVoxels(geometry.xVoxel1, geometry.yVoxel1, geometry.zVoxel2 - iz, geometry.xVoxel2, geometry.yVoxel2, geometry.zVoxel2);

            geometry.setSpace(c, -1);

            c.setVoxels(geometry.xVoxel1 + ixy, geometry.yVoxel1 + ixy, geometry.zVoxel2 - iz, geometry.xVoxel2 - ixy, geometry.yVoxel2 - ixy, geometry.zVoxel2);

            geometry.setSpace(c, 1);

        }

    }

    public boolean initGeometryMitochondriaPreset() {

        double x0, y0, z0;

        GeometryTools.coordinates = new CoordinatesVoxels[4];
        GeometryTools.useCoordinates = true;

        //mito #0 -0.8, -6.85, -0.55, -0.2, 5.2, 0.1
        //xScale=1,yScale=0,zScale=1

        //mito #1 0.2, 0.1, -5.65, 0.8, 0.7, 6.4
        //xScale=1,yScale=1,zScale=0

        //mito #2 -1.05, -6.3, 0.6, -0.45, 5.75, 1.2
        //xScale=1,yScale=0,zScale=1

        //mito #3 -5.65, -0.95, -1.3, 6.4, -0.35, -0.7
        //xScale=0,yScale=1,zScale=1

        GeometryTools.coordinates[0] = new CoordinatesVoxels(project, -0.8, -6.85, -0.55, -0.2, 5.2, 0.1);
        GeometryTools.coordinates[0].xScale = 1;
        GeometryTools.coordinates[0].yScale = 0;
        GeometryTools.coordinates[0].zScale = 1;

        GeometryTools.coordinates[1] = new CoordinatesVoxels(project, 0.2, 0.1, -5.65, 0.8, 0.7, 6.4);
        GeometryTools.coordinates[1].xScale = 1;
        GeometryTools.coordinates[1].yScale = 1;
        GeometryTools.coordinates[1].zScale = 0;

        GeometryTools.coordinates[2] = new CoordinatesVoxels(project, -1.05, -6.3, 0.6, -0.45, 5.75, 1.2);
        GeometryTools.coordinates[2].xScale = 1;
        GeometryTools.coordinates[2].yScale = 0;
        GeometryTools.coordinates[2].zScale = 1;

        GeometryTools.coordinates[3] = new CoordinatesVoxels(project, -5.65, -0.95, -1.3, 6.4, -0.35, -0.7);
        GeometryTools.coordinates[3].xScale = 0;
        GeometryTools.coordinates[3].yScale = 1;
        GeometryTools.coordinates[3].zScale = 1;

        x0 = GeometryTools.coordinates[0].xCenter;
        y0 = GeometryTools.coordinates[0].yCenter;
        z0 = GeometryTools.coordinates[0].zCenter;

        x0 -= 0.15;
        z0 -= 0.15;

        GeometryTools.coordinates[0].setCenter(x0, y0, z0);

        initGeometryMitochondria2();

        for (int i = 0; i < GeometryTools.coordinates.length; i++) {

            if (GeometryTools.coordinates[i] == null) {
                break;
            }

            Master.log("mito #" + i + " " + GeometryTools.coordinates[i].dimensionsToString(true));
            Master.log(GeometryTools.coordinates[i].scalesToString(true));

        }

        return false;

    }

    public Coordinates[] getMitochondriaPreset0() {

        double x0, y0, z0;

        String cstr = "[r=200,g=200,b=200]";

        Coordinates[] mito = new Coordinates[4];

        mito[0] = new Coordinates(project, -0.8, -6.85, -0.55, -0.2, 5.2, 0.1);
        mito[0].setShape("cylindery");
        mito[0].color.setColor(cstr);

        x0 = mito[0].xCenter;
        y0 = mito[0].yCenter;
        z0 = mito[0].zCenter;

        x0 -= 0.15;
        z0 -= 0.10;

        mito[0].setCenter(x0, y0, z0);

        mito[1] = new Coordinates(project, 0.2, 0.1, -5.65, 0.8, 0.7, 6.4);
        mito[1].setShape("cylinderz");
        mito[1].color.setColor(cstr);

        mito[2] = new Coordinates(project, -1.05, -6.3, 0.6, -0.45, 5.75, 1.2);
        mito[2].setShape("cylindery");
        mito[2].color.setColor(cstr);

        mito[3] = new Coordinates(project, -5.65, -0.95, -1.3, 6.4, -0.35, -0.7);
        mito[3].setShape("cylinderx");
        mito[3].color.setColor(cstr);

        return mito;

    }

    public Coordinates[] getMitochondriaPreset1() {

        double xo, yo;

        String cstr = "[r=200,g=200,b=200]";

        Coordinates[] mito = new Coordinates[4];

        mito[0] = new Coordinates(project, -1.0, -0.8, -0.63, 1.0, -0.147, 0.017);
        mito[0].setShape("cylinderx");
        mito[0].color.setColor(cstr);

        mito[1] = new Coordinates(project, -1.0, -0.55464, 0.45621, 1.0, 0.09536, 1.10621);
        mito[1].setShape("cylinderx");
        mito[1].color.setColor(cstr);

        mito[2] = new Coordinates(project, -1.12836, 0.26097, -1.0, -0.47836, 0.91097, 1.0);
        mito[2].setShape("cylinderz");
        mito[2].color.setColor(cstr);

        xo = 0.15;
        yo = 0.1;

        mito[3] = new Coordinates(project, -0.25562 + xo, 0.24019 + yo, -1.0, 0.253 + xo, 0.74882 + yo, 1.0);
        mito[3].setShape("cylinderz");
        mito[3].color.setColor(cstr);

        return mito;

    }

    public Coordinates[] getMitochondriaPreset2() {

        double xo, zo;

        String cstr = "[r=200,g=200,b=200]";

        Coordinates[] mito = new Coordinates[9];

        mito[0] = new Coordinates(project, -1.27927, -0.07558, -1.5, -0.62927, 0.57442, 1.5);
        mito[0].setShape("cylinderz");
        mito[0].color.setColor(cstr);

        mito[1] = new Coordinates(project, -1.5, 0.59841, -1.11563, 1.5, 1.24841, -0.46563);
        mito[1].setShape("cylinderx");
        mito[1].color.setColor(cstr);

        mito[2] = new Coordinates(project, -1.5, -1.19391, -0.78495, 1.5, -0.54391, -0.13495);
        mito[2].setShape("cylinderx");
        mito[2].color.setColor(cstr);

        mito[3] = new Coordinates(project, 0.67073, -1.5, 1.02263, 1.32073, 1.5, 1.67263);
        mito[3].setShape("cylindery");
        mito[3].color.setColor(cstr);

        mito[4] = new Coordinates(project, 0.57183, -1.5, 0.23752, 1.22183, 1.5, 0.88752);
        mito[4].setShape("cylindery");
        mito[4].color.setColor(cstr);

        mito[5] = new Coordinates(project, -0.11824, -1.5, 1.00657, 0.53176, 1.5, 1.65657);
        mito[5].setShape("cylindery");
        mito[5].color.setColor(cstr);

        mito[6] = new Coordinates(project, 0.59867, -1.5, -1.77503, 1.24867, 1.5, -1.12503);
        mito[6].setShape("cylindery");
        mito[6].color.setColor(cstr);

        mito[7] = new Coordinates(project, -0.19043, -1.5, -1.80081, 0.45957, 1.5, -1.15081);
        mito[7].setShape("cylindery");
        mito[7].color.setColor(cstr);

        xo = -0.1;
        zo = 0.25;
        
        mito[8] = new Coordinates(project, -0.41156 + xo, -1.5, -0.02648f + zo, 0.18286 + xo, 1.5, 0.56794 + zo);
        mito[8].setShape("cylindery");
        mito[8].color.setColor(cstr);

        return mito;

    }

    public boolean initGeometryMitochondria() {

        int maxMitoTrials = 500, trialCount = 0;

        boolean error = false;

        if (mitoVolumeFraction <= 0 ) {
            return false;
        }

        Master.log("adding mito...");

        for (int i = 0; i < maxMitoTrials; i++) {

            //Master.log("mito placement trial #" + i);
            error = initGeometryMitochondria2();

            trialCount++;

            if (!error) {
                break;
            }

        }

        if (error) {
            Master.log("failed to add mitochondria, trials = " + trialCount);
            return true;
        }

        Master.log("mitochondria volume fraction = " + mitoVolumeFractionActual + " percent, trials = " + trialCount);

        //for (int i = 0; i < GeometryTools.coordinates.length; i++) {

            //if (GeometryTools.coordinates[i] == null) {
                //break;
            //}

            //Master.log("mito #" + i + " " + GeometryTools.coordinates[i].dimensionsToString(true));
            //Master.log(GeometryTools.coordinates[i].scalesToString(true));

        //}

        return false;

    }

    public boolean initGeometryMitochondria2() {

        double dran, volumeFraction;

        int znob = 2;

        int linkNum = 3;
        int extraVoxelsXYZ = 25;

        RunMonteCarloAZ mc;

        Geometry geometry = project.geometry;

        CoordinatesVoxels c = new CoordinatesVoxels(project, geometry);

        if (nob) {
            c.setVoxels(geometry.xVoxel1, geometry.yVoxel1, geometry.zVoxel1, geometry.xVoxel2, geometry.yVoxel2, geometry.zVoxel2-znob);
        } else {
            c.setVoxels(geometry.xVoxel1, geometry.yVoxel1, geometry.zVoxel1, geometry.xVoxel2, geometry.yVoxel2, geometry.zVoxel2);
        }

        initGeometryMC2();

        if (mito_ijkselect == -2) {

            dran = ran.nextDouble();

            if (dran < 0.3333333333) {
                mito_ijkselect = 0; // xy plane
            } else if (dran < 0.6666666666) {
                mito_ijkselect = 1; // yz plane
            } else {
                mito_ijkselect = 2; // zx plane
            }

        }

        if (mitoVolumeFractionExact) {
            volumeFraction = GeometryTools.addEllipsoids(geometry, c, mitoRadius, mitoAxialRatio, mitoVolumeFraction, -1, mitoVolumeFractionExact, mito_ijkselect, mitoCylinder, linkNum, extraVoxelsXYZ);
        } else {
            volumeFraction = GeometryTools.addEllipsoids(geometry, c, mitoRadius, mitoAxialRatio, mitoVolumeFraction - mitoVolumeFractionTolerance, -1, mitoVolumeFractionExact, mito_ijkselect, mitoCylinder, linkNum, extraVoxelsXYZ);
        }

        if ((volumeFraction < mitoVolumeFraction - mitoVolumeFractionTolerance) || (volumeFraction > mitoVolumeFraction + mitoVolumeFractionTolerance)) {
            return true; // error
        }

        if (project.monteCarlo instanceof RunMonteCarloAZ) {

            mc = (RunMonteCarloAZ) project.monteCarlo;

            double azWidthVoxels = mc.azWidth / project.dx;
            double azHeightVoxels = mc.azHeight / project.dx;
            int extra = 0;//1;

            int i1 = (int) Math.round(geometry.xVoxelCenter - (azWidthVoxels/2) - extra);
            int i2 = (int) Math.round(geometry.xVoxelCenter + (azWidthVoxels/2) + extra);
            int j1 = i1;
            int j2 = i2;
            int k1 = (int) Math.round(geometry.zVoxel2 - azHeightVoxels - extra);
            int k2 = geometry.zVoxel2;

            //Master.log("i1=" + i1);
            //Master.log("i2=" + i2);
            //Master.log("k1=" + k1);
            //Master.log("k2=" + k2);

            for (int k = k1; k <= k2; k++) {
                for (int j = j1; j <= j2; j++) {
                    for (int i = i1; i <= i2; i++) {
                        if (geometry.space[i][j][k] == -1) {
                            //Master.log("AZ non-space");
                            return true; // error, non-space voxel overlaps with active zone
                        }
                    }
                }
            }

        }

        //shape.checkSpace();

        //mitoVolumeFractionActual = (100.0 - 100.0 * shape.spaceVoxels / shape.voxels); // does not work for spherical shape
        mitoVolumeFractionActual = volumeFraction;

        return false;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof InitProject)) {
            return false;
        }

        return false;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof InitProject)) {
            return false;
        }

        String n = o.getName();

        return false;

    }

}
