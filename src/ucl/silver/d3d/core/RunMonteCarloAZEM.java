package ucl.silver.d3d.core;
import java.io.*;
import java.util.Scanner;

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
public class RunMonteCarloAZEM extends RunMonteCarloAZ {
    
    public int EMseries = 1;
    public int EMseriesAZ = 0;

    public double extendAZ = 0; // um
    public double xStep_T = -1; // fraction 0 - 1

    public boolean fixNonSpaceVoxels = true;
    public boolean removeShrinkage = false;

    public boolean importVesicleXYZ = false;
    public String importVesicleXYZ_dr = "";
    //public double removeOverlapLimitD = 0.025;

    public boolean extraReserve = false;
    private int newReserveVoxels = 0;

    public boolean emptyGeometry = false;
    public boolean openGeometry = false;
    private int openGeometryNewVoxels = 0;
    public int openGeometryForwardX = 0;
    public int openGeometryForwardY = 0;
    public int openGeometryForwardZ = 0;
    private int[][][] openGeometrySpace;
    private boolean openGeometrySpaceInit = false;

    public double reserveVesicleDensity = 0.17; // Zoltan
    public double reserveImmobileFraction = 0.25;
    public boolean reserveX1 = true;
    public boolean reserveX2 = true;
    public boolean reserveY1 = true;
    public boolean reserveY2 = true;
    public boolean reserveZ1 = false;
    public boolean reserveZ2 = false;

    public double dockRefractoryPeriod = 0;

    private CoordinatesVoxels EMcoordinates = null;
    private Voxel[][][] voxelSpace = null;

    private static double minD2AZ_x, minD2AZ_y, minD2AZ_z;

    public RunMonteCarloAZEM(Project p) {
        super(p);
    }

    @Override
    public boolean initAll(){

        if (PBC) {
            //Master.exit("PBC not allowed with RunMonteCarloAZEM");
        }

        releaseCounter = 0;

        if (checkVariables()) {
            return true;
        }

        if (initDiffusantVesiclesArray()) {
            return true;
        }

        //if (initDetectors()){
        //     return true;
        //}

        //if (initPSFs()) {
        //    return true;
        //}

        if (initVoxels()) {
            return true;
        }

        if (initActiveZone()) {
            return true;
        }

        if (initGeometry()) {
            return true;
        }

        if (initVesicles()) {
            return true;
        }

        if (checkDX()) {
            return true;
        }

        if (initDT()) {
            return true;
        }

        //if (initVesiclesDocked()){
        //    return true;
        //}

        if (initVesicleRandom && initVesiclesRandom()) {
            return true;
        }

        //minDistanceSlow();

        

        //if (initReserveList()) {
        //    return true;
        //}

        //if (removeVesicleOverlap && removeVesicleOverlap()) {
        //    return true;
        //}

        //if (initVesiclesImmobile()) {
        //    return true;
        //}

        initConnectors();

        vesicleStats();

        replenishCoordinates = new Coordinates(project, geometry);

        if (!freeDiffusion && (vesicleVolumeFraction > 0.44)) {
            replenishRate = 0; // TAKES TOO LONG
            Master.log("TURNED OFF REPLENISHING: replenishRate = 0");
        }

        //initAZ_patch();

        //Master.updatePanel2D();

        if (PBC) {
            voxelSpace = geometry.voxelSpacePBC;
        } else {
            voxelSpace = geometry.voxelSpace;
        }
        
        fillEmptySpace2();
        removeIsolatedNonSpaceVoxels(2);
        //newVoxelSize(0.022); // for EM density analysis
        initVoxels();
        initVesicles();

        if (fixAZvoxels()) {
            initVoxels();
            initVesicles();
        }

        if (removeShrinkage) {
            removeShrinkage();
            initDT();
            vesicleStats();
        }

        testVesicleOverlapEM();

        if (importVesicleXYZ) {
            if (importVesicleXYZCanRead()) {
                removeVesicleOverlap = false;
            } else {
                Master.log("cannot read vesicle positions from file " + importVesicleXYZfilename());
                importVesicleXYZ = false;
            }
        }

        if (removeVesicleOverlap && removeVesicleOverlap()) {
            return true;
        }

        if (extraReserve) {

            addReserveSpace();

            countSpaceVoxels(EMcoordinates);

            if (openGeometry) {
                openGeometry();
            } else if (emptyGeometry) {
                emptyGeometry();
            }

            countSpaceVoxels(EMcoordinates);

            addReserveVesicles();

        } else {

            EMcoordinates = new CoordinatesVoxels(project, geometry.xVoxel1 + 1, geometry.yVoxel1 + 1, geometry.zVoxel1 + 1, geometry.xVoxel2 - 1, geometry.yVoxel2 - 1, geometry.zVoxel2 - 1);

        }

        //percentVoxelSpace(true);

        if (importVesicleXYZ) {
            importVesicleXYZ();
        }

        int overlaps = testVesicleOverlapEM();

        if ((overlaps != 0) && removeVesicleOverlap) {
            Master.exit("encountered vesicle overlaps n=" + overlaps);
        }

        setDistanceToAZ(true);

        double saveImmobilePercent = diffusants[0].setImmobilePercent;

        diffusants[0].setImmobilePercent = reserveImmobileFraction;

        if (initVesiclesImmobile()) {
            return true;
        }

        diffusants[0].setImmobilePercent = saveImmobilePercent;

        if (diffusants[0].setImmobilePercent == 0) {
            removeImmobileVesicles();
        }

        initConnectors();

        vesicleStats();
        initVesiclesDsDff();

        connectVesicles();

        if (initVesiclesDocked()){
            return true;
        }

        if (initDockedList()) {
            return true;
        }

        if (false) {
            vesicleDensity2(true, false, "");
            //vesicleDensityFinal(true, false);
            return false;
        }

        initMSDspatialTimes();

        //EM_vesicleDensity3(0.0044, 0.005, true, false, "");
        //EM_vesicleDensity3(0.03, true, false, "");

        //Master.log("min distance between vesicles " + minDistance());

        if (hydrodynamicsLocalD) {
            localDensityAll(true);
        }

        if (dockingOffAfterInit) {
            dockingOn = false; // TEMPORARY removal of docking
        }

        return false;

    }

    @Override
    public boolean initSimulation() {
        boolean error = super.initSimulation();
        return error;
    }

    @Override
    public void finishSimulation(){
        super.finishSimulation();
        vesicleDensity2(true, true, "finish");
        //vesicleDensityFinal(true, true);
    }

    @Override
    public boolean initActiveZone() {

        double x0, x1, y0, y1, z0, z1;

        if ((geometry == null) || (activeZone == null)) {
            return true;
        }

        x0 = azx1[EMseriesAZ] + extendAZ;
        y0 = azy1[EMseriesAZ] + extendAZ;
        z0 = azz1[EMseriesAZ] + extendAZ;
        x1 = azx2[EMseriesAZ] + extendAZ;
        y1 = azy2[EMseriesAZ] + extendAZ;
        z1 = azz2[EMseriesAZ] + extendAZ;

        activeZone.setShape("cuboid");
        activeZone.setCoordinates(x0, y0, z0, x1, y1, z1);

        if (tetherVesicles) {

            x0 += azTetherLength / 2.0;
            y0 += azTetherLength / 2.0;
            z0 += azTetherLength / 2.0;
            x1 += azTetherLength / 2.0;
            y1 += azTetherLength / 2.0;
            z1 += azTetherLength / 2.0;
            activeZoneTethering = new Coordinates(project, x0, y0, z0, x1, y1, z1);
            activeZoneTethering.setShape("cuboid");

            if (azNumTethers > 0) {
                azDV = new DiffusantVesicle(project, "AZ", 0, 0, 0, 0, 0);
                azDV.connectTo = new DiffusantVesicle[azNumTethers];
                azDV.connectorOffTime = new double[azNumTethers];
            }

        }

        return false;

    } 

    private void initMSDspatialTimes() {

        if (!saveMSDspatial) {
            return;
        }

        project.simTime = MSDspatialTbgn + MSDspatial_win + 2;

        Master.log("simulation time = " + project.simTime);

    }

    @Override
    public void initD20() {

        DiffusantVesicleAZ vaz;
        
        String simTag = Integer.toString((int)time) + "ms";
        
        vesicleDensity2(true, true, simTag);

        releaseRate = 0; // STOP RELEASING VESICLES

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {
                
                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (!vaz.mobile || !vaz.insideGeometry) {
                    continue;
                }

                vaz.x0 = vaz.x;
                vaz.y0 = vaz.y;
                vaz.z0 = vaz.z;

                if (EMcoordinates.isInside(vaz.x, vaz.y, vaz.z)) {
                    vaz.d2AZ0 = minDistanceToAZ(azXYZ, EMseriesAZ, vaz.x, vaz.y, vaz.z);
                } else {
                    vaz.d2AZ0 = -1;
                }

            }

        }

    }

    @Override
    public double MSDspatial(int ibin) {

        double sumSD = 0.0, count = 0.0;
        double x1 = ibin * MSDspatial_binWidth;
        double x2 = x1 + MSDspatial_binWidth;

        DiffusantVesicleAZ vaz;

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {
                
                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (!vaz.mobile || !vaz.insideGeometry) {
                    continue;
                }

                //if (!EMcoordinates.isInside(dv.x, dv.y, dv.z)) {
                //    continue;
                //}

                if (vaz.d2AZ0 < 0) {
                    continue;
                }

                if ((vaz.d2AZ0 >= x1) && (vaz.d2AZ0 < x2)) {
                    sumSD += vaz.squareDisplacement();
                    count += 1.0;
                    //Master.log("" + dv + " " + dv.d2AZ0);
                }

            }

        }

        return (sumSD / count);

    }

    @Override
    public boolean setVesicleLocation(DiffusantVesicle dv, double x, double y, double z, boolean initStartLocation) {

        int xVoxel, yVoxel, zVoxel;

        if (openGeometrySpaceInit) {

            xVoxel = (int) geometry.computeVoxelX(x);
            yVoxel = (int) geometry.computeVoxelY(y);
            zVoxel = (int) geometry.computeVoxelZ(z);

            if (openGeometrySpace[xVoxel][yVoxel][zVoxel] == 0) {
                return false;
            }

        }

        return super.setVesicleLocation(dv, x, y, z, initStartLocation);

    }

    @Override
    public double localDensityVoxels(DiffusantVesicle dv) {

        if (EMcoordinates == null) {
            dv.localStep3 = dv.step3;
            return Double.NaN;
        }

        if (!EMcoordinates.isInside(dv.x, dv.y, dv.z)) {
            dv.localStep3 = dv.step3;
            return Double.NaN;
        }

        return super.localDensityVoxels(dv);

    }

    @Override
    public boolean moveVesicleGauss(DiffusantVesicleAZ vaz) {

        double step3;
        //double stepx, stepy, stepz;
        double dx, dy, dz, d, dlimit;

        if (hydrodynamicsLocalD) {
            step3 = vaz.localStep3;
        } else {
            step3 = vaz.step3;
        }

        if (removingVesicleOverlap) {
            step3 = removeVesicleOverlapStep3;
        }

        //if (rstep <= 0) {
        //    Master.exit("moveVesicleGauss: bad rstep: " + rstep);
        //}

        if (removingVesicleOverlap) {

            if (vaz.d2AZ < 0.2) {
                stepx = ranGauss() * step3 * 0.001;
                stepy = ranGauss() * step3 * 0.001;
                stepz = ranGauss() * step3 * 0.01;
                dlimit = vaz.radius * 1.1; // slightly larger than 1.0 for S04AZ3
            } else {
                stepx = ranGauss() * step3 * 0.01;
                stepy = ranGauss() * step3 * 0.01;
                stepz = ranGauss() * step3;
                dlimit = vaz.radius * 2;
            }

            dx = vaz.x + stepx - vaz.x0;
            dy = vaz.y + stepy - vaz.y0;
            dz = vaz.z + stepz - vaz.z0;

            d = Math.sqrt(dx * dx + dy * dy + dz * dz);

            if (d > dlimit) {
                return false; // too far from original center point
            }

        } else {

            stepx = ranGauss() * step3;
            stepy = ranGauss() * step3;
            stepz = ranGauss() * step3;

        }

        return setVesicleLocation(vaz, vaz.x + stepx, vaz.y + stepy, vaz.z + stepz, false);

    }

    @Override
    public boolean moveVesicleGaussHydroWallz(DiffusantVesicleAZ vaz) {

        double step3, d2az, dlimit, dx, dy, dz, d;
        double b_ll, b_T;

        if (removingVesicleOverlap) {

            step3 = removeVesicleOverlapStep3;

            if (vaz.d2AZ < 0.2) {
                stepx = ranGauss() * step3 * 0.001;
                stepy = ranGauss() * step3 * 0.001;
                stepz = ranGauss() * step3 * 0.01;
                dlimit = vaz.radius * 1.1; // slightly larger than 1.0 for S04AZ3
            } else {
                stepx = ranGauss() * step3 * 0.01;
                stepy = ranGauss() * step3 * 0.01;
                stepz = ranGauss() * step3;
                dlimit = vaz.radius * 2;
            }

            dx = vaz.x + stepx - vaz.x0;
            dy = vaz.y + stepy - vaz.y0;
            dz = vaz.z + stepz - vaz.z0;

            d = Math.sqrt(dx * dx + dy * dy + dz * dz);

            if (d > dlimit) {
                return false; // too far from original center point
            }

        } else {

            step3 = vaz.step3;

            if ((time > 0) && EMcoordinates.isInside(vaz.x, vaz.y, vaz.z)) {

                if (hydrodynamicsLocalD) {
                    step3 = vaz.localStep3;
                }

                d2az = minDistanceToAZ(azXYZ, EMseriesAZ, vaz.x, vaz.y, vaz.z);

                b_ll = hydroWall_ll(vaz.radius, d2az);
                b_T = hydroWall_T(vaz.radius, Math.abs(d2az - vaz.radius));

                //if (d2az < 0.045) {
                //    Master.log("" + d2az + " " + b_T);
                //}

                //if (b_T < 0.001) {
                    //Master.log("" + decimal3.format(Math.abs(d2az - dv.radius)) + ", " + decimal3.format(b_T));
                    //Master.log("" + decimal3.format(hydroWall_T(dv.radius, 0.0001)));
                //}

                b_ll = 1 / (1 + vaz.DsDff * ((1 / b_ll) - 1)); // Michailidou et al. 2009
                b_T = 1 / (1 + vaz.DsDff * ((1 / b_T) - 1)); // Michailidou et al. 2009

                b_ll = Math.sqrt(b_ll);
                b_T = Math.sqrt(b_T);

                dx = ranGauss() * step3;
                dy = ranGauss() * step3;
                dz = ranGauss() * step3;

                stepx = dx * b_T * xStep_T + dx * b_ll * (1 - xStep_T);
                stepy = dy * b_T * (1 - xStep_T) + dy * b_ll * xStep_T;
                stepz = dz * b_ll;

            } else {

                stepx = ranGauss() * step3;
                stepy = ranGauss() * step3;
                stepz = ranGauss() * step3;

            }

        }

        return setVesicleLocation(vaz, vaz.x + stepx, vaz.y + stepy, vaz.z + stepz, false);

    }
    
    @Override
    public boolean isInsideActiveZone(DiffusantVesicleAZ vaz) {
        
        double d;

        boolean inside = activeZone.isInsideCuboid(vaz.x, vaz.y, vaz.z);

        if (inside) { // more refined search

            d = minDistanceToAZ(azXYZ, EMseriesAZ, vaz.x, vaz.y, vaz.z);
            
            inside = (d < vaz.radius);

        }
        
        return inside;

    }

    @Override
    public boolean isInsideActiveZoneTethering(DiffusantVesicleAZ vaz) {

        double d;

        boolean inside = activeZoneTethering.isInside(vaz.x, vaz.y, vaz.z);

        if (inside) { // more refined search

            d = minDistanceToAZ(azXYZ, EMseriesAZ, vaz.x, vaz.y, vaz.z);
            
            inside = (d < vaz.radius + azTetherLength);

        }

        return inside;

    }

    private void setDistanceToAZ(boolean setD2AZ0){

        DiffusantVesicleAZ vaz;

        if (diffusants == null) {
            return;
        }

        if (EMcoordinates == null) {
            return;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {
                
                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (Double.isNaN(vaz.x)) {
                    continue;
                }

                if (EMcoordinates.isInside(vaz.x, vaz.y, vaz.z)) {
                    vaz.d2AZ = minDistanceToAZ(azXYZ, EMseriesAZ, vaz.x, vaz.y, vaz.z);
                    if (setD2AZ0) {
                        vaz.d2AZ0 = vaz.d2AZ;
                    }
                }

            }
        }

    }

    public static double[][] minDistanceToAZ(double[][][] azXYZ, double[][] vesicleXYZ) {

        double vx, vy, vz;

        if ((azXYZ == null) || (vesicleXYZ == null)) {
            return null;
        }

        double[][] vd = new double[azXYZ.length][vesicleXYZ.length];

        for (int iAZ = 0; iAZ < azXYZ.length; iAZ++) {
            for (int iv = 0; iv < vesicleXYZ.length; iv++) {
                vd[iAZ][iv] = Double.NaN;
            }
        }

        for (int iAZ = 0; iAZ < azXYZ.length; iAZ++) {
            for (int iv = 0; iv < vesicleXYZ.length; iv++) {
                vx = vesicleXYZ[iv][0];
                vy = vesicleXYZ[iv][1];
                vz = vesicleXYZ[iv][2];
                vd[iAZ][iv] = minDistanceToAZ(azXYZ, iAZ, vx, vy, vz);
            }
        }

        for (int iAZ = 0; iAZ < vd.length; iAZ++) {

            Master.log("AZ #" + iAZ);

            for (int iv = 0; iv < vd[0].length; iv++) {
                if ((iv == 0) && Double.isInfinite(vd[iAZ][iv])) {
                    break;
                }
                if (Double.isInfinite(vd[iAZ][iv])) {
                    Master.log("found vesicle Inf : " + iv);
                }
                if (Double.isNaN(vd[iAZ][iv])) {
                    Master.log("found vesicle NaN : " + iv);
                }
                Master.log("" + vd[iAZ][iv]);
            }

            Master.log(" ");

        }

        return vd;

    }

    public static double minDistanceToAZ(double[][][] azXYZ, int azSelect, double vx, double vy, double vz) {

        double dx, dy, dz;
        double d, minD = Double.POSITIVE_INFINITY;

        if (azXYZ == null) {
            return Double.NaN;
        }

        if (Double.isNaN(vx * vy * vz)) {
            return Double.NaN;
        }

        for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {

            if (Double.isNaN(azXYZ[azSelect][ipnt][0] * azXYZ[azSelect][ipnt][1] * azXYZ[azSelect][ipnt][2])) {
                continue;
            }

            dx = vx - azXYZ[azSelect][ipnt][0];
            dy = vy - azXYZ[azSelect][ipnt][1];
            dz = vz - azXYZ[azSelect][ipnt][2];

            d = Math.sqrt(dx * dx + dy * dy + dz * dz);

            if (d < minD) {
                minD = d;
                minD2AZ_x = azXYZ[azSelect][ipnt][0];
                minD2AZ_y = azXYZ[azSelect][ipnt][1];
                minD2AZ_z = azXYZ[azSelect][ipnt][2];
            }

        }

        return minD;

    }

    private boolean addReserveSpace() {

        int oldSpaceVoxels = 0, newSpaceVoxels = 0;

        int xVoxelMin = geometry.space.length;
        int yVoxelMin = geometry.space.length;
        int zVoxelMin = geometry.space.length;
        int xVoxelMax = 0;
        int yVoxelMax = 0;
        int zVoxelMax = 0;

        for (int i = 0; i < geometry.space.length; i++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int k = 0; k < geometry.space[0][0].length; k++) {
                    if (geometry.space[i][j][k] == 1) {
                        xVoxelMin = Math.min(xVoxelMin, i);
                        yVoxelMin = Math.min(yVoxelMin, j);
                        zVoxelMin = Math.min(zVoxelMin, k);
                        xVoxelMax = Math.max(xVoxelMax, i);
                        yVoxelMax = Math.max(yVoxelMax, j);
                        zVoxelMax = Math.max(zVoxelMax, k);
                        oldSpaceVoxels++;
                    }
                }
            }
        }

        EMcoordinates = new CoordinatesVoxels(project, xVoxelMin, yVoxelMin, zVoxelMin, xVoxelMax, yVoxelMax, zVoxelMax);
        replenishCoordinatesExclude = new Coordinates(project, EMcoordinates);

        //Master.log(EMcoordinates.voxelsToString(false));

        //Master.log("old space voxels = " + oldSpaceVoxels);
        //Master.log("iVoxelMin = " + xVoxelMin + ", iVoxelMax = " + xVoxelMax);
        //Master.log("jVoxelMin = " + yVoxelMin + ", jVoxelMax = " + yVoxelMax);
        //Master.log("kVoxelMin = " + zVoxelMin + ", kVoxelMax = " + zVoxelMax);

        if (reserveX1) {
            for (int k = 0; k <= geometry.zVoxel2; k++) {
                for (int j = 0; j <= geometry.yVoxel2; j++) {
                    for (int i = 0; i <= xVoxelMin - 1; i++) {
                        geometry.setSpace(i, j, k, 1);
                    }
                }
            }
        }

        if (reserveX2) {
            for (int k = 0; k <= geometry.zVoxel2; k++) {
                for (int j = 0; j <= geometry.yVoxel2; j++) {
                    for (int i = xVoxelMax + 1; i <= geometry.xVoxel2; i++) {
                        geometry.setSpace(i, j, k, 1);
                    }
                }
            }
        }

        if (reserveY1) {
            for (int k = 0; k <= geometry.zVoxel2; k++) {
                for (int j = 0; j <= yVoxelMin - 1; j++) {
                    for (int i = 0; i <= geometry.xVoxel2; i++) {
                        geometry.setSpace(i, j, k, 1);
                    }
                }
            }
        }

        if (reserveY2) {
            for (int k = 0; k <= geometry.zVoxel2; k++) {
                for (int j = yVoxelMax + 1; j <= geometry.yVoxel2; j++) {
                    for (int i = 0; i <= geometry.xVoxel2; i++) {
                        geometry.setSpace(i, j, k, 1);
                    }
                }
            }
        }

        if (reserveZ1) {
            for (int k = 0; k <= zVoxelMin - 1; k++) {
                for (int j = 0; j <= geometry.yVoxel2; j++) {
                    for (int i = 0; i <= geometry.xVoxel2; i++) {
                        geometry.setSpace(i, j, k, 1);
                    }
                }
            }
        }

        if (reserveZ2) {
            for (int k = zVoxelMax + 1; k <= geometry.zVoxel2; k++) {
                for (int j = 0; j <= geometry.yVoxel2; j++) {
                    for (int i = 0; i <= geometry.xVoxel2; i++) {
                        geometry.setSpace(i, j, k, 1);
                    }
                }
            }
        }

        for (int i = 0; i < geometry.space.length; i++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int k = 0; k < geometry.space[0][0].length; k++) {
                    if (geometry.space[i][j][k] == 1) {
                        newSpaceVoxels++;
                    }
                }
            }
        }

        newReserveVoxels = newSpaceVoxels - oldSpaceVoxels;
        Master.log("added reserve space voxels n = " + newReserveVoxels);

        return false;

    }
    
    private String importVesicleXYZfilename() {
        
        String dr, azSeries, fileName;
        
        if (EMseries < 10) {
            azSeries = "S0" + Integer.toString(EMseries) + "AZ" + Integer.toString(EMseriesAZ);
        } else {
            azSeries = "S" + Integer.toString(EMseries) + "AZ" + Integer.toString(EMseriesAZ);
        }

        dr = importVesicleXYZ_dr + azSeries + "/";

        fileName = "XYZS_R_" + azSeries + "_im25_" + Integer.toString(project.batchNum) + ".dat";
        
        return dr + fileName;
        
    }

    private boolean importVesicleXYZCanRead() {

        String filePath = importVesicleXYZfilename();
        
        File f = new File(filePath);

        return f.canRead();

    }

    private boolean importVesicleXYZ() {

        DiffusantVesicleAZ vaz;

        String filePath = importVesicleXYZfilename();

        double[][] vesicleXYZ3 = importVesiclesXYZ(filePath);

        if (vesicleXYZ3 == null) {
            return true;
        }

        if (vesicleXYZ3.length != diffusants[0].vesicles.length) {
            Master.exit("wrong number of vesicles from file " + filePath);
        }

        initVoxels();

        for (int i = 0; i < diffusants[0].vesicles.length; i++) {

            vaz = (DiffusantVesicleAZ) diffusants[0].vesicles[i];

            if (!setVesicleLocation(vaz, vesicleXYZ3[i][0], vesicleXYZ3[i][1], vesicleXYZ3[i][2], true)) {
                return true;
            }

            if (!addToVoxelList(vaz)) {
                return true;
            }

        }

        vesicleStats();

        return false;

    }

    private double[][] importVesiclesXYZ(String filePath) {

        String s, ss;

        int ichar, i, j, k1, k2, istart = 0, counter = 0;
        char c;

        int nVesicles = diffusants[0].vesicles.length;

        double[][] xyz = new double[nVesicles][3];

        for (i = 0; i < xyz.length; i++) {
            xyz[i][0] = Double.NaN;
            xyz[i][1] = Double.NaN;
            xyz[i][2] = Double.NaN;
        }

        try {

            s = new Scanner(new File(filePath)).useDelimiter("\\A").next();

            Master.log("importing vesicle positions from file " + filePath);

            for (int is = 0; is < s.length()-1; is++) {
                if (s.regionMatches(is, "format=Unicode Text", 0, 19)) { // end of header
                    istart = is + 20;
                }
            }

            i = 0;
            j = 0;
            k1 = -1;

            for (int is = istart; is < s.length(); is++) {

                c = s.charAt(is);
                ichar = (int) c;

                //Master.log("" + c + " " + ichar + " " + Character.isDigit(c));

                if (Character.isDigit(c)) {

                    if (k1 == -1) {
                        k1 = is; // start of new number
                    }

                } else {

                    if ((c == '.') || ((c == 'E'))) { // ichar = 46
                        //
                    } else if (c == '-') { // ichar = 45
                        if (k1 == -1) {
                            k1 = is; // start of new number
                        }
                    } else if ((c == '\t') || (ichar == 9) || (ichar == 10) || (ichar == 13)) {

                        k2 = is;
                        ss = s.substring(k1, k2);
                        //Master.log("" + ss);

                        try {
                            xyz[i][j] = Double.parseDouble(ss);
                        } catch (NumberFormatException e) {
                            System.err.println("parse to double error: " + ss);
                            return null;
                        }

                        j++;

                        if (j == 3) {
                            j = 0;
                            i++; // next vesicle
                        }

                        k1 = -1;
                        
                    } else {
                        Master.exit("unrecognized character: " + is + " " + c + " " + ichar + " " + Character.isDigit(c));
                    }

                }

            }

            for (i = 0; i < xyz.length; i++) {
                if (Double.isNaN(xyz[i][0])) {
                    break;
                }
                //Master.log("" + xyz[i][0] + ", " + xyz[i][1] + ", " + xyz[i][2]);
                counter++;
            }

            Master.log("imported non-overlapping vesicle positions n = " + counter);

        } catch (IOException e) {
            //System.err.println("unable to read file: " + WRL_fileName);
            Master.log("unable to read file: " + filePath);
            return null;
        }

        return xyz;

    }

    private boolean addReserveVesicles() {

        DiffusantVesicleAZ vaz;

        int numReserveVesiclesToAdd1, numReserveVesiclesToAdd2, newNumVesicles;

        int oldNumVesicles = numVesicles;
        //double radius = diffusant[0].setMeanRadius;
        //double vesicleVolume = 4.0 * Math.PI * radius * radius * radius / 3.0;
        double addedReserveVolume1 = newReserveVoxels * project.dx * project.dx * project.dx;
        double addedReserveVolume2 = openGeometryNewVoxels * project.dx * project.dx * project.dx;

        Coordinates c = new Coordinates(project, geometry);
        Coordinates emc = new Coordinates(project, EMcoordinates);

        geometry.checkSpace();
        initVoxels();

        numReserveVesiclesToAdd1 = (int) (reserveVesicleDensity * addedReserveVolume1 / vesicleVolume);
        numReserveVesiclesToAdd2 = (int) (reserveVesicleDensity * addedReserveVolume2 / vesicleVolume);

        newNumVesicles = oldNumVesicles + numReserveVesiclesToAdd1 + numReserveVesiclesToAdd2;

        DiffusantVesicleAZ[] vesicleTemp = new DiffusantVesicleAZ[newNumVesicles];

        for (int i = 0; i < oldNumVesicles; i++) {
            vesicleTemp[i] = (DiffusantVesicleAZ) diffusants[0].vesicles[i];
            if (!setVesicleLocation(vesicleTemp[i], vesicleTemp[i].x, vesicleTemp[i].y, vesicleTemp[i].z, true)) {
                return true;
            }
            if (!addToVoxelList(vesicleTemp[i])) {
                return true;
            }
        }

        for (int i = oldNumVesicles; i < vesicleTemp.length; i++) {
            vesicleTemp[i] = new DiffusantVesicleAZ(project, "ready", diffusants[0].setMeanRadius, diffusants[0].D, Double.NaN, Double.NaN, Double.NaN);
        }

        diffusants[0].vesicles = vesicleTemp;
        diffusants[0].numVesicles = diffusants[0].vesicles.length;

        vesicleStats();

        numVesicles = diffusants[0].vesicles.length;

        int ii = oldNumVesicles;

        for (int i = 0; i < numReserveVesiclesToAdd1; i++) {
            if (initVesicleRandom(diffusants[0].vesicles[ii], c, emc, false, true, -1, overlapTrialLimit, absLimit, Double.NaN)) {
                return true;
            }
            ii++;
        }

        Master.log("added reserve vesicles n = " + numReserveVesiclesToAdd1);

        openGeometrySpaceInit = true;

        for (int i = 0; i < numReserveVesiclesToAdd2; i++) {
            if (initVesicleRandom(diffusants[0].vesicles[ii], emc, null, false, true, -1, overlapTrialLimit, absLimit, Double.NaN)) {
                return true;
            }
            vaz = (DiffusantVesicleAZ) diffusants[0].vesicles[ii];
            vaz.d2AZ = minDistanceToAZ(azXYZ, EMseriesAZ, vaz.x, vaz.y, vaz.z);
            //Master.log("" +  dv.d2AZ);
            ii++;
        }

        openGeometrySpaceInit = false;

        Master.log("added open geometry vesicles n = " + numReserveVesiclesToAdd2);

        return false;

    }

    private void removeImmobileVesicles() {

        DiffusantVesicle dv;

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {
            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {
                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                    dv = voxelSpace[i][j][k].firstReady;

                    while (dv != null) {
                        dv.mobile = true;
                        dv = dv.nextReady;
                    }

                }
            }
        }

    }

    private void initVesiclesDsDff() {

        double hydroDsDff;
        //double volumeFraction = 0.17;

        if (hydroWallZ) {

            for (DiffusantVesicles d : diffusants) {

                if ((d == null) || (d.vesicles == null)) {
                    continue;
                }

                //hydroDsDff = DiffusantVesicle.Dratio_short(volumeFraction, diffusant[i].setImmobilePercent);
                //hydroDsDff /= DiffusantVesicle.Dff_short_Banchio(volumeFraction * (1 - diffusant[i].setImmobilePercent));
                hydroDsDff = DiffusantVesicle.Dratio_shortNew(d.mobileVolumeFraction, d.immobileVolumeFraction);
                hydroDsDff /= DiffusantVesicle.Dff_short_Banchio(d.mobileVolumeFraction);

                //Master.log("hydroDsDff = " + hydroDsDff);

                for (DiffusantVesicle v : d.vesicles) {
                    v.DsDff = hydroDsDff;
                }

            }

        }

    }

    @Override
    public boolean initVesiclesDocked() {

        int count = 0;
        double d2AZ, minD2AZ = 99999;
        double minDistanceForDocking = 0.026;
        DiffusantVesicleAZ vaz;

        if (!dockingOn) {
            return false;
        }

        if (diffusants == null) {
            return true;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {
                
                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                d2AZ = minDistanceToAZ(azXYZ, EMseriesAZ, vaz.x, vaz.y, vaz.z);
                minD2AZ = Math.min(d2AZ, minD2AZ);
                
                if ((count < dockedVesiclesMax) && (d2AZ < minDistanceForDocking)) {
                    vaz.setType("docked");
                    count++;
                }

                if (dockRefractoryPeriod > 0) {
                    vaz.dockRefractoryPeriod = dockRefractoryPeriod;
                }

            }

        }

        //Master.log("defined docked vesicles: " + count);

        return false;

    }

    private void fillEmptySpace() {

        boolean empty;
        int extra = 1;

        geometry.simpleCuboid = false;

        for (int i = 0; i < geometry.space.length; i++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int k = 0; k < geometry.space[0][0].length; k++) {
                    if (voxelSpace[i][j][k].firstReady == null) {

                        empty = true;

                        for (int ii = i - extra; ii <= i + extra; ii++) {

                            if ((ii < 0) || (ii >= voxelSpace.length)) {
                                continue;
                            }

                            for (int jj = j - extra; jj <= j + extra; jj++) {

                                if ((jj < 0) || (jj >= voxelSpace[0].length)) {
                                    continue;
                                }

                                for (int kk = k - extra; kk <= k + extra; kk++) {

                                    if ((kk < 0) || (kk >= voxelSpace[0][0].length)) {
                                        continue;
                                    }

                                    if (voxelSpace[ii][jj][kk].firstReady != null) {
                                        empty = false;
                                    }

                                }
                            }
                        }

                        if (empty) {
                            geometry.setSpace(i, j, k, -1);
                            voxelSpace[i][j][k].isSpace = false;
                        }

                    }
                }
            }
        }

        int i1 = (int) geometry.computeVoxelX(activeZone.x1);
        int i2 = (int) geometry.computeVoxelX(activeZone.x2);
        int j1 = (int) geometry.computeVoxelY(activeZone.y1);
        int j2 = (int) geometry.computeVoxelY(activeZone.y2);
        int k1 = (int) geometry.computeVoxelZ(activeZone.z1);
        int k2 = (int) geometry.computeVoxelZ(activeZone.z2);

        for (int i = i1; i <= i2; i++) {
            for (int j = j1; j <= j2; j++) {
                for (int k = k1; k <= k2; k++) {
                    //geometry.setSpace(i, j, k, 1);
                    //voxelSpace[i][j][k].isSpace = true;
                }
            }
        }

    }

    private void fillEmptySpace2() {

        double dlimit;

        DiffusantVesicle dv;
        Voxel v2;

        int i1 = 0;
        int i2 = voxelSpace.length - 1;
        int j1 = 0;
        int j2 = voxelSpace[0].length - 1;
        int k1 = 0;
        int k2 = voxelSpace[0][0].length - 1;

        geometry.simpleCuboid = false;

        for (int i = i1; i <= i2; i++) {
            for (int j = j1; j <= j2; j++) {
                for (int k = k1; k <= k2; k++) {
                    geometry.setSpace(i, j, k, -1);
                }
            }
        }

        for (int i = i1; i <= i2; i++) {
            for (int j = j1; j <= j2; j++) {
                for (int k = k1; k <= k2; k++) {

                    dv = voxelSpace[i][j][k].firstReady;

                    while (dv != null) {

                        geometry.setSpace(i, j, k, 1);
                        //vSpace[i][j][k].isSpace = true;

                        dlimit = project.dx * 0.5 + dv.radius; // um

                        for (int ii = i - 1; ii <= i + 1; ii++) {

                            if ((ii < 0) || (ii >= voxelSpace.length)) {
                                continue;
                            }

                            for (int jj = j - 1; jj <= j + 1; jj++) {

                                if ((jj < 0) || (jj >= voxelSpace[0].length)) {
                                    continue;
                                }

                                for (int kk = k - 1; kk <= k + 1; kk++) {

                                    if ((kk < 0) || (kk >= voxelSpace[0][0].length)) {
                                        continue;
                                    }

                                    v2 = voxelSpace[ii][jj][kk];

                                    if ((dv.x > v2.x - dlimit) && (dv.x < v2.x + dlimit)) {
                                        if ((dv.y > v2.y - dlimit) && (dv.y < v2.y + dlimit)) {
                                            if ((dv.z > v2.z - dlimit) && (dv.z < v2.z + dlimit)) {
                                                geometry.setSpace(ii, jj, kk, 1);
                                                //v2.isSpace = true;
                                            }
                                        }
                                    }

                                }
                            }
                        }

                        dv = dv.nextReady;

                    }

                }
            }
        }
        
    }

     private void countSpaceVoxels(CoordinatesVoxels c) {

        int i1 = c.xVoxel1;
        int i2 = c.xVoxel2;
        int j1 = c.yVoxel1;
        int j2 = c.yVoxel2;
        int k1 = c.zVoxel1;
        int k2 = c.zVoxel2;

        int spaceVoxels = 0;
        int totalVoxels = 0;

        for (int i = i1; i <= i2; i++) {
            for (int j = j1; j <= j2; j++) {
                for (int k = k1; k <= k2; k++) {
                    if (geometry.isSpace(i, j, k)) {
                        spaceVoxels++;
                    }
                    totalVoxels++;
                }
            }
        }

        Master.log("total voxels = " + totalVoxels);
        Master.log("space voxels = " + spaceVoxels);

     }

    private boolean voxelContainsVesicle(int i, int j, int k) {

        Double dlimit;

        if ((i < 0) || (i >= voxelSpace.length)) {
            return false;
        }

        if ((j < 0) || (j >= voxelSpace[0].length)) {
            return false;
        }

        if ((k < 0) || (k >= voxelSpace[0][0].length)) {
            return false;
        }

        Voxel voxel = voxelSpace[i][j][k];

        DiffusantVesicle dv = voxel.firstReady;

        if (dv != null) {
            return true;
        }

        for (int ii = i - 2; ii <= i + 2; ii++) {

            if ((ii < 0) || (ii >= voxelSpace.length)) {
                continue;
            }

            for (int jj = j - 2; jj <= j + 2; jj++) {

                if ((jj < 0) || (jj >= voxelSpace[0].length)) {
                    continue;
                }

                for (int kk = k - 2; kk <= k + 2; kk++) {

                    if ((kk < 0) || (kk >= voxelSpace[0][0].length)) {
                        continue;
                    }

                    dv = voxelSpace[ii][jj][kk].firstReady;

                    while (dv != null) {

                        dlimit = project.dx * 0.5 + dv.radius;

                        if ((dv.x > voxel.x - dlimit) && (dv.x < voxel.x + dlimit)) {
                            if ((dv.y > voxel.y - dlimit) && (dv.y < voxel.y + dlimit)) {
                                if ((dv.z > voxel.z - dlimit) && (dv.z < voxel.z + dlimit)) {
                                    return true;
                                }
                            }
                        }

                        dv = dv.nextReady;

                    }

                }

            }

        }

        return false;

    }

    private void removeIsolatedNonSpaceVoxels(int minNeighbors) {

        int counter, counter2 = 0;
        int i1, i2, j1, j2, k1, k2;

        i1 = 0;
        i2 = voxelSpace.length - 1;
        j1 = 0;
        j2 = voxelSpace[0].length - 1;
        k1 = 0;
        k2 = voxelSpace[0][0].length - 1;

        for (int i = i1 + 1; i <= i2 - 1; i++) {
            for (int j = j1 + 1; j <= j2 - 1; j++) {
                for (int k = k1 + 1; k <= k2 - 1; k++) {

                    if (geometry.space[i][j][k] == 1) {
                        continue;
                    }

                    counter = 0;

                    if (geometry.space[i-1][j][k] < 0) {
                        counter++;
                    }

                    if (geometry.space[i+1][j][k] < 0) {
                        counter++;
                    }

                    if (geometry.space[i][j-1][k] < 0) {
                        counter++;
                    }

                    if (geometry.space[i][j+1][k] < 0) {
                        counter++;
                    }

                    if (geometry.space[i][j][k-1] < 0) {
                        counter++;
                    }

                    if (geometry.space[i][j][k+1] < 0) {
                        counter++;
                    }

                    if (counter <= minNeighbors) {
                        geometry.setSpace(i, j, k, 1);
                        counter2++;
                    }

                }
            }
        }

         Master.log("removed isolated non-space voxels n = " + counter2);

    }

    private void removeSpaceIslands() {

        // this does not make difference in AZ supply rate

        int spaceVoxelsStart = 0;
        int spaceVoxelsEnd = 0;

        for (int k = 0; k < geometry.space[0][0].length; k++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int i = 0; i < geometry.space.length ; i++) {
                    if ((geometry.space[i][j][k] == 1)) {
                        spaceVoxelsStart++;
                    }
                }
            }
        }

        for (int k = 0; k < geometry.space[0][0].length; k++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int i = 0; i < geometry.space.length - 2; i++) {
                    if ((geometry.space[i][j][k] == 1) && (geometry.space[i+1][j][k] != 1) && (geometry.space[i+2][j][k] == 1)) {
                        geometry.setSpace(i+1, j, k, 1);
                    }
                }
            }
        }

        for (int k = 0; k < geometry.space[0][0].length; k++) {
            for (int j = 0; j < geometry.space[0].length - 2; j++) {
                for (int i = 0; i < geometry.space.length; i++) {
                    if ((geometry.space[i][j][k] == 1) && (geometry.space[i][j+1][k] != 1) && (geometry.space[i][j+2][k] == 1)) {
                        geometry.setSpace(i, j+1, k, 1);
                        //Master.log("" + i + ",");
                    }
                }
            }
        }

        for (int k = 0; k < geometry.space[0][0].length - 2; k++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int i = 0; i < geometry.space.length; i++) {
                    if ((geometry.space[i][j][k] == 1) && (geometry.space[i][j][k+1] != 1) && (geometry.space[i][j][k+2] == 1)) {
                        geometry.setSpace(i, j, k+1, 1);
                        //Master.log("" + i + "," + j + "," + k);
                    }
                }
            }
        }

        for (int k = 0; k < geometry.space[0][0].length; k++) {
            for (int j = 0; j < geometry.space[0].length - 1; j++) {
                for (int i = 0; i < geometry.space.length - 1; i++) {
                    if ((geometry.space[i][j][k] != 1) && (geometry.space[i+1][j+1][k] != 1) && (geometry.space[i][j+1][k] == 1) && (geometry.space[i+1][j][k] == 1)) {
                        geometry.setSpace(i, j, k, 1);
                        geometry.setSpace(i+1, j+1, k, 1);
                    }
                    if ((geometry.space[i][j][k] == 1) && (geometry.space[i+1][j+1][k] == 1) && (geometry.space[i][j+1][k] != 1) && (geometry.space[i+1][j][k] != 1)) {
                        geometry.setSpace(i, j+1, k, 1);
                        geometry.setSpace(i+1, j, k, 1);
                    }
                }
            }
        }

        for (int k = 0; k < geometry.space[0][0].length - 1; k++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int i = 0; i < geometry.space.length - 1; i++) {
                    if ((geometry.space[i][j][k] != 1) && (geometry.space[i+1][j][k+1] != 1) && (geometry.space[i][j][k+1] == 1) && (geometry.space[i+1][j][k] == 1)) {
                        geometry.setSpace(i, j, k, 1);
                        geometry.setSpace(i+1, j, k+1, 1);
                    }
                    if ((geometry.space[i][j][k] == 1) && (geometry.space[i+1][j][k+1] == 1) && (geometry.space[i][j][k+1] != 1) && (geometry.space[i+1][j][k] != 1)) {
                        geometry.setSpace(i, j, k+1, 1);
                        geometry.setSpace(i+1, j, k, 1);
                    }
                }
            }
        }

        for (int k = 0; k < geometry.space[0][0].length - 1; k++) {
            for (int j = 0; j < geometry.space[0].length - 1; j++) {
                for (int i = 0; i < geometry.space.length; i++) {
                    if ((geometry.space[i][j][k] != 1) && (geometry.space[i][j+1][k+1] != 1) && (geometry.space[i][j][k+1] == 1) && (geometry.space[i][j+1][k] == 1)) {
                        geometry.setSpace(i, j, k, 1);
                        geometry.setSpace(i, j+1, k+1, 1);
                    }
                    if ((geometry.space[i][j][k] == 1) && (geometry.space[i][j+1][k+1] == 1) && (geometry.space[i][j][k+1] != 1) && (geometry.space[i][j+1][k] != 1)) {
                        geometry.setSpace(i, j, k+1, 1);
                        geometry.setSpace(i, j+1, k, 1);
                    }
                }
            }

        }

        for (int k = 0; k < geometry.space[0][0].length; k++) {
            for (int j = 0; j < geometry.space[0].length; j++) {
                for (int i = 0; i < geometry.space.length ; i++) {
                    if ((geometry.space[i][j][k] == 1)) {
                        spaceVoxelsEnd++;
                    }
                }
            }
        }

        Master.log("% increase in space voxels = " + (100 * (spaceVoxelsEnd-spaceVoxelsStart)/spaceVoxelsStart));

    }

    private void openGeometry() {

        boolean foundSpace;

        int spaceVoxelsStart = 0;
        int spaceVoxelsEnd = 0;

        openGeometrySpace = new int[geometry.xVoxels][geometry.yVoxels][geometry.zVoxels];

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {
            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {
                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {
                    if ((geometry.space[i][j][k] == 1)) {
                        spaceVoxelsStart++;
                    }
                }
            }
        }

        if (openGeometryForwardX == 1) {

            for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

                for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                    foundSpace = false;

                    for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                        if (foundSpace) {
                            if (geometry.space[i][j][k] != 1) {
                                geometry.setSpace(i, j, k, 1);
                                openGeometrySpace[i][j][k] = 1;
                            }

                        } else if (geometry.space[i][j][k] == 1) {
                            foundSpace = true;
                        }

                    }
                }

            }

        } else if (openGeometryForwardX == -1) {

            for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

                for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                    foundSpace = false;

                    for (int i = EMcoordinates.xVoxel2; i >= EMcoordinates.xVoxel1; i--) {

                        if (foundSpace) {
                            if (geometry.space[i][j][k] != 1) {
                                geometry.setSpace(i, j, k, 1);
                                openGeometrySpace[i][j][k] = 1;
                            }

                        } else if (geometry.space[i][j][k] == 1) {
                            foundSpace = true;
                        }

                    }
                }

            }

        }

        if (openGeometryForwardY == 1) {

            for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                    foundSpace = false;

                    for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                        if (foundSpace) {
                            if (geometry.space[i][j][k] != 1) {
                                geometry.setSpace(i, j, k, 1);
                                openGeometrySpace[i][j][k] = 1;
                            }

                        } else if (geometry.space[i][j][k] == 1) {
                            foundSpace = true;
                        }

                    }
                }

            }

        } else if (openGeometryForwardY == -1) {

            for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                    foundSpace = false;

                    for (int j = EMcoordinates.yVoxel2; j >= EMcoordinates.yVoxel1; j--) {

                        if (foundSpace) {
                            if (geometry.space[i][j][k] != 1) {
                                geometry.setSpace(i, j, k, 1);
                                openGeometrySpace[i][j][k] = 1;
                            }

                        } else if (geometry.space[i][j][k] == 1) {
                            foundSpace = true;
                        }

                    }
                }

            }

        }

        if (openGeometryForwardZ == 1) {

            for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                    foundSpace = false;

                    for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

                        if (foundSpace) {
                            if (geometry.space[i][j][k] != 1) {
                                geometry.setSpace(i, j, k, 1);
                                openGeometrySpace[i][j][k] = 1;
                            }

                        } else if (geometry.space[i][j][k] == 1) {
                            foundSpace = true;
                        }

                    }
                }

            }

        } else if (openGeometryForwardZ == -1) {

            for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                    foundSpace = false;

                    for (int k = EMcoordinates.zVoxel2; k >= EMcoordinates.zVoxel1; k--) {

                        if (foundSpace) {
                            if (geometry.space[i][j][k] != 1) {
                                geometry.setSpace(i, j, k, 1);
                                openGeometrySpace[i][j][k] = 1;
                            }

                        } else if (geometry.space[i][j][k] == 1) {
                            foundSpace = true;
                        }

                    }
                }

            }

        }

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {
            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {
                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {
                    if ((geometry.space[i][j][k] == 1)) {
                        spaceVoxelsEnd++;
                    }
                }
            }
        }

        openGeometryNewVoxels = spaceVoxelsEnd-spaceVoxelsStart;

        Master.log("opened geometry, added space voxels = " + openGeometryNewVoxels);
        Master.log("non-space volume fraction = " + (1.0 * openGeometryNewVoxels / (spaceVoxelsEnd * 1.0)));

    }

    private void emptyGeometry() {

        boolean foundSpace;

        int spaceVoxelsStart = 0;
        int spaceVoxelsEnd = 0;

        openGeometrySpace = new int[geometry.xVoxels][geometry.yVoxels][geometry.zVoxels];

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {
            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {
                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {
                    if ((geometry.space[i][j][k] == 1)) {
                        spaceVoxelsStart++;
                    }
                }
            }
        }

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                foundSpace = false;

                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                    if (foundSpace) {
                        if (geometry.space[i][j][k] != 1) {
                            geometry.setSpace(i, j, k, 1);
                            openGeometrySpace[i][j][k] = 1;
                        }

                    } else if (geometry.space[i][j][k] == 1) {
                        foundSpace = true;
                    }

                }
            }

        }

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                foundSpace = false;

                for (int i = EMcoordinates.xVoxel2; i >= EMcoordinates.xVoxel1; i--) {

                    if (foundSpace) {
                        if (geometry.space[i][j][k] != 1) {
                            geometry.setSpace(i, j, k, 1);
                            openGeometrySpace[i][j][k] = 1;
                        }

                    } else if (geometry.space[i][j][k] == 1) {
                        foundSpace = true;
                    }

                }
            }

        }

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

            for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                foundSpace = false;

                for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                    if (foundSpace) {
                        if (geometry.space[i][j][k] != 1) {
                            geometry.setSpace(i, j, k, 1);
                            openGeometrySpace[i][j][k] = 1;
                        }

                    } else if (geometry.space[i][j][k] == 1) {
                        foundSpace = true;
                    }

                }
            }

        }

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

            for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

                foundSpace = false;

                for (int j = EMcoordinates.yVoxel2; j >= EMcoordinates.yVoxel1; j--) {

                    if (foundSpace) {
                        if (geometry.space[i][j][k] != 1) {
                            geometry.setSpace(i, j, k, 1);
                            openGeometrySpace[i][j][k] = 1;
                        }

                    } else if (geometry.space[i][j][k] == 1) {
                        foundSpace = true;
                    }

                }
            }

        }

        for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                foundSpace = false;

                for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {

                    if (foundSpace) {
                        if (geometry.space[i][j][k] != 1) {
                            geometry.setSpace(i, j, k, 1);
                            openGeometrySpace[i][j][k] = 1;
                        }

                    } else if (geometry.space[i][j][k] == 1) {
                        foundSpace = true;
                    }

                }
            }

        }

        for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {

            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {

                foundSpace = false;

                for (int k = EMcoordinates.zVoxel2; k >= EMcoordinates.zVoxel1; k--) {

                    if (foundSpace) {
                        if (geometry.space[i][j][k] != 1) {
                            geometry.setSpace(i, j, k, 1);
                            openGeometrySpace[i][j][k] = 1;
                        }

                    } else if (geometry.space[i][j][k] == 1) {
                        foundSpace = true;
                    }

                }
            }

        }

        for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {
            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {
                for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {
                    if ((geometry.space[i][j][k] == 1)) {
                        spaceVoxelsEnd++;
                    }
                }
            }
        }

        openGeometryNewVoxels = spaceVoxelsEnd-spaceVoxelsStart;

        geometry.setSpace(1);

        Master.log("opened geometry, added space voxels = " + openGeometryNewVoxels);
        Master.log("non-space volume fraction = " + (1.0 * openGeometryNewVoxels / (spaceVoxelsEnd * 1.0)));

    }

    private boolean fixAZvoxels() {

        int i0, j0, k0;
        int iextra = 6;
        int jextra = 6;
        int counter = 0;

        boolean iback = false;
        boolean iforward = false;
        boolean jback = false;
        boolean jforward = false;

        switch(EMseries) {

            case 1:
                iback = true;
                jback = true;
                break;

            case 3:
                iforward = true;
                jback = true;
                break;

            case 4:

                switch (EMseriesAZ) {
                    case 0:
                        iback = true;
                        iextra = 2;
                        jforward = true;
                        jextra = 10;
                        break;
                    case 1:
                        iback = true;
                        jback = true;
                        jextra = 10;
                        break;
                    case 2:
                        iback = true;
                        jforward = true;
                        jextra = 10;
                        break;
                    case 3:
                        iforward = true;
                        iback = true;
                        iextra = 8;
                        jforward = true;
                        break;
                }

                break;

            case 5:
                jforward = true;
                break;

            case 10:

                switch (EMseriesAZ) {
                    case 0:
                        iback = true;
                        iextra = 1;
                        break;
                    case 1:
                        jback = true;
                        jextra = 2;
                        break;
                    case 2:
                        iforward = true;
                        iextra = 2;
                        break;
                    case 3:
                        iback = true;
                        iextra = 3;
                        jback = true;
                        jextra = 3;
                        break;
                    case 4:
                        iback = true;
                        iextra = 3;
                        break;
                }

                break;

            case 16:
                iforward = true;
                iextra = 3;
                jforward = true;
                jextra = 3;
                break;

            case 18:
                iforward = true;
                iextra = 3;
                jforward = true;
                jextra = 3;
                break;

        }

        for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {

            if (Double.isNaN(azXYZ[EMseriesAZ][ipnt][0] * azXYZ[EMseriesAZ][ipnt][1] * azXYZ[EMseriesAZ][ipnt][2])) {
                break;
            }

            i0 = (int) geometry.computeVoxelX(azXYZ[EMseriesAZ][ipnt][0]);
            j0 = (int) geometry.computeVoxelY(azXYZ[EMseriesAZ][ipnt][1]);
            k0 = (int) geometry.computeVoxelZ(azXYZ[EMseriesAZ][ipnt][2]);

            if (iback) {
                for (int i = i0 + 1; i <= i0 + iextra; i++) {
                    if (geometry.isSpace(i, j0, k0) && !voxelContainsVesicle(i, j0, k0)) {
                        geometry.setSpace(i, j0, k0, -1);
                        counter++;
                    }
                }
            }
            if (iforward) {
                for (int i = i0 - 1; i >= i0 - iextra; i--) {
                    if (geometry.isSpace(i, j0, k0) && !voxelContainsVesicle(i, j0, k0)) {
                        geometry.setSpace(i, j0, k0, -1);
                        counter++;
                    }
                }
            }
            if (jback) {
                for (int j = j0 + 1; j <= j0 + jextra; j++) {
                    if (geometry.isSpace(i0, j, k0) && !voxelContainsVesicle(i0, j, k0)) {
                        geometry.setSpace(i0, j, k0, -1);
                        counter++;
                    }
                }
            }
            if (jforward) {
                for (int j = j0 - 1; j >= j0 - jextra; j--) {
                    if (geometry.isSpace(i0, j, k0) && !voxelContainsVesicle(i0, j, k0)) {
                        geometry.setSpace(i0, j, k0, -1);
                        counter++;
                    }
                }
            }

        }

        if (!fixNonSpaceVoxels) {
            return true;
        }

        // use code below for running simulations but not computing density near AZ

        for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {

            if (Double.isNaN(azXYZ[EMseriesAZ][ipnt][0] * azXYZ[EMseriesAZ][ipnt][1] * azXYZ[EMseriesAZ][ipnt][2])) {
                break;
            }

            i0 = (int) geometry.computeVoxelX(azXYZ[EMseriesAZ][ipnt][0]);
            j0 = (int) geometry.computeVoxelY(azXYZ[EMseriesAZ][ipnt][1]);
            k0 = (int) geometry.computeVoxelZ(azXYZ[EMseriesAZ][ipnt][2]);

            if (geometry.space[i0][j0][k0] < 0) {
                Master.log("fixed AZ non-space at " + i0 + "," + j0 + "," + k0);
                geometry.setSpace(i0, j0, k0, 1);
            }

        }

        if (true && (EMseries == 10) && (EMseriesAZ == 0)) {

            for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {

                if (Double.isNaN(azXYZ[EMseriesAZ][ipnt][0] * azXYZ[EMseriesAZ][ipnt][1] * azXYZ[EMseriesAZ][ipnt][2])) {
                    break;
                }

                i0 = (int) geometry.computeVoxelX(azXYZ[EMseriesAZ][ipnt][0]);
                j0 = (int) geometry.computeVoxelY(azXYZ[EMseriesAZ][ipnt][1]);
                k0 = (int) geometry.computeVoxelZ(azXYZ[EMseriesAZ][ipnt][2]);

                for (int i = i0 - 2; i <= i0 + 0; i++) {
                    for (int j = j0 - 0; j <= j0 + 0; j++) {
                        for (int k = k0 - 0; k <= k0 + 2; k++) {
                            if (geometry.space[i][j][k] < 0) {
                                Master.log("fixed AZ non-space at " + i + "," + j + "," + k);
                                geometry.setSpace(i, j, k, 1);
                            }
                        }
                    }
                }

            }

        }

        return true;

    }

    private double[] percentVoxelSpace(boolean print) {

        double percentSpace;
        double vx, vy, vz, dmin;

        DiffusantVesicle dv;

        double binWidth = 0.05; // um
        double maxDistance = 0.5; // um

        int numBins = 1 + (int) (maxDistance / binWidth);
        int ibin;

        double[] sum = new double[numBins];
        double[] count = new double[numBins];
        double[] percentSpaceVoxels = new double[numBins];

        for (int i = 0; i < sum.length; i++) {
            sum[i] = 0;
            count[i] = 0;
        }

        for (int i = 0; i < voxelSpace.length; i++) {
            for (int j = 0; j < voxelSpace[0].length; j++) {
                for (int k = 0; k < voxelSpace[0][0].length; k++) {

                    vx = voxelSpace[i][j][k].x;
                    vy = voxelSpace[i][j][k].y;
                    vz = voxelSpace[i][j][k].z;

                    if (!EMcoordinates.isInside(vx, vy, vz)) {
                        continue;
                    }

                    dv = voxelSpace[i][j][k].firstReady;

                    if (dv == null) {
                        continue;
                    }

                    dmin = minDistanceToAZ(azXYZ, EMseriesAZ, vx, vy, vz);

                    if (Double.isNaN(dmin) || Double.isInfinite(dmin)) {
                        continue;
                    }

                    ibin = (int) Math.floor(dmin / binWidth);
                    
                    percentSpace = 0;

                    for (int ii = i - 1; ii <= i + 1; ii++) {
                        for (int jj = j - 1; jj <= j + 1; jj++) {
                            for (int kk = k - 1; kk <= k + 1; kk++) {

                                if ((ii == i) && (jj == j) && (kk == k)) {
                                    continue;
                                }

                                if (!geometry.checkBounds(ii, jj, kk)) {
                                    continue;
                                }

                                if (voxelSpace[ii][jj][kk].isSpace) {
                                    percentSpace += 1;
                                }

                            }
                        }
                    }

                    percentSpace /= 26.0;

                    if ((ibin >= 0) && (ibin < numBins)) {
                        sum[ibin] += percentSpace;
                        count[ibin] += 1;
                    }

                }
            }

        }

        for (ibin = 0; ibin < numBins; ibin++) {

            percentSpaceVoxels[ibin] = sum[ibin] / count[ibin];

            if (print) {
                Master.log("" + percentSpaceVoxels[ibin]);
            }
            
        }

        return percentSpaceVoxels;

    }

    private void vesicleDensity1() {

        int extra;
        double density1;
        double ax, ay, az, d, dmin;

        //int numBands = 5;
        //int bandVoxels = 2; // 100 nm

        int numBands = 10;
        int bandVoxels = 1; // 50 nm

        DiffusantVesicle dv;

        int i1 = (int) geometry.computeVoxelX(azx1[EMseriesAZ]);
        int i2 = (int) geometry.computeVoxelX(azx2[EMseriesAZ]);
        int j1 = (int) geometry.computeVoxelY(azy1[EMseriesAZ]);
        int j2 = (int) geometry.computeVoxelY(azy2[EMseriesAZ]);
        int k1 = (int) geometry.computeVoxelZ(azz1[EMseriesAZ]);
        int k2 = (int) geometry.computeVoxelZ(azz2[EMseriesAZ]);

        Master.log("" + i1 + "," + j1 + "," + k1);
        Master.log("" + i2 + "," + j2 + "," + k2);

        double[] volume = new double[numBands];
        int[] nVesicles = new int[numBands];
        double[] minDistance = new double[numBands];

        boolean[][][] g = new boolean[voxelSpace.length][voxelSpace[0].length][voxelSpace[0][0].length];
        
        for (boolean[][] i : g) {
            for (boolean[] j : i) {
                for (boolean ijk : j) {
                    ijk = false;
                }
            }
        }

        for (int iband = 0; iband < numBands; iband++) {

            volume[iband] = 0;
            nVesicles[iband] = 0;
            minDistance[iband] = 0;

            extra = iband * bandVoxels;

            for (int i = i1 - extra; i <= i2 + extra; i++) {

                if ((i < 0) || (i >= voxelSpace.length)) {
                    continue;
                }

                for (int j = j1 - extra; j <= j2 + extra; j++) {

                    if ((j < 0) || (j >= voxelSpace[0].length)) {
                        continue;
                    }

                    for (int k = k1 - extra; k <= k2 + extra; k++) {

                        if ((k < 0) || (k >= voxelSpace[0][0].length)) {
                            continue;
                        }

                        if (g[i][j][k]) {
                            continue; // already tested
                        }

                        if (voxelSpace[i][j][k].isSpace) {

                            volume[iband] += voxelVolume;

                            dv = voxelSpace[i][j][k].firstReady;

                            while (dv != null) {

                                dmin = Double.POSITIVE_INFINITY;

                                for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {

                                    if (Double.isNaN(azXYZ[EMseriesAZ][ipnt][0])) {
                                        continue;
                                    }

                                    ax = azXYZ[EMseriesAZ][ipnt][0];
                                    ay = azXYZ[EMseriesAZ][ipnt][1];
                                    az = azXYZ[EMseriesAZ][ipnt][2];

                                    d = Math.sqrt((dv.x - ax) * (dv.x - ax) + (dv.y - ay) * (dv.y - ay) + (dv.z - az) * (dv.z - az));
                                    dmin = Math.min(dmin, d);

                                }

                                nVesicles[iband]++;
                                minDistance[iband] += dmin;

                                dv = dv.nextReady;

                            }

                        }

                        g[i][j][k] = true;

                    }
                }
            }

        }

        double nves = 0;

        for (int iband = 0; iband < volume.length; iband++) {
            //Master.log("" + volume[iband] + "," + nVesicles[iband]);
            //Master.log("" + nVesicles[iband]);
        }

        for (int iband = 0; iband < nVesicles.length; iband++) {
            density1 = nVesicles[iband] / volume[iband];
            //density2 = 4.403E-05 * density1;
            if (nVesicles[iband] > 0) {
                minDistance[iband] /= nVesicles[iband];
            } else {
                minDistance[iband] = 0;
            }
            //Master.log("" + density1 + "," + density2);
            //Master.log("" + density1);
            Master.log("" + minDistance[iband]);
            nves += nVesicles[iband];
        }

        Master.log("# ves " + nves);

    }

    @Override
    public boolean saveVesicleDensity(String saveTag) {
        vesicleDensity2(true, true, saveTag);
        return false;
    }
    
    public void vesicleDensityFinal(boolean print, boolean save) {
        
        String saveTag = "finish";

        double newDX = project.dx * 0.5;

        newVoxelSize(newDX); // for EM density analysis
        initVoxels();
        initVesicles();

        //if (fixAZvoxels()) {
            //initVoxels();
            //initVesicles();
        //}

        vesicleDensity2(print, save, saveTag);

    }

    private double[] vesicleDensity2(boolean print, boolean save, String saveTag) {

        double density1, density2;
        double dmin;

        double binWidth = project.dx; // USE THIS
        //double binWidth = 0.5 * project.dx; // um
        double maxDistance = 0.59; // um

        int numBins = 1 + (int) (maxDistance / binWidth);
        int ibin, numVes, numVox = 0;

        DiffusantVesicle dv;

        double[] volume = new double[numBins];
        double[] density = new double[numBins];
        int[] nVesicles = new int[numBins];

        for (ibin = 0; ibin < numBins; ibin++) {
            volume[ibin] = 0;
            density[ibin] = 0;
            nVesicles[ibin] = 0;
        }

        for (Voxel[][] i : voxelSpace) {
            for (Voxel[] j : i ) {
                for (Voxel v : j) {

                    if (!v.isSpace) {
                        continue;
                    }

                    if (!EMcoordinates.isInside(v.x, v.y, v.z)) {
                        continue;
                    }

                    numVes = 0;

                    dv = v.firstReady;

                    while (dv != null) {
                        numVes++;
                        dv = dv.nextReady;
                    }

                    dmin = minDistanceToAZ(azXYZ, EMseriesAZ, v.x, v.y, v.z);

                    if (Double.isNaN(dmin) || Double.isInfinite(dmin)) {
                        continue;
                    }

                    //dmin -= radius;

                    ibin = (int) Math.floor(dmin / binWidth);

                    if ((ibin >= 0) && (ibin < numBins)) {
                        volume[ibin] += voxelVolume;
                        nVesicles[ibin] += numVes;
                    } else if (dmin < maxDistance) {
                        Master.log("bin out of range: " + ibin + ", dmin = " + dmin);
                    }

                    if (ibin == 0) {
                        numVox++;
                    }

                }
            }

        }

        //Master.log("bin #0 voxels = " + numVox);

        int i1 = (int) geometry.computeVoxelX(azx1[EMseriesAZ]);
        int i2 = (int) geometry.computeVoxelX(azx2[EMseriesAZ]);
        int j1 = (int) geometry.computeVoxelY(azy1[EMseriesAZ]);
        int j2 = (int) geometry.computeVoxelY(azy2[EMseriesAZ]);
        int k1 = (int) geometry.computeVoxelZ(azz1[EMseriesAZ]);
        int k2 = (int) geometry.computeVoxelZ(azz2[EMseriesAZ]);

        double icounter1 = 0;
        double icounter2 = 0;

        for (int i = i1; i <= i2; i++) {
            for (int j = j1; j <= j2; j++) {
                for (int k = k1; k <= k2; k++) {
                    icounter2++;
                    if (voxelSpace[i][j][k].isSpace) {
                        icounter1++;
                    }
                }
            }
        }

        //volume[0] -= azVolume[azSelect] * icounter1 / icounter2;
        //Master.log("% AZ included = " + (icounter1 / icounter2) + ", " + icounter1 + "/" + icounter2);

        for (ibin = 0; ibin < numBins; ibin++) {
            density1 = nVesicles[ibin] / volume[ibin];
            density[ibin] = density1;
            density2 = vesicleVolume * density1;
            if (print) {
                //Master.log("" + density1);
                Master.log("" + density2);
                //Master.log("" + nVesicles[ibin] + ", " + volume[ibin]);
            }
            //Master.log("" + volume[ibin]);
            //Master.log("" + nVesicles[ibin]);
            //nves += nVesicles[ibin];
        }

        //Master.log("# vesicles " + nves);

        if (save) {
            saveArray(density, "VesicleDensity", saveTag);
        }

        return density;

    }

    private static void vesicleDensity3(double[][][] azXYZ, int azSelect, double[][] vesicleXYZ, Geometry g, double radius, double newDX, double binWidth) {

        int ibin;
        double x, y, z;
        double dx, dy, dz, d, dmin;
        double ivoxel, jvoxel, kvoxel;
        double spaceVoxels = 0, totalVoxels = 0, vesicleVoxels = 0;

        double maxDistanceFromAZ = 0.5; // um
        int numBins = 1 + (int) (maxDistanceFromAZ / binWidth);

        double[] vesicle = new double[numBins];
        int[] total = new int[numBins];

        for (ibin = 0; ibin < numBins; ibin++) {
            vesicle[ibin] = 0;
            total[ibin] = 0;
        }

        int xVoxels = (int) (g.xWidth / newDX);
        int yVoxels = (int) (g.yWidth / newDX);
        int zVoxels = (int) (g.zWidth / newDX);

        int xVoxel1 = 0;
        int yVoxel1 = 0;
        int zVoxel1 = 0;

        double xVoxelCenter = xVoxel1 + (xVoxels - 1) / 2.0;
        double yVoxelCenter = yVoxel1 + (yVoxels - 1) / 2.0;
        double zVoxelCenter = zVoxel1 + (zVoxels - 1) / 2.0;

        int[][][] space = new int[xVoxels][yVoxels][zVoxels];
        double[] x0 = new double[xVoxels];
        double[] y0 = new double[yVoxels];
        double[] z0 = new double[zVoxels];

        for (int i = 0; i < g.xVoxels; i++) {
            for (int j = 0; j < g.yVoxels; j++) {
                for (int k = 0; k < g.zVoxels; k++) {
                    if (g.isSpace(i, j, k)) {
                        spaceVoxels++;
                    }
                    totalVoxels++;
                }
            }
        }

        Master.log("old percent space voxels = " + (spaceVoxels / totalVoxels));

        spaceVoxels = 0;
        totalVoxels = 0;

        for (int i = 0; i < xVoxels; i++) {
            x0[i] = (i - xVoxelCenter) * newDX;
            ivoxel = g.computeVoxelX(x0[i]);
            for (int j = 0; j < yVoxels; j++) {
                y0[j] = (j - yVoxelCenter) * newDX;
                jvoxel = g.computeVoxelY(y0[j]);
                for (int k = 0; k < zVoxels; k++) {
                    z0[k] = (k - zVoxelCenter) * newDX;
                    kvoxel = g.computeVoxelZ(z0[k]);
                    if (g.isSpace(ivoxel, jvoxel, kvoxel)) {
                        space[i][j][k] = 1;
                        spaceVoxels++;
                    } else {
                        space[i][j][k] = 0;
                    }
                    totalVoxels++;
                }
            }
        }

        Master.log("new percent space voxels = " + (spaceVoxels / totalVoxels));

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {
                    for (int iv = 0; iv < vesicleXYZ.length; iv++) {

                        dx = vesicleXYZ[iv][0] - x0[i];
                        dy = vesicleXYZ[iv][1] - y0[j];
                        dz = vesicleXYZ[iv][2] - z0[k];

                        d = Math.sqrt(dx * dx + dy * dy + dz * dz);

                        if (d <= radius) {
                            space[i][j][k] = 2;
                            vesicleVoxels += 1;
                            break;
                        }

                    }
                }
            }
        }

        Master.log("percent vesicles = " + (vesicleVoxels / spaceVoxels));

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    if (space[i][j][k] == 0) {
                        continue;
                    }

                    dmin = minDistanceToAZ(azXYZ, azSelect, x0[i], y0[j], z0[k]);

                    if (Double.isNaN(dmin) || Double.isInfinite(dmin)) {
                        continue;
                    }

                    ibin = (int) Math.floor(dmin / binWidth);

                    if ((ibin >= 0) && (ibin < numBins)) {
                        total[ibin]++;
                        if (space[i][j][k] == 2) {
                            vesicle[ibin]++;
                        }
                    }

                }
            }
        }

        for (ibin = 0; ibin < numBins; ibin++) {
            if (total[ibin] > 0) {
                vesicle[ibin] /= total[ibin];
            } else {
                vesicle[ibin] = 0;
            }
            Master.log("" + vesicle[ibin]);
        }

    }

    private static double[] vesicleDensity3(Geometry geometry, double[][][] azXYZ, int azSelect, double[][] vxyz, double radius, double newDX, double binWidth, boolean print, boolean save, String saveTag) {

        int ibin;
        double x, y, z;
        double dx, dy, dz, d, dmin;
        double ivoxel, jvoxel, kvoxel;
        double spaceVoxels = 0, totalVoxels = 0, vesicleVoxels = 0;

        //double[][] vesicleXYZ = diffusant[0].xyz;
        //double radius = diffusant[0].setMeanRadius;

        double maxDistanceFromAZ = 0.5; // um
        int numBins = 1 + (int) (maxDistanceFromAZ / binWidth);

        double[] vesicle = new double[numBins];
        double[] density = new double[numBins];
        int[] total = new int[numBins];

        for (ibin = 0; ibin < numBins; ibin++) {
            vesicle[ibin] = 0;
            total[ibin] = 0;
        }

        int xVoxels = (int) (geometry.xWidth / newDX);
        int yVoxels = (int) (geometry.yWidth / newDX);
        int zVoxels = (int) (geometry.zWidth / newDX);

        int xVoxel1 = 0;
        int yVoxel1 = 0;
        int zVoxel1 = 0;

        double xVoxelCenter = xVoxel1 + (xVoxels - 1) / 2.0;
        double yVoxelCenter = yVoxel1 + (yVoxels - 1) / 2.0;
        double zVoxelCenter = zVoxel1 + (zVoxels - 1) / 2.0;

        int[][][] space = new int[xVoxels][yVoxels][zVoxels];
        double[] x0 = new double[xVoxels];
        double[] y0 = new double[yVoxels];
        double[] z0 = new double[zVoxels];

        for (int i = geometry.xVoxel1; i <= geometry.xVoxel2; i++) {
            for (int j = geometry.yVoxel1; j <= geometry.yVoxel2; j++) {
                for (int k = geometry.zVoxel1; k <= geometry.zVoxel2; k++) {
                    if (geometry.isSpace(i, j, k)) {
                        spaceVoxels++;
                    }
                    totalVoxels++;
                }
            }
        }

        Master.log("old percent space voxels = " + (spaceVoxels / totalVoxels));

        spaceVoxels = 0;
        totalVoxels = 0;

        for (int i = 0; i < xVoxels; i++) {
            x0[i] = (i - xVoxelCenter) * newDX;
            ivoxel = geometry.computeVoxelX(x0[i]);
            for (int j = 0; j < yVoxels; j++) {
                y0[j] = (j - yVoxelCenter) * newDX;
                jvoxel = geometry.computeVoxelY(y0[j]);
                for (int k = 0; k < zVoxels; k++) {
                    z0[k] = (k - zVoxelCenter) * newDX;
                    kvoxel = geometry.computeVoxelZ(z0[k]);
                    if (geometry.isSpace(ivoxel, jvoxel, kvoxel)) {
                        space[i][j][k] = 1;
                        spaceVoxels++;
                    } else {
                        space[i][j][k] = 0;
                    }
                    totalVoxels++;
                }
            }
        }

        Master.log("new percent space voxels = " + (spaceVoxels / totalVoxels));

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {
                    for (int iv = 0; iv < vxyz.length; iv++) {

                        dx = vxyz[iv][0] - x0[i];
                        dy = vxyz[iv][1] - y0[j];
                        dz = vxyz[iv][2] - z0[k];

                        d = Math.sqrt(dx * dx + dy * dy + dz * dz);

                        if (d <= radius) {
                            space[i][j][k] = 2;
                            vesicleVoxels += 1;
                            break;
                        }

                    }
                }
            }
        }

        Master.log("percent vesicles = " + (vesicleVoxels / spaceVoxels));

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    if (space[i][j][k] == 0) {
                        continue;
                    }

                    dmin = minDistanceToAZ(azXYZ, azSelect, x0[i], y0[j], z0[k]);

                    ibin = (int) Math.floor(dmin / binWidth);

                    if ((ibin >= 0) && (ibin < numBins)) {
                        total[ibin]++;
                        if (space[i][j][k] == 2) {
                            vesicle[ibin]++;
                        }
                    }

                }
            }
        }

        for (ibin = 0; ibin < numBins; ibin++) {
            if (total[ibin] > 0) {
                density[ibin] = vesicle[ibin] / total[ibin];
            } else {
                density[ibin] = 0;
            }
            if (print) {
                Master.log("" + density[ibin]);
            }
        }

        if (save) {
            //EM_saveArray(density, "VesicleDensity", saveTag);
        }

        return density;

    }

    private double[] vesicleDensity3(double binWidth, boolean print, boolean save, String saveTag) {

        int ibin;
        double x, y, z;
        double dx, dy, dz, d, dmin;
        double ivoxel, jvoxel, kvoxel;
        double spaceVoxels = 0, totalVoxels = 0, vesicleVoxels = 0;

        DiffusantVesicleAZ vaz;

        if (EMcoordinates == null) {
            return null;
        }

        //double[][] vesicleXYZ = diffusant[0].xyz;
        double radius = diffusants[0].setMeanRadius;

        double maxDistanceFromAZ = 0.5; // um
        int numBins = 1 + (int) (maxDistanceFromAZ / binWidth);

        double[] vesicle = new double[numBins];
        double[] density = new double[numBins];
        int[] total = new int[numBins];

        for (ibin = 0; ibin < numBins; ibin++) {
            vesicle[ibin] = 0;
            total[ibin] = 0;
        }

        int increaseResolutionFactor = 10;

        double newDX = project.dx / (1.0 * increaseResolutionFactor);

        int xVoxels = EMcoordinates.xVoxels * increaseResolutionFactor;
        int yVoxels = EMcoordinates.yVoxels * increaseResolutionFactor;
        int zVoxels = EMcoordinates.zVoxels * increaseResolutionFactor;

        int xVoxel1 = 0;
        int yVoxel1 = 0;
        int zVoxel1 = 0;

        double xVoxelCenter = xVoxel1 + (xVoxels - 1) / 2.0;
        double yVoxelCenter = yVoxel1 + (yVoxels - 1) / 2.0;
        double zVoxelCenter = zVoxel1 + (zVoxels - 1) / 2.0;

        int[][][] space = new int[xVoxels][yVoxels][zVoxels];
        double[] x0 = new double[xVoxels];
        double[] y0 = new double[yVoxels];
        double[] z0 = new double[zVoxels];

        for (int i = EMcoordinates.xVoxel1; i <= EMcoordinates.xVoxel2; i++) {
            for (int j = EMcoordinates.yVoxel1; j <= EMcoordinates.yVoxel2; j++) {
                for (int k = EMcoordinates.zVoxel1; k <= EMcoordinates.zVoxel2; k++) {
                    if (geometry.isSpace(i, j, k)) {
                        spaceVoxels++;
                    }
                    totalVoxels++;
                }
            }
        }

        Master.log("old percent space voxels = " + (spaceVoxels / totalVoxels));

        spaceVoxels = 0;
        totalVoxels = 0;

        double minX = 9999, maxX = -9999;

        for (int i = 0; i < xVoxels; i++) {
            x0[i] = (i - xVoxelCenter) * newDX;
            ivoxel = geometry.computeVoxelX(x0[i]);
            minX = Math.min(x0[i], minX);
            maxX = Math.max(x0[i], maxX);
            for (int j = 0; j < yVoxels; j++) {
                y0[j] = (j - yVoxelCenter) * newDX;
                jvoxel = geometry.computeVoxelY(y0[j]);
                for (int k = 0; k < zVoxels; k++) {
                    z0[k] = (k - zVoxelCenter) * newDX;
                    kvoxel = geometry.computeVoxelZ(z0[k]);
                    if (geometry.isSpace(ivoxel, jvoxel, kvoxel)) {
                        space[i][j][k] = 1;
                        spaceVoxels++;
                    } else {
                        space[i][j][k] = 0;
                    }
                    totalVoxels++;
                }
            }
        }

        Master.log("new percent space voxels = " + (spaceVoxels / totalVoxels));

        Master.log("x min = " + minX);
        Master.log("x max = " + maxX);

        int counter = 0;

        for (DiffusantVesicle v : diffusants[0].vesicles) {

            if (!(v instanceof DiffusantVesicleAZ)) {
                continue;
            }

            vaz = (DiffusantVesicleAZ) v;

            if (EMcoordinates.isInside(vaz.x, vaz.y, vaz.z)) {
                counter++;
            }

        }

        int[] vindex = new int[counter];

        counter = 0;

        for (int iv = 0; iv < diffusants[0].vesicles.length; iv++) {

            vaz = (DiffusantVesicleAZ) diffusants[0].vesicles[iv];

            if (EMcoordinates.isInside(vaz.x, vaz.y, vaz.z)) {
                vindex[counter] = iv;
                counter++;
            }

        }

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {
                    for (int iv = 0; iv < vindex.length; iv++) {

                        counter = vindex[iv];

                        vaz = (DiffusantVesicleAZ) diffusants[0].vesicles[counter];

                        dx = vaz.x - x0[i];
                        dy = vaz.y - y0[j];
                        dz = vaz.z - z0[k];

                        d = Math.sqrt(dx * dx + dy * dy + dz * dz);

                        if (d <= radius) {
                            space[i][j][k] = 2;
                            vesicleVoxels += 1;
                            break;
                        }

                    }
                }
            }
        }

        Master.log("percent vesicles = " + (vesicleVoxels / spaceVoxels));

        for (int i = 0; i < xVoxels; i++) {
            for (int j = 0; j < yVoxels; j++) {
                for (int k = 0; k < zVoxels; k++) {

                    if (space[i][j][k] == 0) {
                        continue;
                    }

                    dmin = minDistanceToAZ(azXYZ, EMseriesAZ, x0[i], y0[j], z0[k]);

                    ibin = (int) Math.floor(dmin / binWidth);

                    if ((ibin >= 0) && (ibin < numBins)) {
                        total[ibin]++;
                        if (space[i][j][k] == 2) {
                            vesicle[ibin]++;
                        }
                    }

                }
            }
        }

        for (ibin = 0; ibin < numBins; ibin++) {
            if (total[ibin] > 0) {
                density[ibin] = vesicle[ibin] / total[ibin];
            } else {
                density[ibin] = 0;
            }
            if (print) {
                Master.log("" + density[ibin]);
            }
        }

        if (save) {
            saveArray(density, "VesicleDensity", saveTag);
        }

        return density;

    }

    @Override
    public void saveVoxelResidenceTime() {

        //double d2AZlimit = 0.022;

        double shiftx = 1; // for reserve
        double extra = project.dx;

        residenceTimeCoordinates = new CoordinatesVoxels[2];

        residenceTimeCoordinates[0] = new CoordinatesVoxels(project, activeZone.x1 - extra, activeZone.y1 - extra, activeZone.z1 - extra, activeZone.x2 + extra, activeZone.y2 + extra, activeZone.z2 + extra);

        extra *= 2;

        residenceTimeCoordinates[1] = new CoordinatesVoxels(project, activeZone.x1 - extra - shiftx, activeZone.y1 - extra, activeZone.z1 - extra, activeZone.x2 + extra - shiftx, activeZone.y2 + extra, activeZone.z2 + extra);

        voxelResidenceTime(residenceTimeCoordinates[0], true); // AZ
        voxelResidenceTime(residenceTimeCoordinates[1], false); // reserve

    }

    private double voxelResidenceTime(CoordinatesVoxels c, boolean az) {

        int numVoxels = 0;
        double sum = 0, count = 0;
        double t = 0, n = 0, tTotal = 0, nTotal = 0;
        double vx, vy, vz, dmin;

        if (c == null) {
            return Double.NaN;
        }

        for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

                    if (voxelSpace != null) {
                        if (voxelSpace[i][j][k].isSpace) {

                            vx = voxelSpace[i][j][k].x;
                            vy = voxelSpace[i][j][k].y;
                            vz = voxelSpace[i][j][k].z;
                            dmin = minDistanceToAZ(azXYZ, EMseriesAZ, vx, vy, vz);

                            t = 0;
                            n = 0;

                            if (!az || (az && (dmin >= 0) && (dmin < 0.05))) {
                                t = voxelSpace[i][j][k].sumResidenceTime;
                                n = voxelSpace[i][j][k].residenceCounter;
                                if (t > 0) {
                                    tTotal += t;
                                    nTotal += n;
                                    numVoxels++;
                                }
                            }

                        }
                    } else {
                        return Double.NaN;
                    }

                    if ((t > 0) && (n > 0)) {
                        sum += t / n;
                        count++;
                    }

                }
            }
        }

        if (az) {
            Master.log("AZ: avg voxel residence time = " + (sum / count) + " ms (" + numVoxels + " voxels, t = " + tTotal + ", n = " + nTotal + ")");
        } else {
            Master.log("Reserve: avg voxel residence time = " + (sum / count) + " ms (" + numVoxels + " voxels, t = " + tTotal + ", n = " + nTotal + ")");
        }

        return (sum / count);

    }

    private int testVesicleOverlapEM(){

        int overlaps, counter = 0;
        double dx, dy, dz, dis, minD = Double.POSITIVE_INFINITY;

        Voxel voxel, kvoxel;
        DiffusantVesicle kvesicle;

        if (diffusants == null) {
            return -1;
        }

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            for (DiffusantVesicle v : d.vesicles) {

                if (v == null) {
                    continue;
                }

                if (Double.isNaN(v.x)) {
                    continue;
                }

                voxel = v.voxel;

                if (voxel == null) {
                    continue;
                }

                for (int k = 0; k < voxel.numNeighbors; k++) {

                    kvoxel = voxel.neighbors[k];
                    kvesicle = kvoxel.firstReady;

                    while (kvesicle != null) {

                        if ((kvesicle != v) && (!Double.isNaN(kvesicle.x))) {

                            dx = v.x - kvesicle.x;
                            dy = v.y - kvesicle.y;
                            dz = v.z - kvesicle.z;

                            dis = Math.sqrt(dx * dx + dy * dy + dz * dz);

                            minD = Math.min(minD, dis);

                            if (dis < v.radius + kvesicle.radius) {
                                //d2AZ = EM_minDistanceToAZ(azXYZ, azEM_Select, dv.x, dv.y, dv.z);
                                //Master.log("" + d2AZ);
                                counter++;
                            }

                        }

                        kvesicle = kvesicle.nextReady;

                    }

                }

            }
        }

        overlaps = (counter/2);

        Master.log("min distance between vesicles " + minD);
        Master.log("overlaps n = " + overlaps);

        return overlaps;

    }

    private void removeShrinkage(){

        double shrinkfactor = 1.0 / 0.89;
        double x, y, z;

        DiffusantVesicleAZ vaz;

        if (diffusants == null) {
            return;
        }

        project.dx *= shrinkfactor;
        voxelVolume = project.dx * project.dx * project.dx;

        geometry.update();

        initVoxels();

        for (int iAZ = 0; iAZ < azXYZ.length; iAZ++) {

            azx1[iAZ] *= shrinkfactor;
            azy1[iAZ] *= shrinkfactor;
            azz1[iAZ] *= shrinkfactor;
            azx2[iAZ] *= shrinkfactor;
            azy2[iAZ] *= shrinkfactor;
            azz2[iAZ] *= shrinkfactor;

            for (int ipnt = 0; ipnt < azXYZ[0].length; ipnt++) {
                if (!Double.isNaN(azXYZ[iAZ][ipnt][0])) {
                    azXYZ[iAZ][ipnt][0] *= shrinkfactor;
                    azXYZ[iAZ][ipnt][1] *= shrinkfactor;
                    azXYZ[iAZ][ipnt][2] *= shrinkfactor;
                }
            }

        }

        extendAZ *= shrinkfactor;
        initActiveZone();

        //x1 = activeZone.x1 * shrinkfactor;
        //y1 = activeZone.y1 * shrinkfactor;
        //z1 = activeZone.z1 * shrinkfactor;
        //x2 = activeZone.x2 * shrinkfactor;
        //y2 = activeZone.y2 * shrinkfactor;
        //z2 = activeZone.z2 * shrinkfactor;

        //activeZone.setCoordinates(x1, y1, z1, x2, y2, z2);

        for (DiffusantVesicles d : diffusants) {

            if ((d == null) || (d.vesicles == null)) {
                continue;
            }

            d.setMeanRadius *= shrinkfactor;

            for (DiffusantVesicle v : d.vesicles) {
                
                if (!(v instanceof DiffusantVesicleAZ)) {
                    continue;
                }

                vaz = (DiffusantVesicleAZ) v;

                if (Double.isNaN(vaz.x)) {
                    continue;
                }

                vaz.radius *= shrinkfactor;

                x = vaz.x * shrinkfactor;
                y = vaz.y * shrinkfactor;
                z = vaz.z * shrinkfactor;

                removeFromVoxelList(vaz, vaz.voxel);
                setVesicleLocation(vaz, x, y, z, true);
                addToVoxelList(vaz);

            }

        }

    }

    private void newVoxelSize(double dx) {

        double x, y, z;
        int vx, vy, vz;

        double dxOld = project.dx;
        int xVoxelsOld = geometry.xVoxels;
        int yVoxelsOld = geometry.yVoxels;
        int zVoxelsOld = geometry.zVoxels;

        int xVoxels = (int) (geometry.xWidth / dx);
        int yVoxels = (int) (geometry.yWidth / dx);
        int zVoxels = (int) (geometry.zWidth / dx);

        int[][][] space = new int[xVoxelsOld][yVoxelsOld][zVoxelsOld];

        for (int i = 0; i < xVoxelsOld; i++) {
            for (int j = 0; j < yVoxelsOld; j++) {
                for (int k = 0; k < zVoxelsOld; k++) {
                    space[i][j][k] = geometry.space[i][j][k];
                }
            }
        }

        project.dx = dx;
        geometry.resizeWithSpace(xVoxels, yVoxels, zVoxels);

        for (int i = 0; i < geometry.space.length; i++) {

            x = geometry.computeX(i);
            vx = (int) ((x / dxOld) + (xVoxelsOld / 2.0));

            for (int j = 0; j < geometry.space[0].length; j++) {

                y = geometry.computeY(j);
                vy = (int) ((y / dxOld) + (yVoxelsOld / 2.0));

                for (int k = 0; k < geometry.space[0][0].length; k++) {

                    z = geometry.computeZ(k);
                    vz = (int) ((z / dxOld) + (zVoxelsOld / 2.0));

                    geometry.setSpace(i, j, k, space[vx][vy][vz]);

                }
            }
        }

    }

}
