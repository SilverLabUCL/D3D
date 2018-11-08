package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.*;

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
final public class GeometryBouton {

    public static CoordinatesVoxels bouton = null; // coordinates of bouton
    public static CoordinatesVoxels axonLeft = null; // coordinates of axon
    public static CoordinatesVoxels axonRight = null; // coordinates of axon

    private GeometryBouton() {
    }

    // create axonal bouton structure
    public static void create(Geometry geometry, double axon_width,
            double bouton_width, double bouton_length,
            double total_length, boolean enpassant) {

        boolean evenNumVoxels = true;

        int iBgn, iEnd, jkBgn, jkEnd, iMax, jkMax;
        int axon_width_voxels, axon_vRadius, ainc = 0, offset = 0;
        int bouton_length_voxels, bouton_slope_len;
        int xVoxel1, yVoxel1, zVoxel1, xVoxel2, yVoxel2, zVoxel2;

        double dx = geometry.project.dx;

        iMax = (int) (total_length / dx);
        jkMax = (int) (bouton_width / dx);

        if (evenNumVoxels) { // use even for quarter symmetry simulations
            iMax = Utility.makeEven(iMax);
            jkMax = Utility.makeEven(jkMax);
            offset = 1;
        }

        geometry.resizeWithSpace(iMax, jkMax, jkMax);
        geometry.fill(); // fill with non-space

        bouton_length_voxels = (int) Math.floor(bouton_length / dx);
        axon_width_voxels = (int) Math.floor(axon_width / dx);
        
        if (evenNumVoxels) {
            bouton_length_voxels = Utility.makeEven(bouton_length_voxels);
            axon_width_voxels = Utility.makeEven(axon_width_voxels);
        }

        axon_vRadius = axon_width_voxels / 2;

        if (enpassant) {
            xVoxel1 = (int) (geometry.xVoxelCenter - (bouton_length / (2*dx))) + offset;
            xVoxel2 = (int) (geometry.xVoxelCenter + (bouton_length / (2*dx)));
        } else {
            xVoxel1 = iMax - 1 - bouton_length_voxels + offset;
            xVoxel2 = iMax - 1;
        }

        yVoxel1 = geometry.yVoxel1;
        yVoxel2 = geometry.yVoxel2;

        zVoxel1 = geometry.zVoxel1;
        zVoxel2 = geometry.zVoxel2;

        bouton = new CoordinatesVoxels(geometry.project, xVoxel1, yVoxel1, zVoxel1, xVoxel2, yVoxel2, zVoxel2);

        bouton_slope_len = (jkMax - axon_width_voxels) / 2;

        // create axon on left side

        iBgn = geometry.xVoxel1;
        iEnd = bouton.xVoxel1;
        
        jkBgn = (int) (geometry.yVoxelCenter - axon_vRadius + offset);
        jkEnd = (int) (geometry.yVoxelCenter + axon_vRadius);

        axonLeft = new CoordinatesVoxels(geometry.project, iBgn, jkBgn, jkBgn, iEnd, jkEnd, jkEnd);
        geometry.setSpace(axonLeft, 1);

        // create left side of bouton (tapered)

        iBgn = bouton.xVoxel1;
        iEnd = bouton.xVoxel1 + bouton_slope_len - 1;

        for (int i = iBgn; i <= iEnd; i++, ainc++) {

            axon_vRadius = (axon_width_voxels / 2) + ainc;

            jkBgn = (int) (geometry.yVoxelCenter - axon_vRadius + offset);
            jkEnd = (int) (geometry.yVoxelCenter + axon_vRadius);

            for (int k = jkBgn; k <= jkEnd; k++) {
                for (int j = jkBgn; j <= jkEnd; j++) {
                    geometry.setSpace(i, j, k, 1);
                }
            }

        }

        // create bouton center (cube)

        if (enpassant) {
            iEnd = bouton.xVoxel2 - bouton_slope_len;
        } else {
            iEnd = bouton.xVoxel2;
        }

        iBgn = bouton.xVoxel1 + bouton_slope_len;

        jkBgn = bouton.yVoxel1;
        jkEnd = bouton.yVoxel2;

        CoordinatesVoxels c = new CoordinatesVoxels(geometry.project, iBgn, jkBgn, jkBgn, iEnd, jkEnd, jkEnd);
        geometry.setSpace(c, 1);

        axonRight = null;

        if (!enpassant) {
            return;
        }

        // create right side of bouton (tapered)

        iBgn = bouton.xVoxel2 - bouton_slope_len;
        iEnd = bouton.xVoxel2;

        for (int i = iBgn; i <= iEnd; i++, ainc--) {

            axon_vRadius = (axon_width_voxels / 2) + ainc;

            jkBgn = (int) (geometry.yVoxelCenter - axon_vRadius + offset);
            jkEnd = (int) (geometry.yVoxelCenter + axon_vRadius);

            for (int k = jkBgn; k <= jkEnd; k++) {
                for (int j = jkBgn; j <= jkEnd; j++) {
                    geometry.setSpace(i, j, k, 1);
                }
            }

        }

        // create right side of axon

        iBgn = bouton.xVoxel2 + 1;
        iEnd = iMax - 1;
        axon_vRadius = axon_width_voxels / 2;
        jkBgn = (int) (geometry.yVoxelCenter - axon_vRadius + offset);
        jkEnd = (int) (geometry.yVoxelCenter + axon_vRadius);

        axonRight = new CoordinatesVoxels(geometry.project, iBgn, jkBgn, jkBgn, iEnd, jkEnd, jkEnd);
        geometry.setSpace(axonRight, 1);

        geometry.checkSpace();

    }

    public static void addVesicles(Geometry geometry, int N, double vwidth, long seed) { // vesicles are cubes
        if (bouton != null) {
            GeometryTools.addCuboids(geometry, bouton, N, vwidth, vwidth, vwidth, seed); // add to bouton
        }
    }

}
