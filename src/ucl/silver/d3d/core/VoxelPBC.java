package ucl.silver.d3d.core;

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
public class VoxelPBC extends Voxel {

    public int numPBCneighbors;

    public transient Voxel[] PBCneighbors; // up to 19

    public int[] PBCi; // translation indexes of PBC neighbors
    public int[] PBCj;
    public int[] PBCk;

    public VoxelPBC(){
        PBCneighbors = new Voxel[19];
        PBCi = new int[19];
        PBCj = new int[19];
        PBCk = new int[19];
    }

}
