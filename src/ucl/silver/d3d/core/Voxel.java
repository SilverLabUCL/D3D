package ucl.silver.d3d.core;

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
public class Voxel extends ParamVector {

    public double x, y, z; // center location in um
    public double PSFi = -1; // illumination PSF
    public double PSFd = -1; // detection PSF
    public double PSFweight = 1; // non-space neighbor weight (8/8 or 4/8 or 2/8 or 1/8 or 0 )

    public boolean isSpace = true;

    public int numNeighbors, numNonSpaceNeighbors;

    public double sumResidenceTime = 0;
    public int residenceCounter = 0;

    public transient Voxel[] neighbors; // 3 x 3 x 3 = 27
    public transient Voxel[] nonSpaceNeighbors; // up to 26

    public DiffusantVesicle firstReady = null;

    public Voxel(){
        super(null);
        neighbors = new Voxel[27];
        nonSpaceNeighbors = new Voxel[26];
    }

}
