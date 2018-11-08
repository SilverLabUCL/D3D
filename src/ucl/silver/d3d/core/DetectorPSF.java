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
public class DetectorPSF extends Detector {

    public DetectorPSF(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PSF PSF) {
        super(p, NAME, DiffusantNum, c);
        psf = PSF;
        color = new ColorD3D(name + "_color", "Heat");
        save.dataPoints = 1;
        createVector(true);
    }

    @Override
    public boolean saveInit() {
        if (save == null) {
            return false;
        }

        int dataPoints = 1;

        return save.init(name, coordinates(), -1, dataPoints);
    }

    @Override
    public void detect(RunFiniteDifference fd, Geometry geometry) {

        int index;
        int xVoxel1, yVoxel1, zVoxel1, xVoxels, yVoxels, zVoxels;
        double psfValue, conc = 0, norm = 0;

        if (save == null) {
            return;
        }

        if (!psf.exists() || (fd.diffus == null)) {
            return;
        }

        if ((diffusantNum < 0) || (diffusantNum >= fd.diffus.length)) {
            return;
        }

        if (fd.diffus[0].length == 1) {
            save.saveData(fd.diffus[diffusantNum][0]);
            return; // single compartment
        }

        xVoxel1 = coordinates().xVoxel1;
        yVoxel1 = coordinates().yVoxel1;
        zVoxel1 = coordinates().zVoxel1;
        xVoxels = coordinates().xVoxels;
        yVoxels = coordinates().yVoxels;
        zVoxels = coordinates().zVoxels;

        for (int k = 0, kk = zVoxel1; k < zVoxels; k++, kk++) {
            for (int j = 0, jj = yVoxel1; j < yVoxels; j++, jj++) {
                for (int i = 0, ii = xVoxel1; i < xVoxels; i++, ii++) {

                    if (geometry.isSpace(ii, jj, kk)) {

                        index = geometry.space[ii][jj][kk];
                        psfValue = psf.getArrayValue(i, j, k);

                        norm += psfValue;
                        conc += fd.diffus[diffusantNum][index] * psfValue;

                    }

                }
            }
        }

        conc /= norm;

        save.saveData(conc);

    }

}
