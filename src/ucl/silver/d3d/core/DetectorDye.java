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
public class DetectorDye extends Detector {

    public String normType = "dF/F0"; // "dF/F0" or "(F-F0)/F0" or "F/F0" or "F/Fmax" or "F"
    // "dF/F0" and "(F-F0)/F0" are the same

    private transient DiffusantDye dye_ab; // [ab]
    private transient boolean usePSF = false;

    public DetectorDye(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PSF PSF) {
        super(p, NAME, DiffusantNum, c);
        psf = PSF;
        usePSF = true;
        color = new ColorD3D(name + "_color", "Heat");
        createVector(true);
    }

    public DetectorDye(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, PSF PSF, String NormType) {
        super(p, NAME, DiffusantNum, c);
        psf = PSF;
        usePSF = true;
        normType = NormType;
        color = new ColorD3D(name + "_color", "Heat");
        createVector(true);
    }

    @Override
    public void init() {

        super.init();

        Diffusant[] diffusant = project.diffusants;

        if (diffusant[diffusantNum] instanceof DiffusantDye) {
            dye_ab = (DiffusantDye) diffusant[diffusantNum];
        } else {
            error("init", "diffusantNum", "selected diffusant is not of DiffusantReactantDye class");
        }

        if (normType.equalsIgnoreCase("dF/F0")) {
            normType = "dF/F0";
        } else if (normType.equalsIgnoreCase("(F-F0)/F0")) {
            normType = "(F-F0)/F0";
        } else if (normType.equalsIgnoreCase("F/F0")) {
            normType = "F/F0";
        } else if (normType.equalsIgnoreCase("F/Fmax")) {
            normType = "F/Fmax";
        } else if (normType.equalsIgnoreCase("F")) {
            normType = "F"; // no normalization
        } else {
            error("init", "normType", "unknown F normalization type");
        }

    }

    @Override
    public boolean saveInit() {

        if (save == null) {
            return false;
        }
        
        save.dataPoints = 2; // [ab] and Fnorm

        return save.init(name, coordinates(), -1, save.dataPoints);
        
    }

    @Override
    public void saveDimensions() {

        if ((save == null) || (!save.autoDimensions)) {
            return;
        }

        save.xdim = project.timeUnits;

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            save.ydim = project.concUnits + ";" + normType + ";";
        } else {
            save.ydim = diffusantName + " (" + project.concUnits + ");" + diffusantName + " " + normType + ";";
        }

    }

    @Override
    public void detect(RunFiniteDifference fd, Geometry geometry) {

        double avg_ab = 0;
        double f;
        double psfValue, psfSum = 0;
        int index, numVoxels = 0;
        int xVoxel1, yVoxel1, zVoxel1, xVoxels, yVoxels, zVoxels;

        double[] values = new double[2];

        if (fd.diffus == null) {
            return;
        }

        if (fd.diffus[0].length == 1) { // single compartment

            avg_ab = fd.diffus[diffusantNum][0];

        } else if (usePSF) {

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

                            psfSum += psfValue;
                            avg_ab += fd.diffus[diffusantNum][index] * psfValue;

                        }

                    }
                }
            }

            avg_ab /= psfSum; // weighted avg

        } else if (coordinates().indexRFD != null) {

            numVoxels = coordinates().indexRFD.length;

            for (int i = 0; i < numVoxels; i++) {
                index = coordinates().indexRFD[i];
                avg_ab += fd.diffus[diffusantNum][index];
            }

            avg_ab /= 1.0 * numVoxels;

        } else {

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
                            avg_ab += fd.diffus[diffusantNum][index];
                            numVoxels++;
                        }

                    }
                }
            }

            avg_ab /= 1.0 * numVoxels;

        }

        f = dye_ab.computeF(avg_ab, normType);

        values[0] = avg_ab;
        values[1] = f;

        save.saveData(values);

    }

}
