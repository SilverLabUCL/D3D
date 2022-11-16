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
public class DetectorAvg extends Detector {
    
    public double Cavg; // mM

    public DetectorAvg(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c) {
        super(p, NAME, DiffusantNum, c);
        createVector(true);
    }
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("cavg")) {
            return project.concUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("cavg")) {
            return false;
        }
        return super.canEdit(name);
    }

    @Override
    public boolean saveInit() {

        if (save == null) {
            return false;
        }
        
        save.dataPoints = 1;

        return save.init(name, coordinates(), -1, save.dataPoints);

    }
    
    public boolean save() {
        return save.saveData(Cavg);
    }

    @Override
    public void detect(RunFiniteDifference fd, Geometry geometry) {
        avg(fd, geometry);
        save();
    }
    
    public void avg(RunFiniteDifference fd, Geometry geometry) {

        int index, numVoxels = 0;
        int xVoxel1, yVoxel1, zVoxel1, xVoxels, yVoxels, zVoxels;
        
        Cavg = Double.NaN;

        if (fd.diffus[0].length == 1) { // single compartment
            
            numVoxels = 1;
            Cavg = fd.diffus[diffusantNum][0];
            
        } else if (coordinates().indexRFD == null) {

            xVoxel1 = coordinates().xVoxel1;
            yVoxel1 = coordinates().yVoxel1;
            zVoxel1 = coordinates().zVoxel1;
            xVoxels = coordinates().xVoxels;
            yVoxels = coordinates().yVoxels;
            zVoxels = coordinates().zVoxels;
            
            Cavg = 0;

            for (int k = 0, kk = zVoxel1; k < zVoxels; k++, kk++) {
                for (int j = 0, jj = yVoxel1; j < yVoxels; j++, jj++) {
                    for (int i = 0, ii = xVoxel1; i < xVoxels; i++, ii++) {
                        if (geometry.isSpace(ii, jj, kk)) {
                            index = geometry.space[ii][jj][kk];
                            Cavg += fd.diffus[diffusantNum][index];
                            numVoxels++;
                        }

                    }
                }
            }

        } else {

            numVoxels = coordinates().indexRFD.length;
            Cavg = 0;

            for (int i : coordinates().indexRFD) {
                Cavg += fd.diffus[diffusantNum][i];
            }

        }

        if (numVoxels > 1) {
            Cavg /= 1.0 * numVoxels;
        }

    }

}
