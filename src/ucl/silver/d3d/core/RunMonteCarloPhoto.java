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
public class RunMonteCarloPhoto extends RunMonteCarlo {
    
    public boolean frapOn = false;
    public boolean saveFluorescence = false;
    
    public CoordinatesVoxels initFluorescenceCoordinates = null;
    public double initF_inside = 1;
    public double initF_outside = 0;
    
    public transient PSF PSFi = null; // illumination PSF
    public transient PSF PSFd = null; // detection PSF
    
    public transient DetectorPSF[] detectors = null;

    private int iPSF = -1;
    private int dPSF = -1;
    
    public RunMonteCarloPhoto(Project p) {

        super(p);

        createVector(true);

    }
    
    public void setPSFSelect(int IPSF, int DPSF) {
        iPSF = IPSF;
        dPSF = DPSF;
    }
    
    @Override
    public boolean initAll() {
        
        if (super.initAll()) {
            return true;
        }
        
        if (initPSFs()) {
            return true;
        }
        
        initFluorescence();
        
        countParticlesWithinPSF();

        if (PSFi != null) {
            PSFi.sum(geometry);
        }

        if (PSFd != null) {
            if (geometry.simpleCuboid) {
                PSFd.sum(geometry);
            } else {
                sumPSF(PSFd, geometry.voxelSpace);
            }
        }

        meanFluorescence(true);
        
        Master.updatePanel2D();
        
        return false;
                
    }
    
    @Override
    public boolean initDetectors() {

        int count = 0;

        if (project.detectors == null) {
            return false; // nothing to do
        }

        for (Detector d : project.detectors) {
            if (d instanceof DetectorPSF) {
                count++;
            }
        }

        if (count == 0) {
            return false; // nothing to do
        }

        detectors = new DetectorPSF[project.detectors.length];

        for (int i = 0; i < project.detectors.length; i++) {

            if (project.detectors[i] instanceof DetectorPSF) {
                detectors[i] = (DetectorPSF) project.detectors[i];
            } else {
                detectors[i] = null;
            }

        }

        return false;

    }
    
    @Override
    public boolean initVoxels() {

        if (PBC) {
            geometry.initVoxelsPBC(iPSF, dPSF);
        } else {
            geometry.initVoxels(iPSF, dPSF);
        }

        voxelVolume = project.dx * project.dx * project.dx;
        
        return false;
        
    }
    
    public void sumPSF(PSF psf, Voxel[][][] voxels) {

        double avg = 0, sum = 0, count = 0;

        if (voxels == null) {
            return;
        }

        if (psf.array == null) {
            return;
        }

        if (psf.array.length != voxels.length) {
            return;
        }

        if (psf.array[0].length != voxels[0].length) {
            return;
        }

        if (psf.array[0][0].length != voxels[0][0].length) {
            return;
        }

        for (int i = 0; i < psf.array.length; i++) {
            for (int j = 0; j < psf.array[0].length; j++) {
                for (int k = 0; k < psf.array[0][0].length; k++) {
                    if (voxels[i][j][k].isSpace) {
                        sum += psf.array[i][j][k] * voxels[i][j][k].PSFweight;
                        count++;
                    }
                }
            }
        }

        if (count > 0) {
            avg = sum / count;
        }

        Master.log("weighted dPSF sum: " + sum);

        psf.sum = sum;
        psf.avg = avg;

        psf.setParamObject("sum", sum);
        psf.setParamObject("avg", avg);

    }

    public boolean initPSFs() {

        int count = 0;
        double e2 = 0.135335;

        if (!frapOn) {
            return false; // nothing to do
        }

        if (geometry.PSFs == null) {
            error("initPSFs", "geometry.PSFs", "no PSFs");
            return true;
        }

        if ((iPSF < 0) || (iPSF >= geometry.PSFs.length)) {
            error("initPSFs", "iPSF", "bad value");
            return true;
        }

        if ((dPSF < 0) || (dPSF >= geometry.PSFs.length)) {
            error("initPSFs", "dPSF", "bad value");
            return true;
        }

        PSFi = geometry.PSFs[iPSF];

        PSFi.useGeometryCoordinates = true;
        PSFi.array = null;
        PSFi.checkExists();

        PSFd = geometry.PSFs[dPSF];

        PSFd.useGeometryCoordinates = true;
        PSFd.array = null;
        PSFd.checkExists();
        
        for (double[][] i : PSFi.array) {
            for (double[] j : i) {
                for (double psf : j) {
                    if (psf > e2) {
                        count++;
                    }
                }
            }
            
        }

        Master.log("PSFi(voxel > e^-2) n = " + count);

        count = 0;
        
        for (double[][] i : PSFd.array) {
            for (double[] j : i) {
                for (double psf : j) {
                    if (psf > e2) {
                        count++;
                    }
                }
            }
            
        }

        Master.log("PSFc(voxel > e^-2) n = " + count);

        return false;

    }
    
    public void initFluorescence() {

        if ((!frapOn) || (diffusants == null)) {
            return;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null) || (geometry == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {
                if (initFluorescenceCoordinates != null) {
                    if (p.isinside(initFluorescenceCoordinates)) {
                        p.fluorescence = initF_inside;
                    } else {
                        p.fluorescence = initF_outside;
                    }
                } else if (p.isinside(geometry)) {
                    p.fluorescence = initF_inside;
                } else {
                    p.fluorescence = initF_outside;
                }
            }

        }

    }
    
    public double meanFluorescence(boolean printResults) {

        double w, avg = 0, count = 0;

        if ((!frapOn) || (diffusants == null)) {
            return 0;
        }

        for (DiffusantParticles d : diffusants) {

            if ((d == null) || (d.particles == null)) {
                continue;
            }

            for (DiffusantParticle p : d.particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                if (freeDiffusion || PBC) {
                    w = 1.0;
                } else {
                    w = p.voxel.PSFweight;
                }

                avg += p.fluorescence * p.voxel.PSFd * w;
                count++;

            }

        }

        if ((PSFd.avg > 0) && (count > 0)) {
            //avg /= (PSFd.avg * count);
            avg /= PSFd.avg;
        } else {
            avg = 0;
        }

        if (printResults) {
            Master.log("mean fluorescence: " + avg);
        }

        return avg;

    }
    
    public int countBleached() {

        int dnum, count = 0;
        double w = 1, f, e2 = 0.135335; // e^-2

        if (!frapOn || (detectors == null) || (dPSF == -1)) {
            return 0; // nothing to do
        }

        for (Detector d : detectors) {

            if (d == null) {
                continue;
            }

            dnum = d.diffusantNum;

            if ((dnum < 0) || (dnum >= diffusants.length)) {
                continue;
            }

            if (diffusants[dnum].particles == null) {
                continue;
            }

            for (DiffusantParticle p : diffusants[dnum].particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                //if (freeDiffusion || PBC) {
                //w = 1.0;
                //} else {
                //w = diffusant[dnum].particle[j].voxel.PSFweight; // extra weighting due to non-space voxels
                //}
                //f = diffusants[dnum].particles[j].fluorescence * diffusants[dnum].particles[j].voxel.PSFd * w;
                f = p.fluorescence;

                if (f < 1.0 - e2) {
                    count++;
                }

            }

        }

        Master.log("bleaching > 86% n = " + count);

        return count;

    }

    public int countParticlesWithinPSF() {

        int dnum;
        int count = 0;
        double e2 = 0.135335;

        if (!frapOn || (detectors == null) || (dPSF == -1)) {
            return 0; // nothing to do
        }

        for (Detector d : detectors) {

            if (d == null) {
                continue;
            }

            dnum = d.diffusantNum;

            if ((dnum < 0) || (dnum >= diffusants.length)) {
                continue;
            }

            if (diffusants[dnum].particles == null) {
                continue;
            }

            for (DiffusantParticle p : diffusants[dnum].particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                if (p.voxel.PSFd > e2) {
                    count++;
                }

            }

        }

        Master.log("particles within cPSF n = " + count);

        return count;

    }
    
    @Override
    public void react() {

        DiffusantParticlesPhoto dp;

        if (frapOn) {

            detectFluorescence();

            for (DiffusantParticles d : diffusants) {
                if (d instanceof DiffusantParticlesPhoto) {
                    dp = (DiffusantParticlesPhoto) d;
                    dp.react(itime);
                }
            }

        }

    }

    public void detectFluorescence() {

        int dnum;
        double w, avg = 0, count = 0;

        if (!frapOn || !saveFluorescence || (detectors == null) || (dPSF == -1)) {
            return; // nothing to do
        }

        for (Detector d : detectors) {

            if (d == null) {
                continue;
            }

            dnum = d.diffusantNum;

            if ((dnum < 0) || (dnum >= diffusants.length)) {
                continue;
            }

            if (diffusants[dnum].particles == null) {
                continue;
            }

            for (DiffusantParticle p : diffusants[dnum].particles) {

                if (!p.insideGeometry) {
                    continue;
                }

                if (freeDiffusion || PBC) {
                    w = 1.0;
                } else {
                    w = p.voxel.PSFweight; // extra weighting due to non-space voxels
                }

                avg += p.fluorescence * p.voxel.PSFd * w;
                count++;

            }

            if ((PSFd.avg > 0) && (count > 0)) {
                //avg /= (PSFd.avg * count);
                avg /= PSFd.avg;
            } else {
                avg = 0;
            }

            d.save.saveData(avg);

        }

    }
    
}
