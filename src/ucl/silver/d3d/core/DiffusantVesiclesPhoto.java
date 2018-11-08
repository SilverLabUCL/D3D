package ucl.silver.d3d.core;

import ucl.silver.d3d.utils.MersenneTwisterFast;

import java.awt.Color;

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
public class DiffusantVesiclesPhoto extends DiffusantVesicles {
    
    public double meanFluorescence;
    public double kPhoto; // k of photolysis reaction (1/ms)
    public double kPhotoSave; // variable for saving kPhoto to output file
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("kPhoto")) {
            return "1/" + project.timeUnits;
        }
        return super.units(name);
    }
    
    public DiffusantVesiclesPhoto(Project p, String NAME, double InitialConcentration, double DiffusionConstant, CoordinatesVoxels c,
            double radius, PulseTimer pt, PSF PSF) {

        super(p, NAME, InitialConcentration, DiffusionConstant, c, radius);

        pulseTimer = pt;
        reaction = true;
        psf = PSF;

        if (save != null) {
            save.xdim = project.timeUnits;
            save.ydim = "kPhoto (1/" + project.timeUnits + ")";
        }

        createVector(true);

    }
    
    @Override
    public void vesicleStats() {

        if (project.monteCarlo == null) {
            return;
        }

        if (vesicles == null) {
            return;
        }
        
        super.vesicleStats();

        meanFluorescence = 0;

        for (DiffusantVesicle v : vesicles) {

            if (v == null) {
                continue;
            }

            if (!v.insideGeometry) {
                continue;
            }
            
            meanFluorescence += v.fluorescence;

        }

        meanFluorescence /= numVesicles;

        setParamObject("meanFluorescence", meanFluorescence);

    }
    
    @Override
    public boolean save() {

        if (save != null) {
            save.saveData(kPhotoSave);
        }
        
        return super.save();

    }
    
    public void react(int it) {

        double photolysis;

        if ((pulseTimer == null) || (it >= pulseTimer.timer.length)) {
            return;
        }

        if (pulseTimer.timer[it] == 0) {
            kPhotoSave = 0;
            return;
        }
        
        kPhotoSave = kPhoto;

        for (DiffusantVesicle v : vesicles) {

            if (v == null) {
                continue;
            }

            if (v.voxel == null) {
                continue;
            }

            if (!v.insideGeometry) {
                continue;
            }

            if (v.voxel.PSFi >= 0) {
                photolysis = v.fluorescence * kPhoto * v.voxel.PSFi * project.dt;
                v.fluorescence -= photolysis;
            }

        }

    }
    
    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantVesicles)) {
            return false;
        }
        
        String n = o.getName();

        
        if (n.equalsIgnoreCase("meanFluorescence")) {
            if (v < 0) {
                return false;
            }
            return setVesicles("fluorescence",v);
        }
        if (n.equalsIgnoreCase("kPhoto")) {
            if (v < 0) {
                return false;
            }
            kPhoto = v;
            return true;
        }
        return super.setMyParams(o, v);
    }
    
}
