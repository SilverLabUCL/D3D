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
public class DiffusantParticlesPhoto extends DiffusantParticles {
    
    public double meanFluorescence;
    
    // Kphoto is PulseTimer amplitude
    
    public DiffusantParticlesPhoto(Project p, String NAME, double InitialConcentration, double DiffusionConstant, CoordinatesVoxels c,
            double radius_mean, double radius_stdv, PulseTimer pt, PSF PSF) {

        super(p, NAME, InitialConcentration, DiffusionConstant, c, radius_mean, radius_stdv);

        pulseTimer = pt;
        reaction = true;
        psf = PSF;

        createVector(true);

    }
    
    @Override
    public void init() {

        super.init();
        
        if (pulseTimer == null) {
            //error("init", "pulseTimer", "no pulse timer for photolysis reaction");
            Master.log("warning: no pulseTimer for photolysis reaction: " + this);
        }
        
    }
    
    @Override
    public void particleStats() {

        if (project.monteCarlo == null) {
            return;
        }

        if (particles == null) {
            return;
        }
        
        super.particleStats();

        meanFluorescence = 0;

        for (DiffusantParticle p : particles) {

            if (p == null) {
                continue;
            }

            if (!p.insideGeometry) {
                continue;
            }
            
            meanFluorescence += p.fluorescence;

        }

        meanFluorescence /= numParticles;

        setParamObject("meanFluorescence", meanFluorescence);

    }
    
    @Override
    public void saveFileName(){
        
        super.saveFileName();

        if (save != null) {
            save.fileName(name, "Kphoto");
        }
        
        if (save_XYZ != null) {
            save_XYZ.fileName(name, "XYZF");
        }

    }
    
    @Override
    public void saveDimensions() {
        
        super.saveDimensions();
        
        if ((save != null) && save.autoDimensions) {
            save.xdim = project.timeUnits;
            save.ydim = "Kphoto (1/" + project.timeUnits + ")";
        }

    }
    
    public void react(int it) {

        double photolysis;

        if ((pulseTimer == null) || (it >= pulseTimer.timer.length)) {
            return;
        }
        
        saveValue = pulseTimer.high(it); // Kphoto = PulseTimer amplitude
        
        if (saveValue == 0) {
            return; // nothing to do
        }

        for (DiffusantParticle p : particles) {

            if (p == null) {
                continue;
            }

            if (p.voxel == null) {
                continue;
            }

            if (!p.insideGeometry) {
                continue;
            }

            if (p.voxel.PSFi >= 0) {
                photolysis = p.fluorescence * saveValue * p.voxel.PSFi * project.dt;
                p.fluorescence -= photolysis;
            }

        }

    }
    
    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantParticlesPhoto)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("meanFluorescence")) {
            if (v < 0) {
                return false;
            }
            return setParticles("fluorescence",v);
        }
        return super.setMyParams(o, v);
    }
    
}
