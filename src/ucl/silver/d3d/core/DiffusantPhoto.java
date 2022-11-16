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
public class DiffusantPhoto extends Diffusant {

    // [a] ---> [b]
    //
    // this Diffusant represents [a]
    // which is converted to [b] at rate Kphoto
    //
    public int aDiffusantNum; // this diffusant
    public String bDiffusantName = null; // diffusant name for [b]
    public int bDiffusantNum = -1; // diffusant number for [b], if -1 then [b] is not saved
    
    // Kphoto is PulseTimer amplitude

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("Kphoto")) {
            return "1/" + project.timeUnits;
        }
        return super.units(name);
    }
    
    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("aDiffusantNum")) {
            return false;
        }
        return super.canEdit(name);
    }

    public DiffusantPhoto(Project p, String NAME, double InitialConcentration, double diffusionConstant, CoordinatesVoxels c,
            int DiffusantNumB, PulseTimer pt, PSF PSF) {

        super(p, NAME, InitialConcentration, diffusionConstant, c);

        bDiffusantNum = DiffusantNumB;
        pulseTimer = pt;
        psf = PSF;
        reaction = true;

        createVector(true);

    }

    @Override
    public void init() {

        super.init();
        
        reaction = true;
        aDiffusantNum = -1;

        for (int i = 0; i < project.diffusants.length; i++) {
            if (project.diffusants[i] == this) {
                aDiffusantNum = i;
                break;
            }
        }

        if ((bDiffusantNum >= 0) && (!project.checkDiffusantNum(bDiffusantNum))) {
            error("init", "bDiffusantNum", "out of range");
            return;
        }
        
        bDiffusantName = project.diffusantName(bDiffusantNum);
        
        if (aDiffusantNum >= 0) {
            if (aDiffusantNum == bDiffusantNum) {
                error("init", "bDiffusantNum", "value conflicts with aDiffusantNum");
            }
        }
        
        if (pulseTimer == null) {
            //error("init", "pulseTimer", "no pulse timer for photolysis reaction");
            Master.log("warning: no pulseTimer for photolysis reaction: " + this);
        }

    }
    
    @Override
    public void saveFileName(){

        if (save != null) {
            save.fileName(name, "Kphoto");
        }

    }
    
    @Override
    public void saveDimensions() {
        
        if ((save != null) && save.autoDimensions) {
            save.xdim = project.timeUnits;
            save.ydim = "Kphoto (1/" + project.timeUnits + ")";
        }

    }

    @Override
    public void react(RunFiniteDifference fd, int aNum) {
        double photolysis;
        int ipnt = fd.thisPnt;
        
        saveValue = 0;

        if (aNum != aDiffusantNum) {
            return;
        }
        
        if (fd.it >= pulseTimer.timer.length) {
            return;
        }

        saveValue = pulseTimer.high(fd.it); // Kphoto = PulseTimer amplitude
        
        if (saveValue == 0) {
            return; // nothing to do
        }

        if (fd.psf == null) {
            photolysis = fd.diffnext[aNum][ipnt] * saveValue * project.dt;
        } else {
            photolysis = fd.diffnext[aNum][ipnt] * saveValue * fd.psf[aNum][ipnt] * project.dt;
        }

        fd.diffnext[aNum][ipnt] -= photolysis;

        if (fd.diffnext[aNum][ipnt] < 0) {
            error("react", "fd.diffnext[aDiffusantNum]", "negative concentration for diffusant #" + aNum);
            fd.stopSimulation();
        }

        if (bDiffusantNum >= 0) {
            fd.diffnext[bDiffusantNum][ipnt] += photolysis;
            if (fd.diffnext[bDiffusantNum][ipnt] < 0) {
                error("react", "fd.diffnext[bDiffusantNum]", "negative concentration for diffusant #" + bDiffusantNum);
                fd.stopSimulation();
            }
        }

    }

    public boolean SetDiffusantNumB(int diffusantNum) {

        Diffusant d = project.getDiffusant(diffusantNum);

        if (d != null) {
            bDiffusantNum = diffusantNum;
            bDiffusantName = d.name;
        } else {
            bDiffusantNum = -1;
            bDiffusantName = null;
        }

        updateVectors();

        return true;

    }

    public boolean SetDiffusantNameB(String diffusantName) {

        int i = project.getDiffusantNum(diffusantName);

        if (i >= 0) {
            bDiffusantName = diffusantName;
            bDiffusantNum = i;
        } else {
            bDiffusantNum = -1;
            bDiffusantName = null;
        }

        updateVectors();

        return true;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantPhoto)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("bDiffusantNum")) {
            return SetDiffusantNumB((int) v);
        }
        return super.setMyParams(o, v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantPhoto)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("bDiffusantName")) {
            return SetDiffusantNameB(s);
        }
        return super.setMyParams(o, s);
    }

}
