package ucl.silver.d3d.init;

import ucl.silver.d3d.core.*;
import ucl.silver.d3d.utils.*;

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
public class InitProject extends ParamVector {
    
    public String directory = "";
    public String folder = "Testing";

    Geometry geometry;
    CoordinatesVoxels coordinates; // general purpose coordinates

    public String[] initFuncList = {"initCube"};

    public InitProject(Project p) {
        super(p);
        geometry = project.geometry;
        coordinates = new CoordinatesVoxels(p);
        createVector(true);
    }

    public int initFunctionNum(String initSelect) {
        return Utility.whichListItem(initFuncList, initSelect, true);
    }

    public boolean initFunction(String initSelect) {

        int i = initFunctionNum(initSelect);

        switch (i) {
            case 0:
                return initCube();
            
            default:
                error("initFunction", "initSelect", "failed to find init function");
                return true; // error
        }

    }
    
    @Override
    public boolean canEdit(String name) {
        return super.canEdit(name);
    }

    @Override
    public String units(String name) {
        return super.units(name);
    }

    @Override
    public String help(String name) {
        return super.help(name);
    }
    
    public void initDirectoryJason() {
        
        String os;
        String dirWin = "/Jason/D3D/Simulations/";
        String dirMac = "/Users/jason/Documents/D3D/Simulations/";

        if (!Master.foundMainStartUpArguments) {
            
            os = System.getProperty("os.name").toLowerCase();

            if (os.contains("win")) {
                directory = dirWin;
            } else if (os.contains("mac")) {
                directory = dirMac;
            }

            project.directory = directory;
            project.folder = folder;

        }

    }

    public boolean initCube() {

        int cubeWidth = 10;

        project.name = "Simple Cube";
        project.simTime = 1.0;
        project.dx = 0.05;

        geometry.resizeWithSpace(cubeWidth, cubeWidth, cubeWidth);
        geometry.clear();

        return false;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof InitProject)) {
            return false;
        }
        
        String n = o.getName();

        return super.setMyParams(o, v);

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof InitProject)) {
            return false;
        }

        String n = o.getName();

        return super.setMyParams(o, s);

    }
    
}
