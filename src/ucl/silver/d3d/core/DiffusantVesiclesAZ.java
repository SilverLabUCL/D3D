package ucl.silver.d3d.core;

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
public class DiffusantVesiclesAZ extends DiffusantVesicles {

    public double meanDockRefractoryPeriod; // ms

    public int setNumVesiclesReady;

    public int numVesiclesReady; // number of ready vesicles (numReady = numMobile + numImmobile)

    public int dockedVesiclesInit = 0; // initial number of docked vesicles
    public double dockedVesiclesInitFraction = 0; // initial number of docked vesicles
    public int dockedVesicles; // initial number of docked vesicles
    public int reserveVesiclesInit = 0; // initial number of reserve vesicles
    public int reserveVesicles; // initial number of reserve vesicles

    //public ColorD3D colorDocked = new ColorD3D( "colorDocked", new Color(204,0,0));
    public ColorD3D colorDocked = new ColorD3D( "colorDocked", new Color(0,136,204));
    public ColorD3D colorReserve = new ColorD3D( "colorReserve", new Color(102,255,255));

    //CoordinatesVoxels initCoordinates = null; // coordinates for initial vesicles placement ( leave null for entire shape )

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("meanDockRefractoryPeriod")) {
            return project.timeUnits;
        }
        return super.units(name);
    }

    public DiffusantVesiclesAZ(Project p, String NAME, double InitialConcentration, double DiffusionConstant, CoordinatesVoxels c,
            double radius) {

        super(p, NAME, InitialConcentration, DiffusionConstant, c, radius);

        createVector(true);

    }

    @Override
    public void vesicleStats() {

        double count = 0.0;

        DiffusantVesicleAZ vaz;

        if (vesicles == null) {
            return;
        }

        super.vesicleStats();

        if (numVesicles <= 0) {
            return;
        }

        meanDockRefractoryPeriod = 0;
        numVesiclesReady = 0;
        dockedVesicles = 0;
        reserveVesicles = 0;

        for (DiffusantVesicle v : vesicles) {

            if (v == null) {
                continue;
            }

            if (!(v instanceof DiffusantVesicleAZ)) {
                continue;
            }

            vaz = (DiffusantVesicleAZ) v;

            count++;

            meanDockRefractoryPeriod += vaz.dockRefractoryPeriod;

            if (vaz.isReady) {
                numVesiclesReady++;
            } else if (vaz.isDocked) {
                dockedVesicles++;
            } else if (vaz.isReserve) {
                reserveVesicles++;
            } else {
                //Master.exit("DiffusantVesiclesAZ: vesicleStats: unknown vesicle type " + vesicles[i].name);
            }

        }

        meanDockRefractoryPeriod /= count;

        setParamObject("meanDockRefractoryPeriod", meanDockRefractoryPeriod);

        setParamObject("numVesiclesReady", numVesiclesReady);
        setParamObject("dockedVesicles", dockedVesicles);
        setParamObject("reserveVesicles", reserveVesicles);

    }

    @Override
    public void initVesicles() {

        RunMonteCarloAZ monteCarlo = ( RunMonteCarloAZ) project.monteCarlo;

        double spaceVolume = monteCarlo.spaceVolume() - monteCarlo.activeZone.geometryVolume();

        if (setNumVesiclesReady > 0) {
            numVesiclesReady = setNumVesiclesReady;
        } else if ((setVolumeFraction > 0) && (spaceVolume > 0)) {
            numVesiclesReady = numVesiclesPossible(spaceVolume, setMeanRadius, setVolumeFraction);
        } else if ((setDensity > 0) && (spaceVolume > 0)) {
            numVesiclesReady = (int) (setDensity * spaceVolume);
        }

        numVesicles = numVesiclesReady + reserveVesiclesInit;

        if (numVesicles <= 0) {
            error("bad number of vesicles: " + numVesicles);
            return;
        }

        vesicles = new DiffusantVesicleAZ[numVesicles];

        for (int i = 0; i < reserveVesiclesInit; i++) {
            vesicles[i] = new DiffusantVesicleAZ(project, "reserve", setMeanRadius, D, Double.NaN, Double.NaN, Double.NaN);
        }

        if ((xyz != null) && (xyz.length == numVesicles)) {
            for (int i = reserveVesiclesInit; i < vesicles.length; i++) {
                vesicles[i] = new DiffusantVesicleAZ(project, "ready", setMeanRadius, D, xyz[i][0], xyz[i][1], xyz[i][2]);
            }
        } else {
            for (int i = reserveVesiclesInit; i < vesicles.length; i++) {
                vesicles[i] = new DiffusantVesicleAZ(project, "ready", setMeanRadius, D, Double.NaN, Double.NaN, Double.NaN);
            }
        }

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (colorDocked != null) {
            colorDocked.addUser(pv);
        }

        if (colorReserve != null) {
            colorReserve.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (colorDocked != null) {
            addBlankParam();
            colorDocked.createVector(true);
            addVector(colorDocked.getVector());
            colorDocked.addUser(this);
        }

        if (colorReserve != null) {
            addBlankParam();
            colorReserve.createVector(true);
            addVector(colorReserve.getVector());
            colorReserve.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (colorDocked != null) {
            colorDocked.updateVector(v);
        }

        if (colorReserve != null) {
            colorReserve.updateVector(v);
        }

    }

    @Override
    public boolean setVesicles(String varName, double value) {

        boolean ok, atLeastOne = false;

        DiffusantVesicleAZ vaz;

        if (vesicles == null) {
            return false;
        }

        for (DiffusantVesicle v : vesicles) {

            if ((v != null) && (v instanceof DiffusantVesicleAZ)) {

                vaz = (DiffusantVesicleAZ) v;

                ok = vaz.set(varName, value);

                if (ok) {
                    atLeastOne = true;
                }

            }
        }

        vesicleStats();

        return atLeastOne;

    }

    @Override
    public double meanSquareDisplacement() {

        double sumSD = 0.0, count = 0.0;

        DiffusantVesicleAZ vaz;

        for (DiffusantVesicle v : vesicles) {

            if ((v != null) && (v instanceof DiffusantVesicleAZ)) {

                vaz = (DiffusantVesicleAZ) v;

                if (vaz.mobile && vaz.insideGeometry && (vaz.d2AZ0 >= 0) && (vaz.d2AZ0 < 0.05)) {
                    sumSD += vaz.squareDisplacement();
                    count += 1.0;
                    //Master.log("" + dv + " " + dv.d2AZ0);
                }

            }


        }

        return sumSD / count;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantVesiclesAZ)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("meanDockRefractoryPeriod")) {
            if (v < 0) {
                return false;
            }
            return setVesicles("dockRefractoryPeriod",v);
        }
        if (n.equalsIgnoreCase("setNumVesiclesReady")) {
            if (v < 0) {
                return false;
            }
            setNumVesiclesReady = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("dockedVesiclesInit")) {
            if (v < 0) {
                return false;
            }
            dockedVesiclesInit = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("reserveVesiclesInit")) {
            if (v < 0) {
                return false;
            }
            reserveVesiclesInit = (int) v;
            return true;
        }
        return super.setMyParams(o, v);
    }

}
