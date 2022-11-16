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
public class DiffusantVesicleAZ extends DiffusantParticle {

    //public String name = "ready"; // "ready" or "docked" or "reserve"

    public boolean isReady = true;
    public boolean isDocked = false;
    public boolean isReserve = false;
    public boolean isBound = false;
    public boolean isTetheredToAZ = false;
    //public boolean isInsideAZ = false;

    public double dockTime = -1;
    public double dockRefractoryPeriod = 0;

    public DiffusantVesicleAZ nextReady = null; // for ready vesicle list
    public DiffusantVesicleAZ nextDocked = null; // for docked vesicle list
    public DiffusantVesicleAZ nextReserve = null; // for reserve vesicle list

    public double d2AZ = -1; // distance to AZ um
    public double d2AZ0 = -1; // distance to AZ um at time = 0

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("dockTime")) {
            return project.timeUnits;
        }
        if (name.equalsIgnoreCase("dockRefractoryPeriod")) {
            return project.timeUnits;
        }
        return super.units(name);

    }

    public DiffusantVesicleAZ(Project p, String TypeName, double RADIUS, double d, double X, double Y, double Z) {

        super(p, TypeName, RADIUS, d, X, Y, Z);

        setType(TypeName);

    }

    public double squareDistanceFromActiveZone(Coordinates activeZone) {

        double dx = activeZone.xCenter - x;
        double dy = activeZone.yCenter - y;
        double dz = activeZone.zCenter - z;

        return dx * dx + dy * dy + dz * dz;

    }

    public boolean copy(DiffusantVesicleAZ dv) {

        if (dv == null) {
            return false;
        }

        super.copy(dv);

        x0 = dv.x0;
        y0 = dv.y0;
        z0 = dv.z0;

        isReady = dv.isReady;
        isDocked = dv.isDocked;
        isReserve = dv.isReserve;
        isBound = dv.isBound;
        isTetheredToAZ = dv.isTetheredToAZ;

        dockTime = dv.dockTime;
        dockRefractoryPeriod = dv.dockRefractoryPeriod;

        nextReady = dv.nextReady;
        nextDocked = dv.nextDocked;
        nextReserve = dv.nextReserve;

        connectTo = dv.connectTo;

        d2AZ = dv.d2AZ;
        d2AZ0 = dv.d2AZ0;

        return true;

    }

    final public boolean setType(String typeName) {

        if (typeName.equalsIgnoreCase("ready")) {
            insideGeometry = true;
            isReady = true;
            isDocked = false;
            isReserve = false;
            name = typeName;
            return true;
        } else if (typeName.equalsIgnoreCase("docked")) {
            insideGeometry = true;
            isReady = false;
            isDocked = true;
            isReserve = false;
            name = typeName;
            return true;
        } else if (typeName.equalsIgnoreCase("reserve")) {
            insideGeometry = false;
            isReady = false;
            isDocked = false;
            isReserve = true;
            name = typeName;
            return true;
        }

        return false;

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantVesicleAZ)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("dockTime")) {
            if (v < 0) {
                return false;
            }
            dockTime = v;
        }
        if (n.equalsIgnoreCase("dockRefractoryPeriod")) {
            if (v < 0) {
                return false;
            }
            dockRefractoryPeriod = v;
        }
        return super.setMyParams(o, v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantVesicleAZ)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            return setType(s);
        }
        
        return super.setMyParams(o, s);
        
    }

}
