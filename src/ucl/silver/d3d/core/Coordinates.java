package ucl.silver.d3d.core;

import java.text.*;
import java.awt.Color;

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
public class Coordinates extends ParamVector {

    public double xWidth, yWidth, zWidth; // dimensions (um)
    public double volume; // um^3
    public double geometryVolume; // um^3, volume inside geometry

    public double x1, y1, z1; // corner 1 location (um)
    public double x2, y2, z2; // corner 2 location (um)
    public double xCenter, yCenter, zCenter; // center location (um)

    public String shape = "cuboid"; // cuboid, ellipsoid, cylinderx, cylindery, cylinderz
    public int shapeNum = 0; // 0, 1, 2, 3, 4

    public ColorD3D color = new ColorD3D( "color", Color.white);

    private final static double FOUR_THIRDS_PI = 4.0 * Math.PI / 3.0;
    
    private final NumberFormat formatter = new DecimalFormat("0.####");

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("x1")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("y1")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("z1")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("x2")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("y2")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("z2")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("xCenter")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("yCenter")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("zCenter")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("xWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("yWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("zWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("volume")) {
            return project.spaceUnits + "^3";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("shape")) {
            return true;
        }
        if (name.equalsIgnoreCase("shapeNum")) {
            return true;
        }
        if (name.equalsIgnoreCase("x1")) {
            return true;
        }
        if (name.equalsIgnoreCase("y1")) {
            return true;
        }
        if (name.equalsIgnoreCase("z1")) {
            return true;
        }
        if (name.equalsIgnoreCase("x2")) {
            return true;
        }
        if (name.equalsIgnoreCase("y2")) {
            return true;
        }
        if (name.equalsIgnoreCase("z2")) {
            return true;
        }
        return false;
    }

    public Coordinates(Project p) {
        super(p);
        createVector(true);
    }

    public Coordinates(Project p, double X1, double Y1, double Z1, double X2, double Y2, double Z2) {
        super(p);
        createVector(true);
        setCoordinates(X1, Y1, Z1, X2, Y2, Z2);
    }

    public Coordinates(Project p, Coordinates c) {
        super(p);
        createVector(true);
        copy(c);
    }

    public Coordinates(Project p, CoordinatesVoxels c) {
        super(p);
        createVector(true);
        setCoordinates(c.x1, c.y1, c.z1, c.x2, c.y2, c.z2);
    }
    
    public boolean setShape(int i) {

        switch (i) {
            case 0:
                shape = "cubiod";
                break;
            case 1:
                shape = "ellipsoid";
                break;
            case 2:
                shape = "cylinderx";
                break;
            case 3:
                shape = "cylindery";
                break;
            case 4:
                shape = "cylinderz";
                break;
            default:
                return false;
        }

        shapeNum = i;

        return true;

    }

    public boolean setShape(String s) {

        if (s.matches("cuboid")) {
            shapeNum = 0;
        } else if (s.matches("ellipsoid")) {
            shapeNum = 1;
        } else if (s.matches("cylinderx")) {
            shapeNum = 2;
        } else if (s.matches("cylindery")) {
            shapeNum = 3;
        } else if (s.matches("cylinderz")) {
            shapeNum = 4;
        } else {
            return false;
        }

        shape = s;

        return true;

    }

    public final void updateDimensions() {

        double X1 = Math.min(x1, x2);
        double Y1 = Math.min(y1, y2);
        double Z1 = Math.min(z1, z2);
        double X2 = Math.max(x1, x2);
        double Y2 = Math.max(y1, y2);
        double Z2 = Math.max(z1, z2);

        x1 = X1;
        y1 = Y1;
        z1 = Z1;
        x2 = X2;
        y2 = Y2;
        z2 = Z2;

        xWidth = x2 - x1;
        yWidth = y2 - y1;
        zWidth = z2 - z1;

        xCenter = (x1 + x2) / 2.0;
        yCenter = (y1 + y2) / 2.0;
        zCenter = (z1 + z2) / 2.0;

        volume = volume();
        geometryVolume = geometryVolume();

        updateVectors();

    }

    public final void copy(Coordinates c) {
        x1 = c.x1;
        y1 = c.y1;
        z1 = c.z1;
        x2 = c.x2;
        y2 = c.y2;
        z2 = c.z2;
        shape = c.shape;
        shapeNum = c.shapeNum;
        name = c.name;
        updateDimensions();
    }

    public final void setCoordinates(double X1, double Y1, double Z1,
            double X2, double Y2, double Z2) {
        x1 = X1;
        y1 = Y1;
        z1 = Z1;
        x2 = X2;
        y2 = Y2;
        z2 = Z2;
        updateDimensions();
    }

    public final void setWidths(double XWIDTH, double YWIDTH, double ZWIDTH) {
        x1 = xCenter - XWIDTH / 2.0;
        y1 = yCenter - YWIDTH / 2.0;
        z1 = zCenter - ZWIDTH / 2.0;
        x2 = x1 + XWIDTH;
        y2 = y1 + YWIDTH;
        z2 = z1 + ZWIDTH;
        updateDimensions();
    }

    public final void setXwidth(double XWIDTH) {
        x1 = xCenter - XWIDTH / 2.0;
        x2 = x1 + XWIDTH;
        updateDimensions();
    }

    public final void setYwidth(double YWIDTH) {
        y1 = yCenter - YWIDTH / 2.0;
        y2 = y1 + YWIDTH;
        updateDimensions();
    }

    public final void setZwidth(double ZWIDTH) {
        z1 = zCenter - ZWIDTH / 2.0;
        z2 = z1 + ZWIDTH;
        updateDimensions();
    }

    public final void setCenter(double x, double y, double z) {

        if (!Double.isNaN(x)) {
            x1 = x - xWidth / 2.0;
            x2 = x1 + xWidth;
        }

        if (!Double.isNaN(y)) {
            y1 = y - yWidth / 2.0;
            y2 = y1 + yWidth;
        }

        if (!Double.isNaN(z)) {
            z1 = z - zWidth / 2.0;
            z2 = z1 + zWidth;
        }

        updateDimensions();

    }

    public final void setXcenter(double x) {
        x1 = x - xWidth / 2.0;
        x2 = x1 + xWidth;
        updateDimensions();
    }

    public final void setYcenter(double y) {
        y1 = y - yWidth / 2.0;
        y2 = y1 + yWidth;
        updateDimensions();
    }

    public final void setZcenter(double z) {
        z1 = z - zWidth / 2.0;
        z2 = z1 + zWidth;
        updateDimensions();
    }

    public double distanceToCenter(double x, double y, double z) {

        double dx = x - xCenter;
        double dy = y - yCenter;
        double dz = z - zCenter;

        return Math.sqrt(dx * dx + dy * dy + dz * dz);

    }

    public boolean isInsideGeometryCompletely() {

        Geometry g = project.geometry;

        if ((x1 >= g.x1) && (y1 >= g.y1) && (z1 >= g.z1)) {
            if ((x2 <= g.x2) && (y2 <= g.y2) && (z2 <= g.z2)) {
                return true;
            }
        }

        return false;

    }

    public boolean isInsideCuboidCompletely(Coordinates c) {
        if ((x1 >= c.x1) && (y1 >= c.y1) && (z1 >= c.z1)) {
            if ((x2 <= c.x2) && (y2 <= c.y2) && (z2 <= c.z2)) {
                return true;
            }
        }
        return false;
    }

    public boolean intersectsCuboid(Coordinates c) {
        return intersectsCuboid(c.x1, c.y1, c.z1, c.x2, c.y2, c.z2);
    }

    public boolean intersectsCuboid(double X1, double Y1, double Z1, double X2, double Y2, double Z2) {
        boolean outside = (x2 < X1) || (x1 > X2) || (y2 < Y1) || (y1 > Y2) || (z2 < Z1) || (z1 > Z2);
        return !outside;
    }

    public boolean isInside(double x, double y, double z) {
        switch(shapeNum) {
            case 0:
                return isInsideCuboid(x, y, z);
            case 1:
                return isInsideEllipsoid(x, y, z);
            case 2:
                return isInsideCylinderX(x, y, z);
            case 3:
                return isInsideCylinderY(x, y, z);
            case 4:
                return isInsideCylinderZ(x, y, z);
        }
        return false;
    }

    public boolean isInsideCuboid(double x, double y, double z) {
        return ((x >= x1) && (y >= y1) && (z >= z1) && (x <= x2) && (y <= y2) && (z <= z2));
    }

    public boolean isInsideEllipsoid(double x, double y, double z) {

        double xv = (x - xCenter) / (xWidth / 2.0);
        double yv = (y - yCenter) / (yWidth / 2.0);
        double zv = (z - zCenter) / (zWidth / 2.0);
        
        return (xv * xv + yv * yv + zv * zv <= 1);

    }

    public boolean isInsideCylinderX(double x, double y, double z) {

        double yv = (y - yCenter) / (yWidth / 2.0);
        double zv = (z - zCenter) / (zWidth / 2.0);

        if ((x >= x1) && (x <= x2)) {
            if ((yv * yv + zv * zv) <= 1) {
                return true;
            }
        }

        return false;

    }

    public boolean isInsideCylinderY(double x, double y, double z) {

        double xv = (x - xCenter) / (xWidth / 2.0);
        double zv = (z - zCenter) / (zWidth / 2.0);

        if ((y >= y1) && (y <= y2)) {
            if ((xv * xv + zv * zv) <= 1) {
                return true;
            }
        }

        return false;

    }

    public boolean isInsideCylinderZ(double x, double y, double z) {

        double xv = (x - xCenter) / (xWidth / 2.0);
        double yv = (y - yCenter) / (yWidth / 2.0);

        if ((z >= z1) && (z <= z2)) {
            if ((xv * xv + yv * yv) <= 1) {
                return true;
            }
        }

        return false;

    }

    public boolean isInsideCylinderZ_Perimeter(double x, double y, double z, double fractionInside) {

        double xv = (x - xCenter) / (xWidth / 2.0);
        double yv = (y - yCenter) / (yWidth / 2.0);
        double rxy2 = xv * xv + yv * yv;

        // fractionInside ( 0 to 1 )

        if ((z >= z1) && (z <= z2)) {
            if ((rxy2 <= 1) && (rxy2 >= fractionInside)) {
                return true;
            }
        }

        return false;

    }

    public boolean isInside(double x, double y, double z, double sphereRadius) {
        switch(shapeNum) {
            case 0:
                return isInsideCuboid(x, y, z, sphereRadius);
            case 1:
                return isInsideEllipsoid(x, y, z, sphereRadius);
            case 2:
                return isInsideCylinderX(x, y, z, sphereRadius);
            case 3:
                return isInsideCylinderY(x, y, z, sphereRadius);
            case 4:
                return isInsideCylinderZ(x, y, z, sphereRadius);
        }
        return false;
    }

    private boolean isInsideCuboid(double x, double y, double z, double sphereRadius) {
        if ((x + sphereRadius >= x1) && (y + sphereRadius >= y1) && (z + sphereRadius >= z1)) {
            if ((x - sphereRadius <= x2) && (y - sphereRadius <= y2) && (z - sphereRadius <= z2)) {
                return true;
            }
        }

        return false;

    }

    public boolean isInsideEllipsoid(double x, double y, double z, double sphereRadius) {

        //double xv = (x - xCenter) / ((0.5 * xWidth - 0.5 * sphereRadius));
        //double yv = (y - yCenter) / ((0.5 * yWidth - 0.5 * sphereRadius));
        //double zv = (z - zCenter) / ((0.5 * zWidth - 0.5 * sphereRadius));

        double xv = (x - xCenter) / ((xWidth + 2 * sphereRadius) / 2.0);
        double yv = (y - yCenter) / ((yWidth + 2 * sphereRadius) / 2.0);
        double zv = (z - zCenter) / ((zWidth + 2 * sphereRadius) / 2.0);
        
        return (xv * xv + yv * yv + zv * zv <= 1);

    }

    public boolean isOutsideEllipsoid(double x, double y, double z, double sphereRadius) {
        double xv = (x - xCenter) / ((0.5 * xWidth + sphereRadius));
        double yv = (y - yCenter) / ((0.5 * yWidth + sphereRadius));
        double zv = (z - zCenter) / ((0.5 * zWidth + sphereRadius));
        return (xv * xv + yv * yv + zv * zv > 1);
    }

    private boolean isInsideCylinderX(double x, double y, double z, double sphereRadius) {

        double yv = (y - yCenter) / ((0.5 * yWidth + sphereRadius));
        double zv = (z - zCenter) / ((0.5 * zWidth + sphereRadius));

        if ((x >= x1 - sphereRadius) && (x <= x2 + sphereRadius)) {
            return ((yv * yv + zv * zv) <= 1);
        }

        return false;

    }

    private boolean isInsideCylinderY(double x, double y, double z, double sphereRadius) {

        double xv = (x - xCenter) / ((0.5 * xWidth + sphereRadius));
        double zv = (z - zCenter) / ((0.5 * zWidth + sphereRadius));

        if ((y >= y1 - sphereRadius) && (y <= y2 + sphereRadius)) {
            return ((xv * xv + zv * zv) <= 1);
        }

        return false;

    }

    private boolean isInsideCylinderZ(double x, double y, double z, double sphereRadius) {

        double xv = (x - xCenter) / ((0.5 * xWidth + sphereRadius));
        double yv = (y - yCenter) / ((0.5 * yWidth + sphereRadius));

        if ((z >= z1 - sphereRadius) && (z <= z2 + sphereRadius)) {
            return ((xv * xv + yv * yv) <= 1);
        }

        return false;

    }

    public double volume() {
        switch (shapeNum) {
            case 0: // cuboid
                return xWidth * yWidth * zWidth;
            case 1: // ellipsoid
                return FOUR_THIRDS_PI * (xWidth / 2.0) * (yWidth / 2.0) * (zWidth / 2.0);
            case 2: // cylinderx
                return Math.PI * xWidth * (yWidth / 2.0) * (zWidth / 2.0);
            case 3: // cylindery
                return Math.PI * (xWidth / 2.0) * yWidth * (zWidth / 2.0);
            case 4: // cylinderz
                return Math.PI * (xWidth / 2.0) * (yWidth / 2.0) * zWidth;
        }
        return 0;
    }

    public double geometryVolume() {

        boolean exact = true;

        if (isInsideGeometryCompletely()) {
            return volume();
        }

        if (!exact) {
            return geometryVolumeCuboid();
        }

        switch (shapeNum) {
            case 0: // cuboid
                return geometryVolumeCuboid();
            case 1: // ellipsoid
                return geometryVolumeEllipsoid();
            case 2: // cylinderx
                return geometryVolumeCylinderX();
            case 3: // cylindery
                return geometryVolumeCylinderY();
            case 4: // cylinderz
                return geometryVolumeCylinderZ();
        }

        return 0;

    }

    public double geometryVolumeCuboid() {

        Geometry g = project.geometry;

        double X1 = Math.max(x1, g.x1);
        double Y1 = Math.max(y1, g.y1);
        double Z1 = Math.max(z1, g.z1);
        double X2 = Math.min(x2, g.x2);
        double Y2 = Math.min(y2, g.y2);
        double Z2 = Math.min(z2, g.z2);

        return (X2 - X1) * (Y2 - Y1) * (Z2 - Z1);

    }

    public double geometryVolumeEllipsoid() {

        int isteps, jsteps;
        double X1, Y1, Z1, X2, Y2, Z2;
        double xa;
        double ya, yb, y2ab, ya2b;
        double za, zb, z2ab, za2b;
        double area = 0;

        double istep = 0.001; // um
        double jstep = 0.001; // um

        boolean clipToGeometry = true;

        Geometry g = project.geometry;

        if (clipToGeometry) {
            X1 = Math.max(x1, g.x1);
            Y1 = Math.max(y1, g.y1);
            Z1 = Math.max(z1, g.z1);
            X2 = Math.min(x2, g.x2);
            Y2 = Math.min(y2, g.y2);
            Z2 = Math.min(z2, g.z2);
        } else {
            X1 = x1;
            Y1 = y1;
            Z1 = z1;
            X2 = x2;
            Y2 = y2;
            Z2 = z2;
        }

        isteps = (int) (1 + Math.ceil((X2 - X1) / istep));
        jsteps = (int) (1 + Math.ceil((Y2 - Y1) / jstep));

        for (int i = 0; i < isteps; i++) {

            xa = X1 + i * istep;

            for (int j = 0; j < jsteps; j++) {

                ya = Y1 + j * jstep;
                yb = Y1 + (j + 1) * jstep;
                y2ab = (2 * ya + yb) / 3;
                ya2b = (ya + 2 * yb) / 3;

                if ((ya < Y1) || (ya > Y2)) {
                    continue;
                }

                if ((yb < Y1) || (yb > Y2)) {
                    continue;
                }

                za = zCircle(xa, ya, false, clipToGeometry);
                zb = zCircle(xa, yb, false, clipToGeometry);
                z2ab = zCircle(xa, y2ab, false, clipToGeometry);
                za2b = zCircle(xa, ya2b, false, clipToGeometry);

                area += (jstep / 8) * ((za - zCenter) + 3 * (z2ab - zCenter) + 3 * (za2b - zCenter) + (zb - zCenter));

                za = zCircle(xa, ya, true, clipToGeometry);
                zb = zCircle(xa, yb, true, clipToGeometry);
                z2ab = zCircle(xa, y2ab, true, clipToGeometry);
                za2b = zCircle(xa, ya2b, true, clipToGeometry);

                area += (jstep / 8) * ((zCenter - za) + 3 * (zCenter - z2ab) + 3 * (zCenter - za2b) + (zCenter - zb));

            }

        }

        area *= istep;

        if (Double.isNaN(area)) {
            printDimensions();
            Master.log("ellipsoid area = " + area + ", ( " + (4 * Math.PI * (xWidth / 2) * (yWidth / 2) * (zWidth / 2) / 3.0) + " )");
        }

        return area;
        
    }

    public double zCircle(double x, double y, boolean negative, boolean clipToGeometry) {

        double X = (x - xCenter) / (xWidth / 2.0);
        double Y = (y - yCenter) / (yWidth / 2.0);
        double oneX2Y2 = 1 - X * X - Y * Y;
        double z = zCenter + (zWidth / 2.0) * Math.sqrt(Math.max(oneX2Y2, 0));

        if (negative) {
            z = -z + 2 * zCenter;
        }

        if (clipToGeometry) {
            z = Math.max(z, project.geometry.z1);
            z = Math.min(z, project.geometry.z2);
        }

        return z;

    }

    public double geometryVolumeCylinderX() {

        int jsteps;
        double X1, Y1, Z1, X2, Y2, Z2;
        double ya, yb, y2ab, ya2b;
        double za, zb, z2ab, za2b;
        double area = 0;

        double step = 0.001; // um

        boolean clipToGeometry = true;
        
        Geometry g = project.geometry;

        if (clipToGeometry) {
            X1 = Math.max(x1, g.x1);
            Y1 = Math.max(y1, g.y1);
            Z1 = Math.max(z1, g.z1);
            X2 = Math.min(x2, g.x2);
            Y2 = Math.min(y2, g.y2);
            Z2 = Math.min(z2, g.z2);
        } else {
            X1 = x1;
            Y1 = y1;
            Z1 = z1;
            X2 = x2;
            Y2 = y2;
            Z2 = z2;
        }

        jsteps = (int) (1 + Math.ceil((Y2 - Y1) / step));

        for (int j = 0; j < jsteps; j++) {

            ya = Y1 + j * step;
            yb = Y1 + (j + 1) * step;
            y2ab = (2 * ya + yb) / 3;
            ya2b = (ya + 2 * yb) / 3;

            if ((ya < Y1) || (ya > Y2)) {
                continue;
            }

            if ((yb < Y1) || (yb > Y2)) {
                continue;
            }

            za = zCircle(ya, false, clipToGeometry);
            zb = zCircle(yb, false, clipToGeometry);
            z2ab = zCircle(y2ab, false, clipToGeometry);
            za2b = zCircle(ya2b, false, clipToGeometry);

            area += (step / 8) * ((za - zCenter) + 3 * (z2ab - zCenter) + 3 * (za2b - zCenter) + (zb - zCenter));

            za = zCircle(ya, true, clipToGeometry);
            zb = zCircle(yb, true, clipToGeometry);
            z2ab = zCircle(y2ab, true, clipToGeometry);
            za2b = zCircle(ya2b, true, clipToGeometry);

            area += (step / 8) * ((zCenter - za) + 3 * (zCenter - z2ab) + 3 * (zCenter - za2b) + (zCenter - zb));

        }

        if (Double.isNaN(area)) {
            printDimensions();
            Master.log("cyl-x : area = " + area + ", ( " + (Math.PI * (zWidth / 2) * (yWidth / 2)) + " )");
        }

        return (X2 - X1) * area;
        
    }

    public double geometryVolumeCylinderY() {

        int ksteps;
        double X1, Y1, Z1, X2, Y2, Z2;
        double xa, xb, x2ab, xa2b;
        double za, zb, z2ab, za2b;
        double area = 0;

        double step = 0.001; // um

        boolean clipToGeometry = true;

        Geometry g = project.geometry;

        if (clipToGeometry) {
            X1 = Math.max(x1, g.x1);
            Y1 = Math.max(y1, g.y1);
            Z1 = Math.max(z1, g.z1);
            X2 = Math.min(x2, g.x2);
            Y2 = Math.min(y2, g.y2);
            Z2 = Math.min(z2, g.z2);
        } else {
            X1 = x1;
            Y1 = y1;
            Z1 = z1;
            X2 = x2;
            Y2 = y2;
            Z2 = z2;
        }

        ksteps = (int) (1 + Math.ceil((Z2 - Z1) / step));

        for (int k = 0; k < ksteps; k++) {

            za = Z1 + k * step;
            zb = Z1 + (k + 1) * step;
            z2ab = (2 * za + zb) / 3;
            za2b = (za + 2 * zb) / 3;

            if ((za < Z1) || (za > Z2)) {
                continue;
            }

            if ((zb < Z1) || (zb > Z2)) {
                continue;
            }

            xa = xCircle(za, false, clipToGeometry);
            xb = xCircle(zb, false, clipToGeometry);
            x2ab = xCircle(z2ab, false, clipToGeometry);
            xa2b = xCircle(za2b, false, clipToGeometry);

            area += (step / 8) * ((xa - xCenter) + 3 * (x2ab - xCenter) + 3 * (xa2b - xCenter) + (xb - xCenter));

            xa = xCircle(za, true, clipToGeometry);
            xb = xCircle(zb, true, clipToGeometry);
            x2ab = xCircle(z2ab, true, clipToGeometry);
            xa2b = xCircle(za2b, true, clipToGeometry);

            area += (step / 8) * ((xCenter - xa) + 3 * (xCenter - x2ab) + 3 * (xCenter - xa2b) + (xCenter - xb));

        }

        if (Double.isNaN(area)) {
            printDimensions();
            Master.log("cyl-y : area = " + area + ", ( " + (Math.PI * (zWidth / 2) * (xWidth / 2)) + " )");
        }

        return (Y2 - Y1) * area;

    }

    public double geometryVolumeCylinderZ() {

        int isteps;
        double X1, Y1, Z1, X2, Y2, Z2;
        double xa, xb, x2ab, xa2b;
        double ya, yb, y2ab, ya2b;
        double area = 0;

        double step = 0.001; // um

        boolean clipToGeometry = true;

        Geometry g = project.geometry;

        if (clipToGeometry) {
            X1 = Math.max(x1, g.x1);
            Y1 = Math.max(y1, g.y1);
            Z1 = Math.max(z1, g.z1);
            X2 = Math.min(x2, g.x2);
            Y2 = Math.min(y2, g.y2);
            Z2 = Math.min(z2, g.z2);
        } else {
            X1 = x1;
            Y1 = y1;
            Z1 = z1;
            X2 = x2;
            Y2 = y2;
            Z2 = z2;
        }

        isteps = (int) (1 + Math.ceil((X2 - X1) / step));

        for (int i = 0; i < isteps; i++) {

            xa = X1 + i * step;
            xb = X1 + (i + 1) * step;
            x2ab = (2 * xa + xb) / 3;
            xa2b = (xa + 2 * xb) / 3;

            if ((xa < X1) || (xa > X2)) {
                continue;
            }

            if ((xb < X1) || (xb > X2)) {
                continue;
            }

            ya = yCircle(xa, false, clipToGeometry);
            yb = yCircle(xb, false, clipToGeometry);
            y2ab = yCircle(x2ab, false, clipToGeometry);
            ya2b = yCircle(xa2b, false, clipToGeometry);

            area += (step / 8) * ((ya - yCenter) + 3 * (y2ab - yCenter) + 3 * (ya2b - yCenter) + (yb - yCenter));

            ya = yCircle(xa, true, clipToGeometry);
            yb = yCircle(xb, true, clipToGeometry);
            y2ab = yCircle(x2ab, true, clipToGeometry);
            ya2b = yCircle(xa2b, true, clipToGeometry);

            area += (step / 8) * ((yCenter - ya) + 3 * (yCenter - y2ab) + 3 * (yCenter - ya2b) + (yCenter - yb));

        }

        if (Double.isNaN(area)) {
            printDimensions();
            Master.log("cyl-z : area = " + area + ", ( " + (Math.PI * (xWidth / 2) * (yWidth / 2)) + " )");
        }

        return (Z2 - Z1) * area;

    }
    
    public double xCircle(double z, boolean negative, boolean clipToGeometry) {

        double Z = (z - zCenter) / (zWidth / 2.0);
        double oneZ2 = 1 - Z * Z;
        double x = xCenter + (xWidth / 2.0) * Math.sqrt(Math.max(oneZ2, 0));

        if (negative) {
            x = -x + 2 * xCenter;
        }

        if (clipToGeometry) {
            x = Math.max(x, project.geometry.x1);
            x = Math.min(x, project.geometry.x2);
        }

        return x;

    }

    public double yCircle(double x, boolean negative, boolean clipToGeometry) {

        double X = (x - xCenter) / (xWidth / 2.0);
        double oneX2 = 1 - X * X;
        double y = yCenter + (yWidth / 2.0) * Math.sqrt(Math.max(oneX2, 0));

        if (negative) {
            y = -y + 2 * yCenter;
        }

        if (clipToGeometry) {
            y = Math.max(y, project.geometry.y1);
            y = Math.min(y, project.geometry.y2);
        }

        return y;

    }

    public double zCircle(double y, boolean negative, boolean clipToGeometry) {

        double Y = (y - yCenter) / (yWidth / 2.0);
        double oneY2 = 1 - Y * Y;
        double z = zCenter + (zWidth / 2.0) * Math.sqrt(Math.max(oneY2, 0));

        if (negative) {
            z = -z + 2 * zCenter;
        }

        if (clipToGeometry) {
            z = Math.max(z, project.geometry.z1);
            z = Math.min(z, project.geometry.z2);
        }

        return z;

    }

    public String printDimensions() {
        String c1 = "x1=" + formatter.format(x1) + ",y1=" + formatter.format(y1) + ",z1=" + formatter.format(z1);
        String c2 = "x2=" + formatter.format(x2) + ",y2=" + formatter.format(y2) + ",z2=" + formatter.format(z2);
        String w = "xwidth=" + formatter.format(xWidth) + ",ywidth=" + formatter.format(yWidth) + ",zwidth=" + formatter.format(zWidth);
        String s = "shape=" + shape;
        String v = "volume=" + formatter.format(volume);
        Master.log(c1);
        Master.log(c2);
        Master.log(w);
        Master.log(s);
        Master.log(v);
        return c1 + "," + c2;
    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (color != null) {
            color.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (color != null) {
            addBlankParam();
            color.createVector(true);
            addVector(color.getVector());
            color.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (color != null) {
            color.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof Coordinates)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("x1")) {
            x1 = v;
            updateDimensions();
            return true;
        }
        if (n.equalsIgnoreCase("y1")) {
            y1 = v;
            updateDimensions();
            return true;
        }
        if (n.equalsIgnoreCase("z1")) {
            z1 = v;
            updateDimensions();
            return true;
        }
        if (n.equalsIgnoreCase("x2")) {
            x2 = v;
            updateDimensions();
            return true;
        }
        if (n.equalsIgnoreCase("y2")) {
            y2 = v;
            updateDimensions();
            return true;
        }
        if (n.equalsIgnoreCase("z2")) {
            z2 = v;
            updateDimensions();
            return true;
        }
        if (n.equalsIgnoreCase("xWidth")) {
            setXwidth(v);
            return true;
        }
        if (n.equalsIgnoreCase("yWidth")) {
            setYwidth(v);
            return true;
        }
        if (n.equalsIgnoreCase("zWidth")) {
            setZwidth(v);
            return true;
        }
        if (n.equalsIgnoreCase("xCenter")) {
            setXcenter(v);
            return true;
        }
        if (n.equalsIgnoreCase("yCenter")) {
            setYcenter(v);
            return true;
        }
        if (n.equalsIgnoreCase("zCenter")) {
            setZcenter(v);
            return true;
        }
        if (n.equalsIgnoreCase("shapeNum")) {
            return setShape((int) v);
        }

        return super.setMyParams(o, v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        String n = o.getName();

        if (!(o.paramVector instanceof Coordinates)) {
            return false;
        }

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        if (n.equalsIgnoreCase("shape")) {
            return setShape(s);
        }
        return super.setMyParams(o, s);
    }

}
