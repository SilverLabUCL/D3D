package ucl.silver.d3d.core;

import java.text.*;

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
public class CoordinatesVoxels extends ParamVector {

    public boolean voxelPriority = true; // false - compute voxels from dimensions (um); true - compute dimensions (um) from voxels

    public double xWidth, yWidth, zWidth; // dimensions (um)
    public double volume; // um^3

    public double x1, y1, z1; // corner 1 location (um)
    public double x2, y2, z2; // corner 2 location (um)
    public double xCenter, yCenter, zCenter; // center location (um)

    public double voxelWidth; // um
    public int xVoxels, yVoxels, zVoxels;
    public int voxels, spaceVoxels;
    public double spaceVolume; // um^3

    public int xVoxel1, yVoxel1, zVoxel1; // voxel corner 1 location
    public int xVoxel2, yVoxel2, zVoxel2; // voxel corner 2 location
    public double xVoxelCenter, yVoxelCenter, zVoxelCenter; // center point
    
    public boolean shell = false;
    public boolean ellipsoid = false;
    public int xScale = 1, yScale = 1, zScale = 1; // used with ellipsoid for circular variation (0) no (1) yes

    public transient int[] indexRFD = null; // index counter for RunFiniteDifference

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
        if (name.equalsIgnoreCase("voxelWidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("volume")) {
            return project.spaceUnits + "^3";
        }
        if (name.equalsIgnoreCase("spaceVolume")) {
            return project.spaceUnits + "^3";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {

        if (name.equalsIgnoreCase("shape")) {
            return true;
        }

        if (voxelPriority) {
            if (name.equalsIgnoreCase("xVoxel1")) {
                return true;
            }
            if (name.equalsIgnoreCase("yVoxel1")) {
                return true;
            }
            if (name.equalsIgnoreCase("zVoxel1")) {
                return true;
            }
            if (name.equalsIgnoreCase("xVoxel2")) {
                return true;
            }
            if (name.equalsIgnoreCase("yVoxel2")) {
                return true;
            }
            if (name.equalsIgnoreCase("zVoxel2")) {
                return true;
            }
        } else {
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
        }
        return false;
    }

    public CoordinatesVoxels(Project p) {

        super(p);

        createVector(true);

        if (project != null) {

            Geometry g = project.geometry;

            if (g != null) {
                xVoxel1 = g.xVoxel1;
                yVoxel1 = g.yVoxel1;
                zVoxel1 = g.zVoxel1;
                xVoxel2 = g.xVoxel2;
                yVoxel2 = g.yVoxel2;
                zVoxel2 = g.zVoxel2;
                voxelPriority = true;
                updateVoxels();
            }

        }

    }

    public CoordinatesVoxels(Project p, int Xvoxel1, int Yvoxel1, int Zvoxel1, int Xvoxel2, int Yvoxel2, int Zvoxel2) {
        super(p);
        createVector(true);
        setVoxels(Xvoxel1, Yvoxel1, Zvoxel1, Xvoxel2, Yvoxel2, Zvoxel2);
    }

    public CoordinatesVoxels(Project p, double X1, double Y1, double Z1, double X2, double Y2, double Z2) {
        super(p);
        createVector(true);
        voxelPriority = false;
        setDimensions(X1, Y1, Z1, X2, Y2, Z2);
    }

    public CoordinatesVoxels(Project p, CoordinatesVoxels c) {
        super(p);
        createVector(true);
        matchVoxels(c);
    }

    public void update() {

        if (voxelPriority) {
            updateVoxels();
        } else {
            updateDimensions();
        }

    }

    private void updateVoxels() { // voxel priority

        int xv1, yv1, zv1, xv2, yv2, zv2;

        if ((project == null) || (project.geometry == null)) {
            return;
        }

        voxelWidth = project.dx;

        Geometry g = project.geometry;

        xv1 = Math.min(xVoxel1, xVoxel2);
        yv1 = Math.min(yVoxel1, yVoxel2);
        zv1 = Math.min(zVoxel1, zVoxel2);
        xv2 = Math.max(xVoxel1, xVoxel2);
        yv2 = Math.max(yVoxel1, yVoxel2);
        zv2 = Math.max(zVoxel1, zVoxel2);

        xVoxel1 = xv1;
        yVoxel1 = yv1;
        zVoxel1 = zv1;
        xVoxel2 = xv2;
        yVoxel2 = yv2;
        zVoxel2 = zv2;

        xVoxels = xVoxel2 - xVoxel1 + 1;
        yVoxels = yVoxel2 - yVoxel1 + 1;
        zVoxels = zVoxel2 - zVoxel1 + 1;

        xVoxelCenter = xVoxel1 + (xVoxels - 1) / 2.0;
        yVoxelCenter = yVoxel1 + (yVoxels - 1) / 2.0;
        zVoxelCenter = zVoxel1 + (zVoxels - 1) / 2.0;

        voxels = xVoxels * yVoxels * zVoxels;

        x1 = g.computeXedge(xVoxel1);
        y1 = g.computeYedge(yVoxel1);
        z1 = g.computeZedge(zVoxel1);

        x2 = g.computeXedge(xVoxel2);
        y2 = g.computeYedge(yVoxel2);
        z2 = g.computeZedge(zVoxel2);

        xWidth = xVoxels * voxelWidth;
        yWidth = yVoxels * voxelWidth;
        zWidth = zVoxels * voxelWidth;

        xCenter = x1 + xWidth / 2.0;
        yCenter = y1 + yWidth / 2.0;
        zCenter = z1 + zWidth / 2.0;

        volume = xWidth * yWidth * zWidth;

        spaceVoxels = g.countSpaceVoxels(this);
        spaceVolume = spaceVoxels * voxelWidth * voxelWidth * voxelWidth;

        updateVectors();

    }

    private void updateDimensions() { // dimension priority

        double X1, Y1, Z1, X2, Y2, Z2;

        if ((project == null) || (project.geometry == null)) {
            return;
        }

        voxelWidth = project.dx;

        Geometry g = project.geometry;

        X1 = Math.min(x1, x2);
        Y1 = Math.min(y1, y2);
        Z1 = Math.min(z1, z2);
        X2 = Math.max(x1, x2);
        Y2 = Math.max(y1, y2);
        Z2 = Math.max(z1, z2);

        x1 = X1;
        y1 = Y1;
        z1 = Z1;
        x2 = X2;
        y2 = Y2;
        z2 = Z2;

        xWidth = x2 - x1;
        yWidth = y2 - y1;
        zWidth = z2 - z1;

        xCenter = x1 + xWidth / 2.0;
        yCenter = y1 + yWidth / 2.0;
        zCenter = z1 + zWidth / 2.0;

        volume = xWidth * yWidth * zWidth;

        xVoxel1 = (int) g.computeVoxelX(x1);
        yVoxel1 = (int) g.computeVoxelY(y1);
        zVoxel1 = (int) g.computeVoxelZ(z1);
        xVoxel2 = (int) g.computeVoxelX(x2);
        yVoxel2 = (int) g.computeVoxelY(y2);
        zVoxel2 = (int) g.computeVoxelZ(z2);

        xVoxels = xVoxel2 - xVoxel1 + 1;
        yVoxels = yVoxel2 - yVoxel1 + 1;
        zVoxels = zVoxel2 - zVoxel1 + 1;

        xVoxelCenter = xVoxel1 + (xVoxels - 1) / 2.0;
        yVoxelCenter = yVoxel1 + (yVoxels - 1) / 2.0;
        zVoxelCenter = zVoxel1 + (zVoxels - 1) / 2.0;

        voxels = xVoxels * yVoxels * zVoxels;

        spaceVoxels = g.countSpaceVoxels(this);
        spaceVolume = spaceVoxels * voxelWidth * voxelWidth * voxelWidth;

        updateVectors();

    }

    public int updateIndexRFD(){

        int count = 0;

        indexRFD = null;

        for (int k = zVoxel1; k <= zVoxel2; k++) {
            for (int j = yVoxel1; j <= yVoxel2; j++) {
                for (int i = xVoxel1; i <= xVoxel2; i++) {
                    if (project.geometry.isSpace(i, j, k) && isInside(i, j, k)) {
                        count++;
                    }
                }
            }
        }

        if (count == 0) {
            return 0;
        } else {
            indexRFD = new int[count];
        }

        count = 0;

        for (int k = zVoxel1; k <= zVoxel2; k++) {
            for (int j = yVoxel1; j <= yVoxel2; j++) {
                for (int i = xVoxel1; i <= xVoxel2; i++) {
                    if (project.geometry.isSpace(i, j, k) && isInside(i, j, k)) {
                        indexRFD[count] = project.geometry.space[i][j][k];
                        count++;
                    }
                }
            }
        }

        return count;

    }

    public final boolean matchVoxels(CoordinatesVoxels c) {

        boolean changed = false;

        if (voxelWidth != c.voxelWidth) {
            voxelWidth = c.voxelWidth;
            changed = true;
        }

        if (xVoxel1 != c.xVoxel1) {
            xVoxel1 = c.xVoxel1;
            changed = true;
        }

        if (yVoxel1 != c.yVoxel1) {
            yVoxel1 = c.yVoxel1;
            changed = true;
        }

        if (zVoxel1 != c.zVoxel1) {
            zVoxel1 = c.zVoxel1;
            changed = true;
        }

        if (xVoxel2 != c.xVoxel2) {
            xVoxel2 = c.xVoxel2;
            changed = true;
        }

        if (yVoxel2 != c.yVoxel2) {
            yVoxel2 = c.yVoxel2;
            changed = true;
        }

        if (zVoxel2 != c.zVoxel2) {
            zVoxel2 = c.zVoxel2;
            changed = true;
        }

        if (xScale != c.xScale) {
            xScale = c.xScale;
            changed = true;
        }

        if (yScale != c.yScale) {
            yScale = c.yScale;
            changed = true;
        }

        if (zScale != c.zScale) {
            zScale = c.zScale;
            changed = true;
        }

        if (changed) {
            updateVoxels();
        }

        return changed;

    }

    public void setCenter(double x, double y, double z) {

        Geometry g = project.geometry;

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

    public final boolean matchDimensions(CoordinatesVoxels c) {

        boolean changed = false;

        if (x1 != c.x1) {
            x1 = c.x1;
            changed = true;
        }

        if (y1 != c.y1) {
            y1 = c.y1;
            changed = true;
        }

        if (z1 != c.z1) {
            z1 = c.z1;
            changed = true;
        }

        if (x2 != c.x2) {
            x2 = c.x2;
            changed = true;
        }

        if (y2 != c.y2) {
            y2 = c.y2;
            changed = true;
        }

        if (z2 != c.z2) {
            z2 = c.z2;
            changed = true;
        }

        if (changed) {
            updateDimensions();
        }

        return changed;

    }

    public final void setInside(CoordinatesVoxels c, int extraVoxels) {

        xVoxel1 = c.xVoxel1 + extraVoxels;
        yVoxel1 = c.yVoxel1 + extraVoxels;
        zVoxel1 = c.zVoxel1 + extraVoxels;
        xVoxel2 = c.xVoxel2 - extraVoxels;
        yVoxel2 = c.yVoxel2 - extraVoxels;
        zVoxel2 = c.zVoxel2 - extraVoxels;
        
        xScale = c.xScale;
        yScale = c.yScale;
        zScale = c.zScale;

        updateVoxels();

    }
    
    public final void setOutside(CoordinatesVoxels c, int extraVoxels) {

        xVoxel1 = c.xVoxel1 - extraVoxels;
        yVoxel1 = c.yVoxel1 - extraVoxels;
        zVoxel1 = c.zVoxel1 - extraVoxels;
        xVoxel2 = c.xVoxel2 + extraVoxels;
        yVoxel2 = c.yVoxel2 + extraVoxels;
        zVoxel2 = c.zVoxel2 + extraVoxels;
        
        xScale = c.xScale;
        yScale = c.yScale;
        zScale = c.zScale;

        updateVoxels();

    }

    public final CoordinatesVoxels setDimensions(double X1, double Y1, double Z1,
            double X2, double Y2, double Z2) {

        x1 = X1;
        y1 = Y1;
        z1 = Z1;
        x2 = X2;
        y2 = Y2;
        z2 = Z2;

        updateDimensions();

        return this;

    }

    public final CoordinatesVoxels setVoxels(double Xvoxel1, double Yvoxel1, double Zvoxel1,
            double Xvoxel2, double Yvoxel2, double Zvoxel2) {
        return setVoxels((int) Xvoxel1, (int) Yvoxel1, (int) Zvoxel1, (int) Xvoxel2, (int) Yvoxel2, (int) Zvoxel2);
    }

    public final CoordinatesVoxels setVoxels(int Xvoxel1, int Yvoxel1, int Zvoxel1,
            int Xvoxel2, int Yvoxel2, int Zvoxel2) {

        xVoxel1 = Xvoxel1;
        yVoxel1 = Yvoxel1;
        zVoxel1 = Zvoxel1;
        xVoxel2 = Xvoxel2;
        yVoxel2 = Yvoxel2;
        zVoxel2 = Zvoxel2;

        updateVoxels();

        return this;

    }

    public final CoordinatesVoxels setPoint(double x, double y, double z) {
        return setDimensions(x, y, z, x, y, z);
    }
    
    public final CoordinatesVoxels setVoxelPoint(double xVoxel, double yVoxel, double zVoxel) {
        return setVoxels(xVoxel, yVoxel, zVoxel, xVoxel, yVoxel, zVoxel); // a single voxel
    }

    public final CoordinatesVoxels setVoxelPoint(int xVoxel, int yVoxel, int zVoxel) {
        return setVoxels(xVoxel, yVoxel, zVoxel, xVoxel, yVoxel, zVoxel); // a single voxel
    }

    public final CoordinatesVoxels xySquare(double x, double y, double z, double xyWidth) {

        double half = Math.abs(xyWidth / 2.0);

        x1 = x - half;
        y1 = y - half;
        z1 = z;
        x2 = x1 + xyWidth;
        y2 = y1 + xyWidth;
        z2 = z;

        updateDimensions();

        return this;
    }

    public final CoordinatesVoxels xyVoxelSquare(double xVoxel, double yVoxel, double zVoxel, double xyVoxelWidth) {

        double half = xyVoxelWidth / 2.0;

        xVoxel1 = (int) Math.ceil(xVoxel - half);
        yVoxel1 = (int) Math.ceil(yVoxel - half);
        zVoxel1 = (int) zVoxel;
        xVoxel2 = (int) (xVoxel1 + xyVoxelWidth - 1);
        yVoxel2 = (int) (yVoxel1 + xyVoxelWidth - 1);
        zVoxel2 = (int) zVoxel;
        
        updateVoxels();

        return this;

    }
    
    public final CoordinatesVoxels yzVoxelSquare(double xVoxel, double yVoxel, double zVoxel, double yzVoxelWidth) {

        double half = yzVoxelWidth / 2.0;

        xVoxel1 = (int) xVoxel;
        yVoxel1 = (int) Math.ceil(yVoxel - half);
        zVoxel1 = (int) Math.ceil(zVoxel - half);
        xVoxel2 = (int) xVoxel;
        yVoxel2 = (int) (yVoxel1 + yzVoxelWidth - 1);
        zVoxel2 = (int) (zVoxel1 + yzVoxelWidth - 1);
        
        updateVoxels();

        return this;

    }
    
    public final CoordinatesVoxels zxVoxelSquare(double xVoxel, double yVoxel, double zVoxel, double zxVoxelWidth) {

        double half = zxVoxelWidth / 2.0;

        xVoxel1 = (int) Math.ceil(xVoxel - half);
        yVoxel1 = (int) yVoxel;
        zVoxel1 = (int) Math.ceil(zVoxel - half);
        xVoxel2 = (int) (xVoxel1 + zxVoxelWidth - 1);
        yVoxel2 = (int) yVoxel;
        zVoxel2 = (int) (zVoxel1 + zxVoxelWidth - 1);
        
        updateVoxels();

        return this;

    }

    public boolean isInside(CoordinatesVoxels c) {
        for (int k = c.zVoxel1; k <= c.zVoxel2; k += 1) {
            for (int j = c.yVoxel1; j <= c.yVoxel2; j += 1) {
                for (int i = c.xVoxel1; i <= c.xVoxel2; i += 1) {
                    if (isInside(i, j, k)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    public boolean isInside(int xVoxel, int yVoxel, int zVoxel) {

        if ((xVoxel >= xVoxel1) && (yVoxel >= yVoxel1) && (zVoxel >= zVoxel1)) {
            if ((xVoxel <= xVoxel2) && (yVoxel <= yVoxel2) && (zVoxel <= zVoxel2)) {
                // return true
            } else {
                return false;
            }
        } else {
            return false;
        }

        if (ellipsoid) {
            return isInsideEllipsoid(xVoxel, yVoxel, zVoxel);
        }
        
        if (shell) {
            if ((xVoxel >= xVoxel1 + 1) && (yVoxel >= yVoxel1 + 1) && (zVoxel >= zVoxel1 + 1)) {
                if ((xVoxel <= xVoxel2 - 1) && (yVoxel <= yVoxel2 - 1) && (zVoxel <= zVoxel2 - 1)) {
                    return false;
                } else {
                    return true;
                }
            }
        }

        return true;

    }
    
    public boolean isInsideEllipsoid(int xVoxel, int yVoxel, int zVoxel) {
        
        double xv, yv, zv;

        double xd = (xVoxel - xVoxelCenter) / (xVoxels / 2.0);
        double yd = (yVoxel - yVoxelCenter) / (yVoxels / 2.0);
        double zd = (zVoxel - zVoxelCenter) / (zVoxels / 2.0);

        if ((xd * xd * xScale + yd * yd * yScale + zd * zd * zScale) <= 1) {
            
            if (shell) {
                
                if (xVoxels > 2) {
                    xv = xVoxels - 2;
                } else {
                    xv = xVoxels;
                }
                
                if (yVoxels > 2) {
                    yv = yVoxels - 2;
                } else {
                    yv = yVoxels;
                }
                
                if (zVoxels > 2) {
                    zv = zVoxels - 2;
                } else {
                    zv = zVoxels;
                }
                
                xd = (xVoxel - xVoxelCenter) / (xv / 2.0);
                yd = (yVoxel - yVoxelCenter) / (yv / 2.0);
                zd = (zVoxel - zVoxelCenter) / (zv / 2.0);
                
                return ((xd * xd * xScale + yd * yd * yScale + zd * zd * zScale) > 1);
                
            }
            
            return true;
            
        } else {
            
            return false;
            
        }

    }
    
    public boolean isInsideEllipsoid(double x, double y, double z) {
        
        double xv, yv, zv, dx2;

        double xd = (x - xCenter) / (xWidth / 2.0);
        double yd = (y - yCenter) / (yWidth / 2.0);
        double zd = (z - zCenter) / (zWidth / 2.0);

        if ((xd * xd * xScale + yd * yd * yScale + zd * zd * zScale) <= 1) {
            
            if (shell) {
                
                dx2 = 2 * project.dx;
                
                if (xWidth > dx2) {
                    xv = xWidth - dx2;
                } else {
                    xv = xWidth;
                }
                
                if (yWidth > 2) {
                    yv = yWidth - dx2;
                } else {
                    yv = yWidth;
                }
                
                if (zWidth > 2) {
                    zv = zWidth - dx2;
                } else {
                    zv = zWidth;
                }
                
                xd = (x - xCenter) / (xv / 2.0);
                yd = (y - yCenter) / (yv / 2.0);
                zd = (z - zCenter) / (zv / 2.0);
                
                return ((xd * xd * xScale + yd * yd * yScale + zd * zd * zScale) > 1);
                
            }
            
            return true;
            
        } else {
            
            return false;
            
        }

    }

    public boolean isInside(double x, double y, double z) {
        
        double dx = project.dx;
        
        if ((x >= x1) && (y >= y1) && (z >= z1) && (x <= x2) && (y <= y2) && (z <= z2)) {
            //return true;
        } else {
            return false;
        }
        
        if (ellipsoid) {
            return isInsideEllipsoid(x, y, z);
        }
        
        if (shell) {
            if ((x >= x1 + dx) && (y >= y1 + dx) && (z >= z1 + dx)) {
                if ((x <= x2 - dx) && (y <= y2 - dx) && (z <= z2 - dx)) {
                    return false;
                } else {
                    return true;
                }
            }
        }
        
        return true;
    }
    
    public boolean isInsideGeometry() {

        Geometry g = project.geometry;

        if ((x1 >= g.x1) && (y1 >= g.y1) && (z1 >= g.z1)) {
            if ((x2 <= g.x2) && (y2 <= g.y2) && (z2 <= g.z2)) {
                return true;
            }
        }

        return false;

    }

    public String voxelsToString(boolean noNames) {

        String c1, c2;

        if (noNames) {
            c1 = xVoxel1 + "," + yVoxel1 + "," + zVoxel1;
            c2 = xVoxel2 + "," + yVoxel2 + "," + zVoxel2;
        } else {
            c1 = "xVoxel1=" + xVoxel1 + ",yVoxel1=" + yVoxel1 + ",zVoxel1=" + zVoxel1;
            c2 = "xVoxel2=" + xVoxel2 + ",yVoxel2=" + yVoxel2 + ",zVoxel2=" + zVoxel2;
        }

        return c1 + "," + c2;

    }

    public String dimensionsToString(boolean noNames) {

        String c1, c2;

        if (noNames) {
            c1 = formatter.format(x1) + "," + formatter.format(y1) + "," + formatter.format(z1);
            c2 = formatter.format(x2) + "," + formatter.format(y2) + "," + formatter.format(z2);
        } else {
            c1 = "x1=" + formatter.format(x1) + ",y1=" + formatter.format(y1) + ",z1=" + formatter.format(z1);
            c2 = "x2=" + formatter.format(x2) + ",y2=" + formatter.format(y2) + ",z2=" + formatter.format(z2);
        }

        return c1 + "," + c2;

    }

    public String scalesToString(boolean noNames) {

        String c1;
        
        if (noNames) {
            c1 = formatter.format(xScale) + "," + formatter.format(yScale) + "," + formatter.format(zScale);
        } else {
            c1 = "xScale=" + formatter.format(xScale) + ",yScale=" + formatter.format(yScale) + ",zScale=" + formatter.format(zScale);
        }

        return c1;
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof CoordinatesVoxels)) {
            return false;
        }

        String n = o.getName();

        if (voxelPriority) {
            if (n.equalsIgnoreCase("xVoxel1")) {
                if (v < 0) {
                    return false;
                }
                xVoxel1 = (int) v;
                updateVoxels();
                return true;
            }
            if (n.equalsIgnoreCase("yVoxel1")) {
                if (v < 0) {
                    return false;
                }
                yVoxel1 = (int) v;
                updateVoxels();
                return true;
            }
            if (n.equalsIgnoreCase("zVoxel1")) {
                if (v < 0) {
                    return false;
                }
                zVoxel1 = (int) v;
                updateVoxels();
                return true;
            }
            if (n.equalsIgnoreCase("xVoxel2")) {
                if (v < 0) {
                    return false;
                }
                xVoxel2 = (int) v;
                updateVoxels();
                return true;
            }
            if (n.equalsIgnoreCase("yVoxel2")) {
                if (v < 0) {
                    return false;
                }
                yVoxel2 = (int) v;
                updateVoxels();
                return true;
            }
            if (n.equalsIgnoreCase("zVoxel2")) {
                if (v < 0) {
                    return false;
                }
                zVoxel2 = (int) v;
                updateVoxels();
                return true;
            }
        } else {
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

        if (!(o.paramVector instanceof CoordinatesVoxels)) {
            return false;
        }

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        return super.setMyParams(o, s);
    }

}
