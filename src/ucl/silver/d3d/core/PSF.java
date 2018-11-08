package ucl.silver.d3d.core;

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
public class PSF extends ParamVector {

    public String fileName = null; // import file of pre-computed PSF

    public boolean arbitraryVoxelCenter = false;
    public double xVoxelCenter, yVoxelCenter, zVoxelCenter; // for arbitrary center point

    public double max; // max value
    public double sum; // sum of PSF values
    public double avg; // sum / numvoxels
    public double dx; // voxel size (should be same as project.dx)

    public int computedValues;

    public boolean useGeometryCoordinates = true;

    public double rotateXY, rotateYZ, rotateZX; // degrees

    public boolean xySymmetric = false;
    public boolean normalize = true;

    //int rotZ; // OBSOLETE, psf rotation (0) none (1) swap yz (-1) flip z
    public int zCount = 50; // output counter

    public double[][][] array = null; // 3D psf

    private double iRotated, jRotated, kRotated;

    private CoordinatesVoxels coordinates = null;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("dx")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("rotateXY")) {
            return "degrees";
        }
        if (name.equalsIgnoreCase("rotateYZ")) {
            return "degrees";
        }
        if (name.equalsIgnoreCase("rotateZX")) {
            return "degrees";
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("max")) {
            return false;
        }
        if (name.equalsIgnoreCase("sum")) {
            return false;
        }
        if (name.equalsIgnoreCase("avg")) {
            return false;
        }
        if (name.equalsIgnoreCase("computedValues")) {
            return false;
        }
        if (name.equalsIgnoreCase("dx")) {
            return false;
        }
        if (name.equalsIgnoreCase("psfFile")) {
            return false;
        }
        if (name.equalsIgnoreCase("useGeometryCoordinates")) {
            return false;
        }
        if (name.equalsIgnoreCase("arbitraryVoxelCenter")) {
            return false;
        }
        if (!arbitraryVoxelCenter) {
            if (name.equalsIgnoreCase("xVoxelCenter")) {
                return false;
            }
            if (name.equalsIgnoreCase("yVoxelCenter")) {
                return false;
            }
            if (name.equalsIgnoreCase("zVoxelCenter")) {
                return false;
            }
        }
        return true;
    }

    public PSF(Project p, CoordinatesVoxels c) {

        super(p);

        dx = project.dx;

        if (c == null) {
            useGeometryCoordinates = true;
            coordinates = null;
        } else {
            useGeometryCoordinates = false;
            coordinates = new CoordinatesVoxels(p);
            coordinates.matchVoxels(c);
        }

        createVector(true);

    }

    public CoordinatesVoxels coordinates() {
        if ((useGeometryCoordinates) || (coordinates == null)) {
            return project.geometry;
        } else {
            return coordinates;
        }
    }

    public void updateCoordinates() {

        if (arbitraryVoxelCenter) {
            return;
        }

        xVoxelCenter = coordinates().xVoxelCenter;
        yVoxelCenter = coordinates().yVoxelCenter;
        zVoxelCenter = coordinates().zVoxelCenter;

    }

    public double xVoxelCenter() {
        if (arbitraryVoxelCenter) {
            return xVoxelCenter;
        } else {
            return coordinates().xVoxelCenter;
        }
    }

    public double yVoxelCenter() {
        if (arbitraryVoxelCenter) {
            return yVoxelCenter;
        } else {
            return coordinates().yVoxelCenter;
        }
    }

    public double zVoxelCenter() {
        if (arbitraryVoxelCenter) {
            return zVoxelCenter;
        } else {
            return coordinates().zVoxelCenter;
        }
    }

    public void setVoxelCenter(double XvoxelCenter, double YvoxelCenter, double ZvoxelCenter) {
        arbitraryVoxelCenter = true;
        xVoxelCenter = XvoxelCenter;
        yVoxelCenter = YvoxelCenter;
        zVoxelCenter = ZvoxelCenter;
    }

    public void init() {

        if (dx != project.dx) {
            max = -1;
            array = null; // force recomputation
        }

        if (max <= 0) {
            max = computeValue(xVoxelCenter(), yVoxelCenter(), zVoxelCenter(), false);
            computedValues = 0;
        }

        updateCoordinates();
        updateVectors();

    }

    public boolean exists() {
        return (array != null);
    }

    public void checkExists() {
        if (array == null) {
            compute();
        }
    }

    public double[][][] getPSF() {
        return array;
    }

    public boolean checkArrayXYZ(int xVoxel, int yVoxel, int zVoxel) {
        if (array == null) {
            return false;
        }
        if ((xVoxel < 0) || (xVoxel >= array.length)) {
            return false;
        }
        if ((yVoxel < 0) || (yVoxel >= array[0].length)) {
            return false;
        }
        if ((zVoxel < 0) || (zVoxel >= array[0][0].length)) {
            return false;
        }
        return true; // OK
    }

    public double getArrayValue(int xVoxel, int yVoxel, int zVoxel) {
        if (checkArrayXYZ(xVoxel, yVoxel, zVoxel)) {
            return array[xVoxel][yVoxel][zVoxel];
        } else {
            return -1;
        }
    }

    public double getValue(int xVoxel, int yVoxel, int zVoxel) {
        if (array != null) {
            return getArrayValue(xVoxel, yVoxel, zVoxel);
        } else {
            return computeValue(xVoxel, yVoxel, zVoxel, normalize);
        }
    }

    public boolean xySymmetric() {

        if (rotateXY != 0) {
            return false;
        }

        if (rotateYZ != 0) {
            return false;
        }

        if (rotateZX != 0) {
            return false;
        }

        return xySymmetric;

    }

    public double symI(double xVoxel) {

        double xvCenter = coordinates().xVoxelCenter;

        return (Math.floor(xvCenter - Math.abs(xVoxel - xvCenter)));

    }

    public double symJ(double yVoxel) {

        double yvCenter = coordinates().yVoxelCenter;

        return (Math.floor(yvCenter - Math.abs(yVoxel - yvCenter)));

    }

    public double iRotated() {
        return iRotated;
    }

    public double jRotated() {
        return jRotated;
    }

    public double kRotated() {
        return kRotated;
    }

    public void rotateAll(double xVoxel, double yVoxel, double zVoxel) {

        iRotated = rotateXYi(xVoxel, yVoxel, zVoxel);
        jRotated = rotateXYj(xVoxel, yVoxel, zVoxel);
        kRotated = zVoxel;

        xVoxel = iRotated;
        yVoxel = rotateYZj(iRotated, jRotated, kRotated);
        zVoxel = rotateYZk(iRotated, jRotated, kRotated);

        iRotated = rotateZXi(xVoxel, yVoxel, zVoxel);
        jRotated = yVoxel;
        kRotated = rotateZXk(xVoxel, yVoxel, zVoxel);

    }

    public double rotateXYi(double xVoxel, double yVoxel, double zVoxel) {

        double r, a;

        if (rotateXY != 0) {

            r = Math.sqrt(Math.pow(xVoxel, 2) + Math.pow(yVoxel, 2));
            a = Math.atan2(yVoxel, xVoxel); // convert to polar coordinates

            a += Math.PI * rotateXY / 180.0; // rotate

            xVoxel = r * Math.cos(a); // convert back to cartesian coordinates
            yVoxel = r * Math.sin(a);

        }

        return xVoxel;

    }

    public double rotateZXi(double xVoxel, double yVoxel, double zVoxel) {

        double r, a;

        if (rotateZX != 0) {

            r = Math.sqrt(Math.pow(xVoxel, 2) + Math.pow(zVoxel, 2));
            a = Math.atan2(zVoxel, xVoxel);

            a += Math.PI * rotateZX / 180.0;

            xVoxel = r * Math.cos(a);
            zVoxel = r * Math.sin(a);

        }

        return xVoxel;

    }

    public double rotateXYj(double xVoxel, double yVoxel, double zVoxel) {

        double r, a;

        if (rotateXY != 0) {

            r = Math.sqrt(Math.pow(xVoxel, 2) + Math.pow(yVoxel, 2));
            a = Math.atan2(yVoxel, xVoxel);

            a += Math.PI * rotateXY / 180.0;

            xVoxel = r * Math.cos(a);
            yVoxel = r * Math.sin(a);

        }

        return yVoxel;

    }

    public double rotateYZj(double xVoxel, double yVoxel, double zVoxel) {

        double r, a;

        if (rotateYZ != 0) {

            r = Math.sqrt(Math.pow(yVoxel, 2) + Math.pow(zVoxel, 2));
            a = Math.atan2(zVoxel, yVoxel);

            a += Math.PI * rotateYZ / 180.0;

            yVoxel = r * Math.cos(a);
            zVoxel = r * Math.sin(a);

        }

        return yVoxel;

    }

    public double rotateYZk(double xVoxel, double yVoxel, double zVoxel) {

        double r, a;

        if (rotateYZ != 0) {

            r = Math.sqrt(Math.pow(yVoxel, 2) + Math.pow(zVoxel, 2));
            a = Math.atan2(zVoxel, yVoxel);

            a += Math.PI * rotateYZ / 180.0;

            yVoxel = r * Math.cos(a);
            zVoxel = r * Math.sin(a);

        }

        return zVoxel;

    }

    public double rotateZXk(double xVoxel, double yVoxel, double zVoxel) {

        double r, a;

        if (rotateZX != 0) {

            r = Math.sqrt(Math.pow(xVoxel, 2) + Math.pow(zVoxel, 2));
            a = Math.atan2(zVoxel, xVoxel);

            a += Math.PI * rotateZX / 180.0;

            xVoxel = r * Math.cos(a);
            zVoxel = r * Math.sin(a);

        }

        return zVoxel;

    }

    public void compute() {

        int ii, jj, kk, count = 0;
        int xVoxels, yVoxels, zVoxels;
        double v;

        init();

        xVoxels = coordinates().xVoxels;
        yVoxels = coordinates().yVoxels;
        zVoxels = coordinates().zVoxels;

        if ((xVoxels < 0) || (yVoxels < 0) || (zVoxels < 0)) {
            return;
        }

        array = new double[xVoxels][yVoxels][zVoxels];

        sum = 0;
        computedValues = 0;

        if (name == null) {
            Master.log("computing PSF...");
        } else {
            Master.log("computing " + name + "...");
        }

        for (int k = 0; k < array[0][0].length; k++) {

            if ((k > 0) && (Math.IEEEremainder(k, zCount)) == 0) {
                Master.log("PSF k: " + k);
            }

            for (int j = 0; j < array[0].length; j++) {
                for (int i = 0; i < array.length; i++) {

                    if (xySymmetric()) {

                        ii = (int) symI(i);
                        jj = (int) symJ(j);
                        kk = k;

                        if ((ii >= 0) && (jj >= 0) && (kk >= 0) && (array[ii][jj][kk] > 0)) {
                            v = array[ii][jj][kk]; // x-y symmetrical
                            //v = 0; // for visual testing
                            array[i][j][k] = v;
                            sum += v;
                            count++;
                            continue;
                        }

                    }

                    v = computeValue(i, j, k, normalize);

                    array[i][j][k] = v;
                    sum += v;
                    count++;

                }
            }

        }

        if (count > 0) {
            avg = sum / count;
        }

        //Master.log("PSF sum: " + sum);
        findMax();

        updateVectors();

    }

    public double computeValue(double xVoxel, double yVoxel, double zVoxel, boolean normalize) {

        double value = computeVoxel(xVoxel, yVoxel, zVoxel);

        if (normalize && (max <= 0)) {
            init();
        }

        if (normalize && (max > 0)) {
            value /= max; // normalize
        }

        computedValues += 1;

        if (value < 0) {
            error("computeValue", "value", "negative value at point (" + xVoxel + "," + yVoxel + "," + zVoxel + ")");
        }

        return value;

    }

    public double computeVoxel(double xVoxel, double yVoxel, double zVoxel) { // this method can be over-rided

        double x, y, z;

        dx = project.dx;

        xVoxel -= coordinates().xVoxelCenter;
        yVoxel -= coordinates().yVoxelCenter;
        zVoxel -= coordinates().zVoxelCenter;

        rotateAll(xVoxel, yVoxel, zVoxel); // results saved in irotated, jrotated, krotated

        x = iRotated() * dx; // um
        y = jRotated() * dx; // um
        z = kRotated() * dx; // um

        if ((Math.abs(x) < 1) && (Math.abs(y) < 1) && (Math.abs(z) < 1)) {
            return 1; // this is a very boring PSF!!!
        }

        return 0;

    }

    public double sum() {

        int count = 0;

        sum = 0;

        if (array == null) {
            return Double.NaN;
        }
        
        for (double[][] i : array) {
            for (double[] j :  i) {
                for (double a : j) {
                    sum += a;
                    count++;
                }
            }
        }

        if (count > 0) {
            avg = sum / count;
        }

        //Master.log("PSF sum: " + sum);
        updateVectors();

        return sum;

    }

    public double sum(Geometry g) {

        int count = 0;

        if (array == null) {
            return Double.NaN;
        }

        if (g == null) {
            return Double.NaN;
        }

        sum = 0;

        for (int k = 0; k < array[0][0].length; k++) {
            for (int j = 0; j < array[0].length; j++) {
                for (int i = 0; i < array.length; i++) {
                    if (g.isSpace(i, j, k)) {
                        sum += array[i][j][k];
                        count++;
                    }
                }
            }
        }

        if (count > 0) {
            avg = sum / count;
        }

        //Master.log("PSF sum: " + sum);
        updateVectors();

        return sum;

    }

    public void findMax() {

        int ii = -1, jj = -1, kk = -1;

        if (array == null) {
            return;
        }

        max = 0;

        for (int k = 0; k < array[0][0].length; k++) {
            for (int j = 0; j < array[0].length; j++) {
                for (int i = 0; i < array.length; i++) {
                    if (array[i][j][k] > max) {
                        max = array[i][j][k];
                        ii = i;
                        jj = j;
                        kk = k;
                    }
                }
            }
        }

        Master.log("PSF max: " + max + " at " + ii + ", " + jj + ", " + kk);

    }

    public void ellipsoid(CoordinatesVoxels c, double value) {
        GeometryTools.ellipsoid(array, c, value);
    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (coordinates != null) {
            coordinates.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (coordinates != null) {
            addBlankParam();
            coordinates.createVector(true);
            addVector(coordinates.getVector());
            coordinates.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (coordinates != null) {
            coordinates.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            if (coordinates.setMyParams(o, v)) {
                array = null;
                updateCoordinates();
                return true;
            }
            return false;
        }

        if (!(o.paramVector instanceof PSF)) {
            return false;
        }

        String n = o.getName();

        if (arbitraryVoxelCenter) {
            if (n.equalsIgnoreCase("xVoxelCenter")) {
                if (v < 0) {
                    return false;
                }
                xVoxelCenter = v;
                array = null;
                return true;
            }
            if (n.equalsIgnoreCase("yVoxelCenter")) {
                if (v < 0) {
                    return false;
                }
                yVoxelCenter = v;
                array = null;
                return true;
            }
            if (n.equalsIgnoreCase("zVoxelCenter")) {
                if (v < 0) {
                    return false;
                }
                zVoxelCenter = v;
                array = null;
                return true;
            }
        }
        if (n.equalsIgnoreCase("rotateXY")) {
            rotateXY = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("rotateYZ")) {
            rotateYZ = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("rotateZX")) {
            rotateZX = v;
            array = null;
            return true;
        }
        if (n.equalsIgnoreCase("zCount")) {
            if (v < 2) {
                return false;
            }
            zCount = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("xySymmetric")) {
            xySymmetric = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("normalize")) {
            normalize = (v == 1);
            return true;
        }
        return super.setMyParams(o, v);
    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, String s) {

        if ((o == null) || (s == null)) {
            return false;
        }

        if ((o.paramVector instanceof CoordinatesVoxels) && (coordinates != null)) {
            return coordinates.setMyParams(o, s);
        }

        if (!(o.paramVector instanceof PSF)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        if (n.equalsIgnoreCase("psfFile")) {
            fileName = s;
            return true;
        }
        return super.setMyParams(o, s);
    }
}
