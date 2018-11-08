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
public class RunMonteCarloPhotoMito extends RunMonteCarloPhoto {
    //
    // here mitochondria are defined by coordinates instead of non-space voxels
    // this allows simulations with drift
    //
    public Coordinates[] mito = null; // mito coordinates
    public boolean mitoCoordinatesOn = false;
    public double mitoRadius = (0.25 / 0.89) / 2.0; // um Palay
    public double mitoAxialRatio = 2.0 / 0.25; // Palay
    public double setMitoVolumeFraction = 0.28; // Zoltan average
    public double mitoVolumeFraction = 0;
    public double mitoVolume = 0; // um^3
    //public String mitoShape = "cuboid";
    public String mitoShape = "cylinder";
    //public String mitoShape = "ellipsoid";
    public int mito_ijkselect = -2;
    
    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("mitoRadius")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("mitoVolume")) {
            return project.spaceUnits + "^3";
        }
        return super.units(name);
    }
    
    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("mitoVolumeFraction")) {
            return false;
        }
        if (name.equalsIgnoreCase("mitoVolume")) {
            return false;
        }
        return super.canEdit(name);
    }
    
    public RunMonteCarloPhotoMito(Project p) {
        super(p);
        createVector(true);
    }
    
    @Override
    public double spaceVolume() {
        return project.geometry.spaceVolume - mitoVolume;
    }
    
    @Override
    public boolean initGeometry() {
        return initMitochondria();
    }
    
    public boolean initMitochondria() {

        int count = 0, iMitoTrials = 20, itrialMax = 100, radiusTrials = 2000;
        double d, x = 0, y = 0, z = 0;
        double radius, rstep, vf, volume;
        boolean overlap, finished = false;

        mitoVolume = 0;

        if (!mitoCoordinatesOn || (setMitoVolumeFraction <= 0)) {
            return false;
        }

        if (mito != null) {

            for (Coordinates m : mito) {
                mitoVolume += m.geometryVolume();
            }

            mitoVolumeFraction = mitoVolume / geometry.volume;

            Master.log("final mito volume fraction = " + mitoVolumeFraction);

            return false;

        }

        Coordinates test = new Coordinates(project);
        Coordinates[] temp = new Coordinates[iMitoTrials];

        radius = mitoRadius;
        rstep = radius / radiusTrials;

        if (mito_ijkselect == -2) {

            d = mt.nextDouble();

            if (d < 0.3333333333) {
                mito_ijkselect = 0; // xy plane
            } else if (d < 0.6666666666) {
                mito_ijkselect = 1; // yz plane
            } else {
                mito_ijkselect = 2; // zx plane
            }

        }

        Master.log("initializing mitochondria...");

        for (int imito = 0; imito < iMitoTrials; imito++) {

            overlap = true;

            for (int itrial = 0; itrial < itrialMax; itrial++) {

                x = geometry.x1 + mt.nextDouble() * (geometry.x2 - geometry.x1);
                y = geometry.y1 + mt.nextDouble() * (geometry.y2 - geometry.y1);
                z = geometry.z1 + mt.nextDouble() * (geometry.z2 - geometry.z1);

                if (mito_ijkselect == -1) {

                    d = mt.nextDouble(); // random orientation

                    if (d < 0.3333333333) {
                        mito_ijkselect = 0;
                    } else if (d < 0.6666666666) {
                        mito_ijkselect = 1;
                    } else {
                        mito_ijkselect = 2;
                    }

                }

                setMitoPosition(test, x, y, z, radius, mitoAxialRatio, mito_ijkselect);

                overlap = false;

                for (int jmito = 0; jmito < imito; jmito++) {
                    if (temp[jmito].intersectsCuboid(test)) {
                        overlap = true;
                        break;
                    }
                }

                if (!overlap) {
                    break;
                }

            }

            if (overlap) {
                Master.exit("MonteCarlo: initMitochondria : failed to place mitochondria #" + imito);
            }

            volume = test.geometryVolume();

            vf = (mitoVolume + volume) / geometry.volume;

            if (vf == setMitoVolumeFraction) {

                finished = true;

            } else if (vf > setMitoVolumeFraction) {

                finished = false;

                for (int i = 1; i < radiusTrials; i++) {

                    radius -= i * rstep;

                    if (radius <= 0) {
                        break;
                    }

                    setMitoPosition(test, x, y, z, radius, mitoAxialRatio, mito_ijkselect);

                    volume = test.geometryVolume();

                    vf = (mitoVolume + volume) / geometry.volume;

                    if (vf <= setMitoVolumeFraction) {
                        finished = true;
                        break;
                    }

                }

                if (!finished) {
                    Master.exit("MonteCarlo: initMitochondria : failed to place mitochondria #" + imito);
                }

            }

            temp[imito] = new Coordinates(project, test);
            temp[imito].updateDimensions();
            temp[imito].name = "mito #" + imito;
            temp[imito].color.setColor("[r=155,g=155,b=155]");

            mitoVolume += volume;
            count++;

            Master.log("mito #" + imito + " radius = " + radius);
            //Master.log("mito shape " + temp[imito].shape);
            //Master.log("mito shapeNum " + temp[imito].shapeNum);
            //Master.log("mitoCoordinatesOn % volume = " + vf * 100.0);

            if (finished) {
                break;
            }

        }

        if (count == 0) {
            return true;
        }

        mito = new Coordinates[count];

        System.arraycopy(temp, 0, mito, 0, count);

        mitoVolumeFraction = mitoVolume / geometry.volume;

        Master.log("final mito volume fraction = " + mitoVolumeFraction);

        return false;

    }
    
    public void setMitoPosition(Coordinates c, double x, double y, double z, double radius, double mitoAxialRatio, int ijkSelect) {

        String xyz = "";

        double halfLength = 0.5 * 2.0 * radius * mitoAxialRatio;

        switch (ijkSelect) {
            case 0:
                c.x1 = x - halfLength; // geometry.x1;
                c.y1 = y - radius;
                c.z1 = z - radius;
                c.x2 = x + halfLength; // geometry.x2;
                c.y2 = y + radius;
                c.z2 = z + radius;
                xyz = "x";
                break;
            case 1:
                c.x1 = x - radius;
                c.y1 = y - halfLength; // geometry.y1;
                c.z1 = z - radius;
                c.x2 = x + radius;
                c.y2 = y + halfLength; // geometry.y2;
                c.z2 = z + radius;
                xyz = "y";
                break;
            case 2: // xy
                c.x1 = x - radius;
                c.y1 = y - radius;
                c.z1 = z - halfLength; // geometry.z1;
                c.x2 = x + radius;
                c.y2 = y + radius;
                c.z2 = z + halfLength; // geometry.z2;
                xyz = "z";
                break;
        }

        if (mitoShape.equalsIgnoreCase("cuboid")) {
            c.setShape("cuboid");
        } else if (mitoShape.equalsIgnoreCase("ellipsoid")) {
            c.setShape("ellipsoid");
        } else if (mitoShape.equalsIgnoreCase("cylinder")) {
            c.setShape("cylinder" + xyz);
        }

        c.updateDimensions();

    }
    
    public double mitoHydrodynamics() {

        double dx, dy, dz, h, radiusMito;
        double b1avg = 0, b2avg = 0, b12avg, count = 0;

        if ((mito == null) || (diffusants == null)) {
            return Double.NaN;
        }

        for (Coordinates m : mito) {

            if (m == null) {
                break;
            }

            for (DiffusantVesicles d : diffusants) {

                if (d == null) {
                    continue;
                }

                for (DiffusantVesicle v : d.vesicles) {

                    if (v == null) {
                        continue;
                    }

                    if (!v.insideGeometry) {
                        continue;
                    }

                    if (m.isInside(v.x, v.y, v.z, v.radius)) {
                        continue;
                    }

                    switch (m.shapeNum) {

                        case 2: // cylinderx
                            dy = v.y - m.yCenter;
                            dz = v.z - m.zCenter;
                            h = Math.sqrt(dy * dy + dz * dz);
                            radiusMito = (m.y2 - m.y1) * 0.5;
                            break;

                        case 3: // cylindery
                            dx = v.x - m.xCenter;
                            dz = v.z - m.zCenter;
                            h = Math.sqrt(dx * dx + dz * dz);
                            radiusMito = (m.x2 - m.x1) * 0.5;
                            break;
                        case 4: // cylinderz
                            dx = v.x - m.xCenter;
                            dy = v.y - m.yCenter;
                            h = Math.sqrt(dx * dx + dy * dy);
                            radiusMito = (m.x2 - m.x1) * 0.5;
                            break;
                        default:
                            return Double.NaN;

                    }

                    h -= radiusMito;

                    b1avg += hydroWall_ll(v.radius, h);
                    b2avg += hydroWall_T(v.radius, Math.abs(h - v.radius));
                    count++;

                }

            }

        }

        b1avg /= count;
        b2avg /= count;
        b12avg = b1avg * b2avg;

        Master.log("avg mito hydrodynamic b1 = " + b1avg);
        Master.log("avg mito hydrodynamic b2 = " + b2avg);
        Master.log("b1 * b2 = " + b12avg);

        return b12avg;

    }
    
    @Override
    public boolean testOverlap(DiffusantVesicle v) {
        
        if (mito == null) {
            return false;
        }

        for (Coordinates m : mito) {
            if (m == null) {
                break;
            }
            if (m.isInside(v.x, v.y, v.z, v.radius)) {
                return true;
            }
        }

        return false;

    }
    
    @Override
    public void driftGeometry(double dx, double dy, double dz) {
        double x, y, z;
        if (mito != null) {
            for (Coordinates m : mito) {
                x = m.xCenter + dx;
                y = m.yCenter + dy;
                z = m.zCenter + dz;
                m.setCenter(x, y, z);
            }
        }
    }
    
    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (mito != null) {

            addBlankParam();

            for (Coordinates m : mito) {
                m.createVector(true);
                addVector(m.getVector());
                m.addUser(this);
            }

        }

        if (close) {
            closeVector();
        }

        return true;

    }
    
}
