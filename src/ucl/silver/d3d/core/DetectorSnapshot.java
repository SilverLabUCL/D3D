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
public class DetectorSnapshot extends Detector {

    public DetectorSnapshot(Project p, String NAME, int DiffusantNum, CoordinatesVoxels c, double time) {

        super(p, NAME, DiffusantNum, c);

        pulseTimer = new PulseTimer(p, time, 0);
        color = new ColorD3D(name + "_color", new Color(153, 204, 51));

        save.dataPoints = 1;

        createVector(true);

    }

    @Override
    public void init() {

        super.init();

        if (save != null) {
            save.samples2save = 0;
            save.skipSamples = 0;
            save.sampleRate = -1;
            save.sampleInterval = 0;
            save.saveWhileComputing = true;
            save.updateVectors();
        }

    }

    public boolean setFileName(int usec) {

        String dname = "";

        Diffusant[] diffusants = project.diffusants;

        if (save == null) {
            return false;
        }

        if ((diffusantNum >= 0) && (diffusantNum < diffusants.length)) {
            dname = diffusants[diffusantNum].name;
        }

        save.fileName(name, dname);

        if (usec >= 0) {
            save.outputFile += "_" + Integer.toString(usec);
        }

        return true;

    }

    @Override
    public void saveDimensions() {

        if ((save == null) || (!save.autoDimensions)) {
            return;
        }

        save.xdim = "ZYX voxels";

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            save.ydim = project.concUnits;
        } else {
            save.ydim = diffusantName + " (" + project.concUnits + ")";
        }

    }

    public boolean addTime(double time) {

        if (pulseTimer == null) {
            return false;
        }

        pulseTimer.add(time, 0);

        return true;

    }

    @Override
    public void detect(RunFiniteDifference fd, Geometry geometry) {

        int xVoxel1, yVoxel1, zVoxel1, xVoxel2, yVoxel2, zVoxel2;
        double a;

        int dataPoints = 1;

        if ((fd.diffus == null) || (fd.it >= fd.itmax)) {
            return;
        }

        if ((pulseTimer == null) || (pulseTimer.timer == null)) {
            return; // requires pulse timer
        }

        if (pulseTimer.timer[fd.it] > 0) {

            if (!setFileName((int) (1000 * fd.time))) {
                return;
            }

            if (!save.init(name, coordinates(), fd.time, dataPoints)) {
                return;
            }

            if (fd.diffus[0].length == 1) {

                save.saveData(fd.diffus[diffusantNum][0]); // single compartment

            } else {

                xVoxel1 = coordinates().xVoxel1;
                yVoxel1 = coordinates().yVoxel1;
                zVoxel1 = coordinates().zVoxel1;
                xVoxel2 = coordinates().xVoxel2;
                yVoxel2 = coordinates().yVoxel2;
                zVoxel2 = coordinates().zVoxel2;

                for (int k = zVoxel1; k <= zVoxel2; k++) {
                    for (int j = yVoxel1; j <= yVoxel2; j++) {
                        for (int i = xVoxel1; i <= xVoxel2; i++) {

                            if (geometry.isSpace(i, j, k)) {
                                a = fd.diffus[diffusantNum][geometry.space[i][j][k]];
                            } else {
                                a = -1;
                            }

                            save.saveData(a);

                        }
                    }
                }

            }

            save.finish(name, coordinates(), fd.time);

            Master.log("Snapshot taken at " + fd.time + " " + project.timeUnits);

        }

    }

    @Override
    public boolean saveInit() {
        return true; // nothing to do
    }

    @Override
    public boolean saveFinish() {
        return true; // nothing to do
    }
}
