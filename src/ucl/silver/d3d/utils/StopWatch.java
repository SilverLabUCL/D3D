package ucl.silver.d3d.utils;

import ucl.silver.d3d.core.*;
import java.text.DecimalFormat;

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
public class StopWatch {
    
    public long milliseconds = 0;

    private long start = 0;

    public boolean isRunning = false;

    double msecRemaining = 0; // msecs

    public long iTimer, iTimerMax;
    public int iTimerFirstNum = 200; // gives reasonable estimation of simulation "time remaining"
    public boolean iTimerFirst;

    public String status;

    public DecimalFormat decimal3 = new DecimalFormat("0.000");

    public StopWatch() {
    }

    public void start() {
        start = System.currentTimeMillis();
        milliseconds = 0;
        isRunning = true;
    }

    public void stop() {
        
        long stop = System.currentTimeMillis();

        if (isRunning) {
            milliseconds = stop - start;
        } else {
            milliseconds = 0;
        }

        isRunning = false;

    }

    public void timer(double time) {

        if (!isRunning) {

            if (Master.project.printRate <= 0) {
                return;
            }

            iTimerMax = (int) (1.0 / (Master.project.printRate * Master.project.dt));
            iTimerMax = Math.max(iTimerMax, 1); // at least one time step

            iTimer = 0;
            iTimerFirst = true;

            start();

            return;

        }

        iTimer++;

        if ((iTimerFirst && (iTimer == iTimerFirstNum)) || (iTimer == iTimerMax)) {

            stop();

            msecRemaining = milliseconds * (Master.project.simTime - time) / (iTimer * Master.project.dt); // seconds

            status = "sim t: " + decimal3.format(time) + " ms,"
                    + " " + iTimer + "*dt: " + decimal3.format(milliseconds / 1000.0) + " sec,"
                    + " time remaining: " + stopWatch((long) msecRemaining);

            //status += ", iTimer = " + iTimer;

            Master.log(status);

            iTimer = 0;
            iTimerFirst = false;

            start();

        }

    }

    static public String stopWatch(long msec) {

        double seconds = (msec * 1.0) / 1000.0;

        int HOURS = (int) Math.floor(seconds / 3600.0);

        seconds = seconds - HOURS * 3600.0;

        int MINUTES = (int) Math.floor(seconds / 60.0);

        seconds = seconds - MINUTES * 60.0;

        int SECONDS = (int) Math.floor(seconds);

        return "" + String.format("%02d", HOURS) + ":" + String.format("%02d", MINUTES) + ":" + String.format("%02d", SECONDS);

    }

    @Override
    public String toString() {
        //return Long.toString(milliseconds) + " msec";
        return stopWatch(milliseconds);
    }

}
