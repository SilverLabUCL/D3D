package ucl.silver.d3d.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.text.DecimalFormat;
import ucl.silver.d3d.core.*;

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
final public class Grid
    extends JPanel implements MouseListener, MouseMotionListener {

    private int mx = -1, my = -1; // mouse xy holders
    private int mx0, my0, mx1, my1; // mouse marquee xy holders
    private double iVoxel = -1, jVoxel = -1; // mouse cursor voxel coordinates in ij-plane
    public int kVoxel = 0; // k-plane number to display

    public int displayMode = 4; // (0) shape space (1) diffusants (2) sources (3) detectors (4) All
    public int displayModeArrayNum = -1;
    public int editMode = 0; // (0) no edit (1) non-space (2) space (3) marquee
    public int axesMode = 0; // (0) xy (1) yz (2) zx
    private int copyPlane = -1;
    
    public int pixelsPerVoxel = 15;
    public int pixelsPerVoxelMin = 5;
    public int pixelsPerVoxelMax = 80;

    public double[][] diffusant = null; // for Preview animation display

    private boolean marquee = false, marqueeEdit = false;
    public boolean ijSwitch = false; // switch i and j axes
    public boolean preview = false;
    public boolean cross = false; // mid-line cross
    public boolean sym4 = false; // quarter symmetry view mode
    public boolean voxelGrid = true;
    public boolean frame = true;
    public boolean xyzMarkers = false;
    public boolean redParticlesPBC = true;

    public boolean coordVoxels = false; // (false) um (true) voxels
    public boolean coordClick = false; // (false) update always (true) update on mouse click

    Cursor defaultCursor = new Cursor(Cursor.DEFAULT_CURSOR);
    Cursor crossCursor = new Cursor(Cursor.CROSSHAIR_CURSOR);

    Color psfCenterColor = Color.RED;

    //private String colorScaleStr = "Heat";

    DecimalFormat decimal3 = new DecimalFormat("0.000");

    private Panel2D panel2D = null;
    
    private double phi_Keiding = 0;
    //private double phi_Keiding = 40;

    // class constructor
    public Grid(int w, int h, Panel2D p2D) {
        initVoxelWidth(w, h);
        setGridSize(w, h);
        addMouseListener(this);
        addMouseMotionListener(this);
        panel2D = p2D;
    }

    public void adjustGridSize() {

        int w = getGridWidth();
        int h = getGridHeight();

        if (sym4) {
            setGridSize(w * 2, h * 2);
        } else {
            setGridSize(w, h);
        }

    }

    public void setGridSize(int w, int h) {
        setPreferredSize(new Dimension(w, h));
        revalidate();
        repaint();
    }

    // switch 2D voxelGrid axes
    public void switchAxesToggle() {
        ijSwitch = !ijSwitch;
        adjustGridSize();
    }

    // set square voxel width in pixels
    public void initVoxelWidth(int w, int h) {

        if (w < 0) {
            //w = Math.min(getWidth(), 500);
            w = Math.min(getPreferredSize().width, 500);
        }

        if (h < 0) {
            //h = Math.min(getHeight(), 500);
            h = Math.min(getPreferredSize().height, 500);
        }

        pixelsPerVoxel = Math.min(h / iVoxels(), h / jVoxels());
        pixelsPerVoxel = Math.min(pixelsPerVoxel, h / kVoxels());

        pixelsPerVoxel = Math.min(pixelsPerVoxel, pixelsPerVoxelMax);
        pixelsPerVoxel = Math.max(pixelsPerVoxel, pixelsPerVoxelMin);

    }

    // set square voxel width in pixels
    public void setVoxelWidth(int w) {

        if ((w == pixelsPerVoxel) || (w < 2) || (w > 200)) {
            return;
        }

        pixelsPerVoxel = w;
        adjustGridSize();

    }

    public void setCross(boolean on) {
        cross = on;
        repaint();
    }

    public void setVoxelGrid(boolean on) {
        voxelGrid = on;
        repaint();
    }

    public void setSym4(boolean on) {
        sym4 = on;
        adjustGridSize();
        repaint();
    }

    // set new 2D voxelGrid k-kVoxel
    public int setGridPlane(int k) {

        if ((k != kVoxel) && (k >= 0) && (k < kVoxels())) {
            kVoxel = k;
            repaint();
        }

        return kVoxel;

    }

    // set new 2D display mode (0) xy plane (1) yz plane (2) zx plane
    public int setAxesMode(int iMode) {

        if ((iMode != axesMode) && (iMode >= 0) && (iMode <= 2)) {
            axesMode = iMode;
            marquee = false;
            copyPlane = -1;
            adjustGridSize();
        }

        return axesMode;

    }

    // set edit drawing mode (0) no edit (1) non-space (2) space (3) marquee
    public int setEditMode(int iMode) {

        if ((iMode != editMode) && (iMode >= 0) && (iMode <= 3)) {
            editMode = iMode;
            marquee = false;
            repaint();
        }

        return editMode;

    }

    // set view mode (0) shape space (1) diffusants (2) sources (3) detectors (4) All
    public int setViewMode(int iMode) {

        if ((iMode != displayMode) && (iMode >= 0) && (iMode <= 4)) {
            displayMode = iMode;
            marquee = false;
            repaint();
        }

        return displayMode;

    }

    public int setViewArrayNum(int i) {

        switch (displayMode) {

            case 0: // shape space
                displayModeArrayNum = 0;
                break;

            case 1: // diffusants
                if ((Master.project.diffusants != null) && (i >= -1) &&
                        (i < Master.project.diffusants.length)) {
                    displayModeArrayNum = i;
                }
                break;

            case 2: // sources
                if ((Master.project.sources != null) && (i >= -1) &&
                        (i < Master.project.sources.length)) {
                    displayModeArrayNum = i;
                }
                break;

            case 3: // detectors
                if ((Master.project.detectors != null) && (i >= -1) &&
                        (i < Master.project.detectors.length)) {
                    displayModeArrayNum = i;
                }
                break;

            default:
                displayModeArrayNum = -1;

        }

        return displayModeArrayNum;

    }

    private int xyPixelOffset() {
        return 1 * pixelsPerVoxel;
    }

    public int getGridWidth() {

        int offset = xyPixelOffset();

        if (ijSwitch) {
            return offset + jVoxels() * pixelsPerVoxel;
        } else {
            return offset + iVoxels() * pixelsPerVoxel;
        }

    }

    public int getGridHeight() {

        int offset = xyPixelOffset();

        if (ijSwitch) {
            return offset + iVoxels() * pixelsPerVoxel;
        } else {
            return offset + jVoxels() * pixelsPerVoxel;
        }

    }

    // get max i-dim voxels
    public int iVoxels() {
        switch (axesMode) {
            case 0: // xy
                return Master.project.geometry.xVoxels;
            case 1: // yz
                return Master.project.geometry.yVoxels;
            case 2: // zx
                return Master.project.geometry.zVoxels;
            default:
                return -1;
        }
    }

    // get max j-dim voxels
    public int jVoxels() {
        switch (axesMode) {
            case 0: // xy
                return Master.project.geometry.yVoxels;
            case 1: // yz
                return Master.project.geometry.zVoxels;
            case 2: // zx
                return Master.project.geometry.xVoxels;
            default:
                return -1;
        }
    }

    // get max k-dim voxels
    public int kVoxels() {
        switch (axesMode) {
            case 0: // xy
                return Master.project.geometry.zVoxels;
            case 1: // yz
                return Master.project.geometry.xVoxels;
            case 2: // zx
                return Master.project.geometry.yVoxels;
            default:
                return -1;
        }
    }

    // convert pixels to voxelGrid voxel index
    public double pixels2Voxel(int pixels) {
        return ((pixels - xyPixelOffset()) * 1.0) / (pixelsPerVoxel * 1.0);
    }

    public double iPixels2Voxel(int iPixels, int jPixels) {
        if (ijSwitch) {
            return pixels2Voxel(jPixels);
        } else {
            return pixels2Voxel(iPixels);
        }
    }

    public double jPixels2Voxel(int iPixels, int jPixels) {
        if (ijSwitch) {
            return pixels2Voxel(iPixels);
        } else {
            return pixels2Voxel(jPixels);
        }
    }

    // convert voxelGrid voxel index to number of pixels
    public double voxel2Pixels(double voxel) {
        return voxel * pixelsPerVoxel + xyPixelOffset();
    }

    public double iVoxel2Pixels(double iVoxel, double jVoxel) {
        if (ijSwitch) {
            return voxel2Pixels(jVoxel);
        } else {
            return voxel2Pixels(iVoxel);
        }
    }

    public double jVoxel2Pixels(double iVoxel, double jVoxel) {
        if (ijSwitch) {
            return voxel2Pixels(iVoxel);
        } else {
            return voxel2Pixels(jVoxel);
        }
    }

    // determine if i-j points lie within 2D voxelGrid
    private boolean isGrid(double i, double j) {
        return ((i >= 0) && (i < iVoxels()) && (j >= 0) && (j < jVoxels()) &&
                (kVoxel >= 0) && (kVoxel < kVoxels()));
    }

    // determine if i-j points lie within 2D voxelGrid
    private boolean isInside(double i, double j) {

        if ((i >= 0) && (i < iVoxels()) && (j >= 0) && (j < jVoxels())) {

            switch (displayMode) {
                default:
                    if ((kVoxel > 0) && (kVoxel < kVoxels() - 1)) {
                        return true;
                    }
                case 1: // diffusants
                    if ((kVoxel >= 0) && (kVoxel <= kVoxels() - 1)) {
                        return true;
                    }
            }

        }

        return false;

    }

    // clear current kVoxel
    public void clearPlane() {
        
        int ivoxels = iVoxels();
        int jvoxels = jVoxels();

        for (int i = 0; i < ivoxels; i++) {
            for (int j = 0; j < jvoxels; j++) {
                switch (displayMode) {
                    case 0: // shape space
                        setSpace(i, j, kVoxel, 1);
                        break;
                }
            }
        }

        repaint();

    }

    // copy the current kVoxel
    public void copyPlane() {
        copyPlane = kVoxel;
    }

    // paste to the current kVoxel
    public void pastePlane() {

        if (kVoxel == copyPlane) {
            return;
        }
        
        int ivoxels = iVoxels();
        int jvoxels = jVoxels();

        for (int i = 0; i < ivoxels; i++) {
            for (int j = 0; j < jvoxels; j++) {
                setSpace(i, j, kVoxel, getSpace(i, j, copyPlane));
            }
        }

        repaint();

    }

    // save value to space matrix
    private void setSpace(int i, int j, int k, double value) {
        switch (axesMode) {
            case 0: // xy
                Master.project.geometry.setSpace(i, j, k, (byte) value);
                break;
            case 1: // yz
                Master.project.geometry.setSpace(k, i, j, (byte) value);
                break;
            case 2: // zx
                Master.project.geometry.setSpace(j, k, i, (byte) value);
                break;
        }
    }

    // get value from space matrix
    private int getSpace(int i, int j, int k) {
        switch (axesMode) {
            case 0: // xy
                return Master.project.geometry.getSpace(i, j, k);
            case 1: // yz
                return Master.project.geometry.getSpace(k, i, j);
            case 2: // zx
                return Master.project.geometry.getSpace(j, k, i);
        }
        return -1;
    }

    // get value from space matrix
    private double getSpacePSFweight(int i, int j, int k) {
        switch (axesMode) {
            case 0: // xy
                return Master.project.geometry.getSpacePSFweight(i, j, k);
            case 1: // yz
                return Master.project.geometry.getSpacePSFweight(k, i, j);
            case 2: // zx
                return Master.project.geometry.getSpacePSFweight(j, k, i);
        }
        return -1;
    }

    private double previewValue(int i, int j, int k) {

        int I;

        if (!isInside(i, j) || (diffusant == null) || (displayModeArrayNum < 0) || (displayModeArrayNum >= diffusant.length)) {
            return -1;
        }

        switch (axesMode) {
            case 0: // xy
                I = Master.project.geometry.space[i][j][k];
                return diffusant[displayModeArrayNum][I];
            case 1: // yz
                I = Master.project.geometry.space[k][i][j];
                return diffusant[displayModeArrayNum][I];
            case 2: // zx
                I = Master.project.geometry.space[j][k][i];
                return diffusant[displayModeArrayNum][I];
        }

        return -1;

    }

    private Color displayColorDiffusant(int i, int j, int k) {
        
        double concentration = -999;

        if (preview) {
            concentration = previewValue(i, j, k);
        }

        switch (axesMode) {
            case 0: // xy
                return displayColorDiffusant2(i, j, k, concentration);
            case 1: // yz
                return displayColorDiffusant2(k, i, j, concentration);
            case 2: // zx
                return displayColorDiffusant2(j, k, i, concentration);
        }

        return null;

    }

    private Color displayColorDiffusant2(int xVoxel, int yVoxel, int zVoxel, double concentration) {

        int icount = 0, iredAvg = 0, igreenAvg = 0, iblueAvg = 0;
        int ired, igreen, iblue;
        int idiffusant = displayModeArrayNum;

        Color c;

        Diffusant[] diffusants = Master.project.diffusants;

        if ((diffusants != null) && (idiffusant < diffusants.length)) {

            if (idiffusant < 0) {

                for (Diffusant d : diffusants) {

                    if (d == null) {
                        continue;
                    }

                    c = d.displayColor(xVoxel, yVoxel, zVoxel);

                    if (c == null) {
                        continue;
                    }

                    ired = c.getRed();
                    igreen = c.getGreen();
                    iblue = c.getBlue();

                    if ((ired < 0) || (igreen < 0) || (iblue < 0)) {
                        continue;
                    }

                    iredAvg += ired;
                    igreenAvg += igreen;
                    iblueAvg += iblue;
                    icount++;

                }

                if ((iredAvg == 0) && (igreenAvg == 0) && (iblueAvg == 0)) {

                    return null;

                } else {

                    iredAvg /= icount;
                    igreenAvg /= icount;
                    iblueAvg /= icount;

                    return new Color(iredAvg, igreenAvg, iblueAvg);

                }

            } else {
                if (concentration == -999) {
                    return diffusants[idiffusant].displayColor(xVoxel, yVoxel, zVoxel);
                } else {
                    return diffusants[idiffusant].displayColor(xVoxel, yVoxel, zVoxel, concentration); // preview
                }
            }
        }

        return null;

    }

    // get detectors number (-1, not within detectors)
    private Color displayColorSource(int i, int j, int k, double min, double max) {
        switch (axesMode) {
            case 0: // xy
                return displayColorSource2(i, j, k, min, max);
            case 1: // yz
                return displayColorSource2(k, i, j, min, max);
            case 2: // zx
                return displayColorSource2(j, k, i, min, max);
        }
        return null;
    }

    private Color displayColorSource2(int xVoxel, int yVoxel, int zVoxel, double min, double max) {

        int icount = 0, iredAvg = 0, igreenAvg = 0, iblueAvg = 0;
        int ired, igreen, iblue;
        int isource = displayModeArrayNum;

        Color c;

        Source[] sources = Master.project.sources;

        if ((sources == null) || (isource >= sources.length)) {
            return null;
        }

        if (isource < 0) {

            for (Source s : sources) {

                if ((s == null) || (s.color == null)) {
                    continue;
                }
                
                if (s.isInside(xVoxel, yVoxel, zVoxel)) {

                    //if (!Master.project.insideAnySource(xVoxel, yVoxel, zVoxel, is)) {
                        c = s.color.getColor(s.Ctotal, min, max);
                    //}

                } else {

                    continue;

                }

                if (c == null) {
                    continue;
                }

                ired = c.getRed();
                igreen = c.getGreen();
                iblue = c.getBlue();

                if ((ired < 0) || (igreen < 0) || (iblue < 0)) {
                    continue;
                }

                iredAvg += ired;
                igreenAvg += igreen;
                iblueAvg += iblue;
                icount++;

            }

            if ((iredAvg == 0) && (igreenAvg == 0) && (iblueAvg == 0)) {

                return null;

            } else {

                iredAvg /= icount;
                igreenAvg /= icount;
                iblueAvg /= icount;

                return new Color(iredAvg, igreenAvg, iblueAvg);

            }

        } else {
            
            if ((sources[isource] == null) || (sources[isource].color == null)) {
                return null;
            }
            
            if (sources[isource].isInside(xVoxel, yVoxel, zVoxel)) {

                //if (!Master.project.insideAnySource(xVoxel, yVoxel, zVoxel, isource)) {
                    return sources[isource].color.getColor(sources[isource].Ctotal, min, max);
                //}

            }

            return null;

        }

    }
    
    // get detectors number (-1, not within detectors)
    private Color displayColorDetector(int i, int j, int k, double min, double max) {
        switch (axesMode) {
            case 0: // xy
                return displayColorDetector2(i, j, k, min, max);
            case 1: // yz
                return displayColorDetector2(k, i, j, min, max);
            case 2: // zx
                return displayColorDetector2(j, k, i, min, max);
        }
        return null;
    }

    private Color displayColorDetector2(int xVoxel, int yVoxel, int zVoxel, double min, double max) {

        int icount = 0, iredAvg = 0, igreenAvg = 0, iblueAvg = 0;
        int ired, igreen, iblue;
        int idetector = displayModeArrayNum;

        Color c;

        Detector[] detectors = Master.project.detectors;

        if ((detectors != null) && (idetector < detectors.length)) {

            if (idetector < 0) {

                for (Detector d : detectors) {

                    if (d == null) {
                        continue;
                    }

                    c = d.displayColor(xVoxel, yVoxel, zVoxel, min, max);

                    if (c == null) {
                        continue;
                    }

                    ired = c.getRed();
                    igreen = c.getGreen();
                    iblue = c.getBlue();

                    if ((ired < 0) || (igreen < 0) || (iblue < 0)) {
                        continue;
                    }

                    iredAvg += ired;
                    igreenAvg += igreen;
                    iblueAvg += iblue;
                    icount++;

                }

                if ((iredAvg == 0) && (igreenAvg == 0) && (iblueAvg == 0)) {

                    return null;
                    
                } else {

                    iredAvg /= icount;
                    igreenAvg /= icount;
                    iblueAvg /= icount;

                    return new Color(iredAvg, igreenAvg, iblueAvg);
                    
                }

            } else {
                return detectors[idetector].displayColor(xVoxel, yVoxel, zVoxel, min, max);
            }
        }

        return null;

    }
    
    private double displayValue(int i, int j, int k) {
        switch (axesMode) {
            case 0: // xy
                return displayValue2(i, j, k);
            case 1: // yz
                return displayValue2(k, i, j);
            case 2: // zx
                return displayValue2(j, k, i);
        }
        return -1;
    }

    private double displayValue2(int xVoxel, int yVoxel, int zVoxel) {

        Diffusant[] diffusants = Master.project.diffusants;
        Source[] sources = Master.project.sources;
        Detector[] detectors = Master.project.detectors;

        switch(displayMode) {
            case 0: // shape
                return Master.project.geometry.displayValue(xVoxel, yVoxel, zVoxel);
            case 1: // diffusants
                if (displayModeArrayNum < 0) {
                    break;
                }
                if ((diffusants != null) && (displayModeArrayNum < diffusants.length)) {
                    return diffusants[displayModeArrayNum].displayValue(xVoxel, yVoxel, zVoxel);
                }
                break;
            case 2: // sources
                if (displayModeArrayNum < 0) {
                    break;
                }
                if ((sources != null) && (displayModeArrayNum < sources.length)) {
                    return sources[displayModeArrayNum].displayValue(xVoxel, yVoxel, zVoxel);
                }
                break;
            case 3: // detectors
                if (displayModeArrayNum < 0) {
                    break;
                }
                if ((detectors != null) && (displayModeArrayNum < detectors.length)) {
                    return detectors[displayModeArrayNum].displayValue(xVoxel, yVoxel, zVoxel);
                }
                break;
        }

        return Double.NaN;

    }

    public String getCurrentColor() {

        Diffusant[] diffusants = Master.project.diffusants;
        Source[] sources = Master.project.sources;
        Detector[] detectors = Master.project.detectors;

        switch(displayMode) {
            case 0: // shape
                switch (editMode) {
                    case 1: // non-space
                        return Master.project.geometry.colorNonSpace.colorRGB;
                    case 2: // space
                        return Master.project.geometry.colorSpace.colorRGB;
                }
                break;
            case 1: // diffusants
                if (displayModeArrayNum < 0) {
                    break;
                }
                if ((diffusants != null) && (displayModeArrayNum < diffusants.length)) {
                    return diffusants[displayModeArrayNum].color.colorRGB;
                }
                break;
            case 2: // sources
                if (displayModeArrayNum < 0) {
                    break;
                }
                if ((sources != null) && (displayModeArrayNum < sources.length)) {
                    return sources[displayModeArrayNum].color.colorRGB;
                }
                break;
            case 3: // detectors
                if (displayModeArrayNum < 0) {
                    break;
                }
                if ((detectors != null) && (displayModeArrayNum < detectors.length)) {
                    return detectors[displayModeArrayNum].color.colorRGB;
                }
                break;
        }

        return null;

    }

    public void setCurrentColor(String colorStr) {

        Diffusant[] diffusants = Master.project.diffusants;
        Source[] sources = Master.project.sources;
        Detector[] detectors = Master.project.detectors;

        if (colorStr == null) {
            return;
        }

        switch(displayMode) {
            case 0: // shape
                if (editMode == 1) {
                    Master.project.geometry.colorNonSpace.setColor(colorStr);
                } else if (editMode == 2) {
                    Master.project.geometry.colorSpace.setColor(colorStr);
                }
                break;
            case 1: // diffusants
                if ((diffusants != null) && (displayModeArrayNum < diffusants.length)) {
                    diffusants[displayModeArrayNum].color.setColor(colorStr);
                }
                break;
            case 2: // sources
                if ((sources != null) && (displayModeArrayNum < sources.length)) {
                    sources[displayModeArrayNum].color.setColor(colorStr);
                }
                break;
            case 3: // detectors
                if ((detectors != null) && (displayModeArrayNum < detectors.length)) {
                    detectors[displayModeArrayNum].color.setColor(colorStr);
                }
                break;
        }

    }

    private double getPSFI0() {
        switch (axesMode) {
            case 0: // xy
                return getPSFvoxelCenter(0);
            case 1: // yz
                return getPSFvoxelCenter(1);
            case 2: // zx
                return getPSFvoxelCenter(2);
        }
        return -1;
    }

    private double getPSFJ0() {
        switch (axesMode) {
            case 0: // xy
                return getPSFvoxelCenter(1);
            case 1: // yz
                return getPSFvoxelCenter(2);
            case 2: // zx
                return getPSFvoxelCenter(0);
        }
        return -1;
    }

    private double getPSFK0() {
        switch (axesMode) {
            case 0: // xy
                return getPSFvoxelCenter(2);
            case 1: // yz
                return getPSFvoxelCenter(0);
            case 2: // zx
                return getPSFvoxelCenter(1);
        }
        return -1;
    }

    private double getPSFvoxelCenter(int xyzSelect) {

        Diffusant[] diffusants = Master.project.diffusants;
        Detector[] detectors = Master.project.detectors;

        int i = displayModeArrayNum;

        switch(displayMode) {

            case 0: // shape
                break;

            case 1: // diffusants

                if (diffusants == null) {
                    return -1;
                }

                if ((i >= 0) && (i < diffusants.length) && (diffusants[i] != null)) {
                    
                    if (diffusants[i].psf != null) {
                        switch (xyzSelect) {
                            case 0:
                                return diffusants[i].psf.xVoxelCenter();
                            case 1:
                                return diffusants[i].psf.yVoxelCenter();
                            case 2:
                                return diffusants[i].psf.zVoxelCenter();
                        }
                    }
                }

                break;

            case 2: // sources
                break;

            case 3: // detectors

                if (detectors == null) {
                    return -1;
                }

                if ((i >= 0) && (i < detectors.length) && (detectors[i] != null)) {
                    if (detectors[i].psf != null) {
                        switch (xyzSelect) {
                            case 0:
                                return detectors[i].psf.xVoxelCenter();
                            case 1:
                                return detectors[i].psf.yVoxelCenter();
                            case 2:
                                return detectors[i].psf.zVoxelCenter();
                        }
                    }
                }

                break;

        }

        return -1;

    }

    private double getI(double x, double y, double z) {
        switch (axesMode) {
            case 0: // xy
                return x;
            case 1: // yz
                return y;
            case 2: // zx
                return z;
        }
        return -1;
    }

    private double getJ(double x, double y, double z) {
        switch (axesMode) {
            case 0: // xy
                return y;
            case 1: // yz
                return z;
            case 2: // zx
                return x;
        }
        return -1;
    }

    private double getK(double x, double y, double z) {
        switch (axesMode) {
            case 0: // xy
                return z;
            case 1: // yz
                return x;
            case 2: // zx
                return y;
        }
        return -1;
    }

    // get voxel value from mouse xy pixel position
    private double getValue(int xPixels, int yPixels) {

        int i = (int) iPixels2Voxel(xPixels, yPixels);
        int j = (int) jPixels2Voxel(xPixels, yPixels);

        if (!isGrid(i, j)) {
            return -1;
        }

        switch (displayMode) {
            case 0: // shape space
                return getSpace(i, j, kVoxel);
            default:
                return -1;
        }

    }

    // save value to appropriate matrix
    private void setValue(double i, double j, double value) {

        if (!isInside(i, j)) {
            return;
        }

        switch (displayMode) {
            case 0: // shape space
                if (editMode > 0) {
                    setSpace((int) i, (int) j, kVoxel, value);
                }
                break;
            case 1: // diffusants
                break;
            case 2: // sources
                break;
            case 3: // detectors
                break;
        }

    }

    // set matrix values around center point with radius R
    private void setValueRadius(int xPixels, int yPixels, double value) {

        int R = 5; // pixel radius
        int i0 = (int) iPixels2Voxel(xPixels - R, yPixels - R);
        int j0 = (int) jPixels2Voxel(xPixels - R, yPixels - R);
        int i1 = (int) iPixels2Voxel(xPixels + R, yPixels + R);
        int j1 = (int) jPixels2Voxel(xPixels + R, yPixels + R);

        for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
                setValue(i, j, value);
            }
        }

    }

    private void setVoxelSelect(int xPixels, int yPixels) {

        if ((xPixels >= 0) && (yPixels >= 0)) {
            iVoxel = iPixels2Voxel(xPixels, yPixels);
            jVoxel = jPixels2Voxel(xPixels, yPixels);
        } else {
            iVoxel = -1;
            jVoxel = -1;
        }

        repaint();

    }

    // sets square marquee values based on mouse click state
    private void setMarquee(double value) {

        int xmin = Math.min(mx0, mx1);
        int ymin = Math.min(my0, my1);
        int xmax = Math.max(mx0, mx1);
        int ymax = Math.max(my0, my1);

        int i0 = (int) iPixels2Voxel(xmin, ymin);
        int j0 = (int) jPixels2Voxel(xmin, ymin);
        int i1 = (int) iPixels2Voxel(xmax, ymax);
        int j1 = (int) jPixels2Voxel(xmax, ymax);

        i0 = Math.max(i0, 0);
        j0 = Math.max(j0, 0);
        i1 = Math.min(i1, iVoxels() - 1);
        j1 = Math.min(j1, jVoxels() - 1);

        if (displayMode != 1) {
            for (int i = i0; i <= i1; i++) {
                for (int j = j0; j <= j1; j++) {
                    setValue(i, j, value);
                }
            }
        }

        repaint();

    }

    // get the appropriate painting tool value
    private double getEditValue() {

        switch (editMode) {

            case 1: // non-space
            case 3: // marquee

                switch (displayMode) {
                    case 0: // shape
                        return -1;
                    case 1: // diffusants
                        return 1;
                    case 2: // sources
                        return 0; // get spinner conc value here
                    case 3: // detectors
                    case 4: // All
                }

                return 0;

            case 2: // space

                switch (displayMode) {
                    case 0: // shape
                        return 0;
                    case 1: // diffusants
                    case 2: // sources
                    case 3: // detectors
                        break;
                }

                return 0;

        }

        return 0;

    }

    public void updateCursorStats() {

        double iPosition, jPosition, kPosition;
        double voxelIntensity;

        double dxHalf = Master.project.dx / 2.0;

        String xstr = "", ystr = "", zstr = "", vstr;
        String istr, jstr, kstr;

        if (!isGrid(iVoxel, jVoxel)) {
            //return;
        }

        if (coordVoxels) {

            istr = Integer.toString((int) iVoxel);
            jstr = Integer.toString((int) jVoxel);
            kstr = Integer.toString(kVoxel);

            switch (axesMode) {
                case 0:
                    xstr = "x = " + istr;
                    ystr = "y = " + jstr;
                    zstr = "z = " + kstr;
                    break;
                case 1:
                    xstr = "x = " + kstr;
                    ystr = "y = " + istr;
                    zstr = "z = " + jstr;
                    break;
                case 2:
                    xstr = "x = " + jstr;
                    ystr = "y = " + kstr;
                    zstr = "z = " + istr;
                    break;
            }

        } else {

            switch (axesMode) {
                case 0:

                    iPosition = Master.project.geometry.computeX(iVoxel);
                    jPosition = Master.project.geometry.computeY(jVoxel);
                    kPosition = Master.project.geometry.computeZ(kVoxel);

                    xstr = "x = " + decimal3.format(iPosition - dxHalf);
                    ystr = "y = " + decimal3.format(jPosition - dxHalf);
                    zstr = "z = " + decimal3.format(kPosition);

                    break;

                case 1:

                    iPosition = Master.project.geometry.computeY(iVoxel);
                    jPosition = Master.project.geometry.computeZ(jVoxel);
                    kPosition = Master.project.geometry.computeX(kVoxel);

                    xstr = "x = " + decimal3.format(kPosition);
                    ystr = "y = " + decimal3.format(iPosition - dxHalf);
                    zstr = "z = " + decimal3.format(jPosition - dxHalf);

                    break;

                case 2:

                    iPosition = Master.project.geometry.computeZ(iVoxel);
                    jPosition = Master.project.geometry.computeX(jVoxel);
                    kPosition = Master.project.geometry.computeY(kVoxel);

                    xstr = "x = " + decimal3.format(jPosition - dxHalf);
                    ystr = "y = " + decimal3.format(kPosition);
                    zstr = "z = " + decimal3.format(iPosition - dxHalf);

                    break;

            }

        }

        voxelIntensity = displayValue((int) iVoxel, (int) jVoxel, kVoxel);

        if (Double.isNaN(voxelIntensity)) {
            vstr = "";
        } else {
            vstr = "v = " + decimal3.format(voxelIntensity);
        }

        panel2D.coordinatesSet(xstr, ystr, zstr, vstr);

    }

    @Override
    public void paintComponent(Graphics g) {

        Color nonSpaceColor = Master.project.geometry.colorNonSpace.color;

        g.setColor(nonSpaceColor);
        g.fillRect(0, 0, getWidth(), getHeight());

        if ((kVoxel < 0) || (kVoxel >= kVoxels())) {
            g.setColor(Color.red);
            g.drawString("Plane Off Scale", 50, 50);
            return;
        }

        paintVoxelGrid(g, false, false);

        if (sym4) {
            paintVoxelGrid(g, true, false);
            paintVoxelGrid(g, false, true);
            paintVoxelGrid(g, true, true);
        }

        if (!preview) {
            paintVoxelNumbers(g);
            paintPSFcenter(g);
            paintVoxelMarquee(g);
            paintVoxelCrossAndFrame(g);
        }

        paintMonteCarlo(g);

    }

    private void paintVoxelGrid(Graphics g, boolean isym, boolean jsym) {

        int x0, y0;
        int itemp, jtemp;
        int ivoxels = iVoxels();
        int jvoxels = jVoxels();
        int space;
        
        double maxC0 = 0;
        
        boolean error;

        Color spaceColor = Master.project.geometry.colorSpace.color;
        Color nonSpaceColor = Master.project.geometry.colorNonSpace.color;
        Color tempColor1 = null, tempColor2;

        Source[] sources = Master.project.sources;
        Detector[] detectors = Master.project.detectors;
        
        if (Master.project.finiteDifference != null) {
            maxC0 = Master.project.finiteDifference.maxC0;
        }

        if (maxC0 == 0) {
            maxC0 = 1;
        }

        for (int i = 0; i < ivoxels; i++) {
            for (int j = 0; j < jvoxels; j++) {

                error = false;

                x0 = (int) iVoxel2Pixels(i, j);
                y0 = (int) jVoxel2Pixels(i, j);

                if (isym) {

                    itemp = Math.abs(i - (ivoxels - 1)) + ivoxels;
                    jtemp = Math.abs(j - (jvoxels - 1)) + jvoxels;

                    x0 = (int) iVoxel2Pixels(itemp, jtemp);

                }

                if (jsym) {

                    itemp = Math.abs(i - (ivoxels - 1)) + ivoxels;
                    jtemp = Math.abs(j - (jvoxels - 1)) + jvoxels;

                    y0 = (int) jVoxel2Pixels(itemp, jtemp);

                }

                space = getSpace(i, j, kVoxel);
                //intensity = getSpacePSFweight(i, j, kVoxel);

                if (space < 0) {
                    g.setColor(nonSpaceColor);
                } else {
                    g.setColor(spaceColor);
                    //geometry.setColor(color.getColor(intensity, 0, 1));
                }

                if (preview) {

                    if (space < 0) {
                        g.setColor(Color.white);
                    } else {
                        g.setColor(displayColorDiffusant(i, j, kVoxel));
                    }

                } else {

                    if ((space >= 0) && ((displayMode == 1) || (displayMode == 4))) { // diffusants

                        tempColor1 = displayColorDiffusant(i, j, kVoxel);

                        if (tempColor1 != null) {
                            g.setColor(tempColor1);
                        }

                    }

                    //if ((space >= 0) && (sources != null) && ((displayMode == 2) || (displayMode == 4))) { // sources
                    if ((sources != null) && ((displayMode == 2) || (displayMode == 4))) { // sources
                        
                        tempColor2 = displayColorSource(i, j, kVoxel, 0, maxC0);

                        if (tempColor2 != null) {

                            if (space < 0) {
                                if (!sym4) {
                                    error = true;
                                }
                            } else {
                                
                                if ((displayMode == 4) && (tempColor1 != null)) {
                                    tempColor1 = avgColors(tempColor1, tempColor2);
                                } else {
                                    tempColor1 = tempColor2;
                                }
                                
                                g.setColor(tempColor1);

                            }

                        }

                    }

                    //if ((space >= 0) && (detectors != null) && ((displayMode == 3) || (displayMode == 4))) { // detectorsif ((space >= 0) && (detectors != null) && ((displayMode == 3) || (displayMode == 4))) { // detectorsif ((space >= 0) && (detectors != null) && ((displayMode == 3) || (displayMode == 4))) { // detectors
                    if ((detectors != null) && ((displayMode == 3) || (displayMode == 4))) { // detectors    
                        
                        tempColor2 = displayColorDetector(i, j, kVoxel, 0, 1);

                        if (tempColor2 != null) {

                            if ((displayMode == 4) && (tempColor1 != null)) {
                                tempColor1 = avgColors(tempColor1, tempColor2);
                            } else {
                                tempColor1 = tempColor2;
                            }

                            g.setColor(tempColor1);

                        }

                    }

                }

                if (error) {
                    g.setColor(Color.red);
                }

                if (preview) {
                    g.fillRect(x0, y0, pixelsPerVoxel, pixelsPerVoxel);
                } else {
                    if (voxelGrid) {
                        g.fillRect(x0 + 1, y0 + 1, pixelsPerVoxel - 1, pixelsPerVoxel - 1);
                        //g.fillRect(x0 + 2, y0 + 2, pixelsPerVoxel - 2, pixelsPerVoxel - 2);
                    } else {
                        g.fillRect(x0 + 1, y0 + 1, pixelsPerVoxel, pixelsPerVoxel);
                    }
                }

            }
        }

    }

    private Color avgColors(Color color1, Color color2) {

        int ired = (int) ((color1.getRed() + color2.getRed()) / 2.0);
        int igreen = (int) ((color1.getGreen() + color2.getGreen()) / 2.0);
        int iblue = (int) ((color1.getBlue() + color2.getBlue()) / 2.0);

        return new Color(ired, igreen, iblue);

    }

    private void paintVoxelNumbers(Graphics g) {

        int id, x0, y0, ivox, jvox;

        if (!xyzMarkers) {
            return;
        }

        if (pixelsPerVoxel < 10) {
            return;
        }

        String dstra = "", dstrb = "", dtemp;

        g.setColor(Color.black);

        switch (axesMode) {
            case 0:
                dstra = "x";
                dstrb = "y";
                break;
            case 1:
                dstra = "y";
                dstrb = "z";
                break;
            case 2:
                dstra = "z";
                dstrb = "x";
                break;
        }

        if (ijSwitch) {
            dtemp = dstra;
            dstra = dstrb;
            dstrb = dtemp;
            jvox = iVoxels();
            ivox = jVoxels();
        } else {
            ivox = iVoxels();
            jvox = jVoxels();
        }

        //if (pixelsPerVoxel < 15) {
        //    dstra = "";
        //    dstrb = "";
        //}

        // horizontal

        y0 = Math.max(0, pixelsPerVoxel - 2);

        for (int i = 0; i < ivox; i++) {

            id = (int) Math.IEEEremainder(i, 10);

            if (id < 0) {
                id += 10;
            }

            x0 = pixelsPerVoxel * i + xyPixelOffset() + pixelsPerVoxel / 3;
            
            //geometry.drawString(Integer.toString(id), i0, j0);
            g.drawString(dstra, x0, y0);

        }
        
        //x0 = pixelsPerVoxel * ivox + xyPixelOffset() + pixelsPerVoxel / 3;

        //geometry.drawString(dstra, i0, j0);

        // vertical

        x0 = Math.max( 0, pixelsPerVoxel - 13);

        for (int j = 0; j < jvox; j++) {

            id = (int) Math.IEEEremainder(j, 10);

            if (id < 0) {
                id += 10;
            }

            y0 = pixelsPerVoxel * j + xyPixelOffset() + pixelsPerVoxel * 2 / 3;

            //geometry.drawString(Integer.toString(id), i0, j0);
            g.drawString(dstrb, x0, y0);

        }
        
        //y0 = pixelsPerVoxel * jvox + xyPixelOffset() + pixelsPerVoxel * 2 / 3;

        //geometry.drawString(dstrb, i0, j0);

    }

    private void paintVoxelMarquee(Graphics g) {

        int x0, y0, w, h;

        if (marquee) {

            x0 = Math.min(mx0, mx1);
            y0 = Math.min(my0, my1);
            w = Math.abs(mx1 - mx0);
            h = Math.abs(my1 - my0);

            g.setColor(Color.blue);
            g.drawRect(x0, y0, w, h);

        }

    }

    private void paintVoxelCrossAndFrame(Graphics g) {

        int x0, y0;
        double fi, fj;
        double xc, yc, zc, offset = 0.5;

        if (preview) {
            return;
        }

        int gridWidth = getGridWidth();
        int gridHeight = getGridHeight();

        Geometry geometry = Master.project.geometry;

        if (cross) {

            xc = geometry.xVoxelCenter + offset;
            yc = geometry.yVoxelCenter + offset;
            zc = geometry.zVoxelCenter + offset;

            fi = getI(xc, yc, zc);
            fj = getJ(xc, yc, zc);

            x0 = (int) iVoxel2Pixels(fi, fj);
            y0 = (int) jVoxel2Pixels(fi, fj);

            g.setColor(Color.BLACK);

            if (sym4) {
                g.drawLine(0, gridHeight, gridWidth * 2, gridHeight);
                g.drawLine(gridWidth, 0, gridWidth, gridHeight * 2);
            } else {
                g.drawLine(0, y0, gridWidth + pixelsPerVoxel, y0);
                g.drawLine(x0, 0, x0, gridHeight + pixelsPerVoxel);
            }

        }

        if (frame) {

            g.setColor(Color.BLACK);

            x0 = pixelsPerVoxel;
            y0 = pixelsPerVoxel;

            if (sym4) {
                g.drawRect(x0, y0, 2 * iVoxels() * pixelsPerVoxel, 2 * jVoxels() * pixelsPerVoxel);
            } else {
                g.drawRect(x0, y0, iVoxels() * pixelsPerVoxel, jVoxels() * pixelsPerVoxel);
            }

        }

    }

    private void paintPSFcenter(Graphics g) {

        int x0, y0, spot;
        double fi, fj, fk;

        if (preview) {
            return;
        }

        fk = getPSFK0();

        if (((displayMode == 1) || (displayMode == 3)) && ((kVoxel >= fk - 0.5) && (kVoxel <= fk))) {

            fi = getPSFI0();
            fj = getPSFJ0();

            x0 = (int) iVoxel2Pixels(fi, fj);
            y0 = (int) jVoxel2Pixels(fi, fj);

            spot = pixelsPerVoxel;

            x0 += 1 + (spot / 4);
            y0 += 1 + (spot / 4);

            //geometry.setColor(psfCenterColor);
            //geometry.fillOval(i0, j0, spot / 2, spot / 2);

        }

    }

    private void paintMonteCarlo(Graphics g) {
        
        if (Master.project.monteCarlo == null) {
            return;
        }

        if (Master.project.monteCarlo instanceof RunMonteCarloAZEM) {
            paintMonteCarloTheActiveZoneEM(g);
        } else if (Master.project.monteCarlo instanceof RunMonteCarloAZ) {
            paintMonteCarloTheActiveZone(g);
        }
        
        paintMonteCarloDiffusants(g);
        paintMonteCarloMito(g);

    }

    private void paintMonteCarloDiffusants(Graphics g) {
        
        double radius, area = 0;
        int icount = 0;

        DiffusantParticles d;

        Color color;

        Diffusant[] diffusants = Master.project.diffusants;

        RunMonteCarlo monteCarlo = Master.project.monteCarlo;

        if (diffusants == null) {
            return;
        }

        if ((displayMode != 1) && (displayMode != 4)) {
            return;
        }

        for (int i = 0; i < diffusants.length; i++) {

            if ((displayModeArrayNum >= 0) && (i != displayModeArrayNum)) {
                continue;
            }

            if (diffusants[i] instanceof DiffusantParticles) {

                d = (DiffusantParticles) diffusants[i];

                if (d.particles == null) {
                    break;
                }

                for (DiffusantParticle p : d.particles) {

                    if (p == null) {
                        continue;
                    }

                    if (d instanceof DiffusantVesiclesAZ) {
                        color = getMonteCarloAZColor((DiffusantVesiclesAZ) d, (DiffusantVesicleAZ) p);
                    } else {
                        color = getMonteCarloColor(d, p);
                    }

                    radius = paintMonteCarloDiffusantParticle(g, p, color, monteCarlo.PBC);
                    
                    //if (radius > 0) {
                        //Master.log("" + radius);
                        //area += Math.PI * Math.pow(radius, 2.0);
                        //icount++;
                    //}

                }

            }
        }
        
       //Master.log("particles = " + (ilarge) + " (" + ((ilarge + ismall)) + ")");
       //Master.log("max diameter = " + diameterMax);
       //Master.log("particles area = " + area + " um^2");
       //Master.log("" + icount + "\t" + area);
    }

    private Color getMonteCarloColor(DiffusantParticles d, DiffusantParticle p) {

        Color color = Color.YELLOW;

        if (p.insideGeometry) {
            
            if (Master.project.monteCarlo instanceof RunMonteCarloPhoto) {

                RunMonteCarloPhoto mc = (RunMonteCarloPhoto) Master.project.monteCarlo;

                if (mc.frapOn) {
                    color = d.colorReady.getColor(p.fluorescence, 0, 1);
                } else {
                    color = d.colorReady.color;
                }

            } else {
                color = d.colorReady.color;
            }



            //if (!preview && !p.mobile) {
            //    color = dvs.colorImmobile.color;
            //}

            if (p.isConnected()) {
                color = d.colorConnected.color;
            }
            
            if (!p.mobile) {
                color = d.colorImmobile.color;
            }

        }

        return color;

    }

    private Color getMonteCarloAZColor(DiffusantVesiclesAZ d, DiffusantVesicleAZ v) {

        Color color = Color.YELLOW;

        RunMonteCarloAZ montecarlo = (RunMonteCarloAZ) Master.project.monteCarlo;

        if (v.isReady) {

            color = d.colorReady.color;

            //if (!preview && !p.mobile) {
            //    color = d.colorImmobile.color;
            //}

            if (!v.mobile) {
                color = d.colorImmobile.color;
            }

            if (v.isConnected()) {
                color = d.colorConnected.color;
            }

            if ((montecarlo.azDV != null) && (montecarlo.azDV.connectTo != null)) {
                for (DiffusantParticle v2 : montecarlo.azDV.connectTo) {
                    if (v == v2) {
                        color = Color.PINK;
                    }
                }
            }

        }

        if (v.isDocked) {
            color = d.colorDocked.color;
        }

        if (v.isReserve) {
            color = d.colorReserve.color;
        }

        return color;

    }

    private double paintMonteCarloDiffusantParticle(Graphics g, DiffusantParticle p, Color color, boolean PBC) {

        double xvoxel, yvoxel, zvoxel;
        double ivoxel, jvoxel, kvoxel, halfvoxelwidth;
        double addx, addy, addz, addi, addj;
        double radius;

        int i0, j0, particlePixels;

        boolean beyondLimits;

        Geometry geometry = Master.project.geometry;

        if (p == null) {
            return -1;
        }

        xvoxel = geometry.computeVoxelX(p.x);
        yvoxel = geometry.computeVoxelY(p.y);
        zvoxel = geometry.computeVoxelZ(p.z);

        kvoxel = getK(xvoxel, yvoxel, zvoxel);

        if ((kvoxel < kVoxel) || (kvoxel >= kVoxel + 1)) {
            return -1;
        }

        ivoxel = getI(xvoxel, yvoxel, zvoxel);
        jvoxel = getJ(xvoxel, yvoxel, zvoxel);

        i0 = (int) iVoxel2Pixels(ivoxel, jvoxel);
        j0 = (int) jVoxel2Pixels(ivoxel, jvoxel);

        halfvoxelwidth = getHalfVoxelWidthK(kvoxel, p.radius / Master.project.dx);
        //halfvoxelwidth = p.radius / Master.project.dx;
        
        if (phi_Keiding > 0) {
            double radius_phi = p.radius * Math.sin(phi_Keiding * Math.PI / 180);
            if (halfvoxelwidth < radius_phi) {
                return -1;
            }
        }

        particlePixels = (int) (halfvoxelwidth * pixelsPerVoxel);

        particlePixels = Math.max(particlePixels, 1);

        i0 -= particlePixels;
        j0 -= particlePixels;

        g.setColor(color);

        g.fillOval(i0, j0, particlePixels * 2, particlePixels * 2);
        
        //radius = Master.project.dx * ( particlePixels * 1.0) / ( pixelsPerVoxel * 1.0);

        if (!PBC || !redParticlesPBC) {
            //return radius;
            return (particlePixels * 2.0);
        }

        beyondLimits = false;

        addx = 0;
        addy = 0;
        addz = 0;

        if (p.x < geometry.x1 + p.radius) {
            beyondLimits = true;
            addx = 1;
        } else if (p.x > geometry.x2 - p.radius) {
            beyondLimits = true;
            addx = -1;
        }

        if (p.y < geometry.y1 + p.radius) {
            beyondLimits = true;
            addy = 1;
        } else if (p.y > geometry.y2 - p.radius) {
            beyondLimits = true;
            addy = -1;
        }

        if (p.z < geometry.z1 + p.radius) {
            beyondLimits = true;
            addz = 1;
        } else if (p.z > geometry.z2 - p.radius) {
            beyondLimits = true;
            addz = -1;
        }

        if (!beyondLimits) {
            return -1;
        }

        addi = getI(addx, addy, addz);
        addj = getJ(addx, addy, addz);

        if (addi + addj == 0) {
            return -1;
        }

        if (ijSwitch) {
            i0 += addj * jVoxels() * pixelsPerVoxel;
            j0 += addi * iVoxels() * pixelsPerVoxel;
        } else {
            i0 += addi * iVoxels() * pixelsPerVoxel;
            j0 += addj * jVoxels() * pixelsPerVoxel;
        }

        g.setColor(Color.RED);
        g.fillOval(i0, j0, particlePixels * 2, particlePixels * 2);

        return (particlePixels * 2.0);

    }

    private double getHalfVoxelWidthK(double kvoxelCenter, double ijkVoxelRadius) {

        double kdelta = Math.abs(kvoxelCenter - kVoxel - 0.5);
        //return ijkVoxelRadius;
        return Math.sqrt(ijkVoxelRadius * ijkVoxelRadius - kdelta * kdelta); // sphere

    }

    private double getHalfWidthKplane(double kvoxelCenter, double kVoxelRadius, double ijVoxelRadius) {

        double kdelta = Math.abs(kvoxelCenter - kVoxel - 0.5);
        double k2 = kdelta * kdelta / (kVoxelRadius * kVoxelRadius);

        return ijVoxelRadius * Math.sqrt(1 - k2); // ellipsoid

    }

    private void paintMonteCarloMito(Graphics g) {
        
        if (!(Master.project.monteCarlo instanceof RunMonteCarloPhotoMito)) {
            return;
        }
        
        RunMonteCarloPhotoMito mc = (RunMonteCarloPhotoMito) Master.project.monteCarlo;

        Coordinates[] mito = mc.mito;

        //boolean fill = false;
        boolean fill = true;

        if (mito == null) {
            return;
        }

        if ((displayMode != 1) && (displayMode != 4)) {
            return;
        }

        for (Coordinates c : mito) {

            if (c == null) {
                continue;
            }

            paintCoordinates(g, c, c.color.color, fill);

        }

    }

    private void paintCoordinates(Graphics g, Coordinates c, Color color, boolean fill) {

        if (c == null) {
            return;
        }

        switch (c.shapeNum) {
            case 0:
                paintCoordinatesCuboid(g, c, color, fill);
                break;
            case 1:
                paintCoordinatesEllipsoid(g, c, color, fill);
                break;
            case 2:
                paintCoordinatesCylinderX(g, c, color, fill);
                break;
            case 3:
                paintCoordinatesCylinderY(g, c, color, fill);
                break;
            case 4:
                paintCoordinatesCylinderZ(g, c, color, fill);
                break;
        }
        
    }

    private void paintCoordinatesCuboid(Graphics g, Coordinates c, Color color, boolean fill) {

        int ipixelwidth, jpixelwidth;

        Geometry geometry = Master.project.geometry;

        if (c == null) {
            return;
        }

        double xvoxel = geometry.computeVoxelX(c.x1);
        double yvoxel = geometry.computeVoxelY(c.y1);
        double zvoxel = geometry.computeVoxelZ(c.z1);

        double ivoxel1 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel1 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel1 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel1 >= kVoxel + 1) {
            return;
        }

        int ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
        int jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);

        xvoxel = geometry.computeVoxelX(c.x2);
        yvoxel = geometry.computeVoxelY(c.y2);
        zvoxel = geometry.computeVoxelZ(c.z2);

        double ivoxel2 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel2 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel2 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel2 < kVoxel) {
            return;
        }

        int ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
        int jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

        ipixelwidth = ipixel2 - ipixel1;
        jpixelwidth = jpixel2 - jpixel1;

        g.setColor(color);

        if (fill) {
            g.fillRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        } else {
            g.drawRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        }

    }

    private void paintCoordinatesEllipsoid(Graphics g, Coordinates c, Color color, boolean fill) {

        int ipixelwidth, jpixelwidth;

        Geometry geometry = Master.project.geometry;

        if (c == null) {
            return;
        }

        double xvoxel = geometry.computeVoxelX(c.x1);
        double yvoxel = geometry.computeVoxelY(c.y1);
        double zvoxel = geometry.computeVoxelZ(c.z1);

        double ivoxel1 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel1 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel1 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel1 >= kVoxel + 1) {
            return;
        }

        int ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
        int jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);

        xvoxel = geometry.computeVoxelX(c.x2);
        yvoxel = geometry.computeVoxelY(c.y2);
        zvoxel = geometry.computeVoxelZ(c.z2);

        double ivoxel2 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel2 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel2 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel2 < kVoxel) {
            return;
        }

        int ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
        int jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

        ipixelwidth = ipixel2 - ipixel1;
        jpixelwidth = jpixel2 - jpixel1;

        g.setColor(color);

        // UNFINISHED CODE HERE >>

        if (fill) {
            g.fillOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        } else {
            g.drawOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        }

    }

    private void paintCoordinatesCylinderX(Graphics g, Coordinates c, Color color, boolean fill) {

        int ipixel1, jpixel1, ipixel2, jpixel2;
        int ipixelwidth, jpixelwidth;

        Geometry geometry = Master.project.geometry;

        if (c == null) {
            return;
        }

        double xvoxel = geometry.computeVoxelX(c.x1);
        double yvoxel = geometry.computeVoxelY(c.y1);
        double zvoxel = geometry.computeVoxelZ(c.z1);

        double ivoxel1 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel1 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel1 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel1 >= kVoxel + 1) {
            return;
        }

        xvoxel = geometry.computeVoxelX(c.x2);
        yvoxel = geometry.computeVoxelY(c.y2);
        zvoxel = geometry.computeVoxelZ(c.z2);

        double ivoxel2 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel2 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel2 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel2 < kVoxel) {
            return;
        }

        g.setColor(color);

        if (axesMode == 1) { // yz

            ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
            jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);
            ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
            jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

            ipixelwidth = ipixel2 - ipixel1;
            jpixelwidth = jpixel2 - jpixel1;

            if (fill) {
                g.fillOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
            } else {
                g.drawOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
            }
            
            return;

        }

        double ivoxelwidth = ivoxel2 - ivoxel1;
        double jvoxelwidth = jvoxel2 - jvoxel1;
        double kvoxelwidth = kvoxel2 - kvoxel1;

        xvoxel = geometry.computeVoxelX(c.xCenter);
        yvoxel = geometry.computeVoxelY(c.yCenter);
        zvoxel = geometry.computeVoxelZ(c.zCenter);

        double ivoxel0 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel0 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel0 = getK(xvoxel, yvoxel, zvoxel);

        if (axesMode == 0) { // xy

            //ivoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, ivoxelwidth / 2.0);
            jvoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, jvoxelwidth / 2.0);

            //ivoxel1 = ivoxel0 - ivoxelwidth;
            //ivoxel2 = ivoxel0 + ivoxelwidth;
            jvoxel1 = jvoxel0 - jvoxelwidth;
            jvoxel2 = jvoxel0 + jvoxelwidth;

        } else if (axesMode == 2) { // zx

            ivoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, ivoxelwidth / 2.0);
            //jvoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, jvoxelwidth / 2.0);

            ivoxel1 = ivoxel0 - ivoxelwidth;
            ivoxel2 = ivoxel0 + ivoxelwidth;
            //jvoxel1 = jvoxel0 - jvoxelwidth;
            //jvoxel2 = jvoxel0 + jvoxelwidth;

        }

        ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
        jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);
        ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
        jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

        ipixelwidth = ipixel2 - ipixel1;
        jpixelwidth = jpixel2 - jpixel1;

        if ((ipixelwidth <= 0) || (jpixelwidth <= 0)) {
            return;
        }

        if (fill) {
            g.fillRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        } else {
            g.drawRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        }

    }

    private void paintCoordinatesCylinderY(Graphics g, Coordinates c, Color color, boolean fill) {

        int ipixel1, jpixel1, ipixel2, jpixel2;
        int ipixelwidth, jpixelwidth;

        Geometry geometry = Master.project.geometry;

        if (c == null) {
            return;
        }

        double xvoxel = geometry.computeVoxelX(c.x1);
        double yvoxel = geometry.computeVoxelY(c.y1);
        double zvoxel = geometry.computeVoxelZ(c.z1);

        double ivoxel1 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel1 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel1 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel1 >= kVoxel + 1) {
            return;
        }

        xvoxel = geometry.computeVoxelX(c.x2);
        yvoxel = geometry.computeVoxelY(c.y2);
        zvoxel = geometry.computeVoxelZ(c.z2);

        double ivoxel2 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel2 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel2 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel2 < kVoxel) {
            return;
        }

        g.setColor(color);

        if (axesMode == 2) { // zx

            ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
            jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);
            ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
            jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

            ipixelwidth = ipixel2 - ipixel1;
            jpixelwidth = jpixel2 - jpixel1;

            if (fill) {
                g.fillOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
            } else {
                g.drawOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
            }
            
            return;

        }

        double ivoxelwidth = ivoxel2 - ivoxel1;
        double jvoxelwidth = jvoxel2 - jvoxel1;
        double kvoxelwidth = kvoxel2 - kvoxel1;

        xvoxel = geometry.computeVoxelX(c.xCenter);
        yvoxel = geometry.computeVoxelY(c.yCenter);
        zvoxel = geometry.computeVoxelZ(c.zCenter);

        double ivoxel0 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel0 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel0 = getK(xvoxel, yvoxel, zvoxel);
        
        if (axesMode == 1) { // yz

            //ivoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, ivoxelwidth / 2.0);
            jvoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, jvoxelwidth / 2.0);

            //ivoxel1 = ivoxel0 - ivoxelwidth;
            //ivoxel2 = ivoxel0 + ivoxelwidth;
            jvoxel1 = jvoxel0 - jvoxelwidth;
            jvoxel2 = jvoxel0 + jvoxelwidth;

        } else if (axesMode == 0) { // xy

            ivoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, ivoxelwidth / 2.0);
            //jvoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, jvoxelwidth / 2.0);

            ivoxel1 = ivoxel0 - ivoxelwidth;
            ivoxel2 = ivoxel0 + ivoxelwidth;
            //jvoxel1 = jvoxel0 - jvoxelwidth;
            //jvoxel2 = jvoxel0 + jvoxelwidth;

        }

        ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
        jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);
        ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
        jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

        ipixelwidth = ipixel2 - ipixel1;
        jpixelwidth = jpixel2 - jpixel1;

        if ((ipixelwidth <= 0) || (jpixelwidth <= 0)) {
            return;
        }

        if (fill) {
            g.fillRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        } else {
            g.drawRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        }

    }

    private void paintCoordinatesCylinderZ(Graphics g, Coordinates c, Color color, boolean fill) {

        int ipixel1, jpixel1, ipixel2, jpixel2;
        int ipixelwidth, jpixelwidth;

        Geometry geometry = Master.project.geometry;

        if (c == null) {
            return;
        }

        double xvoxel = geometry.computeVoxelX(c.x1);
        double yvoxel = geometry.computeVoxelY(c.y1);
        double zvoxel = geometry.computeVoxelZ(c.z1);

        double ivoxel1 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel1 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel1 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel1 >= kVoxel + 1) {
            return;
        }

        xvoxel = geometry.computeVoxelX(c.x2);
        yvoxel = geometry.computeVoxelY(c.y2);
        zvoxel = geometry.computeVoxelZ(c.z2);

        double ivoxel2 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel2 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel2 = getK(xvoxel, yvoxel, zvoxel);

        if (kvoxel2 < kVoxel) {
            return;
        }

        g.setColor(color);

        if (axesMode == 0) { // xy

            ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
            jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);
            ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
            jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

            ipixelwidth = ipixel2 - ipixel1;
            jpixelwidth = jpixel2 - jpixel1;

            if (fill) {
                g.fillOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
            } else {
                g.drawOval(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
            }

            return;

        }

        double ivoxelwidth = ivoxel2 - ivoxel1;
        double jvoxelwidth = jvoxel2 - jvoxel1;
        double kvoxelwidth = kvoxel2 - kvoxel1;

        xvoxel = geometry.computeVoxelX(c.xCenter);
        yvoxel = geometry.computeVoxelY(c.yCenter);
        zvoxel = geometry.computeVoxelZ(c.zCenter);

        double ivoxel0 = getI(xvoxel, yvoxel, zvoxel);
        double jvoxel0 = getJ(xvoxel, yvoxel, zvoxel);
        double kvoxel0 = getK(xvoxel, yvoxel, zvoxel);

        if (axesMode == 2) { // zx

            //ivoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, ivoxelwidth / 2.0);
            jvoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, jvoxelwidth / 2.0);

            //ivoxel1 = ivoxel0 - ivoxelwidth;
            //ivoxel2 = ivoxel0 + ivoxelwidth;
            jvoxel1 = jvoxel0 - jvoxelwidth;
            jvoxel2 = jvoxel0 + jvoxelwidth;

        } else if (axesMode == 1) { // yz

            ivoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, ivoxelwidth / 2.0);
            //jvoxelwidth = getHalfWidthKplane(kvoxel0, kvoxelwidth / 2.0, jvoxelwidth / 2.0);

            ivoxel1 = ivoxel0 - ivoxelwidth;
            ivoxel2 = ivoxel0 + ivoxelwidth;
            //jvoxel1 = jvoxel0 - jvoxelwidth;
            //jvoxel2 = jvoxel0 + jvoxelwidth;

        }

        ipixel1 = (int) iVoxel2Pixels(ivoxel1, jvoxel1);
        jpixel1 = (int) jVoxel2Pixels(ivoxel1, jvoxel1);
        ipixel2 = (int) iVoxel2Pixels(ivoxel2, jvoxel2);
        jpixel2 = (int) jVoxel2Pixels(ivoxel2, jvoxel2);

        ipixelwidth = ipixel2 - ipixel1;
        jpixelwidth = jpixel2 - jpixel1;

        if ((ipixelwidth <= 0) || (jpixelwidth <= 0)) {
            return;
        }

        if (fill) {
            g.fillRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        } else {
            g.drawRect(ipixel1, jpixel1, ipixelwidth, jpixelwidth);
        }

    }

    private void paintMonteCarloTheActiveZone(Graphics g) {

        Graphics2D g2 = (Graphics2D) g;
        Coordinates center;
        double ax1, ay1, az1;
        double dxHalf = Master.project.dx / 2.0;
        double d = 0.005;

        //boolean fill = false;
        boolean fill = true;

        //ColorD3D azColor = new ColorD3D( "azColor", new Color(150,0,150));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(0,100,255));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(64,110,255));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(190,190,190));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(102,204,0));

        RunMonteCarloAZ mc = (RunMonteCarloAZ) Master.project.monteCarlo;

        if (mc.activeZoneTethering != null) {

            g2.setStroke(new BasicStroke(3));
            paintCoordinates(g2, mc.activeZoneTethering, Color.darkGray, fill);


        } else if (mc.activeZone != null) {

            g2.setStroke(new BasicStroke(3));

            paintCoordinates(g2, mc.activeZone, mc.activeZone.color.color, fill);

            //Coordinates activeZone2 = new Coordinates(Master.project, montecarlo.activeZone);
            //activeZone2.x1 -= 0.025;
            //activeZone2.x2 += 0.025;
            //activeZone2.y1 -= 0.025;
            //activeZone2.y2 += 0.025;
            //activeZone2.z1 -= 0.025;
            //activeZone2.z2 += 0.025;
            //paintCoordinates(g2, activeZone2, azColor.color, fill);

        }

    }

    private void paintMonteCarloTheActiveZoneEM(Graphics g) {

        Graphics2D g2 = (Graphics2D) g;
        Coordinates center;
        double ax1, ay1, az1;
        double dxHalf = Master.project.dx / 2.0;
        double d = 0.005;
        
        //boolean allPoints = false;
        boolean allPoints = true;

        boolean fill = false;
        //boolean fill = true;

        //ColorD3D azColor = new ColorD3D( "azColor", new Color(150,0,150));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(0,100,255));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(64,110,255));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(190,190,190));
        //ColorD3D azColor = new ColorD3D( "azColor", new Color(102,204,0));

        RunMonteCarloAZEM mc = (RunMonteCarloAZEM) Master.project.monteCarlo;

        if (false && (mc.activeZoneTethering != null)) {

            g2.setStroke(new BasicStroke(3));
            paintCoordinates(g2, mc.activeZoneTethering, Color.darkGray, fill);

        } else if (mc.activeZone != null) {

            g2.setStroke(new BasicStroke(3));

            if (allPoints) {

                g2.setStroke(new BasicStroke(6));

                fill = false;

                if (mc.azXYZ != null) {
                    for (int ipnt = 0; ipnt < mc.azXYZ[mc.EMseriesAZ].length; ipnt++) {

                        ax1 = mc.azXYZ[mc.EMseriesAZ][ipnt][0];
                        ay1 = mc.azXYZ[mc.EMseriesAZ][ipnt][1];
                        az1 = mc.azXYZ[mc.EMseriesAZ][ipnt][2];

                        center = new Coordinates(Master.project, ax1-d, ay1-d, az1-dxHalf, ax1+d, ay1+d, az1+dxHalf);

                        paintCoordinates(g2, center, mc.activeZone.color.color, true);

                    }
                }

            } else {

                paintCoordinates(g2, mc.activeZone, mc.activeZone.color.color, fill);

            }

            //Coordinates activeZone2 = new Coordinates(Master.project, montecarlo.activeZone);
            //activeZone2.x1 -= 0.025;
            //activeZone2.x2 += 0.025;
            //activeZone2.y1 -= 0.025;
            //activeZone2.y2 += 0.025;
            //activeZone2.z1 -= 0.025;
            //activeZone2.z2 += 0.025;
            //paintCoordinates(g2, activeZone2, azColor.color, fill);

        }

    }

    // MouseListener and MouseMotionListener functions
    @Override
    public void mouseClicked(MouseEvent e) {

        double value, ivox, jvox;

        mx = e.getX();
        my = e.getY();

        if (displayMode == 0) {

            ivox = iPixels2Voxel(mx, my);
            jvox = jPixels2Voxel(mx, my);

            if (marquee) {
                value = getValue(mx, my);
                setMarquee(value);
            } else {
                setValue(ivox, jvox, getEditValue());
                repaint();
            }

            marquee = false;
            marqueeEdit = false;

        }

    }

    @Override
    public void mousePressed(MouseEvent e) {
        if ((displayMode == 0) && (editMode == 3) && (!marquee)) {
            // finish marquee
            mx0 = e.getX();
            my0 = e.getY();
            mx1 = mx0;
            my1 = my0;
            marquee = true;
            marqueeEdit = true;
        }
    }

    @Override
    public void mouseReleased(MouseEvent e) {

        if ((displayMode == 0) && (editMode == 3)) {
            marqueeEdit = false; // reset marquee
        }

        if (coordClick) {
            updateCursorStats();
        }

    }

    @Override
    public void mouseEntered(MouseEvent e) {
        if ((displayMode == 0) && (editMode == 3)) { // marquee
            setCursor(crossCursor);
        } else {
            setCursor(defaultCursor);
        }
    }

    @Override
    public void mouseExited(MouseEvent e) {
        iVoxel = -1;
        jVoxel = -1;
    }

    @Override
    public void mouseDragged(MouseEvent e) {

        mx = e.getX();
        my = e.getY();

        setVoxelSelect(mx, my);

        if (displayMode == 0) {

            switch (editMode) {

                case 1: // non-space
                case 2: // space
                    setValueRadius(mx, my, getEditValue());
                    break;

                case 3: // marquee
                    if (marqueeEdit) {
                        mx1 = mx;
                        my1 = my;
                    }
                    break;
            }

            repaint();

        }

    }

    @Override
    public void mouseMoved(MouseEvent e) {

        mx = e.getX();
        my = e.getY();

        setVoxelSelect(mx, my);

        if (!coordClick) {
            updateCursorStats();
        }

    }

}
