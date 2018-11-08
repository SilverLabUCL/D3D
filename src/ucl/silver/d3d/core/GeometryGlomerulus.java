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
public final class GeometryGlomerulus {

  public static int kcleft;
  public static int denX, denY;

  public static CoordinatesVoxels coordinates = null;
  public static CoordinatesVoxels dendrites = null; // coordinates of all dendrites

  private GeometryGlomerulus() {
  }

  public static int createTomGeom(Geometry geometry) {

    double denHeight = 2; // dendritic height um
    //double denWidth = 0.42; // dendritic width um
    double denWidth = 0.62; // dendritic width um
    double cleft = 0.02; // cleft width um

    int den = 12;

    Project project = geometry.project;
    double dx = project.dx;

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);

    create(geometry, den, den, dWidth, dHeight, dSpace, dSpace);
    addDetector(geometry, false, den / 2);

    return 1;

  }

  // create GeometryGlomerulus synaptic structure
  public static void create(Geometry geometry, int dNumX, int dNumY, int dWidth,
                            int dHeight, int dSpace, int cSpace) {
      
    int extra = 0;

    int iMax = dWidth * dNumX + (dNumX + 1) * dSpace + extra;
    int jMax = dWidth * dNumY + (dNumY + 1) * dSpace + extra;
    int kMax = dHeight + cSpace;

    kcleft = 0;
    denX = dNumX;
    denY = dNumY;

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    dendrites = new CoordinatesVoxels(geometry.project, 0, 0, kcleft, iMax - 1, iMax - 1, kMax - 1);
    addDendrites(geometry, dendrites, dWidth, dSpace, cSpace);

    geometry.checkSpace();
    geometry.set("name", "Glomerulus " + dNumX + "x" + dNumY);

  }

  public static void addDendrites(Geometry geometry, CoordinatesVoxels c, int dWidth,
                                  int dSpace, int cSpace) {
    int ii, jj, ui, uj;
    int unit = dSpace + dWidth;

    for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
      for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
        for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

          if ( (k >= c.zVoxel1) && (k <= c.zVoxel1 + cSpace - 1)) {
            geometry.setSpace(i, j, k, 1); // cleft
            continue;
          }

          if ( (k >= c.zVoxel2 - dSpace + 1) && (k <= c.zVoxel2)) {
            geometry.setSpace(i, j, k, 1); // top ends are in space
            continue;
          }

          ii = i - (c.xVoxel1 - 1);
          jj = j - (c.yVoxel1 - 1);

          ui = (ii - 1) / unit;
          uj = (jj - 1) / unit;

          if ( ( (ii - ui * unit >= 1) && (ii - ui * unit <= dSpace)) ||
              ( (jj - uj * unit >= 1) && (jj - uj * unit <= dSpace))) {
            geometry.setSpace(i, j, k, 1);
          }
          else {
            geometry.setSpace(i, j, k, -1);
          }

        }
      }

    }

  }

  public static void addCellBoarder(Geometry geometry, CoordinatesVoxels c, int imid,
                                    int jmid, int cWidth, boolean quarter) {

    for (int k = c.zVoxel1; k <= c.zVoxel2; k++) {
      for (int j = c.yVoxel1; j <= c.yVoxel2; j++) {
        for (int i = c.xVoxel1; i <= c.xVoxel2; i++) {

          if (k == c.zVoxel2) {
            //geom.setSpace(i, j, k, 1);
          }

          if ( (i == 1) || (j == 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if (!quarter && ( (i == c.xVoxel2) || (j == c.yVoxel2))) {
            geometry.setSpace(i, j, k, 1);
          }

          // lines radiating from middle

          if ( (i <= imid) && (j == jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          // cells border lines on outside

          if ( (i <= imid) && (j == jmid - cWidth - 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1 - cWidth - 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i <= imid) && (j == jmid + cWidth + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1 + cWidth + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1 - cWidth - 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid - cWidth - 1) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1 + cWidth + 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + cWidth + 1) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

  }

  public static void addDetector(Geometry geometry, boolean quarter, int claws) {
    int x0, y0, z0;
    int x1, y1, x2, y2, i, j, iend, jend;

    Project project = geometry.project;
    double dx = project.dx;

    CoordinatesVoxels c = new CoordinatesVoxels(project);

    double psd = 0.2; // detector width um
    int dPSD = (int) (psd / dx);

    double denWidth = 0.6;
    int dw = 1 + (int) (denWidth / dx);

    if (quarter) {
      x0 = (geometry.xVoxels / 1) - (dPSD / 2) - 1;
      y0 = (geometry.yVoxels / 1) - (dPSD / 2) - 1;
    }
    else {
      x0 = (geometry.xVoxels / 2) - (dPSD / 2);
      y0 = (geometry.yVoxels / 2) - (dPSD / 2);
    }

    z0 = kcleft;

    iend = 1 + (claws / 2);
    jend = 1; // first row only
    //jend = iend;

    for (j = 0; j < jend; j++) {

      y1 = y0 - j * dw;
      y2 = y1 + dPSD - 1;

      for (i = 0; i < iend; i++) {

        x1 = x0 - i * dw;
        x2 = x1 + dPSD - 1;

        if (i >= j) {
          Master.addDetectorAvg(0, c.setVoxels(x1, y1, z0, x2, y2, z0));
          if ( (i == 0) && (j == 0)) {
            //mast.addDetectorAvg(1, c.set(x1, y1, z0, x2, y2, z0)); // record MNI
          }
        }

      }

    }
  }

  public static void createCrissCross(Geometry geometry) {
    int ui, uj;
    int ibgn, iend, jbgn, jend, kbgn, kend;
    int dNumX = 9;
    int dNumY = dNumX;
    double denHeight = 0.62; // dendritic height um
    double denWidth = 0.62; // dendritic width um
    double cleft = 0.02; // cleft width um

    Project project = geometry.project;
    double dx = project.dx;

    kcleft = 1;

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);
    int cSpace = dSpace;

    int iMax = dWidth * dNumX + (dNumX + 1) * dSpace + 2;
    int jMax = dWidth * dNumY + (dNumY + 1) * dSpace + 2;
    int kMax = cSpace + dHeight + dWidth * 2 + dSpace * 3 + 2;
    int unit = dSpace + dWidth;

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    kMax = cSpace + dHeight + 3;

    ibgn = 1;
    jbgn = 1;
    kbgn = 1;

    iend = iMax - 2;
    jend = jMax - 2;
    kend = kMax - 2;

    dendrites = new CoordinatesVoxels(geometry.project, ibgn, jbgn, kbgn, iend, jend, kend);
    addDendrites(geometry, dendrites, dWidth, dSpace, cSpace);

    kbgn = kend + 1;
    kend = kbgn + dWidth;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if ( (i == ibgn) || (i == iend) || (k == kend)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          uj = j / unit;

          if ( (j - uj * unit >= 1) && (j - uj * unit <= dSpace)) {
            geometry.setSpace(i, j, k, 1);
          }
          else {
            geometry.setSpace(i, j, k, -1);
          }

        }
      }

    }

    kbgn = kend + 1;
    kend = kbgn + dWidth;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if ( (j == ibgn) || (j == iend) || (k == kend)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          ui = i / unit;

          if ( (i - ui * unit >= 1) && (i - ui * unit <= dSpace)) {
            geometry.setSpace(i, j, k, 1);
          }
          else {
            geometry.setSpace(i, j, k, -1);
          }

        }
      }

    }

    addDetector(geometry, false, dNumX);

    geometry.checkSpace();
    geometry.set("name", "Glomerulus CrissCross " + dNumX + "x" + dNumY);

  }

  public static void createClaw(Geometry geometry) {
    int ui, uj;
    int ibgn, iend, jbgn, jend, kbgn, kend;
    int dNumX = 9;
    int dNumY = dNumX;
    double denHeight = 0.62; // dendritic height um
    double denWidth = 0.62; // dendritic width um
    double cleft = 0.02; // cleft width um
    double remain;

    Project project = geometry.project;
    double dx = project.dx;

    kcleft = 1;

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);
    int cSpace = dSpace;

    int iMax = dWidth * dNumX + (dNumX + 1) * dSpace + 2;
    int jMax = dWidth * dNumY + (dNumY + 1) * dSpace + 2;
    int kMax = cSpace + dHeight + dWidth * 2 + dSpace * 3 + 2;
    int unit = dSpace + dWidth;

    boolean space;

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    kMax = cSpace + dHeight + 3;

    ibgn = 1;
    iend = iMax - 2;
    jbgn = 1;
    jend = jMax - 2;
    kbgn = 1;
    kend = kMax - 2;

    //this.addDendrites(ibgn, jbgn, kbgn, iend, jend, kend, dWidth, dSpace, cSpace);

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if ( (k >= 1) && (k <= cSpace)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          ui = (i - 1) / unit;
          uj = (j - 1) / unit;

          space = ( ( (i - ui * unit >= 1) && (i - ui * unit <= dSpace)) ||
                   ( (j - uj * unit >= 1) && (j - uj * unit <= dSpace)));

          remain = Math.abs(Math.IEEEremainder(ui + uj, 2));

          if (space) {
            geometry.setSpace(i, j, k, 1);
          }
          else {

            if (k == kend) {

              if (remain > 0) {
                geometry.setSpace(i, j, k, -1);
              }
              else {
                geometry.setSpace(i, j, k, 1);
              }

            }
            else {
              geometry.setSpace(i, j, k, -1);
            }

          }

        }
      }

    }

    kbgn = kend + 1;
    kend = kbgn + dWidth;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if ( (i == ibgn) || (i == iend) || (k == kend)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          uj = j / unit;

          if ( (j - uj * unit >= 1) && (j - uj * unit <= dSpace)) {
            geometry.setSpace(i, j, k, 1);
          }
          else {
            geometry.setSpace(i, j, k, -1);
          }

        }
      }

    }

    kbgn = kend + 1;
    kend = kbgn + dWidth;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if ( (j == ibgn) || (j == iend) || (k == kend)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          ui = i / unit;

          if ( (i - ui * unit >= 1) && (i - ui * unit <= dSpace)) {
            geometry.setSpace(i, j, k, 1);
          }
          else {
            geometry.setSpace(i, j, k, -1);
          }

        }
      }

    }

    addDetector(geometry, false, dNumX);

    geometry.checkSpace();

  }

  public static void createSmallSpot(Geometry geometry, boolean quarter, boolean large) {

    int ibgn, iend, jbgn, jend, kbgn, kend;
    int dNum = 11;

    //double denHeight = 2.00; // dendritic height um
    //double denWidth = 0.62; // dendritic width um (small spot sims)

    double denHeight = 1.24; // dendritic height um
    double denWidth = 0.60; // dendritic width um (small spot sims)
    double cleft = 0.02; // cleft width um
    double cellHeight = 3;
    double cellWidth = 7;

    double dx = geometry.project.dx;

    kcleft = 1;

    if (large) {
      dNum = 11;
      cellWidth = 7.0;
    }
    else {
      dNum = 5;
      cellWidth = 3;
    }

    int cHeight = (int) (cellHeight / dx);

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);

    int dSpace = (int) (cleft / dx);
    int cSpace = dSpace;
    int cWidth = (int) (cellWidth / dx);

    int iMax = dWidth * dNum + (dNum + 1) * dSpace + 2;
    int jMax = dWidth * dNum + (dNum + 1) * dSpace + 2;
    int kMax = cSpace + dHeight + cHeight + dSpace + 2;

    int imid = (int) ( (iMax - 1.0) / 2.0);
    int jmid = (int) ( (jMax - 1.0) / 2.0);

    if (quarter) {
      iMax = (iMax / 2) + 1;
      jMax = (jMax / 2) + 1;
    }

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    // create packed dendrites

    ibgn = 1;
    jbgn = ibgn;

    if (quarter) {
      iend = iMax - 2;
      jend = jMax - 2;
    }
    else {
      iend = iMax - 2;
      jend = jMax - 2;
    }

    kbgn = 1;
    kend = cSpace + dHeight + 1;

    dendrites = new CoordinatesVoxels(geometry.project, ibgn, jbgn, kbgn, iend, jend, kend);
    addDendrites(geometry, dendrites, dWidth, dSpace, cSpace);

    // cell border on top

    kbgn = kend + 1;
    kend = kbgn + cHeight - 1;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if ( (i == imid) || (j == jmid)) {
            //geom.setSpace(i, j, k, 1);
          }

          if ( (i <= imid) && (j == jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    addDetector(geometry, quarter, dNum);

    geometry.set("name", "Glomerulus Small Spot");

    geometry.checkSpace();

  }

  public static void createBigSpot(Geometry geometry, boolean quarter) {
    int ibgn, iend, jbgn, jend, kbgn, kend;
    int dNumX = 11;
    int dNumY = dNumX;

    double denHeight = 1.25; // dendritic height um
    double denWidth = 0.60; // dendritic width um
    double cleft = 0.02; // cleft width um
    double cellHeight = 3;
    double cellBoarder = 14.04;

    double dx = geometry.project.dx;

    kcleft = 1;

    int cHeight = (int) (cellHeight / dx);

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);
    int cSpace = dSpace;
    int cWidth = (int) (cellBoarder / dx);

    int iMax = 2 * cWidth + dWidth * dNumX + (dNumX + 1) * dSpace + 2;
    int jMax = 2 * cWidth + dWidth * dNumY + (dNumY + 1) * dSpace + 2;
    int kMax = cSpace + dHeight + cHeight + dSpace + 2;

    int imid = (int) ( (iMax - 1.0) / 2.0);
    int jmid = (int) ( (jMax - 1.0) / 2.0);

    if (quarter) {
      iMax = (iMax / 2) + 1;
      jMax = (jMax / 2) + 1;
    }

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    // create packed dendrites

    kMax = cSpace + dHeight + 3;

    ibgn = cWidth + 1;
    jbgn = cWidth + 1;

    kbgn = 1;
    kend = kMax - 3;

    if (quarter) {
      iend = iMax - 2;
      jend = jMax - 2;
    }
    else {
      iend = iMax - cWidth - 2;
      jend = jMax - cWidth - 2;
    }

    dendrites = new CoordinatesVoxels(geometry.project, ibgn, jbgn, kbgn, iend, jend, kend);
    addDendrites(geometry, dendrites, dWidth, dSpace, cSpace);

    // cell boarder on sides

    ibgn = 1;
    iend = iMax - 2;
    jbgn = 1;
    jend = jMax - 2;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if (k == kend) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (i == 1) || (j == 1)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (i == cWidth + 1) || (j == cWidth + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == cWidth / 2 + 1) || (j == cWidth / 2 + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if (!quarter && ( (i == iMax - cWidth - 2) || (j == jMax - cWidth - 2))) {
            geometry.setSpace(i, j, k, 1);
          }

          if (!quarter &&
              ( (i == iMax - cWidth / 2 - 2) || (j == jMax - cWidth / 2 - 2))) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    // cell border on top

    kbgn = kend + 1;
    kend = kbgn + cHeight - 1;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if ( (i == 1) || (j == 1)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (i == imid / 2) || (j == jmid / 2)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == 3 * imid / 2) || (j == 3 * jmid / 2)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i <= imid) && (j == jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    addDetector(geometry, quarter, dNumX);

    geometry.set("name", "Glomerulus Big Spot");

    geometry.checkSpace();

  }

  public static void createBigSpot2(Geometry geometry, boolean quarter, boolean large) {
    int ui, uj, ii, jj;
    int ibgn, iend, jbgn, jend, kbgn, kend;
    int idbgn, idend, jdbgn, jdend;
    int dNum = 11;

    double denHeight = 1.24; // dendritic height um
    double denWidth = 0.60; // dendritic width um
    double cleft = 0.02; // cleft width um
    double cellWidth = 7;

    double dx = geometry.project.dx;

    if (large) {
      dNum = 11;
      cellWidth = 7;
    }
    else {
      dNum = 5;
      cellWidth = 3;
    }

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);
    int cWidth = (int) (cellWidth / dx);

    int iMax = (2 * cWidth) + (dWidth * dNum) + ( (dNum + 3) * dSpace) + 2;
    int jMax = iMax;
    int kMax = (4 * dSpace) + (2 * cWidth) + dHeight + 2;
    int unit = dSpace + dWidth;

    int imid = (int) ( (iMax - 1.0) / 2.0);
    int jmid = (int) ( (jMax - 1.0) / 2.0);

    if (quarter) {
      iMax = (iMax / 2) + 1;
      jMax = (jMax / 2) + 1;
    }

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    // create packed dendrites

    idbgn = dSpace + cWidth + 1;
    jdbgn = idbgn;

    kbgn = cWidth + dSpace + 1;
    kend = kbgn + dHeight + 1;

    if (quarter) {
      idend = iMax - 2;
    }
    else {
      idend = iMax - cWidth - dSpace - 2;
    }

    jdend = idend;

    kcleft = kbgn;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jdbgn; j <= jdend; j++) {
        for (int i = idbgn; i <= idend; i++) {

          if ( (k >= kbgn) && (k <= kbgn + dSpace - 1)) {
            geometry.setSpace(i, j, k, 1); // cleft
            continue;
          }

          ii = i - (cWidth + dSpace);
          jj = j - (cWidth + dSpace);

          ui = (ii - 1) / unit;
          uj = (jj - 1) / unit;

          if ( ( (ii - ui * unit >= 1) && (ii - ui * unit <= dSpace)) ||
              ( (jj - uj * unit >= 1) && (jj - uj * unit <= dSpace))) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    // cell boarder on sides

    kbgn = 1;

    ibgn = 1;
    iend = iMax - 2;
    jbgn = 1;
    jend = jMax - 2;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if (k == kend) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (k == 1) && (i > idbgn) && (j > jdbgn) && (i <= idend) &&
              (j <= jdend)) {
            //geom.setSpace(i, j, k, 1);
          }

          if (k == kend - cWidth - 1) {
            if ( (i < idbgn) || (j < jdbgn) || (i > idend) || (j > jdend)) {
              geometry.setSpace(i, j, k, 1);
            }
          }

          if ( (i == 1) || (j == 1)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if (!quarter && ( (i == iend) || (j == iend))) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (i == cWidth + dSpace + 1) || (j == cWidth + dSpace + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if (!quarter &&
              ( (i == iMax - cWidth - dSpace - 2) ||
               (j == jMax - cWidth - dSpace - 2))) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    // cell border on top

    kbgn = kend + 1;
    kend = kbgn + cWidth;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if (k == kend) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          // lines radiating from middle

          if ( (i <= imid) && (j == jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          // cells border lines on outside

          if ( (i <= imid) && (j == jmid - cWidth - 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1 - cWidth - 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i <= imid) && (j == jmid + cWidth + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i > imid) && (j == jmid + 1 + cWidth + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1 - cWidth - 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid - cWidth - 1) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + 1 + cWidth + 1) && (j <= jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == imid + cWidth + 1) && (j > jmid)) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    addDetector(geometry, quarter, dNum);

    geometry.set("name", "Glomerulus Big Spot 2");

    geometry.checkSpace();

  }

  public static void createBigSpot3(Geometry geometry, boolean quarter, boolean large) {
    int ui, uj, ii, jj;
    int ibgn, iend, jbgn, jend, kbgn, kend;
    int idbgn, idend, jdbgn, jdend;
    int dNum = 11;

    double denHeight = 1.24; // dendritic height um
    double denWidth = 0.60; // dendritic width um
    double cleft = 0.02; // cleft width um
    double cellWidth = 7;

    Project project = geometry.project;
    double dx = project.dx;

    if (large) {
      dNum = 11;
      cellWidth = 7;
    }
    else {
      dNum = 5;
      cellWidth = 3;
    }

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);
    int cWidth = (int) (cellWidth / dx);

    int iMax = (2 * cWidth) + (dWidth * dNum) + ( (dNum + 3) * dSpace) + 2;
    int jMax = iMax;
    int kMax = (4 * dSpace) + (2 * cWidth) + 1;
    int unit = dSpace + dWidth;

    int imid = (int) ( (iMax - 1.0) / 2.0);
    int jmid = (int) ( (jMax - 1.0) / 2.0);

    if (quarter) {
      iMax = (iMax / 2) + 1;
      jMax = (jMax / 2) + 1;
    }

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    // create packed dendrites

    idbgn = dSpace + cWidth + 1;
    jdbgn = idbgn;

    kbgn = (cWidth - dHeight) + dSpace;
    kend = kbgn + dHeight + 1;

    kcleft = kbgn;

    if (quarter) {
      idend = iMax - 2;
    }
    else {
      idend = iMax - cWidth - dSpace - 2;
    }

    jdend = idend;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jdbgn; j <= jdend; j++) {
        for (int i = idbgn; i <= idend; i++) {

          if ( (k >= kbgn) && (k <= kbgn + dSpace - 1)) {
            geometry.setSpace(i, j, k, 1); // cleft
            continue;
          }

          ii = i - (cWidth + dSpace);
          jj = j - (cWidth + dSpace);

          ui = (ii - 1) / unit;
          uj = (jj - 1) / unit;

          if ( ( (ii - ui * unit >= 1) && (ii - ui * unit <= dSpace)) ||
              ( (jj - uj * unit >= 1) && (jj - uj * unit <= dSpace))) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    // cell boarder on sides

    kbgn = 1;

    ibgn = 1;
    iend = iMax - 2;
    jbgn = 1;
    jend = jMax - 2;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if (k == kend) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (k == 1) && (i > idbgn) && (j > jdbgn) && (i <= idend) &&
              (j <= jdend)) {
            //geom.setSpace(i, j, k, 1);
          }

          if (k == kend - cWidth - 1) {
            if ( (i < idbgn) || (j < jdbgn) || (i > idend) || (j > jdend)) {
              geometry.setSpace(i, j, k, 1);
            }
          }

          if ( (i == 1) || (j == 1)) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if (!quarter && ( (i == iend) || (j == iend))) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (i == cWidth + dSpace + 1) || (j == cWidth + dSpace + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if (!quarter &&
              ( (i == iMax - cWidth - dSpace - 2) ||
               (j == jMax - cWidth - dSpace - 2))) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    // cell border on top

    kbgn = kend + 1;
    kend = kbgn + cWidth;

    coordinates = new CoordinatesVoxels(project);
    coordinates.setVoxels(ibgn, jbgn, kbgn, iend, jend, kend);

    addCellBoarder(geometry, coordinates, imid, jmid, cWidth, quarter);

    addDetector(geometry, quarter, dNum);

    geometry.set("name", "Glomerulus Big Spot 3");

    geometry.checkSpace();

  }

  public static void createBigSpot4(Geometry geometry, boolean quarter, boolean large) {
    int ibgn, iend, jbgn, jend, kbgn, kend;
    int idbgn, idend, jdbgn, jdend;

    double denHeight = 1.24; // value used in 2007 paper
    //double denHeight = 2.00; // dendritic height um
    double denWidth = 0.60; // dendritic width um
    double cleft = 0.02; // cleft width um
    double cellWidth = 7;
    double zBoarder = 0.5;
    //double zBoarder = 1.5;

    Master.log("Creating Glomerulus BigSpot4");

    Project project = geometry.project;
    double dx = project.dx;

    if (large) {
      denX = 11;
      denY = 11;
      cellWidth = 7;
    }
    else {
      denX = 5;
      denY = 5;
      cellWidth = 3;
    }

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);
    int cWidth = (int) (cellWidth / dx);
    int zSpace = (int) (zBoarder / dx);

    int iMax = (4 * cWidth) + (5 * dSpace) + 2 - 1;
    int jMax = iMax;
    int gWidth = (dWidth * denX) + ( (denX + 1) * dSpace);
    int kMax = dHeight + (2 * dSpace) + (2 * zSpace) + 2;

    int imid = (int) ( (iMax - 1.0) / 2.0);
    int jmid = (int) ( (jMax - 1.0) / 2.0);

    kcleft = zSpace + dSpace;

    if (quarter) {
      iMax = (iMax / 2) + 1;
      jMax = (jMax / 2) + 1;
    }

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    // create packed dendrites

    kbgn = kcleft;
    kend = kbgn + dHeight + 1;

    if (quarter) {
      idbgn = (int) (iMax - (gWidth / 2) - 1);
      idend = iMax - 1;
    }
    else {
      idbgn = (int) ( (iMax - gWidth) / 2);
      idend = idbgn + gWidth - 1;
    }

    jdbgn = idbgn;
    jdend = idend;

    dendrites = new CoordinatesVoxels(project, idbgn, jdbgn, kbgn, idend, jdend, kend);
    addDendrites(geometry, dendrites, dWidth, dSpace, dSpace);

    // cell boarder on sides

    ibgn = 1;
    iend = iMax - 2;
    jbgn = 1;
    jend = jMax - 2;

    kbgn = 1;

    for (int k = kbgn; k <= kend; k++) {
      for (int j = jbgn; j <= jend; j++) {
        for (int i = ibgn; i <= iend; i++) {

          if (k == 1) {
            //geom.setSpace(i, j, k, 1);
            //continue;
          }

          if (k == kend) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( (i == idbgn) || (j == jdbgn)) {
            geometry.setSpace(i, j, k, 1);
          }

          if (!quarter && ( (i == idend) || (j == jdend))) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == idbgn - cWidth - 1) || (j == jdbgn - cWidth - 1)) {
            geometry.setSpace(i, j, k, 1);
          }

          if ( (i == idend + cWidth + 1) || (j == jdend + cWidth + 1)) {
            geometry.setSpace(i, j, k, 1);
          }

        }
      }

    }

    // cell border on top

    kbgn = kend + 1;
    kend = kMax - 2;

    coordinates = new CoordinatesVoxels(project);
    coordinates.setVoxels(ibgn, jbgn, kbgn, iend, jend, kend);

    addCellBoarder(geometry, coordinates, imid, jmid, cWidth, quarter);

    geometry.set("name", "Glomerulus Big Spot 4");

    geometry.checkSpace();

  }

  public static void createQuanta(Geometry geometry, boolean quarter, boolean large) {
    int imax, jmax, dNum;

    double denHeight = 2.00; // dendritic height um
    double denWidth = 0.60; // dendritic width um
    double cleft = 0.02; // cleft width um

    Project project = geometry.project;
    double dx = project.dx;

    kcleft = 1;

    if (large) {
      dNum = 11;
    }
    else {
      dNum = 5;
    }

    int dHeight = (int) (denHeight / dx);
    int dWidth = (int) (denWidth / dx);
    int dSpace = (int) (cleft / dx);

    create(geometry, dNum, dNum, dWidth, dHeight, dSpace, dSpace);

    CoordinatesVoxels c = new CoordinatesVoxels(project, 1, 1, geometry.zVoxels - 2, geometry.xVoxels - 2,
                                    geometry.yVoxels - 2, geometry.zVoxels - 2);

    geometry.setSpace(c, 1);

    if (quarter) {
      imax = (geometry.xVoxels / 2) + 1;
      jmax = (geometry.yVoxels / 2) + 1;
      geometry.resizeWithSpace(imax, jmax, geometry.zVoxels);
    }

    geometry.checkSpace();

    //addDetectorAvg(mast, quarter, dNum);

  }

  // create GeometryGlomerulus synaptic structure
  public static void createGlomerulusOffset(Geometry geometry, int dNumX, int dNumY,
                                            int dWidth, int dHeight, int dSpace, int cSpace) {

    int ui, uj, ioffset, koffset, ii, jj;
    int iMax = dWidth * dNumX + (dNumX + 1) * dSpace + 2;
    int jMax = dWidth * dNumY + (dNumY + 1) * dSpace + 2;
    int kMax = dHeight + cSpace + 2;
    int unit = dSpace + dWidth;

    boolean ieven, jeven;

    geometry.resizeAndFill(iMax, jMax, kMax, Geometry.NONSPACEVALUE);

    ioffset = 0;
    koffset = 1;

    for (int k = 1; k < kMax - 1; k++) {
      for (int j = 1; j < jMax - 1; j++) {
        for (int i = 1; i < iMax - 1; i++) {

          ui = (i - 1) / unit;
          uj = (j - 1) / unit;

          ieven = (Math.IEEEremainder(ui, 2) == 0);
          jeven = (Math.IEEEremainder(uj, 2) == 0);

          if (jeven) {
            ioffset = 0;
          }
          else {
            //ioffset = dWidth / 2;
            ioffset = 0;
          }

          if ( (ieven && jeven) || (!ieven && !jeven)) {
            koffset = 1;
          }
          else {
            koffset = 9;
          }

          ii = i - ioffset - ui * unit;
          jj = j - uj * unit;

          if ( ( (ii >= 1) && (ii <= dWidth + 1)) &&
              ( (jj >= 1) && (jj <= dWidth + 1)) &&
              ( (k >= koffset) && (k < koffset + cSpace))) {
            geometry.setSpace(i, j, k, 1);
            continue;
          }

          if ( ( (ii >= 1) && (ii <= dSpace)) || ( (jj >= 1) && (jj <= dSpace))) {
            geometry.setSpace(i, j, k, 1);
          }
          else {
            geometry.setSpace(i, j, k, -1);
          }

        }
      }
    }

    geometry.checkSpace();

  }

}
