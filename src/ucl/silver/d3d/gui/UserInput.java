package ucl.silver.d3d.gui;

import ucl.silver.d3d.core.*;
import java.awt.event.ActionEvent;
import java.beans.PropertyChangeEvent;
import ucl.silver.d3d.utils.Utility;

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
public class UserInput
        implements UserInputAction {

    private static final int fieldWidth = 5; // prompt field display width

    private UserInputPrompt prompt = null;
    private UserInputField fieldArray[] = null;

    public double getFieldValue(int index) {
        Number n = (Number) fieldArray[index].getUserValue();
        return n.doubleValue();
    }

    public String getFieldName(int index) {
        return fieldArray[index].vName;
    }

    @Override
    public void actionField(PropertyChangeEvent e, String name, Number n) {
        //NOTHING TO DO, but must be declared
    }

    @Override
    public void actionButton(ActionEvent e) {
        String s, id;

        if (prompt == null) {
            return;
        }

        s = e.getActionCommand();

        if (s.equalsIgnoreCase("Cancel")) {
            prompt.exit();
            return;
        }

        id = prompt.callerID;

        if (id.equalsIgnoreCase("createGlomerulus")) {
            createGlomerulusPrompt2();
        } else if (id.equalsIgnoreCase("createBouton")) {
            createBoutonPrompt2();
        } else if (id.equalsIgnoreCase("KillBatch")) {
            killBatchPrompt2();
        } else if (id.equalsIgnoreCase("AddDiffusant")) {
            addDiffusantPrompt2();
        } else if (id.equalsIgnoreCase("AddDiffusantBuffer")) {
            addDiffusantBufferPrompt2();
        } else if (id.equalsIgnoreCase("KillDiffusant")) {
            killDiffusantPrompt2();
        } else if (id.equalsIgnoreCase("AddSourceImpulse")) {
            addSourceImpulsePrompt2();
        } else if (id.equalsIgnoreCase("AddSourceGauss")) {
            addSourceGaussPrompt2();
        } else if (id.equalsIgnoreCase("KillSource")) {
            killSourcePrompt2();
        } else if (id.equalsIgnoreCase("AddDetector")) {
            addDetectorPrompt2();
        } else if (id.equalsIgnoreCase("AddDetectorWeighted")) {
            addDetectorWeightedPrompt2();
        } else if (id.equalsIgnoreCase("AddDetectorSnap")) {
            addDetectorSnapPrompt2();
        } else if (id.equalsIgnoreCase("KillDetector")) {
            killDetectorPrompt2();
        } else if (id.equalsIgnoreCase("addDiffusantPhotolysis")) {
            addDiffusantPhotolysisPrompt2();
        } else if (id.equalsIgnoreCase("writePSF")) {
            writePSF2();
        }

        prompt.exit();
        Master.project.init();
        Master.mainframe.initTabs();

    }

    public void createGlomerulusPrompt() {

        int dNumX = 3, dNumY = 3, dWidth = 31, dHeight = 10, dSpace = 1, cSpace = 1;

        fieldArray = new UserInputField[6];

        fieldArray[0] = new UserInputField("dNumX", "dendrites in x-dir",
                new Double(dNumX), fieldWidth, this);
        fieldArray[1] = new UserInputField("dNumY", "dendrites in y-dir",
                new Double(dNumY), fieldWidth, this);
        fieldArray[2] = new UserInputField("dWidth", "dendritic width (voxels)",
                new Double(dWidth), fieldWidth, this);
        fieldArray[3] = new UserInputField("dHeight", "dendritic height (voxels)",
                new Double(dHeight), fieldWidth, this);
        fieldArray[4] = new UserInputField("dSpace",
                "inter-dendrite space (voxels)",
                new Double(dSpace), fieldWidth, this);
        fieldArray[5] = new UserInputField("cSpace", "cleft space (voxels)",
                new Double(cSpace), fieldWidth, this);

        prompt = new UserInputPrompt("CreateGlomerulus", "Create Glomerulus",
                fieldArray, this);

    }

    public void createGlomerulusPrompt2() {

        if (fieldArray.length != 6) {
            return;
        }

        int dNumX = (int) getFieldValue(0);
        int dNumY = (int) getFieldValue(1);
        int dWidth = (int) getFieldValue(2);
        int dHeight = (int) getFieldValue(3);
        int dSpace = (int) getFieldValue(4);
        int cSpace = (int) getFieldValue(5);

        GeometryGlomerulus.create(Master.project.geometry, dNumX, dNumY, dWidth, dHeight, dSpace, cSpace);
        Master.mainframe.panel2D.resetVoxelWidth();
        Master.resetAllParamVectors();

    }

    public void createBoutonPrompt() {

        double axon_width = 0.2;
        double bouton_width = 0.8;
        double bouton_length = 1.2;
        double total_length = 5.0;
        double enpassant = 1;

        fieldArray = new UserInputField[5];

        fieldArray[0] = new UserInputField("axon_width", "axon width (um)",
                axon_width, fieldWidth, this);
        fieldArray[1] = new UserInputField("bouton_width", "bouton width (um)",
                bouton_width, fieldWidth, this);
        fieldArray[2] = new UserInputField("bouton_length", "bouton length (um)",
                bouton_length, fieldWidth, this);
        fieldArray[3] = new UserInputField("total_length", "total length (um)",
                total_length, fieldWidth, this);
        fieldArray[4] = new UserInputField("enpassant", "enpassant (0) no (1) yes",
                enpassant, fieldWidth, this);

        prompt = new UserInputPrompt("CreateBouton", "Create Axonal Bouton",
                fieldArray, this);

    }

    public void createBoutonPrompt2() {

        if (fieldArray.length != 5) {
            return;
        }

        double axon_width = getFieldValue(0);
        double bouton_width = getFieldValue(1);
        double bouton_length = getFieldValue(2);
        double total_length = getFieldValue(3);
        double enpassant = getFieldValue(4);
        boolean enpass = true;

        if (enpassant == 0) {
            enpass = false;
        }

        GeometryBouton.create(Master.project.geometry, axon_width, bouton_width, bouton_length,
                total_length, enpass);

        Master.mainframe.panel2D.resetVoxelWidth();
        Master.resetAllParamVectors();

    }

    public void addDiffusantPrompt() {

        double c0 = 0, dd = 0.5;

        fieldArray = new UserInputField[2];

        fieldArray[0] = new UserInputField("c0", "initial concentration (mM)",
                c0, fieldWidth, this);
        fieldArray[1] = new UserInputField("d", "diffusion coefficient",
                dd, fieldWidth, this);

        prompt = new UserInputPrompt("AddDiffusant", "Add Simple Diffusant",
                fieldArray, this);

    }

    public void addDiffusantPrompt2() {
        if (fieldArray.length != 2) {
            return;
        }
        double c0 = getFieldValue(0);
        double dd = getFieldValue(1);
        Master.addDiffusant(c0, dd);
    }

    public void addDiffusantBufferPrompt() {

        double cTotal = 1, dd = 0.5;
        int dnum = 0;
        double Kon = 100, Koff = 10;

        int plimit = Master.project.diffusants.length;

        fieldArray = new UserInputField[5];

        fieldArray[0] = new UserInputField("cTotal", "total concentration (mM)",
                cTotal, fieldWidth, this);
        fieldArray[1] = new UserInputField("d", "diffusion coefficient",
                dd, fieldWidth, this);
        fieldArray[2] = new UserInputField("dnum",
                "reacts with diffusant number ( < " +
                plimit + " )",
                new Double(dnum), fieldWidth, this);
        fieldArray[3] = new UserInputField("Kon", "on rate",
                Kon, fieldWidth, this);
        fieldArray[4] = new UserInputField("Koff", "off rate",
                Koff, fieldWidth, this);

        prompt = new UserInputPrompt("AddDiffusantBuffer",
                "Add Reaction Diffusant (i.e. Buffer)",
                fieldArray, this);

    }

    public void addDiffusantBufferPrompt2() {

        if (fieldArray.length != 5) {
            return;
        }

        double cTotal = getFieldValue(0);
        double dd = getFieldValue(1);
        int dnum = (int) getFieldValue(2);
        double Kon = getFieldValue(3);
        double Koff = getFieldValue(4);

        Master.addDiffusantReactant("Buffer", cTotal, dd, dnum, Kon, Koff);

    }

    public void killDiffusantPrompt(int arrayNum) {

        Diffusant[] diffusant = Master.project.diffusants;

        if (diffusant == null) {
            return;
        }

        if (arrayNum >= diffusant.length) {
            arrayNum = 0;
        }

        fieldArray = new UserInputField[1];

        fieldArray[0] = new UserInputField("pnum", "kill diffusant",
                new Double(arrayNum), fieldWidth, this);

        prompt = new UserInputPrompt("KillDiffusant", "Kill Diffusant", fieldArray, this);

    }

    public boolean killDiffusantPrompt2() {

        if (fieldArray.length != 1) {
            return false;
        }

        int pnum = (int) getFieldValue(0);

        boolean dd = Master.project.killDiffusant(pnum);

        Master.updatePanel2D();

        return dd;

    }

    public void addSourceImpulsePrompt() {

        int pnum = 0, x = 1, y = 1, z = 1, molecules = 0;
        double c = 1, t = 0;
        int plimit = Master.project.diffusants.length;

        fieldArray = new UserInputField[7];

        fieldArray[0] = new UserInputField("pnum",
                "add source to diffusant number ( < " +
                plimit +
                ")",
                new Double(pnum), fieldWidth, this);
        fieldArray[1] = new UserInputField("x", "x voxel location",
                new Double(x), fieldWidth, this);
        fieldArray[2] = new UserInputField("y", "y voxel location",
                new Double(y), fieldWidth, this);
        fieldArray[3] = new UserInputField("z", "z voxel location",
                new Double(z), fieldWidth, this);
        fieldArray[4] = new UserInputField("t", "time of release",
                new Double(t), fieldWidth, this);
        fieldArray[5] = new UserInputField("c", "concentration (mM)",
                new Double(c), fieldWidth, this);
        fieldArray[6] = new UserInputField("molecules",
                "OR, specify number of molecules",
                new Double(molecules), fieldWidth, this);

        prompt = new UserInputPrompt("AddSourceImpulse", "Add Quanta Source",
                fieldArray, this);

    }

    public void addSourceImpulsePrompt2() {

        if (fieldArray.length != 7) {
            return;
        }

        int pnum = (int) getFieldValue(0);

        CoordinatesVoxels coordinates = new CoordinatesVoxels(Master.project);

        coordinates.set("xVoxel1", (int) getFieldValue(1));
        coordinates.set("yVoxel1", (int) getFieldValue(2));
        coordinates.set("zVoxel1", (int) getFieldValue(3));

        double t = getFieldValue(4);
        double ctotal = getFieldValue(5);
        int molecules = (int) getFieldValue(6);

        if (molecules > 0) {
            ctotal = Utility.N2C(molecules, Master.project.dx, coordinates.spaceVoxels);
        }

        Master.addSourceImpulse(pnum, coordinates, t, ctotal, "C");

    }

    public void addSourceGaussPrompt() {

        int pnum = 0, x0 = 1, y0 = 1, z0 = 1, x1 = 1, y1 = 1, z1 = 1;
        double tpeak = 1, stdv = 1; // ms
        double cTotal = 1; // mM
        int plimit = Master.project.diffusants.length;

        fieldArray = new UserInputField[10];

        fieldArray[0] = new UserInputField("pnum",
                "add source to diffusant number ( < " +
                plimit +
                ")",
                new Double(pnum), fieldWidth, this);
        fieldArray[1] = new UserInputField("x0", "corner 1, x0 voxel",
                new Double(x0), fieldWidth, this);
        fieldArray[2] = new UserInputField("y0", "corner 1, y0 voxel",
                new Double(y0), fieldWidth, this);
        fieldArray[3] = new UserInputField("z0", "corner 1, z0 voxel",
                new Double(z0), fieldWidth, this);
        fieldArray[4] = new UserInputField("x1", "corner 2, x1 voxel",
                new Double(x1), fieldWidth, this);
        fieldArray[5] = new UserInputField("y1", "corner 2, y1 voxel",
                new Double(y1), fieldWidth, this);
        fieldArray[6] = new UserInputField("z1", "corner 2, z1 voxel",
                new Double(z1), fieldWidth, this);
        fieldArray[7] = new UserInputField("time2peak", "time of peak (ms)",
                new Double(tpeak), fieldWidth, this);
        fieldArray[8] = new UserInputField("width", "Gaussian stdv (ms)",
                new Double(stdv), fieldWidth, this);
        fieldArray[9] = new UserInputField("release rate", "Ctotal (mM)",
                new Double(cTotal), fieldWidth, this);

        prompt = new UserInputPrompt("AddSourceGauss", "Add Gaussian Source",
                fieldArray, this);

    }

    public void addSourceGaussPrompt2() {

        if (fieldArray.length != 10) {
            return;
        }

        int pnum = (int) getFieldValue(0);

        CoordinatesVoxels coordinates = new CoordinatesVoxels(Master.project);

        coordinates.set("xVoxel1", (int) getFieldValue(1));
        coordinates.set("yVoxel1", (int) getFieldValue(2));
        coordinates.set("zVoxel1", (int) getFieldValue(3));
        coordinates.set("xVoxel2", (int) getFieldValue(4));
        coordinates.set("yVoxel2", (int) getFieldValue(5));
        coordinates.set("zVoxel2", (int) getFieldValue(6));

        double tpeak = getFieldValue(7);
        double stdv = getFieldValue(8);
        double cTotal = getFieldValue(9);

        Master.addSourceGauss(pnum, coordinates, tpeak, stdv, cTotal);

    }

    public void killSourcePrompt(int arrayNum) {

        Source[] source = Master.project.sources;

        if (source == null) {
            return;
        }

        fieldArray = new UserInputField[1];

        if (arrayNum >= source.length) {
            arrayNum = 0;
        }

        fieldArray[0] = new UserInputField("num", "kill source num ",
                new Double(arrayNum), fieldWidth, this);

        prompt = new UserInputPrompt("KillSource", "Kill Source", fieldArray, this);

    }

    public boolean killSourcePrompt2() {

        if (fieldArray.length != 1) {
            return false;
        }

        int pnum = (int) getFieldValue(0);

        boolean s = Master.project.killSource(pnum);

        Master.updatePanel2D();

        return s;

    }

    public void addDetectorPrompt() {

        int pnum = 0, x0 = 1, y0 = 1, z0 = 1, x1 = 2, y1 = 2, z1 = 2;
        int plimit = Master.project.diffusants.length;

        fieldArray = new UserInputField[7];

        fieldArray[0] = new UserInputField("pnum",
                "detect diffusant number ( < " + plimit + ")", new Double(pnum), fieldWidth, this);
        fieldArray[1] = new UserInputField("x0", "corner 1, x0 voxel",
                new Double(x0), fieldWidth, this);
        fieldArray[2] = new UserInputField("y0", "corner 1, y0 voxel",
                new Double(y0), fieldWidth, this);
        fieldArray[3] = new UserInputField("z0", "corner 1, z0 voxel",
                new Double(z0), fieldWidth, this);
        fieldArray[4] = new UserInputField("x1", "corner 2, x1 voxel",
                new Double(x1), fieldWidth, this);
        fieldArray[5] = new UserInputField("y1", "corner 2, y1 voxel",
                new Double(y1), fieldWidth, this);
        fieldArray[6] = new UserInputField("z1", "corner 2, z1 voxel",
                new Double(z1), fieldWidth, this);

        prompt = new UserInputPrompt("AddDetector", "Add Detector", fieldArray, this);

    }

    public void addDetectorPrompt2() {

        if (fieldArray.length != 7) {
            return;
        }

        int pnum = (int) getFieldValue(0);

        CoordinatesVoxels coordinates = new CoordinatesVoxels(Master.project);

        coordinates.set("xVoxel1", (int) getFieldValue(1));
        coordinates.set("yVoxel1", (int) getFieldValue(2));
        coordinates.set("zVoxel1", (int) getFieldValue(3));
        coordinates.set("xVoxel2", (int) getFieldValue(4));
        coordinates.set("yVoxel2", (int) getFieldValue(5));
        coordinates.set("zVoxel2", (int) getFieldValue(6));

        Master.addDetectorAvg("", pnum, coordinates);

    }

    public void addDetectorWeightedPrompt() {

        int pnum = 0, x0 = 0, y0 = 0, z0 = 0;
        int x1 = Master.project.geometry.xVoxels - 1;
        int y1 = Master.project.geometry.yVoxels - 1;
        int z1 = Master.project.geometry.zVoxels - 1;
        double gaussSDxy = 0.5, gaussSDz = 0.5;
        int plimit = Master.project.diffusants.length;

        fieldArray = new UserInputField[9];

        fieldArray[0] = new UserInputField("pnum",
                "detect diffusant number ( < " + plimit +
                ")",
                new Double(pnum), fieldWidth, this);
        fieldArray[1] = new UserInputField("x0", "corner 1, x0 voxel",
                new Double(x0), fieldWidth, this);
        fieldArray[2] = new UserInputField("y0", "corner 1, y0 voxel",
                new Double(y0), fieldWidth, this);
        fieldArray[3] = new UserInputField("z0", "corner 1, z0 voxel",
                new Double(z0), fieldWidth, this);
        fieldArray[4] = new UserInputField("x1", "corner 2, x1 voxel",
                new Double(x1), fieldWidth, this);
        fieldArray[5] = new UserInputField("y1", "corner 2, y1 voxel",
                new Double(y1), fieldWidth, this);
        fieldArray[6] = new UserInputField("z1", "corner 2, z1 voxel",
                new Double(z1), fieldWidth, this);
        fieldArray[7] = new UserInputField("SDxy", "Gaussian STDV xy plane",
                new Double(gaussSDxy), fieldWidth, this);
        fieldArray[8] = new UserInputField("SDz", "Gaussian STDV z plane",
                new Double(gaussSDz), fieldWidth, this);

        prompt = new UserInputPrompt("AddDetectorWeighted",
                "Add Detector Gaussian Weighted", fieldArray, this);

    }

    public void addDetectorWeightedPrompt2() {

        if (fieldArray.length != 9) {
            return;
        }

        int pnum = (int) getFieldValue(0);

        CoordinatesVoxels coordinates = new CoordinatesVoxels(Master.project);

        coordinates.set("xVoxel1", (int) getFieldValue(1));
        coordinates.set("yVoxel1", (int) getFieldValue(2));
        coordinates.set("zVoxel1", (int) getFieldValue(3));
        coordinates.set("xVoxel2", (int) getFieldValue(4));
        coordinates.set("yVoxel2", (int) getFieldValue(5));
        coordinates.set("zVoxel2", (int) getFieldValue(6));

        double SDxy = getFieldValue(7);
        double SDz = getFieldValue(8);

        Master.addDetectorPSF_Gauss("", pnum, coordinates, SDxy, SDxy, SDz);

    }

    public void addDetectorSnapPrompt() {
        int pnum = 0, x0 = 0, y0 = 0, z0 = 0, z1 = 0;
        int x1 = Master.project.geometry.xVoxels - 1;
        int y1 = Master.project.geometry.yVoxels - 1;
        double t = 0;
        int plimit = Master.project.diffusants.length;

        fieldArray = new UserInputField[8];

        fieldArray[0] = new UserInputField("pnum",
                "detect diffusant number ( < " + plimit +
                ")", new Double(pnum), fieldWidth, this);
        fieldArray[1] = new UserInputField("x0", "corner 1, x0 voxel",
                new Double(x0), fieldWidth, this);
        fieldArray[2] = new UserInputField("y0", "corner 1, y0 voxel",
                new Double(y0), fieldWidth, this);
        fieldArray[3] = new UserInputField("z0", "corner 1, z0 voxel",
                new Double(z0), fieldWidth, this);
        fieldArray[4] = new UserInputField("x1", "corner 2, x1 voxel",
                new Double(x1), fieldWidth, this);
        fieldArray[5] = new UserInputField("y1", "corner 2, y1 voxel",
                new Double(y1), fieldWidth, this);
        fieldArray[6] = new UserInputField("z1", "corner 2, z1 voxel",
                new Double(z1), fieldWidth, this);
        fieldArray[7] = new UserInputField("t", "take snapshot at time (ms)",
                new Double(t), fieldWidth, this);

        prompt = new UserInputPrompt("AddDetectorSnap", "Add Snapshot Detector",
                fieldArray, this);

    }

    public void addDetectorSnapPrompt2() {

        if (fieldArray.length != 8) {
            return;
        }

        int pnum = (int) getFieldValue(0);

        CoordinatesVoxels coordinates = new CoordinatesVoxels(Master.project);

        coordinates.set("xVoxel1", (int) getFieldValue(1));
        coordinates.set("yVoxel1", (int) getFieldValue(2));
        coordinates.set("zVoxel1", (int) getFieldValue(3));
        coordinates.set("xVoxel2", (int) getFieldValue(4));
        coordinates.set("yVoxel2", (int) getFieldValue(5));
        coordinates.set("zVoxel2", (int) getFieldValue(6));

        double t = getFieldValue(7);

        Master.addDetectorSnapshot("", pnum, coordinates, t);

    }

    public void killDetectorPrompt(int arrayNum) {

        Detector[] detector = Master.project.detectors;

        if (detector == null) {
            return;
        }

        if (arrayNum >= detector.length) {
            arrayNum = 0;
        }

        fieldArray = new UserInputField[1];

        fieldArray[0] = new UserInputField("num", "kill detector num ",
                new Double(arrayNum), fieldWidth, this);

        prompt = new UserInputPrompt("KillDetector", "Kill Detector", fieldArray, this);

    }

    public boolean killDetectorPrompt2() {

        if (fieldArray.length != 1) {
            return false;
        }

        int pnum = (int) getFieldValue(0);

        boolean dd = Master.project.killDetector(pnum);

        Master.updatePanel2D();

        return dd;

    }

    public void addDiffusantPhotolysisPrompt() {

        int dnum = 0, compound = 0;
        double c = 8, dd = 0.25;
        int plimit = Master.project.diffusants.length;

        fieldArray = new UserInputField[4];

        fieldArray[0] = new UserInputField("dnum",
                "photolysis compound reacts with diffusant # ( < " +
                plimit + ")", new Double(dnum), fieldWidth, this);
        fieldArray[1] = new UserInputField("c", "concentration (mM)",
                new Double(c), fieldWidth, this);
        fieldArray[2] = new UserInputField("d", "diffusion coefficient",
                new Double(dd), fieldWidth, this);
        fieldArray[3] = new UserInputField("psf",
                "(0) Gaussian PSF (1) Torok PSF",
                new Double(compound), fieldWidth, this);

        prompt = new UserInputPrompt("addDiffusantPhotolysis",
                "Add Photolysis Diffusant", fieldArray, this);

    }

    public void addDiffusantPhotolysisPrompt2() {

        PSF p = null;
        PulseTimer pt = new PulseTimer(Master.project, 0.1, 0.1, 1.0, "1/ms");

        if (fieldArray.length != 4) {
            return;
        }

        int dnum = (int) getFieldValue(0);
        double c = (int) getFieldValue(1);
        double dd = getFieldValue(2);
        double sdxy = 0.1, sdz = 0.1;
        int ptype = (int) getFieldValue(3);

        if (ptype == 1) {
            p = new PSFgauss(Master.project, null, sdxy, sdxy, sdz);
        } else {
            p = new PSFtorok(Master.project, null);
        }

        Master.project.addDiffusant(new DiffusantPhoto(Master.project, "LightReactant", c, dd, null, dnum, pt, p));

        Master.updatePanel2D();
        Master.resetAllParamVectors();

    }

    public void killBatchPrompt(int arrayNum) {

        Batch[] batch = Master.project.batches;

        if (batch == null) {
            return;
        }

        if (arrayNum >= batch.length) {
            arrayNum = 0;
        }

        fieldArray = new UserInputField[1];

        fieldArray[0] = new UserInputField("pnum", "kill batch", new Double(arrayNum), fieldWidth, this);

        prompt = new UserInputPrompt("KillBatch", "Kill Batch", fieldArray, this);

    }

    public boolean killBatchPrompt2() {

        if (fieldArray.length != 1) {
            return false;
        }

        int pnum = (int) getFieldValue(0);

        return Master.project.killBatch(pnum);

    }

    public void writePSF(String fName) {

        int imax = Master.project.geometry.xVoxels;
        int jmax = Master.project.geometry.yVoxels;
        int kmax = Master.project.geometry.zVoxels;

        if (!Master.photoExists() || (Master.getPhoto().psf == null)) {
            Master.log("photolysis diffusant PSF does not appear to exist.");
            return;
        }

        fieldArray = new UserInputField[4];

        fieldArray[0] = new UserInputField("iMax", "x-voxels", new Double(imax),
                fieldWidth, this);
        fieldArray[1] = new UserInputField("jMax", "y-voxels", new Double(jmax),
                fieldWidth, this);
        fieldArray[2] = new UserInputField("kMax", "z-voxels", new Double(kmax),
                fieldWidth, this);
        fieldArray[3] = new UserInputField(fName, "", new Double(Double.NaN),
                0, this);

        prompt = new UserInputPrompt("writePSF", "Write Random Access PSF",
                fieldArray, this);

    }

    public void writePSF2() {

        if (fieldArray.length != 4) {
            return;
        }

        int imax = (int) getFieldValue(0);
        int jmax = (int) getFieldValue(1);
        int kmax = (int) getFieldValue(2);

        String fName = getFieldName(3);

        Master.writePSF(fName, imax, jmax, kmax);

    }

}
