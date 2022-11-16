package ucl.silver.d3d.core;

import java.io.*;
import ucl.silver.d3d.init.*;
import ucl.silver.d3d.utils.*;

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
public class Project extends ParamVector {

    public String initClassAndFunction = null; // (e.g. "InitProject.initCube")
    public String file = null; // file name, if this Project was opened from previously saved class
    public String psfFile = null; // psf file to open
    public String geomFile = null; // geometry file to open

    public double simTime = 1; // simulation time (ms)
    public double dt = 0.1; //Double.NaN; // time step (ms), computed by RunFiniteDifference or RunMonteCarlo
    public double dx = 1; // voxel width (um)
    
    public double simTemp = Double.NaN; // simulation temperature // can be used for Q10 scaling

    public double saveRate = 100; // sample rate of output files (kHz)
    public double saveDT = 1 / saveRate; // sample rate of output files (ms)

    public double printRate = 1; // rate simulation results are printed (kHz)
    public double printDT = 1 / printRate; // rate simulation results are printed (ms)

    public int batchNum = -1; // file batch number, use negative number (-1) for no batch number
    private int batchCounter = -1; // Batch array counter

    public String seedMT = ""; // seed for Mersenne Twister random number generator

    public String timeUnits = "ms";
    public String freqUnits = "kHz"; // inverse time
    public String spaceUnits = "μm";
    public String concUnits = "mM";
    public String currentUnits = "pA"; // nS * mV = pA
    public String voltUnits = "mV";
    public String conductanceUnits = "nS";
    public String resistanceUnits = "GΩ"; // 1 / nS
    public String tempUnits = "°C"; // temperature

    public String directory = ""; // diretory where folder/files are saved
    public String folder = "D3Doutput";
    private String subfolder = ""; // subdirectory for batches

    public String date = "";

    public InitProject initProject = null;
    public Geometry geometry = null;
    public Diffusant[] diffusants = null;
    public Source[] sources = null;
    public Detector[] detectors = null;
    public Batch[] batches = null;
    private Batch saveBatch = null;
    public D3Derror[] errors = null;

    public double[] sourceArray = null;

    public RunFiniteDifference finiteDifference = null;
    public RunMonteCarlo monteCarlo = null;

    public transient StopWatch timer = new StopWatch();

    public int EM_iseries = 1; // MC AZEM sims
    //public int EM_azSelect = 0; // MC AZEM sims

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("simTime")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("simTemp")) {
            return tempUnits;
        }
        if (name.equalsIgnoreCase("dt")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("dx")) {
            return spaceUnits;
        }
        if (name.equalsIgnoreCase("saveRate")) {
            return freqUnits;
        }
        if (name.equalsIgnoreCase("saveDT")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("printRate")) {
            return freqUnits;
        }
        if (name.equalsIgnoreCase("printDT")) {
            return timeUnits;
        }
        if (name.equalsIgnoreCase("directory")) {
            return "DIR";
        }
        return super.units(name);
    }

    @Override
    public String help(String name) {
        if (name.equalsIgnoreCase("name")) {
            return "project name";
        }
        if (name.equalsIgnoreCase("initClassAndFunction")) {
            return "initialization Class and function";
        }
        if (name.equalsIgnoreCase("simTime")) {
            return "simulation time";
        }
        if (name.equalsIgnoreCase("simTemp")) {
            return "simulation temperature";
        }
        if (name.equalsIgnoreCase("dt")) {
            return "time step";
        }
        if (name.equalsIgnoreCase("dx")) {
            return "voxel cube dimensions";
        }
        if (name.equalsIgnoreCase("timeUnits")) {
            return "time units";
        }
        if (name.equalsIgnoreCase("freqUnits")) {
            return "inverse time units";
        }
        if (name.equalsIgnoreCase("spaceUnits")) {
            return "space dimension units";
        }
        if (name.equalsIgnoreCase("concUnits")) {
            return "concentration units";
        }
        if (name.equalsIgnoreCase("saveRate")) {
            return "default rate for saving data";
        }
        if (name.equalsIgnoreCase("saveDT")) {
            return "default time step for saving data";
        }
        if (name.equalsIgnoreCase("printRate")) {
            return "rate of simulation progress display";
        }
        if (name.equalsIgnoreCase("printDT")) {
            return "time step of simulation progress display";
        }
        if (name.equalsIgnoreCase("seedMT")) {
            return "seed for Mersenne Twister random number generator";
        }
        return super.help(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("file")) {
            return false;
        }
        if (name.equalsIgnoreCase("initClassAndFunction")) {
            return false;
        }
        if (name.equalsIgnoreCase("psfFile")) {
            return false;
        }
        if (name.equalsIgnoreCase("geomFile")) {
            return false;
        }
        if (name.equalsIgnoreCase("dt")) {
            return false;
        }
        return super.canEdit(name);
    }

    public Project() {
        super(null);
        geometry = new Geometry(this, 10, 10, 10);
        finiteDifference = new RunFiniteDifference(this);
        name = "D3D Finite Difference";
        if (directory.isEmpty()) {
            directory = System.getProperty("user.home");
        }
        createVector(true);
    }

    public Project(String type) {

        super(null);

        geometry = new Geometry(this, 10, 10, 10);

        if (type.equalsIgnoreCase("FD")) {
            newFiniteDifference();
        } else if (type.equalsIgnoreCase("MC")) {
            newMonteCarlo();
        } else if (type.equalsIgnoreCase("MCAZ")) {
            newMonteCarloAZ();
        }

        createVector(true);

    }

    public RunFiniteDifference newFiniteDifference() {
        monteCarlo = null;
        finiteDifference = new RunFiniteDifference(this);
        name = "D3D Finite Difference";
        Master.log("started Finite Difference project");
        return finiteDifference;
    }

    public RunMonteCarlo newMonteCarlo() {
        finiteDifference = null;
        monteCarlo = new RunMonteCarlo(this);
        name = "D3D Monte Carlo";
        Master.log("started Monte Carlo project");
        return monteCarlo;
    }
    
    public RunMonteCarloPhoto newMonteCarloPhoto() {
        finiteDifference = null;
        RunMonteCarloPhoto mc = new RunMonteCarloPhoto(this);
        monteCarlo = mc;
        name = "D3D Monte Carlo Photolysis";
        Master.log("started Monte Carlo Photolysis project");
        return mc;
    }
    
    public RunMonteCarloPhotoMito newMonteCarloPhotoMito() {
        finiteDifference = null;
        RunMonteCarloPhotoMito mc = new RunMonteCarloPhotoMito(this);
        monteCarlo = mc;
        name = "D3D Monte Carlo with Mitochondria";
        Master.log("started Monte Carlo Mito project");
        return mc;
    }

    public RunMonteCarloAZ newMonteCarloAZ() {
        finiteDifference = null;
        RunMonteCarloAZ mc = new RunMonteCarloAZ(this);
        monteCarlo = mc;
        name = "D3D Monte Carlo Active Zone";
        Master.log("started Monte Carlo Active Zone project");
        return mc;
    }

    public RunMonteCarloAZEM newMonteCarloAZEM() {
        finiteDifference = null;
        RunMonteCarloAZEM mc = new RunMonteCarloAZEM(this);
        monteCarlo = mc;
        name = "D3D Monte Carlo Active Zone EM";
        Master.log("started Monte Carlo Active Zone EM project");
        return mc;
    }

    public void setDate(){
        date = Master.currentDate();
        updateVectors();
    }

    public void nullArrays() {
        diffusants = null;
        sources = null;
        detectors = null;
        batches = null;
        errors = null;
    }
    
    public void setDirectoryJason() {
        
        String os = System.getProperty("os.name");
        
        if (os.startsWith("Mac")) {
            directory = "/Users/jason/Documents/D3D/Simulations/";
        } else {
            directory = "/Jason/D3D/Simulations/";
        }
        
    }

    public void init() {

        setDate();
        
        setParamError("dt", null);
        
        geometry.init();

        if ((monteCarlo == null) && (finiteDifference == null)) {
            return;
        }

        if (monteCarlo != null) {
            monteCarlo.init();
        }
        
        if (finiteDifference != null) {
            finiteDifference.init();
        }
        
        if ((dt <= 0 ) || (Double.isNaN(dt)) || (Double.isInfinite(dt))) {
            //error("Project.init", "dt", "bad value");
        }

        initDiffusants(-1);
        initSources(-1);
        initDetectors(-1);
        initBatches();

        updateVectors();

        Master.updateMainFrameTitle();

        //Master.log("initialized project: " + name);

    }

    public void setOutputRate(double newRate) {

        if (newRate > 0) {

            saveRate = newRate;
            saveDT = 1.0 / newRate;

            if (diffusants != null) {
                for (int i = 0; i < diffusants.length; i += 1) {
                    diffusants[i].save.setOutputRate(newRate);
                }
            }

            if (detectors != null) {
                for (int i = 0; i < detectors.length; i += 1) {
                    detectors[i].save.setOutputRate(newRate);
                }
            }

            if (sources != null) {
                for (int i = 0; i < sources.length; i += 1) {
                    sources[i].save.setOutputRate(newRate);
                }
            }

            if (monteCarlo != null) {
                monteCarlo.setOutputRate(newRate);
            }

            updateVectors();

        }

    }

    public void setPrintRate(double newRate){

        if (newRate > 0) {
            printRate = newRate;
            printDT = 1.0 / newRate;
            updateVectors();
        }

    }

    public void setSaveRate(double newRate){

        if (newRate > 0) {
            saveRate = newRate;
            saveDT = 1.0 / newRate;
            updateVectors();
        }

    }
    
    public int simPoints() {
        return (int) (simTime / dt) + 1;
    }

    public String fullDirectory() {

        String dir = directory;

        if (folder.length() > 0) {
            dir += folder + "/";
        }

        if (subfolder.length() > 0) {
            dir += subfolder + "/";
        }

        return dir;

    }

    public boolean checkDirectory() {

        File f = new File(fullDirectory());

        if (!f.isDirectory()) {
            f.mkdirs();
        }

        return f.isDirectory();

    }

    //
    // DIFFUSANT FUNCTIONS
    //

    public boolean reactionExists() {
        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                if ((d.reaction)) {
                    return true;
                }
            }
        }
        return false;
    }

    public int numDiffusants() {
        if (diffusants == null) {
            return 0;
        }
        return diffusants.length;
    }

    public boolean checkDiffusantNum(int arrayNum) {
        return ((diffusants != null) && (arrayNum >= 0) && (arrayNum < diffusants.length));
    }

    public String diffusantName(int arrayNum) {
        if ((diffusants != null) && (arrayNum >= 0) && (arrayNum < diffusants.length)) {
            return diffusants[arrayNum].name;
        }
        return null;
    }

    public void diffusantsUpdateVectors() {
        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                d.updateVectors();
            }
        }
    }

    public void checkDiffusants() {
        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                d.check();
            }
        }
    }

    public void initDiffusantsPulseTimer() {
        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                d.initPulseTimer();
            }
        }
    }

    public void openDiffusants() {

        checkDirectory();

        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                d.saveInit();
            }
        }

    }

    public void closeDiffusants() {
        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                d.saveFinish();
            }
        }
    }

    public void initDiffusants(int dNum) {

        int ibgn = dNum;
        int iend = dNum;

        if (diffusants == null) {
            return;
        }
        
        if ((ibgn < 0) || (iend >= diffusants.length)) {
            ibgn = 0;
            iend = diffusants.length - 1;
        }

        for (int i = ibgn; i <= iend; i++) {

            if (diffusants[i].name.length() == 0) {
                diffusants[i].name = "Diffusant" + Integer.toString(i);
            }

            diffusants[i].init();

        }

    }

    public Diffusant getDiffusant(int diffusantNum) {

        if (diffusants == null) {
            return null;
        }

        if ((diffusantNum >= 0) && (diffusantNum < diffusants.length)) {
            return diffusants[diffusantNum];
        }

        return null;

    }

    public Diffusant getDiffusant(String diffusantName) {

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            return null;
        }

        if (diffusants == null) {
            return null;
        }

        for (Diffusant d : diffusants) {

            if (d == null) {
                continue;
            }

            if (d.name.equalsIgnoreCase(diffusantName)) {
                return d;
            }

        }

        return null;

    }

    public int getDiffusantNum(String diffusantName) {

        if ((diffusantName == null) || (diffusantName.length() == 0)) {
            return -1;
        }

        if (diffusants == null) {
            return -1;
        }

        for (int i = 0; i < diffusants.length; i++) {

            if (diffusants[i] == null) {
                continue;
            }

            if (diffusants[i].name.equalsIgnoreCase(diffusantName)) {
                return i;
            }

        }

        return -1;

    }

    public int addDiffusant(Diffusant newDiffusant) {

        int i = 0;

        if (diffusants != null) {
            i = diffusants.length;
        }

        Diffusant[] newArray = new Diffusant[i+1];

        if (i > 0) {
            System.arraycopy(diffusants, 0, newArray, 0, i);
        }

        newArray[i] = newDiffusant;

        diffusants = newArray; // replace old array with new one

        Master.log("added Diffusant #" + i + " : " + newDiffusant.name);

        initDiffusants(i);

        return i;

    }
    
    public void killDiffusantsAll() {
        int i = diffusants.length;
        diffusants = null;
        if (i == 1) {
            Master.log("killed Diffusant #0");
        } else if (i > 1) {
            Master.log("killed Diffusants #0-" + i);
        }
    }

    public boolean killDiffusant(int i) {

        int k = 0;

        if (!checkDiffusantNum(i)) {
            return false;
        }

        if (diffusants.length == 1) {
            diffusants = null;
            Master.log("killed Diffusant #0");
            return true;
        }

        Diffusant[] part = new Diffusant[diffusants.length - 1]; // new array

        for (int j = 0; j < diffusants.length; j++) {

            if (j == i) {
                continue;
            }

            part[k] = diffusants[j];
            k++;

        }

        diffusants = part; // replace old array with new one

        Master.log("killed Diffusant #" + i);

        return true;

    }

    public boolean setDiffusantParam(int arrayNum, String varName, double v) {

        if ((diffusants == null) || (arrayNum < 0) || (arrayNum >= diffusants.length)) {
            return false;
        }

        return diffusants[arrayNum].set(varName, v);

    }

    public boolean setDiffusantParam(int arrayNum, String varName, String s) {

        if ((diffusants == null) || (arrayNum < 0) || (arrayNum >= diffusants.length)) {
            return false;
        }

        return diffusants[arrayNum].set(varName, s);

    }

    //
    // SOURCE FUNCTIONS
    //

    public boolean sourceExists() {
        return (sources != null);
    }

    public int numSources() {

        if (sources == null) {
            return 0;
        }

        return sources.length;

    }

    public boolean checkSourceNum(int i) {
        return ((sources != null) && (i >= 0) && (i < sources.length));
    }

    public void sourcesUpdateVectors() {
        if (sources != null) {
            for (Source s : sources) {
                s.updateVectors();
            }
        }
    }

    public void checkSources() {
        if (sources != null) {
            for (Source s : sources) {
                s.check();
            }
        }
    }

    public void initSourcesPulseTimer() {
        if (sources != null) {
            for (Source s : sources) {
                s.initPulseTimer();
            }
        }
    }

    public void openSources() {

        checkDirectory();

        if (sources != null) {
            for (Source s : sources) {
                s.saveInit();
            }
        }

    }

    public void closeSources() {
        if (sources != null) {
            for (Source s : sources) {
                s.saveFinish();
            }
        }
    }

    public void initSources(int sNum) {
        
        int ibgn = sNum;
        int iend = sNum;

        if (sources == null) {
            return;
        }

        if ((sNum < 0) || (sNum >= sources.length)) {
            ibgn = 0;
            iend = sources.length - 1;
        }

        for (int i = ibgn; i <= iend; i++) {

            if (sources[i].name.length() == 0) {
                sources[i].name = "Source" + Integer.toString(i);
            }

            sources[i].init();

        }

    }

    public Source getSource(int sourceNum) {

        if (sources == null) {
            return null;
        }

        if ((sourceNum >= 0) && (sourceNum < sources.length)) {
            return sources[sourceNum];
        }

        return null;

    }

    public Source getSource(String sourceName) {

        if ((sourceName == null) || (sourceName.length() == 0)) {
            return null;
        }

        if (sources == null) {
            return null;
        }

        for (Source s : sources) {

            if (s == null) {
                continue;
            }

            if (s.name.equalsIgnoreCase(sourceName)) {
                return s;
            }

        }

        return null;

    }

    public int getSourceNum(String sourceName) {

        if ((sourceName == null) || (sourceName.length() == 0)) {
            return -1;
        }

        if (sources == null) {
            return -1;
        }

        for (int i = 0; i < sources.length; i++) {

            if (sources[i] == null) {
                continue;
            }

            if (sources[i].name.equalsIgnoreCase(sourceName)) {
                return i;
            }

        }

        return -1;

    }

    public int addSource(Source newSource) {

        int i = 0;

        if (sources != null) {
            i = sources.length;
        }

        Source[] newArray = new Source[i+1];
        
        if (i > 0) {
            System.arraycopy(sources, 0, newArray, 0, i);
        }

        newArray[i] = newSource;

        sources = newArray; // replace old array with new one

        Master.log("added Source #" + i + " : " + newSource.name);

        initSources(i);

        return i;

    }
    
    public void killSourcesAll() {
        int i = sources.length;
        sources = null;
        if (i == 1) {
            Master.log("killed Source #0");
        } else if (i > 1) {
            Master.log("killed Sources #0-" + i);
        }
    }

    public boolean killSource(int i) {
        int k = 0;

        if (!checkSourceNum(i)) {
            return false;
        }

        if (sources.length == 1) {
            sources = null;
            Master.log("killed Source #0");
            return true;
        }

        Source[] s = new Source[sources.length - 1]; // new array

        for (int j = 0; j < sources.length; j++) {
            if (j == i) {
                continue;
            }
            s[k] = sources[j];
            k++;
        }

        sources = s; // replace old array with new one

        Master.log("killed Source #" + i);

        return true;

    }
    
    public boolean insideAnySource(int xVoxel, int yVoxel, int zVoxel, int skip) {
        
        if (sources == null) {
            return false;
        }
        
        for (int is = 0; is < sources.length; is++) {

            if ((is == skip) || (sources[is] == null)) {
                continue;
            }

            if (sources[is].isInside(xVoxel, yVoxel, zVoxel)) {
                return true;
            }
        }
        
        return false;
        
    }

    public boolean setSourceParam(int arrayNum, String varName, double v) {
        if ((sources == null) || (arrayNum < 0) || (arrayNum >= sources.length)) {
            return false;
        }
        return sources[arrayNum].set(varName, v);
    }

    public boolean setSourceParam(int arrayNum, String varName, String s) {
        if ((sources == null) || (arrayNum < 0) || (arrayNum >= sources.length)) {
            return false;
        }
        return sources[arrayNum].set(varName, s);
    }

    //
    // DETECTOR FUNCTIONS
    //
    public boolean detectorExists() {
        return (detectors != null);
    }

    public int numDetectors() {
        if (detectors == null) {
            return 0;
        }
        return detectors.length;
    }

    public boolean checkDetectorNum(int i) {
        return ((detectors != null) && (i >= 0) && (i < detectors.length));
    }

    public void detectorsUpdateVectors() {
        if (detectors != null) {
            for (Detector d : detectors) {
                d.updateVectors();
            }
        }
    }

    public void checkDetectors() {
        if (detectors != null) {
            for (Detector d : detectors) {
                d.check();
            }
        }
    }

    public void initDetectorsPulseTimer() {
        if (detectors != null) {
            for (Detector d : detectors) {
                d.initPulseTimer();
            }
        }
    }

    public void openDetectors() {

        checkDirectory();

        if (detectors != null) {
            for (Detector d : detectors) {
                d.saveInit();
            }
        }

    }

    public void closeDetectors() {
        if (detectors != null) {
            for (Detector d : detectors) {
                d.saveFinish();
            }
        }
    }

    public void initDetectors(int dNum) {

        int ibgn = dNum;
        int iend = dNum;

        if (detectors == null) {
            return;
        }

        if ((dNum < 0) || (dNum >= detectors.length)) {
            ibgn = 0;
            iend = detectors.length - 1;
        }

        for (int i = ibgn; i <= iend; i++) {

            if (detectors[i].name.length() == 0) {
                detectors[i].name = "Detector" + Integer.toString(i);
            }

            detectors[i].init();

        }

    }

    public Detector getDetector(int detectorNum) {

        if (detectors == null) {
            return null;
        }

        if ((detectorNum >= 0) && (detectorNum < detectors.length)) {
            return detectors[detectorNum];
        }

        return null;

    }

    public Detector getDetector(String detectorName) {

        if ((detectorName == null) || (detectorName.length() == 0)) {
            return null;
        }

        if (detectors == null) {
            return null;
        }

        for (Detector d : detectors) {

            if (d == null) {
                continue;
            }

            if (d.name.equalsIgnoreCase(detectorName)) {
                return d;
            }

        }

        return null;

    }

    public int getDetectorNum(String detectorName) {

        if ((detectorName == null) || (detectorName.length() == 0)) {
            return -1;
        }

        if (detectors == null) {
            return -1;
        }

        for (int i = 0; i < detectors.length; i++) {

            if (detectors[i] == null) {
                continue;
            }

            if (detectors[i].name.equalsIgnoreCase(detectorName)) {
                return i;
            }

        }

        return -1;

    }

    public int addDetector(Detector newDetector) {

        int i = 0;

        if (detectors != null) {
            i = detectors.length;
        }

        Detector[] newArray = new Detector[i+1];
        
        if (i > 0) {
            System.arraycopy(detectors, 0, newArray, 0, i);
        }

        newArray[i] = newDetector;

        detectors = newArray; // replace old array with new one

        Master.log("added detector #" + i + " : " + newDetector.name);

        initDetectors(i);

        return i;

    }
    
    public void killDetectorsAll() {
        int i = detectors.length;
        detectors = null;
        if (i == 1) {
            Master.log("killed Detector #0");
        } else if (i > 1) {
            Master.log("killed Detectors #0-" + i);
        }
    }

    public boolean killDetector(int i) {
        int k = 0;

        if (!checkDetectorNum(i)) {
            return false;
        }

        if (detectors.length == 1) {
            detectors = null;
            Master.log("killed Detector #0");
            return true;
        }

        Detector[] detect = new Detector[detectors.length - 1]; // new array

        for (int j = 0; j < detectors.length; j++) {
            if (j == i) {
                continue;
            }
            detect[k] = detectors[j];
            k++;
        }

        detectors = detect; // replace old array with new one

        Master.log("killed Detector #" + i);

        return true;

    }

    public boolean setDetectorParam(int arrayNum, String varName, double v) {
        if ((detectors == null) || (arrayNum < 0) || (arrayNum >= detectors.length)) {
            return false;
        }
        return detectors[arrayNum].set(varName, v);
    }

    public boolean setDetectorParam(int arrayNum, String varName, String s) {
        if ((detectors == null) || (arrayNum < 0) || (arrayNum >= detectors.length)) {
            return false;
        }
        return detectors[arrayNum].set(varName, s);
    }

    //
    // BATCH FUNCTIONS
    //

    public boolean setBatchNum(int i) {

        if (i >= 0) {
            batchNum = i;
        } else {
            batchNum = -1;
        }

        init();

        return true;
    }

    public int numBatches() {
        if (batches == null) {
            return 0;
        }
        return batches.length;
    }

    public boolean checkBatchNum(int i) {
        return ((batches != null) && (i >= 0) && (i < batches.length));
    }

    public void initBatches() {

        String parameter, strValue, fname;

        if (batches != null) {

            for (Batch b : batches) {

                parameter = b.parameter;

                if (b.isString()) {
                    strValue = "";
                } else {
                    strValue = Float.toString((float) b.value);
                }

                strValue = strValue.replace('.', '_');
                strValue = strValue.replace(';', '_');
                strValue = strValue.replace(',', '_');

                b.init();

                if (!paramVectorVariableExists(parameter)) {
                    b.error("Project.initBatches", parameter, "does not exist");
                }

                fname = "Batch" + Integer.toString(b.batchNum) + "_" + parameter;

                if (strValue.length() > 0) {
                    fname += "_" + strValue;
                }

                fname = fname.replace('.', '_');

                if (b.execute) {
                    if (b.folder.length() == 0) {
                        b.folder = fname;
                    }
                }

            }

        }

    }

    public boolean nextBatch() {

        double saveVar;
        String saveStr, printStr;

        if (batches == null) {
            return false;
        }

        if (saveBatch != null) { // restore original parameter value if it exists

            if (saveBatch.isString()) {
                setObjectParam(saveBatch.parameter, saveBatch.strValue, false, false);
            } else {
                setObjectParam(saveBatch.parameter, saveBatch.value);
            }

        }

        for (int i = 0; i < batches.length; i += 1) {

            if (!batches[i].finished) {

                batchCounter = i;
                batchNum = batches[i].batchNum;
                subfolder = batches[i].folder;

                if (batches[i].isString()) {
                    saveStr = getObjectParamStr(batches[i].parameter);
                    saveBatch = new Batch(batchNum, batches[i].parameter, saveStr, batches[i].execute);
                    setObjectParam(batches[i].parameter, batches[i].strValue, false, false);
                    printStr = "New Batch: " + batches[i].parameter + ", " + batches[i].strValue;
                } else {
                    saveVar = getObjectParamVar(batches[i].parameter);
                    saveBatch = new Batch(batchNum, batches[i].parameter, saveVar, batches[i].execute);
                    setObjectParam(batches[i].parameter, batches[i].value);
                    printStr = "New Batch: " + batches[i].parameter + ", " + batches[i].value;
                }

                if (batches[i].execute) {
                    Master.log(printStr);
                    return true;
                }

            }

        }

        return false;

    }

    public void closeBatch() {
        if (checkBatchNum(batchCounter)) {
            batches[batchCounter].finished = true;
            subfolder = "";
        }
    }

    public int addBatch(Batch newBatch) {

        int i = 0;

        if (batches != null) {
            i = batches.length;
        }

        Batch[] newArray = new Batch[i+1];

        if (i > 0) {
            System.arraycopy(batches, 0, newArray, 0, i);
        }

        newArray[i] = newBatch;

        batches = newArray; // replace old array with new one

        //Master.log("added batch #" + i);

        return i;

    }
    
    public void killBatchesAll() {
        int i = batches.length;
        batches = null;
        if (i == 1) {
            Master.log("killed Batch #0");
        } else if (i > 1) {
            Master.log("killed Batches #0-" + i);
        }
    }

    public boolean killBatch(int i) {

        int k = 0;

        if (!checkBatchNum(i)) {
            return false;
        }

        if (batches.length == 1) {
            batches = null;
            Master.log("killed Batch #0");
            return true;
        }

        Batch[] b = new Batch[batches.length - 1]; // new array

        for (int j = 0; j < batches.length; j++) {
            if (j == i) {
                continue;
            }
            b[k] = batches[j];
            k++;
        }

        batches = b; // replace old array with new one

        Master.log("killed Batch #" + i);

        return true;

    }

    //
    // ERROR FUNCTIONS
    //

    public int addError(String errorStr) {
        return addError(new D3Derror(toString(), "", "", errorStr));
    }

    public int addError(D3Derror newError) {

        int i = 0;

        if (errors != null) {

            for (i = 0; i < errors.length; i++) {

                if (newError.object.equalsIgnoreCase(errors[i].object)) {
                    if (newError.function.equalsIgnoreCase(errors[i].function)) {
                        if (newError.parameter.equalsIgnoreCase(errors[i].parameter)) {
                            if (newError.error.equalsIgnoreCase(errors[i].error)) {
                                return -1; // already exists
                            }
                        }
                    }
                }

            }

            i = errors.length;

        }

        D3Derror[] newArray = new D3Derror[i+1];

        if (i > 0) {
            System.arraycopy(errors, 0, newArray, 0, i);
        }

        newArray[i] = newError;

        errors = newArray; // replace old array with new one

        return i;

    }

    public int numErrors() {

        if (errors == null) {
            return 0;
        }

        return errors.length;

    }
    
    public boolean checkErrorNum(int i) {
        return ((errors != null) && (i >= 0) && (i < errors.length));
    }
    
    public void killErrorsAll() {
        int i = errors.length;
        errors = null;
        if (i == 1) {
            Master.log("killed Error #0");
        } else if (i > 1) {
            Master.log("killed Errors #0-" + i);
        }
    }

    public boolean killError(int i) {
        int k = 0;

        if (!checkErrorNum(i)) {
            return false;
        }

        if (errors.length == 1) {
            errors = null;
            Master.log("killed Error #0");
            return true;
        }

        D3Derror[] e = new D3Derror[errors.length - 1]; // new array

        for (int j = 0; j < errors.length; j++) {
            if (j == i) {
                continue;
            }
            e[k] = errors[j];
            k++;
        }

        errors = e; // replace old array with new one

        Master.log("killed Error #" + i);

        return true;

    }

    //
    // Finite Difference and Monte Carlo Simulation Functions
    //

    public boolean simulationInit(boolean preview) {

        Master.writeToFiles = false;

        if (batches != null) {
            if (!nextBatch()) {
                return true; // finished batches
            }
        }

        if (monteCarlo != null) {
            if (monteCarlo.initSimulation()){
                error("ABORTED D3D init simulation");
                return true;
            }
        }

        init();

        if (numErrors() > 0) {
            error("ABORTED D3D init simulation");
            return true;
        }

        if (preview) {
            Master.log("preview D3D simulation, " + Master.currentDate());
        } else {
            Master.writeToFiles = true;
            Master.log("run D3D simulation, " + Master.currentDate());
        }

        if (monteCarlo != null) {
            monteCarlo.initSave();
        }

        geometry.checkSpace();
        geometry.checkPSFs();

        checkDiffusants();
        checkSources();
        checkDetectors();

        initDiffusantsPulseTimer();
        initDetectorsPulseTimer();
        initSourcesPulseTimer();

        openDiffusants();
        openSources();
        openDetectors();

        writeParamFile();

        return false;

    }

    public void simulationStart(boolean preview) {

        timer.start();

        init();

        if (numErrors() > 0) {
            Master.log("cannot start simulation due to registered error(s)");
            return;
        }

        if (monteCarlo != null) {

            monteCarlo.startSimulation(preview, true);

        } else {

            if ((finiteDifference == null) || finiteDifference.finished) {
                finiteDifference = new RunFiniteDifference(this);
            }

            finiteDifference.startSimulation(preview, true);

        }

    }

    public void simulationPause() {

        if (numErrors() > 0) {
            Master.log("cannot pause simulation due to registered error(s)");
            return;
        }

        if (monteCarlo != null) {
            monteCarlo.pauseSimulation();
        }

        if (finiteDifference != null) {
            finiteDifference.pauseSimulation();
        }

    }

    public void simulationFinish() {

        if (monteCarlo != null) {
            monteCarlo.finishSimulation();
        }

        closeDiffusants();
        closeSources();
        closeDetectors();
        closeBatch();

        timer.stop();

        Master.log(Master.currentDate());
        Master.log("total simulation session time = " + timer.toString());

    }

    public void simulationCancel() {

        if (numErrors() > 0) {
            Master.log("cannot cancel simulation due to registered error(s)");
            return;
        }

        if (monteCarlo != null) {
            monteCarlo.cancelSimulation();
        }

        if (finiteDifference != null) {
            finiteDifference.stopSimulation();
        }

    }

    //
    // Log File Functions
    //

    public String logFileName(){

        checkDirectory();

        String logFile = "D3Dlog";

        if (batchNum >= 0) {
            logFile += "_" + Integer.toString(batchNum);
        }

        logFile += ".dat";

        return logFile;

    }

    //
    // PARAMETER VECTOR FUNCTIONS
    //
    public void writeParamFile() {

        if (!Master.writeToFiles) {
            return;
        }

        String fileName = "D3Dconfig";
        
        if (batchNum >= 0) {
            fileName += "_" + Integer.toString(batchNum);
        }
        
        checkDirectory();

        geometry.updateVectors();
        diffusantsUpdateVectors();
        sourcesUpdateVectors();
        detectorsUpdateVectors();
        initProject.updateVectors();

        writeParamVectors(fileName);
        writeParamVectorsXML(fileName);

    }

    private boolean writeParamVectors(String fileName) {

        if (!Master.writeParamVectorOpen(fullDirectory(), fileName)) {
            return false;
        }
        
        Master.writeParamVector(this);
        Master.writeParamVector(geometry);
        Master.writeParamVectorArray(diffusants);
        Master.writeParamVectorArray(sources);
        Master.writeParamVectorArray(detectors);
        Master.writeParamVector(initProject);
        Master.writeParamVectorClose();

        return true;

    }

    private boolean writeParamVectorsXML(String fileName) {

        if (!Master.writeParamVectorXMLopen(fullDirectory(), fileName)) {
            return false;
        }
        
        Master.writeParamVectorXML(this);
        Master.writeParamVectorXML(geometry);
        Master.writeParamVectorArrayXML(diffusants);
        Master.writeParamVectorArrayXML(sources);
        Master.writeParamVectorArrayXML(detectors);
        Master.writeParamVectorXML(initProject);
        Master.writeParamVectorXMLclose();

        return true;

    }

    public boolean paramVectorVariableExists(String longName) {
        return (getObjectParamStr(longName) != null);
    }

    public double getObjectParamVar(String longName) {

        String parent = Utility.parentClass(longName);
        String child = Utility.childClass(longName);
        String parameter = Utility.parameterName(longName);
        
        //Master.log(parent + ", " + child + ", " + parameter);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.getVar(parameter);
        }

        if (parent.equalsIgnoreCase("project")) {
            return getVar(parameter);
        }

        if (parent.equalsIgnoreCase("geometry")) {
            if (geometry == null) {
                return Double.NaN;
            }
            return geometry.getVar(parameter);
        }

        if (parent.equalsIgnoreCase("finiteDifferences")) {
            if (finiteDifference == null) {
                return Double.NaN;
            }
            return finiteDifference.getVar(parameter);
        }

        if (parent.equalsIgnoreCase("monteCarlo")) {

            if (monteCarlo == null) {
                return Double.NaN;
            }

            if (child == null) {
                return monteCarlo.getVar(parameter);
            }

            return Double.NaN;

        }

        if (parent.equalsIgnoreCase("monteCarloActiveZone")) {

            if (monteCarlo == null) {
                return Double.NaN;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.getVar(parameter);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return Double.NaN;
                }
                return mcaz.activeZone.getVar(parameter);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return Double.NaN;
                }
                return mcaz.activeZone.getVar(parameter);
            }

            if (mcaz.save_NumDocked != null) {
                if (mcaz.save_NumDocked.name.equalsIgnoreCase(child)) {
                    return mcaz.save_NumDocked.getVar(parameter);
                }
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.getVar(parameter);
                }
            }

            return Double.NaN;

        }

        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return Double.NaN;
                        }
                        return d.save.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return Double.NaN;
                        }
                        return d.pulseTimer.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null) {
                            return Double.NaN;
                        }
                        return d.psf.getVar(parameter);
                    }
                    return Double.NaN;
                }
            }
        }

        if (sources != null) {
            for (Source s : sources) {
                if (s == null) {
                    continue;
                }
                if (s.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return s.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (s.save == null) {
                            return Double.NaN;
                        }
                        return s.save.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (s.pulseTimer == null) {
                            return Double.NaN;
                        }
                        return s.pulseTimer.getVar(parameter);
                    }
                    return Double.NaN;
                }
            }
        }

        if (detectors != null) {
            for (Detector d : detectors) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return Double.NaN;
                        }
                        return d.save.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return Double.NaN;
                        }
                        return d.pulseTimer.getVar(parameter);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null) {
                            return Double.NaN;
                        }
                        return d.psf.getVar(parameter);
                    }
                    return Double.NaN;
                }
            }
        }

        return Double.NaN;

    }

    public String getObjectParamStr(String longName) {
        int i;

        String parent = Utility.parentClass(longName);
        String child = Utility.childClass(longName);
        String parameter = Utility.parameterName(longName);

        //Master.log(parent + ", " + child + ", " + parameter);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.getStr(parameter);
        }

        if (parent.equalsIgnoreCase("project")) {
            return getStr(parameter);
        }

        if (parent.equalsIgnoreCase("geometry")) {
            if (geometry == null) {
                return null;
            }
            return geometry.getStr(parameter);
        }

        if (parent.equalsIgnoreCase("FiniteDifferences")) {
            if (finiteDifference == null) {
                return null;
            }
            return finiteDifference.getStr(parameter);
        }

        if (parent.equalsIgnoreCase("MonteCarlo")) {

            if (monteCarlo == null) {
                return null;
            }

            if (child == null) {
                return monteCarlo.getStr(parameter);
            }

            return null;

        }

        if (parent.equalsIgnoreCase("MonteCarloActiveZone")) {

            if (monteCarlo == null) {
                return null;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.getStr(parameter);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return null;
                }
                return mcaz.activeZone.getStr(parameter);
            }

            if (mcaz.save_NumDocked != null) {
                if (mcaz.save_NumDocked.name.equalsIgnoreCase(child)) {
                    return mcaz.save_NumDocked.getStr(parameter);
                }
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.getStr(parameter);
                }
            }

            return null;

        }

        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return null;
                        }
                        return d.save.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return null;
                        }
                        return d.pulseTimer.getStr(parameter);
                    }
                    if (child.startsWith("Pulse")) {
                        
                        if ((d.pulseTimer == null) || (d.pulseTimer.pulses == null)) {
                            return null;
                        }
                        
                        i = Pulse.pulseNum(child);

                        if ((i >= 0) && (i < d.pulseTimer.pulses.length)) {
                            return d.pulseTimer.pulses[i].getStr(parameter);
                        }
                        
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null) {
                            return null;
                        }
                        return d.psf.getStr(parameter);
                    }
                    return null;
                }
            }
        }

        if (sources != null) {
            for (Source s : sources) {
                if (s == null) {
                    continue;
                }
                if (s.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return s.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (s.save == null) {
                            return null;
                        }
                        return s.save.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (s.pulseTimer == null) {
                            return null;
                        }
                        return s.pulseTimer.getStr(parameter);
                    }
                    if (child.startsWith("Pulse")) {
                        
                        if ((s.pulseTimer == null) || (s.pulseTimer.pulses == null)) {
                            return null;
                        }
                        
                        i = Pulse.pulseNum(child);

                        if ((i >= 0) && (i < s.pulseTimer.pulses.length)) {
                            return s.pulseTimer.pulses[i].getStr(parameter);
                        }
                        
                    }
                    return null;
                }
            }
        }

        if (detectors != null) {
            for (Detector d : detectors) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return null;
                        }
                        return d.save.getStr(parameter);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return null;
                        }
                        return d.pulseTimer.getStr(parameter);
                    }
                    if (child.startsWith("Pulse")) {
                        
                        if ((d.pulseTimer == null) || (d.pulseTimer.pulses == null)) {
                            return null;
                        }
                        
                        i = Pulse.pulseNum(child);

                        if ((i >= 0) && (i < d.pulseTimer.pulses.length)) {
                            return d.pulseTimer.pulses[i].getStr(parameter);
                        }
                        
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null) {
                            return null;
                        }
                        return d.psf.getStr(parameter);
                    }
                    return null;
                }
            }
        }

        return null;

    }

    public boolean setObjectParam(String longName, double value) {
        int i;

        String parent = Utility.parentClass(longName);
        String child = Utility.childClass(longName);
        String parameter = Utility.parameterName(longName);
        
        //Master.log(parent + ", " + child + ", " + parameter);

        if (parent.equalsIgnoreCase("initProject")) {
            if (initProject == null) {
                return false;
            }
            return initProject.set(parameter, value);
        }

        if (parent.equalsIgnoreCase("project")) {
            return set(parameter, value);
        }

        if (parent.equalsIgnoreCase("geometry")) {
            if (geometry == null) {
                return false;
            }
            return geometry.set(parameter, value);
        }

        if (parent.equalsIgnoreCase("FiniteDifferences")) {
            if (finiteDifference == null) {
                return false;
            }
            return finiteDifference.set(parameter, value);
        }

        if (parent.equalsIgnoreCase("MonteCarlo")) {

            if (monteCarlo == null) {
                return false;
            }

            if (child == null) {
                return monteCarlo.set(parameter, value);
            }

            return false;

        }

        if (parent.equalsIgnoreCase("MonteCarloActiveZone")) {

            if (monteCarlo == null) {
                return false;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.set(parameter, value);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return false;
                }
                return mcaz.activeZone.set(parameter, value);
            }

            if (mcaz.save_NumDocked != null) {
                if (mcaz.save_NumDocked.name.equalsIgnoreCase(child)) {
                    return mcaz.save_NumDocked.set(parameter, value);
                }
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.set(parameter, value);
                }
            }

            return false;

        }

        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return false;
                        }
                        return d.save.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return false;
                        }
                        return d.pulseTimer.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null) {
                            return false;
                        }
                        return d.psf.set(parameter, value);
                    }
                    return false;
                }
            }
        }

        if (sources != null) {
            for (Source s : sources) {
                if (s == null) {
                    continue;
                }
                if (s.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return s.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (s.save == null) {
                            return false;
                        }
                        return s.save.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (s.pulseTimer == null) {
                            return false;
                        }
                        return s.pulseTimer.set(parameter, value);
                    }
                    if (child.startsWith("Pulse")) {
                        
                        if ((s.pulseTimer == null) || (s.pulseTimer.pulses == null)) {
                            return false;
                        }
                        
                        i = Pulse.pulseNum(child);

                        if ((i >= 0) && (i < s.pulseTimer.pulses.length)) {
                            return s.pulseTimer.pulses[i].set(parameter, value);
                        }
                        
                    }
                    return false;
                }
            }
        }

        if (detectors != null) {
            for (Detector d : detectors) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return false;
                        }
                        return d.save.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return false;
                        }
                        return d.pulseTimer.set(parameter, value);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null){
                            return false;
                        }
                        return d.psf.set(parameter, value);
                    }
                    return false;
                }
            }
        }

        return false;

    }
    
    public boolean setObjectParamProject(String longName, String strValue) {

        String parent = Utility.parentClass(longName);
        String child = Utility.childClass(longName);
        String parameter = Utility.parameterName(longName);

        //Master.log(parent + "," + child + "," + parameter + "," + strValue);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.set(parameter, strValue);
        }

        return false;

    }

    public boolean setObjectParamInitProject(String longName, String strValue) {

        String parent = Utility.parentClass(longName);
        String child = Utility.childClass(longName);
        String parameter = Utility.parameterName(longName);

        //Master.log(parent + "," + child + "," + parameter + "," + strValue);

        if (parent.equalsIgnoreCase("initProject")) {
            return initProject.set(parameter, strValue);
        }

        return false;

    }

    public boolean setObjectParam(String longName, String strValue, boolean onlyProject, boolean onlyInitProject) {

        boolean set;

        String parent = Utility.parentClass(longName);
        String child = Utility.childClass(longName);
        String parameter = Utility.parameterName(longName);
        
        //Master.log(parent + ", " + child + ", " + parameter);
        
        if (onlyProject) {

            if (parent.equalsIgnoreCase("project")) {
                return set(parameter, strValue);
            }

            return false; // nothing

        }

        if (onlyInitProject) {

            if (parent.equalsIgnoreCase("initProject")) {
                return initProject.set(parameter, strValue);
            }

            return false; // nothing

        }

        if (parent.equalsIgnoreCase("geometry")) {

            if (geometry == null) {
                return false; // error
            }

            return geometry.set(parameter, strValue);

        }

        if (parent.equalsIgnoreCase("finitedifferences")) {
            if (finiteDifference == null) {
                return false;
            }
            return finiteDifference.set(parameter, strValue);
        }

        if (parent.equalsIgnoreCase("MonteCarlo")) {

            if (monteCarlo == null) {
                return false;
            }

            if (child == null) {
                return monteCarlo.set(parameter, strValue);
            }

            return false;

        }

        if (parent.equalsIgnoreCase("MonteCarloActiveZone")) {

            if (monteCarlo == null) {
                return false;
            }

            RunMonteCarloAZ mcaz = (RunMonteCarloAZ) monteCarlo;

            if (child == null) {
                return mcaz.set(parameter, strValue);
            }

            if (child.equalsIgnoreCase("ActiveZone")) {
                if (mcaz.activeZone == null) {
                    return false;
                }
                return mcaz.activeZone.set(parameter, strValue);
            }

            if (mcaz.save_NumDocked != null) {
                if (mcaz.save_NumDocked.name.equalsIgnoreCase(child)) {
                    return mcaz.save_NumDocked.set(parameter, strValue);
                }
            }

            if (mcaz.save_Release != null) {
                if (mcaz.save_Release.name.equalsIgnoreCase(child)) {
                    return mcaz.save_Release.set(parameter, strValue);
                }
            }

            return false;

        }

        if (diffusants != null) {
            for (Diffusant d : diffusants) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return false;
                        }
                        return d.save.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return false;
                        }
                        return d.pulseTimer.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null) {
                            return false;
                        }
                        return d.psf.set(parameter, strValue);
                    }
                    return false;
                }
            }
        }

        if (sources != null) {
            for (Source s : sources) {
                if (s == null) {
                    continue;
                }
                if (s.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return s.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (s.save == null) {
                            return false;
                        }
                        return s.save.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (s.pulseTimer == null) {
                            return false;
                        }
                        return s.pulseTimer.set(parameter, strValue);
                    }
                    return false;
                }
            }
        }

        if (detectors != null) {
            for (Detector d : detectors) {
                if (d == null) {
                    continue;
                }
                if (d.name.equalsIgnoreCase(parent)) {
                    if (child == null) {
                        return d.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("Save")) {
                        if (d.save == null) {
                            return false;
                        }
                        return d.save.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PulseTimer")) {
                        if (d.pulseTimer == null) {
                            return false;
                        }
                        return d.pulseTimer.set(parameter, strValue);
                    }
                    if (child.equalsIgnoreCase("PSF")) {
                        if (d.psf == null) {
                            return false;
                        }
                        return d.psf.set(parameter, strValue);
                    }
                    return false;
                }
            }
        }

        return false;

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (finiteDifference != null){
            finiteDifference.addUser(pv);
        }

        if (monteCarlo != null){
            monteCarlo.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (finiteDifference != null){
            addBlankParam();
            finiteDifference.createVector(true);
            addVector(finiteDifference.getVector());
            finiteDifference.addUser(this);
        }

        if (monteCarlo != null){
            addBlankParam();
            monteCarlo.createVector(true);
            addVector(monteCarlo.getVector());
            monteCarlo.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (finiteDifference != null){
            finiteDifference.updateVector(v);
        }

        if (monteCarlo != null){
            monteCarlo.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {
        
        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof Project)) {
            return false;
        }
        
        String n = o.getName();

        if (n.equalsIgnoreCase("dx")) {
            if (v <= 0) {
                return false;
            }
            dx = v;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("dt")) {
            if (v <= 0) {
                return false;
            }
            init();
            return true;
        }
        if (n.equalsIgnoreCase("simTime")) {
            if (v <= 0) {
                return false;
            }
            simTime = v;
            return true;
        }
        if (n.equalsIgnoreCase("simTemp")) {
            if (v < 0) {
                return false;
            }
            simTemp = v;
            return true;
        }
        if (n.equalsIgnoreCase("maxT")) {
            if (v <= 0) {
                return false;
            }
            simTime = v;
            return true;
        }
        if (n.equalsIgnoreCase("saveRate")) {
            if (v <= 0) {
                return false;
            }
            setOutputRate(v);
            return true;
        }
        if (n.equalsIgnoreCase("saveDT")) {
            if (v <= 0) {
                return false;
            }
            setOutputRate(1/v);
            return true;
        }
        if (n.equalsIgnoreCase("printRate")) {
            if (v < 0) {
                return false;
            }
            setPrintRate(v);
            return true;
        }
        if (n.equalsIgnoreCase("printDT")) {
            if (v < 0) {
                return false;
            }
            setPrintRate(1/v);
            return true;
        }
        if (n.equalsIgnoreCase("batchNum")) {
            return setBatchNum((int) v);
        }
        if (n.equalsIgnoreCase("EM_iseries")) {
            if (v < 0) {
                return false;
            }
            EM_iseries = (int) v;
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

        if (!(o.paramVector instanceof Project)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            Master.updateMainFrameTitle();
            return true;
        }
        if (n.equalsIgnoreCase("directory")) {
            directory = s;
            init(); // reset output file names
            return true;
        }
        if (n.equalsIgnoreCase("folder")) {
            folder = s;
            init(); // reset output file names
            return true;
        }
        if (n.equalsIgnoreCase("timeUnits")) {
            timeUnits = s;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("spaceUnits")) {
            spaceUnits = s;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("concUnits")) {
            concUnits = s;
            init();
            return true;
        }
        if (n.equalsIgnoreCase("date")) {
            if (s.length() == 0) {
                setDate();
            } else {
                date = s;
            }
            return true;
        }
        if (n.equalsIgnoreCase("seedMT")) {
            try {
                if (s.length() > 0) {
                    Long.parseLong(s); // throws NumberFormatException
                    seedMT = s;
                    Master.initMersenneTwister();
                } else {
                    seedMT = "";
                }
            } catch (NumberFormatException e) {
                // do nothing
            }
            return true;
        }
        return super.setMyParams(o, s);
    }
}
