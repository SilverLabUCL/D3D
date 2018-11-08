package ucl.silver.d3d.core;

import java.io.*;

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
public class Save extends ParamVector {

    public String outputFile = "";
    private String file = "";

    public int dataPoints = 1; // number of data elements to save (e.g. number of voxels)

    public double sampleRate = 100; // binaryWriter rate (kHz)
    public double sampleInterval = 1 / sampleRate; // binaryWriter time step (ms)
    public int samples2save, skipSamples; // time samples

    public boolean save2TextFile = false; // save data to external Unicode text file
    public boolean save2BinaryFile = false; // save data to external binary file (binary Double)
    public boolean saveWhileComputing = true; // save data while computing, otherwise once at the end of simulation
    public boolean saveWhileComputingAppend = false; // save while computing, closing file after each save
    public boolean save2array = false; // save data to an array

    private boolean save2file = false; // save2TextFile || save2BinaryFile
    private boolean saveWhileComputing2 = false; // quick boolean flag used during simulations

    public String xdim = ""; // x dimensions
    public String ydim = ""; // y dimensions
    public boolean autoDimensions = true; // automatically set dimension names

    public int skipCounter = 0, dataCounter = 0;

    public double[][] data = null; // saved results

    private BufferedWriter textWriter = null;
    private DataOutputStream binaryWriter = null;

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("sampleRate")) {
            return project.freqUnits;
        }
        if (name.equalsIgnoreCase("sampleInterval")) {
            return project.timeUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("samples2save")) {
            return false;
        }
        if (name.equalsIgnoreCase("dataPoints")) {
            return false;
        }
        if (name.equalsIgnoreCase("save2array")) {
            return false;
        }
        if (outputFile.length() == 0) {
            return false;
        }
        if (sampleRate < 0) {
            return false;
        }
        return true;
    }

    public Save(Project p) {
        super(p);
    }

    public void init() {
        setOutputRate(project.saveRate);
    }

    public void setOutputRate(double newRate) {

        double maxRate = 0;
        double projectDT, simTime;

        skipSamples = 0;
        sampleInterval = 0;
        sampleRate = 0;
        samples2save = 0;

        if (project != null) {

            projectDT = project.dt;
            simTime = project.simTime;

            if (projectDT > 0) {
                maxRate = 1 / projectDT;
            }

            if (newRate > maxRate) {
                sampleInterval = projectDT;
                skipSamples = 0;
                sampleRate = maxRate;
            } else {
                skipSamples = (int) (maxRate / newRate) - 1;
                sampleInterval = projectDT * (skipSamples + 1);
                sampleRate = 1 / sampleInterval;
            }

            samples2save = (int) (simTime / sampleInterval) + 1;

        }

        updateVectors();

    }

    public void setOutputRate(int newSkipSamples) {

        int maxSamples;
        double projectDT, simTime;

        skipSamples = 0;
        sampleInterval = 0;
        sampleRate = 0;
        samples2save = 0;

        if (project != null) {

            projectDT = project.dt;
            simTime = project.simTime;

            if ((projectDT > 0) && (simTime > 0)) {

                maxSamples = (int) (((simTime / projectDT) + 1.0) / 2.0);

                if (newSkipSamples >= maxSamples) {
                    skipSamples = maxSamples;
                } else {
                    skipSamples = newSkipSamples;
                }

            }

            sampleInterval = projectDT * (skipSamples + 1);
            sampleRate = 1 / sampleInterval;
            samples2save = (int) (simTime / sampleInterval) + 1;

        }

        updateVectors();

    }

    public void fileName(String fxn, String diffusantName) {

        if (fxn.length() == 0) {
            return;
        }

        String batch = "";

        if (project.batchNum >= 0) {
            batch = "_" + Integer.toString(project.batchNum);
        }

        if (diffusantName.length() > 0) {
            outputFile = fxn + "_" + diffusantName + batch;
        } else {
            outputFile = fxn + batch;
        }

    }

    public String fullPathFileName(int type) {

        String dir = project.fullDirectory();

        switch (type) {
            case 0: // Unicode text
                return dir + outputFile + ".dat";
            case 1: // binary double
                return dir + outputFile + ".bin";
        }

        return ".d3d"; // unknown

    }

    public boolean init(String name, CoordinatesVoxels c, double time, int dPoints) {

        skipCounter = 0;
        dataCounter = 0;

        if (dPoints < 1) {
            dPoints = 1;
        }

        dataPoints = dPoints;

        if (Master.writeToFiles) {
            save2file = save2TextFile || save2BinaryFile;
        } else {
            save2file = false;
        }

        if (save2file && !saveWhileComputing) {
            save2array = true; // must save to array
        }

        if (save2array && (dataPoints > 0) && (samples2save > 0)) {
            data = new double[dataPoints][samples2save];
        }

        saveWhileComputing2 = (save2file && saveWhileComputing);

        if (save2file && saveWhileComputing && (outputFile.length() > 0)) {
            return openFile(name, c, time);
        }

        return true;

    }

    private boolean openFile(String name, CoordinatesVoxels c, double time) {

        boolean text = false, binary = false;

        if (save2TextFile) {
            text = openFileText(name, c, time);
        }

        if (save2BinaryFile) {
            binary = openFileBinary(name, c, time);
        }

        if (saveWhileComputingAppend) {
            closeFile();
        }

        return text || binary;

    }

    private boolean openFileText(String name, CoordinatesVoxels c, double time) {

        file = fullPathFileName(0);

        try {
            textWriter = new BufferedWriter(new FileWriter(file)); // open new file
            writeHeader(textWriter, name, c, time, 0);
            Master.log("created text output file " + file);
            return true;
        } catch (IOException e) {
            error("failed to create text output file " + file);
            return false;
        }

    }

    private boolean openFileText() { // for appending

        try {
            textWriter = new BufferedWriter(new FileWriter(file, true));
            return true;
        } catch (IOException e) {
            error("failed to open text output file " + file);
            return false;
        }

    }

    private boolean openFileBinary(String name, CoordinatesVoxels c, double time) {

        file = fullPathFileName(1);

        BufferedWriter bw;

        try {

            bw = new BufferedWriter(new FileWriter(file));
            writeHeader(bw, name, c, time, 1); // first write text header
            bw.close();
            Master.log("created binary output file " + file);

            binaryWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true))); // open and append

            return true;

        } catch (IOException e) {
            error("failed to create binary output file " + file);
            return false;
        }

    }

    private boolean openFileBinary() { // for appending

        try {
            binaryWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true))); // open and append
            return true;
        } catch (IOException e) {
            error("failed to open binary output file " + file);
            return false;
        }

    }

    private boolean writeHeader(BufferedWriter bw, String name, CoordinatesVoxels c, double time, int dataType) {

        double dx = project.dx;

        String typeStr = "Unknown";

        if (bw == null) {
            return false;
        }

        switch (dataType) {
            case 0:
                typeStr = "Unicode Text";
                break;
            case 1:
                typeStr = "Binary Double";
                break;
        }

        try {

            bw.write("D3D=" + Master.D3D_VERSION);
            bw.newLine();

            if (name.length() > 0) {
                bw.write("name=" + name);
                bw.newLine();
            }

            if (c != null) {
                bw.write("x1=" + Integer.toString(c.xVoxel1));
                bw.newLine();
                bw.write("y1=" + Integer.toString(c.yVoxel1));
                bw.newLine();
                bw.write("z1=" + Integer.toString(c.zVoxel1));
                bw.newLine();
                bw.write("x2=" + Integer.toString(c.xVoxel2));
                bw.newLine();
                bw.write("y2=" + Integer.toString(c.yVoxel2));
                bw.newLine();
                bw.write("z2=" + Integer.toString(c.zVoxel2));
                bw.newLine();
            }

            bw.write("pointsPerSample=" + dataPoints);
            bw.newLine();

            bw.write("samples2save=" + samples2save);
            bw.newLine();

            if (dx > 0) {
                bw.write("dx=" + Double.toString(dx));
                bw.newLine();
            }

            if (sampleInterval > 0) {
                bw.write("dt=" + Double.toString(sampleInterval));
                bw.newLine();
            }

            if (time >= 0) {
                bw.write("t=" + Double.toString(time));
                bw.newLine();
            }

            bw.write("xdim=" + xdim);
            bw.newLine();

            bw.write("ydim=" + ydim);
            bw.newLine();

            bw.write("format=" + typeStr);
            bw.newLine();

            return true;

        } catch (IOException e) {
            error("unable to write to text output file " + bw.toString());
            return false;
        }

    }

    public boolean saveData(double value) {

        boolean ok = false;

        if (skipCounter == 0) {

            if (saveWhileComputing2) {
                ok = writeData(value);
            }

            if (save2array) {

                data[0][dataCounter] = value;
                dataCounter++;

                if (dataCounter >= data[0].length) {
                    save2array = false;
                    saveWhileComputing2 = false;
                }

            }

        }

        skipCounter++;

        if (skipCounter > skipSamples) {
            skipCounter = 0;
        }

        return ok;

    }

    public boolean saveData(double[] values) {

        boolean ok = false;

        if (skipCounter == 0) {

            if (saveWhileComputing2) {
                ok = writeData(values);
            }

            if (save2array && (values.length == data.length)) {

                for (int i = 0; i < values.length; i++) {
                    data[i][dataCounter] = values[i];
                }

                dataCounter++;

                if (dataCounter >= data[0].length) {
                    save2array = false;
                    saveWhileComputing2 = false;
                }

            }

        }

        skipCounter++;

        if (skipCounter > skipSamples) {
            skipCounter = 0;
        }

        return ok;

    }

    public boolean saveData(double[][] values, int dataPoint, int[] indexes) {

        int numValues, index;
        boolean ok = false;

        if (skipCounter == 0) {

            if (saveWhileComputing2) {
                ok = writeData(values, dataPoint, indexes);
            }

            if (save2array) {

                numValues = indexes.length;

                for (int i = 0; i < numValues; i++) {
                    index = indexes[i];
                    data[i][dataCounter] = values[dataPoint][index];
                }

                dataCounter++;

                if (dataCounter >= data[0].length) {
                    save2array = false;
                    saveWhileComputing2 = false;
                }

            }

        }

        skipCounter++;

        if (skipCounter > skipSamples) {
            skipCounter = 0;
        }

        return ok;

    }

    public boolean writeDataArray() {

        if (data == null) {
            return false;
        }

        for (int i = 0; i < dataCounter; i++) {
            for (int j = 0; j < data.length; j++) {
                writeData(data[j][i]);
            }
        }

        return true;

    }

    public boolean writeData(double value) {

        boolean text = false, binary = false;

        if (save2TextFile) {
            text = writeDataText(value);
        }

        if (save2BinaryFile) {
            binary = writeDataBinary(value);
        }

        return text || binary;

    }

    public boolean writeData(double[] values) {

        boolean text = false, binary = false;

        if (save2TextFile) {
            text = writeDataText(values);
        }

        if (save2BinaryFile) {
            binary = writeDataBinary(values);
        }

        return text || binary;

    }

    public boolean writeData(double[][] values, int dataPoint, int[] indexes) {

        boolean text = false, binary = false;

        if (save2TextFile) {
            text = writeDataText(values, dataPoint, indexes);
        }

        if (save2BinaryFile) {
            binary = writeDataBinary(values, dataPoint, indexes);
        }

        return text || binary;

    }

    private boolean writeDataText(double value) {

        if (saveWhileComputingAppend) {
            try {
                textWriter = new BufferedWriter(new FileWriter(file, true));
                textWriter.write(Double.toString(value));
                textWriter.newLine();
                textWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to text output file " + file);
                return false;
            }
        }

        try {
            textWriter.write(Double.toString(value));
            textWriter.newLine();
            return true;
        } catch (IOException e) {
            error("failed to write to text output file " + file);
            return false;
        }

    }

    private boolean writeDataText(double[] values) {

        if (saveWhileComputingAppend) {
            try {
                textWriter = new BufferedWriter(new FileWriter(file, true));
                for (int i = 0; i < values.length; i++) {
                    textWriter.write(Double.toString(values[i]));
                    textWriter.newLine();
                }
                textWriter.newLine();
                textWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to text output file " + file);
                return false;
            }
        }

        try {

            for (int i = 0; i < values.length; i++) {
                textWriter.write(Double.toString(values[i]));
                textWriter.newLine();
            }

            return true;

        } catch (IOException e) {
            error("failed to write to text output file " + file);
            return false;
        }

    }

    private boolean writeDataText(double[][] values, int dataPoint, int[] indexes) {

        int index;

        int numVoxels = indexes.length;

        if (saveWhileComputingAppend) {
            try {
                textWriter = new BufferedWriter(new FileWriter(file, true));
                for (int i = 0; i < numVoxels; i++) {
                    index = indexes[i];
                    textWriter.write(Double.toString(values[dataPoint][index]));
                    textWriter.newLine();
                }
                textWriter.newLine();
                textWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to text output file " + file);
                return false;
            }
        }

        try {

            for (int i = 0; i < numVoxels; i++) {
                index = indexes[i];
                textWriter.write(Double.toString(values[dataPoint][index]));
                textWriter.newLine();
            }

            return true;

        } catch (IOException e) {
            error("failed to write to text output file " + file);
            return false;
        }

    }

    public boolean writeString(String s) {

        boolean text = false, binary = false;

        if (save2TextFile) {
            text = writeStringText(s);
        }

        if (save2BinaryFile) {
            binary = writeStringBinary(s);
        }

        return text || binary;

    }

    private boolean writeStringText(String s) {

        if (saveWhileComputingAppend) {
            try {
                textWriter = new BufferedWriter(new FileWriter(file, true));
                textWriter.write(s);
                textWriter.newLine();
                textWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to text output file " + file);
                return false;
            }
        }

        try {
            textWriter.write(s);
            textWriter.newLine();
            return true;
        } catch (IOException e) {
            error("unable to write to text output file " + fullPathFileName(0));
            return false;
        }

    }

    private boolean writeStringBinary(String s) {

        if (saveWhileComputingAppend) {
            try {
                binaryWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true)));
                binaryWriter.writeUTF(s);
                binaryWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to binary output file " + file);
                return false;
            }
        }

        try {
            binaryWriter.writeUTF(s);
            return true;
        } catch (IOException e) {
            error("unable to write to binary output file " + fullPathFileName(1));
            return false;
        }

    }

    private boolean writeDataBinary(double value) {

        if (saveWhileComputingAppend) {
            try {
                binaryWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true)));
                binaryWriter.writeDouble(value);
                binaryWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to binary output file " + file);
                return false;
            }
        }

        try {
            binaryWriter.writeDouble(value);
            return true;
        } catch (IOException e) {
            error("failed to write to binary output file " + file);
            return false;
        }

    }

    private boolean writeDataBinary(double[] values) {

        if (saveWhileComputingAppend) {
            try {
                binaryWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true)));
                for (int i = 0; i < values.length; i++) {
                    binaryWriter.writeDouble(values[i]);
                }
                binaryWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to binary output file " + file);
                return false;
            }
        }

        try {

            for (int i = 0; i < values.length; i++) {
                binaryWriter.writeDouble(values[i]);
            }

            return true;

        } catch (IOException e) {
            error("failed to write to binary output file " + file);
            return false;
        }

    }

    private boolean writeDataBinary(double[][] values, int dataPoint, int[] indexes) {

        int index;

        if (saveWhileComputingAppend) {
            try {
                binaryWriter = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file, true)));
                for (int i = 0; i < indexes.length; i++) {
                    index = indexes[i];
                    binaryWriter.writeDouble(values[dataPoint][index]);
                }
                binaryWriter.close();
                return true;
            } catch (IOException e) {
                error("failed to write to binary output file " + file);
                return false;
            }
        }

        try {

            for (int i = 0; i < indexes.length; i++) {
                index = indexes[i];
                binaryWriter.writeDouble(values[dataPoint][index]);
            }

            return true;

        } catch (IOException e) {
            error("failed to write to binary output file " + file);
            return false;
        }

    }

    public boolean finish(String name, CoordinatesVoxels c, double time) {

        if (!save2file) {
            return false; // nothing to do
        }

        if (saveWhileComputing) {

            return closeFile();

        } else {

            if (!openFile(name, c, time)) {
                return false;
            }

            writeDataArray();

            return closeFile();

        }

    }

    private boolean closeFile() {

        boolean text = false, binary = false;

        if (save2TextFile) {
            text = closeFileText();
        }

        if (save2BinaryFile) {
            binary = closeFileBinary();
        }

        return text || binary;

    }

    private boolean closeFileText() {

        String file = fullPathFileName(0);

        if (textWriter == null) {
            return false;
        }

        try {
            textWriter.close();
            Master.log("closed text output file " + file);
            textWriter = null;
            return true;
        } catch (IOException e) {
            error("unable to close text output file " + file);
            return false;
        }

    }

    private boolean closeFileBinary() {

        String file = fullPathFileName(1);

        if (binaryWriter == null) {
            return false;
        }

        try {
            binaryWriter.close();
            Master.log("closed binary output file " + file);
            binaryWriter = null;
            return true;
        } catch (IOException e) {
            error("unable to close binary output file " + file);
            return false;
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof Save)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("sampleRate")) {
            if (v <= 0) {
                return false;
            }

            setOutputRate(v);
            return true;
        }
        if (n.equalsIgnoreCase("sampleInterval")) {
            if (v <= 0) {
                return false;
            }
            setOutputRate(1 / v);
            return true;
        }
        if (n.equalsIgnoreCase("skipSamples")) {
            if (v < 0) {
                return false;
            }
            setOutputRate((int) v);
            return true;
        }
        if (n.equalsIgnoreCase("samples2save")) {
            if (v < 0) {
                return false;
            }
            samples2save = (int) v;
            return true;
        }
        if (n.equalsIgnoreCase("save2TextFile")) {
            save2TextFile = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("save2BinaryFile")) {
            save2BinaryFile = (v == 1);
            return true;
        }
        if (n.equalsIgnoreCase("saveWhileComputing")) {
            if (v == 1) {
                saveWhileComputing = true;
            } else {
                saveWhileComputing = false;
                save2array = true;
            }
            return true;
        }
        if (n.equalsIgnoreCase("autoDimensions")) {
            autoDimensions = (v == 1);
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

        if (!(o.paramVector instanceof Save)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        if (n.equalsIgnoreCase("outputFile")) {
            outputFile = s;
            return true;
        }
        if (n.equalsIgnoreCase("xdim")) {
            ydim = s;
            autoDimensions = false;
            return true;
        }
        if (n.equalsIgnoreCase("ydim")) {
            ydim = s;
            autoDimensions = false;
            return true;
        }
        return super.setMyParams(o, s);
    }
}
