package ucl.silver.d3d.init;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import ucl.silver.d3d.core.*;
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
 * @version 2.1.projection.demo
 */
public class InitMC_Projection_Demo extends InitProject {

    public double geometry_xywidth = 1.0; // um
    public double geometry_zwidth = -1; // auto computeProjection
    public boolean geometry_PBC = true;

    public double particle_D_mean = 46 / 1000.0; // ET10 // um
    public double particle_D_stdv = 4 / 1000.0;

    public double set_particle_VF = 0.40;
    public int set_num_particles = -1;
    public double avg_num_particles = 0; // average

    public int num_sections = 1;
    public int num_sections_extra = 4; // for creating geometry
    public double section_thickness = 0; // um

    public String phiList = ""; // list of mean cap angles // degrees
    public double phiCV = 0; // phi CV // degrees

    public double save_diam_scalefactor = 1000; // convert um to nm
    public boolean save_histo_xbins = true;

    public boolean save_to_disk = false;
    //public boolean save_to_disk = true;

    public String MT_seed = "";
    //public String MT_seed = "1603737488075"; // for simulations in manuscript

    private String title = "Projection Demo";

    DiffusantParticlesProjectionKeiding projection;

    private final NumberFormat numformat = new DecimalFormat("0.######");

    public InitMC_Projection_Demo(Project p) {

        super(p);

        projection = new DiffusantParticlesProjectionKeiding(p, "", null);

        projection.print_space_units = "nm";
        directory = "/Users/jason/Documents/FrapVGlut1_hMFTs/Manuscript Stereology/D3D/Demo/";
        projection.directory = directory;

        String[] flist = {"initPrompt", "initFig4", "initFig5", "initETzstack", "initDisector"};
        initFuncList = flist;

        createVector(true);

    }

    @Override
    public boolean initFunction(String initSelect) {

        int i = initFunctionNum(initSelect);

        if (!Master.foundMainStartUpArguments) {
            //i = 3;
        }

        switch (i) {
            case 0:
                return init_prompt();
            case 1:
                //return initFig4_call();
                return initFig4_prompt();
            case 2:
                //return initFig5_all();
                return initFig5_prompt();
            case 3:
                //return init_ET_zstack_call();
                return init_ET_zstack_prompt();
            case 4:
                //return initDisectorCall();
                return initDisector_prompt();
            default:
                error("initFunction", "initSelect", "failed to find init function");
                return true; // error
        }

    }
    
    private boolean init_prompt() {
        
        String Tud_List[] = {"Fig4", "Fig5", "ET z-stack", "Disector"};
        String pstr = Master.promptForInput(Tud_List, "choose projection simulation", title, "Fig4");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }
        
        if (pstr.equalsIgnoreCase("Fig4")) {
            return initFig4_prompt();
        } else if (pstr.equalsIgnoreCase("Fig5")) {
            return initFig5_prompt();
        } if (pstr.equalsIgnoreCase("ET z-stack")) {
            return init_ET_zstack_prompt();
        } if (pstr.equalsIgnoreCase("Disector")) {
            return initDisector_prompt();
        }
        
        return true;
        
    }

    private boolean initFig4_prompt() {

        double D_SD;
        String pstr;

        title += " Fig 4";

        String Tud_List[] = {"0", "1"};
        pstr = Master.promptForInput(Tud_List, "choose section thickness (unit diameters)", title, "0");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        double Tud = Double.parseDouble(pstr);

        String N_List[] = {"200", "300", "500", "2000"};
        pstr = Master.promptForInput(N_List, "choose number of particles per projection", title, "500");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        int numParticles = Integer.parseInt(pstr);

        String DCV_List[] = {"0.04", "0.09", "0.13", "0.17"};
        pstr = Master.promptForInput(DCV_List, "choose CV of 3D diameters", title, "0.09");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        if (pstr.equals("0.04")) {
            D_SD = 2.0 / 1000.0;
        } else if (pstr.equals("0.09")) {
            D_SD = 4.0 / 1000.0;
        } else if (pstr.equals("0.13")) {
            D_SD = 6.0 / 1000.0;
        } else if (pstr.equals("0.17")) {
            D_SD = 8.0 / 1000.0;
        } else {
            return true;
        }

        String Phi_List[] = {"10", "20", "30", "40", "45", "50", "55", "60", "65", "70", "75", "80"};
        pstr = Master.promptForInput(Phi_List, "choose ϕ°", title, "20");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        int phi[] = new int[1];
        phi[0] = Integer.parseInt(pstr);

        String PhiCV_List[] = {"0", "0.2"};
        pstr = Master.promptForInput(PhiCV_List, "choose CV of ϕ", title, "0");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        double phi_CV = Double.parseDouble(pstr);

        if (setDirectoryPrompt(title)) {
            return true;
        }

        return initFig4(Tud, numParticles, D_SD, phi, phi_CV);

    }

    private boolean initFig4_call() {

        double Tud = 0; // ud
        //double Tud = 1; // ud

        //int num_particles = 200;
        //int num_particles = 300;
        int num_particles = 500;
        //int num_particles = 2000;

        //double D_SD = 2.0 / 1000.0; // um
        double D_SD = 4.0 / 1000.0; // um
        //double D_SD = 6.0 / 1000.0; // um
        //double D_SD = 8.0 / 1000.0; // um

        //int[] phi = phiList(10, 10, 4);
        //int[] phi = phiList(45, 5, 8);
        int[] phi = phiList(20, 10, 1);

        double phi_CV = 0;
        //double phi_CV = 0.2;

        return initFig4(Tud, num_particles, D_SD, phi, phi_CV);

    }

    private boolean initFig4(double Tud, int numParticles, double D_SD, int[] phi, double phi_CV) {

        String filename;
        String Tstr, Nstr, DCVstr, phiStr, phiCVstr, IZstr;

        boolean exact_num_particles = false;
        //boolean exact_num_particles = true; // do not use if one needs to compute density // use for phi-cutoff analysis

        //num_sections = 10; // testing
        num_sections = 100;
        avg_num_particles = numParticles;

        if (exact_num_particles) {
            projection.num_particles = numParticles;
            set_num_particles = numParticles;
        }

        phiCV = phi_CV;
        phiCVstr = "CV" + Integer.toString((int) (phi_CV * 10));

        projection.name = "Fig4 Projection";

        particle_D_stdv = D_SD;

        int iDSD = (int) (particle_D_stdv * 1000);
        //String DSDstr = "_DSD" + Integer.toString(iDSD);
        //String VFstr = "_VF" + Integer.toString((int)(set_particle_VF * 100));

        if (D_SD < 5) {
            DCVstr = "_DCV0" + Integer.toString((int) Math.round(100 * particle_D_stdv / particle_D_mean));
        } else {
            DCVstr = "_DCV" + Integer.toString((int) Math.round(100 * particle_D_stdv / particle_D_mean));
        }

        projection.directory = directory;
        Master.log("directory = " + projection.directory);

        Tstr = "_T" + Integer.toString((int) Tud);
        Nstr = "_N" + Integer.toString(numParticles);
        //IZstr = IZmax_str(projection.Iz_max);

        section_thickness = particle_D_mean * Tud;

        projection.compute_Gd = true;

        switch (iDSD) {
            case 2:
                projection.histo_d_max = 65 / 1000.0;
                //projection.histo_d_max = 70 / 1000.0;
                projection.histo_d_bwidth = 1 / 1000.0;
                break;
            case 4:
                projection.histo_d_max = 70 / 1000.0;
                projection.histo_d_bwidth = 2 / 1000.0;
                //if (numParticles == 2000) {
                //projection.histo_d_bwidth = 1 / 1000.0;
                //}
                break;
            case 6:
                projection.histo_d_max = 75 / 1000.0;
                projection.histo_d_bwidth = 3 / 1000.0;
                //if (numParticles == 2000) {
                //projection.histo_d_bwidth = 2 / 1000.0;
                //}
                break;
            case 8:
                projection.histo_d_max = 80 / 1000.0;
                projection.histo_d_bwidth = 4 / 1000.0;
                //if (numParticles == 2000) {
                //projection.histo_d_bwidth = 2 / 1000.0;
                //}
                break;
        }

        //save_histo_xbins = true;
        Master.log("G(d) bin width = " + projection.histo_d_bwidth);

        projection.sections_stats_save_phi = true;

        phiList = "";

        for (int i = 0; i < phi.length; i++) {

            projection.phi_true = phi[i];
            projection.phi_CV = phi_CV;
            phiStr = Integer.toString(phi[i]);
            phiList += phiStr + ";";

            geometry_xywidth = get_geometry_xywidth(numParticles, phi[i]);

            if (exact_num_particles) {
                geometry_xywidth *= 1.2; // increase geometry size
            }

            Master.log("geometry xy width = " + geometry_xywidth + " " + project.spaceUnits);

            if (initProjection()) {
                return true;
            }

            Master.log("Phi = " + Integer.toString(phi[i]) + "°, Imax = " + Double.toString(projection.Iz_max));
            //Master.log("3D density = " + projection.dvs.density);
            //Master.log("geometry xy-area = " + (geometry.xWidth * geometry.yWidth));

            //projection.roi_xy_width = geometry.xWidth * 0.9; // TEST ROI // gives similar 2D density
            projection.computeProjection();
            projection.print_section_stats_avg();

            if (!save_to_disk) {
                continue;
            }

            //filename = "D3D_Fig4" + Tstr + Nstr + DSDstr + "_P" + phiStr + phiCVstr;
            filename = "D3D_Fig4" + Tstr + Nstr + DCVstr + "_P" + phiStr + phiCVstr;

            //projection.save_Fd_total(filename + "_F" , save_diam_scalefactor);
            projection.save_Gd(filename + "_G", save_diam_scalefactor);
            projection.save_section_stats(filename + "_stats", save_diam_scalefactor);

        }

        return false;

    }

    private boolean initFig5_prompt() {

        String pstr;

        double VF = 0.30;
        double phi = 0;

        title += " Fig 5";

        String Tud_List[] = {"A", "B", "C", "D"};
        pstr = Master.promptForInput(Tud_List, "choose panel simulations to compute", title, "B");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        String panel = pstr;

        String VF_List[] = {"0.15", "0.30", "0.45"};
        String phi_List[] = {"0", "70"};

        switch (panel) {

            case "A":
            case "B":

                pstr = Master.promptForInput(VF_List, "choose particle volume fraction", title, "0.30");

                if ((pstr == null) || (pstr.length() == 0)) {
                    return true; // cancel
                }

                VF = Double.parseDouble(pstr);

                break;

            case "C":
            case "D":

                pstr = Master.promptForInput(phi_List, "choose ϕ°", title, "0");

                if ((pstr == null) || (pstr.length() == 0)) {
                    return true; // cancel
                }

                phi = Double.parseDouble(pstr);

                break;

            default:
                return true;
        }

        if (setDirectoryPrompt(title)) {
            return true;
        }

        if (panel.equalsIgnoreCase("B")) {
            initFig5("BO", VF, phi); // opaque sims
        }

        return initFig5(panel, VF, phi);

    }

    private boolean initFig5_all() {

        double phi_0 = 0;
        String fileSuffix = "";

        if (initFig5("A", 0.15, phi_0)) {
            return true;
        }

        if (initFig5("A", 0.30, phi_0)) {
            return true;
        }

        if (initFig5("A", 0.45, phi_0)) {
            return true;
        }

        if (initFig5("B", 0.15, phi_0)) {
            return true;
        }

        if (initFig5("B", 0.30, phi_0)) {
            return true;
        }

        if (initFig5("B", 0.45, phi_0)) {
            return true;
        }

        if (initFig5("C", 0.30, phi_0)) {
            return true;
        }

        if (initFig5("C", 0.30, 70)) {
            return true;
        }

        if (initFig5("D", 0.30, phi_0)) {
            return true;
        }

        if (initFig5("D", 0.30, 70)) {
            return true;
        }

        return false;

    }

    private boolean initFig5(String panel, double volumeFraction, double phi) {

        int num_particles;
        int Tud;
        String filename, Tstr, Nstr, VFstr, pstr;

        boolean SEM = false;

        double[] izmax = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        boolean[] overlap_connect = {false, false};

        projection.histo_z_min = -180 / 1000.0; // um
        projection.histo_z_max = 50 / 1000.0; // um
        projection.histo_z_bwidth = 5 / 1000.0; // um
        projection.histo_z_norm = 1; // density, only for B

        //projection.sections_stats_save_phi = false;
        //num_sections = 1; // testing
        num_sections = 20;

        if (volumeFraction > 0) {
            set_particle_VF = volumeFraction;
        }

        switch (panel) {
            case "A":
                Tud = 2;
                projection.compute_Iz = true;
                izmax[1] = -1; // do not run
                break;
            case "B":
                Tud = 2;
                projection.compute_Cz = true;
                izmax[1] = 0.25;
                break;
            case "BO": // opaque
                Tud = 2;
                projection.compute_Cz = true;
                overlap_connect[1] = true;
                break;
            case "C":
                Tud = 1;
                //projection.phi_true = 0;
                //projection.phi_true = 70;
                projection.compute_Gd = true;
                izmax[1] = 0.25;
                break;
            case "D":
                Tud = 1;
                //projection.phi_true = 0;
                //projection.phi_true = 70;
                projection.compute_Gd = true;
                projection.histo_d_max = 85 / 1000.0;
                //imax[0] = 0.25;
                //imax[1] = 0.25;
                overlap_connect[1] = true;
                break;
            default:
                Master.log("error: unknown panel : " + panel);
                return true;
        }

        projection.phi_true = phi;

        section_thickness = particle_D_mean * Tud;

        //num_particles = (int) (500 * Tud);
        num_particles = 500;
        avg_num_particles = num_particles;

        geometry_xywidth = get_geometry_xywidth(num_particles, projection.phi_true);

        if (initProjection()) {
            return true;
        }

        Tstr = "_T" + Integer.toString(Tud);
        VFstr = "_VF" + Integer.toString((int) (100 * set_particle_VF));

        projection.name = "Fig5" + panel + " Particle Intersections";
        //projection.directory += "D3D_Fig5/";

        pstr = "_P" + Integer.toString((int) projection.phi_true);
        Nstr = "_N" + Integer.toString(num_particles);

        for (int i = 0; i < izmax.length; i++) {

            if (izmax[i] <= 0) {
                continue;
            }

            projection.Iz_max = izmax[i];
            projection.connect_overlapping_particles = overlap_connect[i];
            projection.computeProjection();
            projection.print_section_stats_avg();

            if (!save_to_disk) {
                continue;
            }

            filename = "D3D_Fig5" + panel + Tstr + Nstr + VFstr + pstr + IZmax_str(projection.Iz_max);

            if (overlap_connect[i]) {
                filename += "_Opaque";
            }

            if (projection.compute_Iz) {
                projection.save_Iz_stats(filename + "_Iz", SEM);
                projection.save_section_stats(filename + "_stats", save_diam_scalefactor);
            }

            if (projection.compute_Cz) {
                projection.save_Cz_stats(filename + "_Cz", SEM);
                projection.save_section_stats(filename + "_stats", save_diam_scalefactor);
            }

            if (projection.compute_Gd) {
                //projection.save_Gd(filename + "_G", save_diam_scalefactor);
                projection.save_Gd_stats(filename + "_G", save_diam_scalefactor, SEM);
                projection.save_section_stats(filename + "_stats", save_diam_scalefactor);
            }

        }

        return false;

    }

    private boolean init_ET_zstack_prompt() {

        String pstr, ET, analysis;

        title += " ET z-stack";

        String ET_List[] = {"ET10 - size", "ET11 - density"};
        pstr = Master.promptForInput(ET_List, "choose z-stack to simulate", title, "ET10 - size");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }
        
        if (pstr.startsWith("ET10")) {
            ET = "ET10";
            analysis = "size";
        } else if (pstr.startsWith("ET11")) {
            ET = "ET11";
            analysis = "density";
        } else {
            return true;
        }

        if (setDirectoryPrompt(title)) {
            return true;
        }

        return init_ET_zstack(ET, analysis);

    }

    private boolean init_ET_zstack_call() {

        String ET = "ET10";
        //String ET = "ET11";
        String analysis = "size";
        //String analysis = "density";

        return init_ET_zstack(ET, analysis);

    }

    private boolean init_ET_zstack(String ET, String analysis) {

        int phi_pnts, pscale = 1;
        int num_zscans, vesicles_per_zstack_limit = -1, diameters_limit = -1;

        String fileprefix, cvstr, tstr = "T0", folder;

        boolean particles_completely_within_section = false;

        //int num_zstacks = 1; // testing
        int num_zstacks = 100;

        int num_vesicles;

        double dz_scan;

        double zPSF = 0;

        double[][] phi;

        //set_particle_VF = 0.30;
        set_particle_VF = 0.45;
        MT_seed = "1635862732625"; // this gives 0.45 VF

        switch (ET) {
            case "ET10":
                //MT_seed = "1640375385243";
                num_vesicles = 142;
                particle_D_mean = 46.0 / 1000.0; // um
                particle_D_stdv = 4.0 / 1000.0; // um
                dz_scan = 0.6 / 1000.0; // um // from global fit
                break;
            case "ET11":
                num_vesicles = 271;
                particle_D_mean = 42.7 / 1000.0; // um
                particle_D_stdv = 3.4 / 1000.0; // um
                dz_scan = 0.5 / 1000.0; // um // from global fit
                break;
            default:
                return true;
        }

        if (ET.equals("ET10") && analysis.equals("density")) {

            phi = new double[2][2];

            phi[0][0] = 41;
            phi[0][1] = 0; // CV

            phi[1][0] = 41;
            phi[1][1] = 0.21; // CV

            //num_zscans = 197; // subsection, old value
            num_zscans = 181; // subsection, new value 
            //vesicles_per_zstack_limit = -1; // allows density analysis
            //geometry_xywidth = 0.3; // 6 x 6 voxels

        } else if (ET.equals("ET11") && analysis.equals("density")) {

            phi = new double[2][2];

            phi[0][0] = 42;
            phi[0][1] = 0; // CV

            phi[1][0] = 42;
            phi[1][1] = 0.16; // CV

            //num_zscans = 145; // subsection
            num_zscans = 261; // full z-stack
            //vesicles_per_zstack_limit = -1; // allows density analysis
            //geometry_xywidth = 0.3; // 6 x 6 voxels

        } else if (ET.equals("ET10") && analysis.equals("size")) {

            pscale = 10; // allows names with half degree

            //double[] phi_all = {30, 31, 32, 33}; // measured
            double[] phi_all = {30, 31, 32, 33, 34, 35, 36, 37, 37.5, 38, 38.5, 39, 39.5, 40, 40.5, 41, 41.5, 42, 42.5, 43, 43.5, 44, 44.5, 45, 45.5, 46, 46.5, 47, 47.5, 48, 49, 50};

            phi_pnts = phi_all.length;

            phi = new double[phi_pnts * 2][2];

            for (int i = 0, j = 0; i < phi_all.length; i += 1, j += 2) {
                phi[j][0] = phi_all[i];
                phi[j][1] = 0; // CV
                phi[j + 1][0] = phi_all[i];
                phi[j + 1][1] = 0.21; // CV
            }

            num_zscans = 291;
            //vesicles_per_zstack_limit = 141;
            //diameters_limit = 6704;
            //geometry_xywidth = 0.4;
            //set_particle_VF = 0.45; // lower density is OK

        } else {
            return true; // bad parameters
        }

        num_sections = num_zstacks;
        num_sections_extra = 10;

        section_thickness = (num_zscans - 1) * dz_scan;

        //int vesicles_per_zstack_extra = 230; // adjusted to get enough vesicles per section
        geometry_xywidth = get_geometry_xywidth(num_vesicles, projection.phi_true);

        Master.log("geometry xy width = " + geometry_xywidth + " " + project.spaceUnits);

        if (initProjection()) {
            return true;
        }

        if (projection.d.particles.length < set_num_particles) {
            Master.log("error: not enough vesicles : " + projection.d.particles.length + " of " + set_num_particles);
            return true;
        }

        projection.name = ET + " z-stack";

        //projection.directory += folder + "/";
        projection.compute_Fd = true;
        projection.compute_Gd = true;
        projection.compute_Pa = true;

        for (int i = 0; i < phi.length; i++) {

            if (!Double.isFinite(phi[i][0] * phi[i][1])) {
                continue;
            }

            projection.phi_true = phi[i][0];
            projection.phi_CV = phi[i][1];

            projection.computeZstack(dz_scan, zPSF, vesicles_per_zstack_limit, diameters_limit, particles_completely_within_section);
            projection.print_section_stats_avg();

            if (!save_to_disk) {
                continue;
            }

            fileprefix = "D3D_" + ET + "_" + tstr + "_";
            fileprefix += "P" + Integer.toString((int) Math.round(phi[i][0] * pscale));

            if (phi[i][1] == 0) {
                cvstr = "CV00";
            } else {
                cvstr = "CV" + Integer.toString((int) Math.round((100 * phi[i][1])));
            }

            fileprefix += cvstr;
            //fileprefix = "STDV0_" + fileprefix;

            //projection.save_Fd(fileprefix + "_F", save_diam_scalefactor);
            projection.save_Gd(fileprefix + "_G", save_diam_scalefactor);
            //projection.save_Pa(fileprefix + "_P");

            projection.save_section_stats(fileprefix + "_stats", save_diam_scalefactor);

        }

        //Master.log("true 3D density = " + projection.dvs.density + " particles/um^3");
        return false;

    }

    private boolean initDisector_prompt() {

        String pstr;

        title += " Disector";

        String N_List[] = {"300", "500", "2000"};
        pstr = Master.promptForInput(N_List, "choose number of particles per projection", title, "500");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        int numParticles = Integer.parseInt(pstr);

        String phi_lookup_list[] = {"10", "20", "30", "40", "10-40"};
        pstr = Master.promptForInput(phi_lookup_list, "choose lookup-section ϕ°", title, "20");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }
        
        int phi_lookup[];
        
        if (pstr.equalsIgnoreCase("10-40")) {
            phi_lookup = new int[4];
            phi_lookup[0] = 10;
            phi_lookup[1] = 20;
            phi_lookup[2] = 30;
            phi_lookup[3] = 40;
        } else {
            phi_lookup = new int[1];
            phi_lookup[0] = Integer.parseInt(pstr);
        }

        String phi_ref_bias_list[] = {"0", "5", "10", "15", "20", "0-20"};
        pstr = Master.promptForInput(phi_ref_bias_list, "choose reference-section bias ϕ°", title, "20");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }
        
        int phi_ref_bias[];

        if (pstr.equalsIgnoreCase("0-20")) {
            phi_ref_bias = new int[5];
            phi_ref_bias[0] = 0;
            phi_ref_bias[1] = 5;
            phi_ref_bias[2] = 10;
            phi_ref_bias[3] = 15;
            phi_ref_bias[4] = 20;
        } else {
            phi_ref_bias = new int[1];
            phi_ref_bias[0] = Integer.parseInt(pstr);
        }

        if (setDirectoryPrompt(title)) {
            return true;
        }

        return initDisector(numParticles, phi_lookup, phi_ref_bias);

    }

    private boolean initDisectorCall() {

        //int numParticles = 200;
        //int numParticles = 300;
        int numParticles = 500;
        //int numParticles = 2000;

        //int[] phi_lookup = {10, 20, 30, 40}; // degrees
        int[] phi_lookup = {40};

        //int[] phi_ref_bias = {0, 5, 10, 15, 20}; // degrees
        int[] phi_ref_bias = {0}; // no bias

        return initDisector(numParticles, phi_lookup, phi_ref_bias);

    }

    private boolean initDisector(int numParticles, int[] phi_lookup, int[] phi_ref_bias) {

        String filename, tstr, pstr, bstr;

        num_sections = 100;

        section_thickness = 0.3 * particle_D_mean;
        
        avg_num_particles = numParticles;

        tstr = "T" + (Math.round(10 * section_thickness / particle_D_mean) / 10.0);
        tstr = tstr.replace('.', 'p');
        
        set_particle_VF = 0.40;

        projection.name = "Disector";
        //projection.directory += "D3D_Disector_N" +Integer.toString(numParticles) + "/";

        //if (initProjection()) {
        //    return true;
        //}

        for (int i = 0; i < phi_lookup.length; i++) {

            projection.phi_true = phi_lookup[i];
            pstr = phi_str(phi_lookup[i]);

            geometry_xywidth = get_geometry_xywidth(numParticles, phi_lookup[i]);
            Master.log("geometry xy width = " + geometry_xywidth + " " + project.spaceUnits);

            if (initProjection()) {
                return true;
            }

            for (int j = 0; j < phi_ref_bias.length; j++) {

                projection.phi_disector_ref_bias = phi_ref_bias[j];
                bstr = "B" + Integer.toString((int) phi_ref_bias[j]);

                projection.computeDisector();
                //projection.print_section_stats_avg();

                if (!save_to_disk) {
                    continue;
                }

                filename = "D3D_Disector_" + tstr + "_" + pstr + "_" + bstr;

                projection.save_section_stats(filename + "_stats", save_diam_scalefactor);

            }

        }

        return false;

    }

    private boolean setDirectoryPrompt(String title) {

        String save_List[] = {"no", "yes"};

        String pstr = Master.promptForInput(save_List, "save results to disk?", title, "yes");

        if ((pstr == null) || (pstr.length() == 0)) {
            return true; // cancel
        }

        if (pstr.equalsIgnoreCase("yes")) {

            pstr = Master.getDirectory("");

            if ((pstr == null) || (pstr.length() == 0)) {
                return true; // cancel
            }

            directory = pstr;
            projection.directory = pstr;

            Master.log("projection directory = " + projection.directory);

            save_to_disk = true;

        } else {

            save_to_disk = false;

        }

        return false;

    }

    private double get_geometry_xywidth(int num_particles, double phi) {

        double radius = 0.5 * particle_D_mean;
        double particle_volume = 4.0 * Math.PI * Math.pow(radius, 3.0) / 3.0;
        double density3D = set_particle_VF / particle_volume; // particles/um^3

        double zeta = DiffusantParticlesProjectionKeiding.zeta(section_thickness, particle_D_mean, phi);

        double density2D = density3D * zeta; // particles/um^2

        double square_area = num_particles / density2D; // um^2

        return Math.sqrt(square_area); // um

    }

    private boolean initProjection() {

        double dx, radius_mean, radius_stdv, zwidth;
        String str;

        initDirectoryJason();

        DiffusantParticles dvs;

        RunMonteCarlo mc = project.newMonteCarlo();

        mc.PBC = geometry_PBC;

        mc.minParticleStep = 0.001;

        //mc.removeParticleOverlap = true;
        mc.initParticleRandomTrialLimit = (long) 1e8;

        project.name = "Monte Carlo 2D Projection";

        if (MT_seed.length() > 0) {
            project.seedMT = MT_seed;
            Master.initMersenneTwister();
        }

        radius_mean = particle_D_mean / 2.0;
        radius_stdv = particle_D_stdv / 2.0;

        str = "set 3D particle diameter = ";
        str += numformat.format(2 * 1000 * radius_mean);
        str += " ± " + numformat.format(2 * 1000 * radius_stdv);
        str += " nm";
        Master.log(str);

        dx = particle_D_mean;

        if (geometry_zwidth > 0) {
            zwidth = geometry_zwidth;
        } else {
            zwidth = (section_thickness + particle_D_mean) * (num_sections + num_sections_extra);
            zwidth = Math.max(zwidth, 1);
        }

        project.dx = dx;
        project.simTime = 1000;

        geometry.forceEvenVoxels = true;
        geometry.resizeWithSpace(geometry_xywidth, geometry_xywidth, zwidth);

        //int xvoxels = geometry.xVoxels;
        //int yvoxels = geometry.yVoxels;
        //int zvoxels = geometry.zVoxels;
        double d = 0.060e-3;

        dvs = new DiffusantParticles(project, "Particles", 0, d, null, radius_mean, radius_stdv);

        if (project.numDiffusants() > 0) {
            project.killDiffusantsAll();
        }

        project.addDiffusant(dvs);

        dvs.save.save2TextFile = false;

        dvs.colorReady.setColor("[r=148,g=183,b=80]");
        dvs.color.setColor("[r=0,g=0,b=153]");
        dvs.color.setGradientColor("[r=255,g=255,b=255]");

        if (set_num_particles > 0) {
            dvs.setNumParticles = set_num_particles;
        } else if (set_particle_VF > 0) {
            dvs.setVolumeFraction = set_particle_VF;
        } else {
            Master.exit("unknown particle density");
        }

        dvs.setImmobilePercent = 0;

        if (set_particle_VF > 0.30) {
            mc.sortParticlesByRadius = true;
        }

        //if (particle_diam_min > 0) {
        //    dvs.setRadiusMin = 0.5 * particle_diam_min;
        //}
        //if (particle_diam_max > 0) {
        //    dvs.setRadiusMax = 0.5 * particle_diam_max;
        //}
        mc.initAll();

        if (num_sections <= 0) {
            return true;
        }

        //projection.directory = directory;
        projection.num_sections = num_sections;
        projection.section_thickness = section_thickness;
        projection.histo_save_xbins = save_histo_xbins;

        projection.d = dvs;

        return false;

    }

    private String phi_str(double phi_degrees) {
        if (phi_degrees < 10) {
            return "P0" + Integer.toString((int) phi_degrees);
        }
        return "P" + Integer.toString((int) phi_degrees);
    }

    private int[] phiList(int ibgn, int idelta, int ipnts) {
        int[] phi = new int[ipnts];
        for (int i = 0; i < ipnts; i++) {
            phi[i] = ibgn + i * idelta;
        }
        return phi;
    }

    private String IZmax_str(double imax_fraction) {
        if (Double.isInfinite(imax_fraction)) {
            return "_I99";
        }
        int i = (int) (imax_fraction * 100);
        if (i < 10) {
            return "_I0" + Integer.toString((int) (i));
        }
        return "_I" + Integer.toString((int) (i));
    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (projection != null) {
            projection.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (projection != null) {
            addBlankParam();
            projection.createVector(true);
            addVector(projection.getVector());
            projection.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (projection != null) {
            projection.updateVector(v);
        }

    }

}
