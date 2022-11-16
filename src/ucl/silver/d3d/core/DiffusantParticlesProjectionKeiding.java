package ucl.silver.d3d.core;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

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
 * @version 2.1.projection.demo
 */
public class DiffusantParticlesProjectionKeiding extends ParamVector {

    public String directory = "";

    public DiffusantParticles d = null; // the spherical particles

    public int num_sections = 1; // number of sections
    public double section_thickness = 0; // um
    
    public int num_particles = -1; // number of particles:
                                        // (-1) grab all particles
                                        // (> 0) set number limit

    public double roi_xy_width = -1; // width of square ROI
                                     // -1 to use full xy-width of geometry

    public boolean above = true; // count particles above section
    public boolean within = true; // count particles within section
    public boolean below = true; // count particles below section

    public double phi_true = 0; // cap angle limit // degrees
    public double phi_CV = 0; // coefficient of variation of phi
    public double phi_disector_ref_bias = 0; // degress

    public double Iz_max = Double.POSITIVE_INFINITY; // max limit for the sum of 2D intersections
    public boolean connect_overlapping_particles = false; // connect overlapping particles

    public String print_space_units = ""; // log outputs (analysis is in project.spaceUnits)

    public boolean compute_Fd = false; // Fd(d) // 3D diam histo, only particles in projection
    public boolean compute_Gd = false; // G(d) // 2D diam histo
    public boolean compute_Gd_Goldsmith_transform = false; // Goldsmith_G_to_F_Transform
                                            // compute F from G via Goldsmith transform 

    public boolean compute_Pa = false; // histo of phi as a function of degrees (angle)

    public boolean compute_Cz = false; // density (concentration) // histo as a function of z-axis
    public boolean compute_Iz = false; // sum of 2D intersections  // histo as a function of z-axis

    public double histo_d_min = 0; // diameter histograms
    public double histo_d_max = 70 / 1000.0;
    public double histo_d_bwidth = 1 / 1000.0;
    public int histo_d_norm = 1; // (0) count (1) probability (2) frequency

    public double histo_a_min = 0; // angle histograms
    public double histo_a_max = 90;
    public double histo_a_bwidth = 1;
    public int histo_a_norm = 1; // (0) count (1) probability (2) frequency

    public double histo_z_min = -150 / 1000.0; // z-distance histograms
    public double histo_z_max = 50 / 1000.0;
    public double histo_z_bwidth = 2 / 1000.0;
    public int histo_z_norm = 1; // (0) count (1) density

    public boolean histo_save_xbins = true;
    public boolean sections_stats_save_phi = true;

    public double D_3D_mean_true = 0;
    public double D_3D_stdv_true = 0;
    public double den3D_true = 0;

    public SectionStats[] particle_stats;

    public double[][] Fd = null; // j-sections, j-bins
    public double[][] Gd = null; // j-sections, j-bins
    public double[][] Gd_transform = null; // j-sections, j-bins
    public double[][] Pa = null; // j-sections, j-bins
    public double[][] Cz = null; // z-density histo (concentration) // j-sections, j-bins
    public double[][] Iz = null; // z-intersection histo // j-sections, j-bins

    public Save save = null;

    private Geometry geometry = null;

    private BufferedWriter bw = null;

    private final NumberFormat formatter = new DecimalFormat("0.######");

    private double z1, z2; // bottom/top of current section

    @Override
    public String units(String name) {
        if (name.equalsIgnoreCase("section_thickness")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("roi_xy_width")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("phi_degrees")) {
            return "°";
        }
        if (name.equalsIgnoreCase("phi_stdv_degrees")) {
            return "°";
        }
        if (name.equalsIgnoreCase("histo_d_min")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("histo_d_max")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("histo_d_bwidth")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("histo_z_min")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("histo_z_max")) {
            return project.spaceUnits;
        }
        if (name.equalsIgnoreCase("histo_z_bwidth")) {
            return project.spaceUnits;
        }
        return super.units(name);
    }

    @Override
    public boolean canEdit(String name) {
        return super.canEdit(name);
    }

    private double spaceScale(String spaceUnits) {
        if (project.spaceUnits.equalsIgnoreCase("μm")) {
            if (spaceUnits.equalsIgnoreCase("μm")) {
                return 1;
            }
            if (spaceUnits.equalsIgnoreCase("nm")) {
                return 1000;
            }
        }
        return 1;
    }

    public DiffusantParticlesProjectionKeiding(Project p, String NAME, DiffusantParticles D) {

        super(p);

        if (NAME.length() > 0) {
            name = NAME;
        } else {
            name = "Projection";
        }

        geometry = p.geometry;

        d = D;

        print_space_units = project.spaceUnits;

        save = new Save(p);

        createVector(true);

    }

    private double[] update_all_D_3D_stats() {

        if (d == null) {
            return null;
        }

        double scale = spaceScale(print_space_units);

        double[] D_3D = d.getParameter("diameter");
        
        D_3D_mean_true = 2 * d.setRadiusMean;
        D_3D_stdv_true = 2 * d.setRadiusStdv;
        den3D_true = d.density;

        String str = "true 3D diameters = ";
        str += formatter.format(D_3D_mean_true * scale);
        str += " ± " + formatter.format(D_3D_stdv_true * scale);
        str += " " + print_space_units;
        str += " (n = " + D_3D.length + ")";
        
        Master.log(str);
        
        str = "true 3D density = " + den3D_true + " particles/um^3";
        
        Master.log(str);

        return D_3D;

    }

    public double getPhiRadians() {

        double p = phi_true;

        if (phi_CV > 0) {
            p += Master.randomGauss() * phi_CV * p;
        }

        p = Math.max(p, 0);
        p = Math.min(p, 90);

        return p * Math.PI / 180.0;

    }

    public double section_xwidth() {
        if (roi_xy_width > geometry.xWidth) {
            Master.exit("roi_xy_width is larger than geometry");
        }
        if (roi_xy_width > 0) {
            return roi_xy_width;
        } else {
            return geometry.xWidth;
        }
    }

    public double section_ywidth() {
        if (roi_xy_width > geometry.yWidth) {
            Master.exit("roi_xy_width is larger than geometry");
        }
        if (roi_xy_width > 0) {
            return roi_xy_width;
        } else {
            return geometry.yWidth;
        }
    }

    public double section_area() {
        return section_xwidth() * section_ywidth();
    }

    public double zeta_true() {
        return zeta(section_thickness, D_3D_mean_true, phi_true);
    }

    public static double zeta(double SECTION_THICKNESS, double DIAM_3D_MEAN, double PHI_DEGREES) {
        double phi_radians = PHI_DEGREES * Math.PI / 180.0;
        return (SECTION_THICKNESS + DIAM_3D_MEAN * Math.cos(phi_radians));
    }
    
    public static double kv_Weibel(double SECTION_THICKNESS, double DIAM_3D_MEAN, double DIAM_3D_STDV, double PHI_DEGREES) {
        double g = SECTION_THICKNESS / DIAM_3D_MEAN;
        double phi_radians = PHI_DEGREES * Math.PI / 180.0;
        double x = 1 - Math.cos(phi_radians);
        double dm2 = DIAM_3D_MEAN * DIAM_3D_MEAN;
        double dm3 = DIAM_3D_MEAN * DIAM_3D_MEAN * DIAM_3D_MEAN;
        double ds2 = DIAM_3D_STDV * DIAM_3D_STDV;
        double m2 = (dm2 + ds2) / dm2;
        double m3 = DIAM_3D_MEAN * (dm2 + 3 * ds2) / dm3;
        double kv = 2 * m3 / (2 * m3 + 3 * g * m2 - 3 * x * x + x * x * x);
        return kv;
    }
    
    public double phi_eq2(double den2D, double sectionT) {
        
        double zeta = den2D / den3D_true; // Eq3
        double cos_phi = (zeta - sectionT) / D_3D_mean_true; // Eq2
        double phi_radians; // phi computed from measured den2D
        
        if ((cos_phi >= 0) && (cos_phi <= 1)) {
            phi_radians = Math.acos(cos_phi);
            return (phi_radians * 180.0 / Math.PI);
        } else {
            return Double.NaN;
        }
        
    }

    private int histo_d_nbins() {
        return 1 + (int) ((histo_d_max - histo_d_min) / histo_d_bwidth);
    }

    private int histo_a_nbins() {
        return 1 + (int) ((histo_a_max - histo_a_min) / histo_a_bwidth);
    }

    private int histo_z_nbins() {
        return 1 + (int) ((histo_z_max - histo_z_min) / histo_z_bwidth);
    }

    public boolean computeDisector() {

        int i_particle;
        int iRef, iRefAndLookup, iDisectorRef;
        int iLostCapRef, iLostCapLookup;
        double avgRef = 0, avgRefAndLookup = 0, avgDisectorRef = 0;
        double avgLostCapRef = 0, avgLostCapLookup = 0;
        double zPlaneRef, zPlaneLookup, z_offset, z_step;
        double dv_xyz, norm_scale;
        double phiRef, phiLookup;
        
        int histo_z_nbins = histo_z_nbins();

        double scale = spaceScale(print_space_units);
        double area_ROI = section_area();
        
        double[] dz_top;
        double[][] histo;
        
        Particle[] p_ref; // particles in reference section
        Particle[] p_lookup; // particles in lookup section

        if (num_sections <= 1) {
            return true;
        }
        
        if (update_all_D_3D_stats() == null) {
            return true;
        }

        Master.log("computing 2D projections... ");
        Master.log("# sections = " + Integer.toString(num_sections));
        Master.log("x width = " + formatter.format(section_xwidth()) + " " + project.spaceUnits);
        Master.log("y width = " + formatter.format(section_ywidth()) + " " + project.spaceUnits);
        Master.log("area = " + formatter.format(area_ROI) + " " + project.spaceUnits + "^2");
        Master.log("T = " + formatter.format(section_thickness) + " " + project.spaceUnits);
        
        particle_stats = new SectionStats[num_sections];

        for (int i = 0; i < particle_stats.length; i++) {
            particle_stats[i] = new SectionStats();
        }
        
        if (compute_Cz) {
            Cz = new double[num_sections + 1][histo_z_nbins];
        }

        z_step = geometry.zWidth / (num_sections + 4);
        z_offset = 2 * z_step;
        
        Master.log("z-distance between sections = " + formatter.format(z_step * scale) + " " + print_space_units);
        
        if (phi_disector_ref_bias > 0) {
            phiRef = phi_true + phi_disector_ref_bias;
        } else {
            phiRef = phi_true;
        }
        
        phiLookup = phi_true;
        
        Master.log("lookup phi = " + formatter.format(phiLookup) + " ± " + formatter.format(phi_CV * phiLookup) + "°");
        Master.log("reference phi = " + formatter.format(phiRef) + " ± " + formatter.format(phi_CV * phiRef) + "° (bias = +" + formatter.format(phi_disector_ref_bias) + "°)");

        for (int isection = 0; isection < num_sections; isection++) {

            zPlaneRef = geometry.z1 + z_offset + isection * z_step;
            zPlaneLookup = zPlaneRef + section_thickness;

            phi_true = phiRef;
            p_ref = getSectionParticles(zPlaneRef); // grab all particles for zeta = T + D

            phi_true = phiLookup;
            p_lookup = getSectionParticles(zPlaneLookup);
            
            //Master.log("z-plane = " + zPlaneRef);
            //Master.log("z-plane lookup = " + zPlaneLookup);

            if (p_ref == null) {
                return true;
            }
            
            for (Particle prf : p_ref) {

                if (!prf.within && (prf.phi_radians > 0)) {

                    if (prf.r2D < prf.r2D_min) {
                        prf.exclude = true;
                        particle_stats[isection].exclude_dmin++;
                    }

                }

            }

            for (Particle plu : p_lookup) {

                if (!plu.within && (plu.phi_radians > 0)) {

                    if (plu.r2D < plu.r2D_min) {
                        plu.exclude = true;
                        //particle_stats[isection].exclude_dmin++;
                    }

                }

            }
            
            iRef = 0;
            iRefAndLookup = 0;
            iDisectorRef = 0;
            
            iLostCapRef = 0;
            iLostCapLookup = 0;

            for (Particle prf : p_ref) {
                
                iRef++;

                if (prf.exclude) {
                    iLostCapRef++; // cannot see particle in reference section, do not count
                    continue;
                }

                for (Particle plu : p_lookup) {
                    if (prf.id == plu.id) {
                        iRefAndLookup++; // particle is in both sections
                        if (plu.exclude) {
                            iLostCapLookup++;
                            prf.exclude = false;  // cannot see reference particle in lookup section, so reference particle is counted
                        } else {
                            prf.exclude = true;  // can see reference particle in lookup section, so reference particle is NOT counted
                        }
                        break;
                    }
                }

                if (!prf.exclude) {
                    iDisectorRef++;
                }

            }
            
            // compute stats of disector particles
            
            dz_top = new double[p_ref.length];

            i_particle = 0;

            for (Particle prf : p_ref) {
                
                if (prf.exclude) {
                    dz_top[i_particle] = Double.NaN;
                    i_particle++;
                    continue;
                }
                
                dz_top[i_particle] = prf.dz_top;

                if (prf.above) {
                    particle_stats[isection].above++;
                } else if (prf.within) {
                    particle_stats[isection].within++;
                } else if (prf.below) {
                    particle_stats[isection].below++;
                }

                i_particle++;

            }

            particle_stats[isection].N3D = (int) iDisectorRef;
            particle_stats[isection].N2D = (int) iDisectorRef;
            particle_stats[isection].den3D = den3D_true;
            particle_stats[isection].den2D = 1.0 * iDisectorRef / area_ROI;
            particle_stats[isection].exclude_dmin *= 100.0 / (double) iRef;
            particle_stats[isection].exclude_imax *= 100.0 / (double) iRef;
            particle_stats[isection].above *= 100.0 / (double) iDisectorRef;
            particle_stats[isection].within *= 100.0 / (double) iDisectorRef;
            particle_stats[isection].below *= 100.0 / (double) iDisectorRef;
            
            if (Cz != null) {

                histo = Utility.histogram(dz_top, histo_z_min, histo_z_max, histo_z_bwidth, 0);

                dv_xyz = area_ROI * histo_z_bwidth;

                switch (histo_z_norm) {
                    case 0: // count
                        norm_scale = 1;
                        break;
                    case 1: // density, particles/um^3
                        norm_scale = 1 / dv_xyz;
                        break;
                    //case 2: // volume fraction // using dvs.volumeMean is only approximate
                        //norm_scale = dvs.volumeMean / dv_xyz;
                        //break;
                    default:
                        return true;
                }

                for (int ih = 0; ih < histo.length; ih++) {
                    Cz[0][ih] = histo[ih][0];
                    Cz[isection + 1][ih] = histo[ih][1] * norm_scale;
                }

            }
            
            avgRef += iRef;
            avgRefAndLookup += iRefAndLookup;
            avgDisectorRef += iDisectorRef;
            
            avgLostCapRef += iLostCapRef;
            avgLostCapLookup += iLostCapLookup;

        }
        
        avgRef /= 1.0 * num_sections;
        avgRefAndLookup /= 1.0 * num_sections;
        avgDisectorRef /= 1.0 * num_sections;
        
        avgLostCapRef /= 1.0 * num_sections;
        avgLostCapLookup /= 1.0 * num_sections;
        
        double avgDen3D = avgDisectorRef / (area_ROI * section_thickness);
        double deltaDen3D = 100 * (avgDen3D - den3D_true) / den3D_true;
        
        Master.log("Disector averages for reference (ref) and lookup sections (n=" + num_sections + "):");
        Master.log("(1) # particles in ref = " + avgRef);
        Master.log("(2) # lost caps in ref = " + avgLostCapRef);
        Master.log("(3) # adjusted particles in ref (1)-(2) = " + (avgRef - avgLostCapRef));
        Master.log("(4) # particles in ref & lookup = " + avgRefAndLookup);
        Master.log("(5) # lost caps in lookup = " + avgLostCapLookup);
        Master.log("(6) # adjusted particles in ref & lookup (4)-(5) = " + (avgRefAndLookup - avgLostCapLookup));
        Master.log("(7) # final disector particles in ref (3)-(6) = " + avgDisectorRef);
        Master.log("(8) estimated 3D density = " + avgDen3D + " particles/um^3 [" + deltaDen3D + "%]");
        Master.log("(9) true 3D density = " + den3D_true + " particles/um^3");

        return false;

    }

    public boolean computeProjection() {

        double dp_xyz;
        double z_plane, z_offset, z_step = 0;
        double z_midpoint = Double.NaN;
        double z_min, z_max, zeta_min_max, den2D;
        double i_total;
        double norm_scale;

        double[][] histo, histo_transform;

        Particle[] particles;

        double scale = spaceScale(print_space_units);

        int histo_d_nbins = histo_d_nbins();
        int histo_a_nbins = histo_a_nbins();
        int histo_z_nbins = histo_z_nbins();

        double[] diam2D, diam3D, dz_top, phi_degrees;
        double[][] zintersections;

        if (!project.monteCarlo.PBC) {
            Master.exit("error: projections should be computed with PBC to get correct density measurements");
        }

        if ((section_thickness == 0) && !Double.isInfinite(Iz_max)) {
            Master.exit("it does not make sense to use Iz_max when section_thickness = 0");
        }

        if (update_all_D_3D_stats() == null) {
            return true;
        }

        double areaxy = section_area();

        Master.log("computing 2D projections... ");
        Master.log("# sections = " + Integer.toString(num_sections));
        Master.log("x width = " + formatter.format(section_xwidth()) + " " + project.spaceUnits );
        Master.log("y width = " + formatter.format(section_ywidth()) + " " + project.spaceUnits );
        Master.log("section area = " + formatter.format(areaxy) + " " + project.spaceUnits + "^2");
        Master.log("T = " + formatter.format(section_thickness * scale) + " " + print_space_units);
        Master.log("zeta = " + formatter.format(zeta_true() * scale) + " " + print_space_units);

        Master.log("phi = " + formatter.format(phi_true) + " ± " + formatter.format(phi_CV * phi_true) + "°");

        Master.log("zintersect_max = " + formatter.format(Iz_max));

        if (num_sections > 1) {
            z_step = geometry.zWidth / (num_sections + 4);
            Master.log("z-distance between sections = " + formatter.format(z_step * scale) + " " + print_space_units);
        }

        particle_stats = new SectionStats[num_sections];

        for (int i = 0; i < particle_stats.length; i++) {
            particle_stats[i] = new SectionStats();
        }

        if (compute_Fd) {
            Fd = new double[num_sections + 1][histo_d_nbins];
        }

        if (compute_Gd) {
            Gd = new double[num_sections + 1][histo_d_nbins];
        }

        if (compute_Gd_Goldsmith_transform) {
            Gd_transform = new double[num_sections + 1][histo_d_nbins];
        }

        if (compute_Pa) {
            Pa = new double[num_sections + 1][histo_a_nbins];
        }

        if (compute_Cz) {
            Cz = new double[num_sections + 1][histo_z_nbins];
        }

        if (compute_Iz) {
            Iz = new double[num_sections + 1][histo_z_nbins];
        }

        if (num_sections > 1) {
            //z_offset = 1 * (section_thickness + D_3D_mean);
            z_offset = 2 * z_step; // avoid z-borders
        } else {
            z_offset = 0;
        }

        for (int isection = 0; isection < num_sections; isection++) {

            if (num_sections > 1) {
                z_plane = geometry.z1 + z_offset + isection * z_step;
            } else if (Double.isFinite(z_midpoint) && (z_midpoint > geometry.z1) && (z_midpoint < geometry.z2)) {
                z_plane = z_midpoint;
            } else {
                z_plane = geometry.zCenter;
            }

            //Master.log("z-plane = " + zPlaneRef);

            particles = getSectionParticles(z_plane); // grab all particles for zeta = T + D

            if (particles == null) {
                return true;
            }

            i_total = 0;

            for (Particle p : particles) { // exclude particles

                if (!above && p.above) {
                    p.exclude = true;
                    continue;
                }

                if (!within && p.within) {
                    p.exclude = true;
                    continue;
                }

                if (!below && p.below) {
                    p.exclude = true;
                    continue;
                }

                if (Double.isFinite(Iz_max) && (p.areaXY_intersect_norm() > Iz_max)) {
                    p.exclude = true;
                    particle_stats[isection].exclude_imax++;
                    continue;
                }

                if (!p.within && (p.phi_radians > 0)) { // Keiding cap-angle limit

                    if (p.r2D < p.r2D_min) {
                        p.exclude = true;
                        particle_stats[isection].exclude_dmin++;
                        continue;
                    }

                }
                
                if ((num_particles > 0) && (i_total >= num_particles)) {
                    p.exclude = true;
                }

                if (!p.exclude) {
                    i_total++;
                }

            }
            
            if ((num_particles > 0) && (i_total < num_particles)) {
                Master.exit("error: number of particles is less than limit: " + num_particles);
            }

            if (i_total == 0) {
                break;
            }
            
            if (connect_overlapping_particles) {
                particles = connectOverlappingParticles(particles);
            }

            diam2D = new double[particles.length];
            diam3D = new double[particles.length];
            dz_top = new double[particles.length];
            phi_degrees = new double[particles.length];
            zintersections = new double[particles.length][2];

            i_total = 0;
            z_min = Double.POSITIVE_INFINITY;
            z_max = Double.NEGATIVE_INFINITY;
            
            for (int i = 0; i < particles.length; i++) {
                
                if (particles[i] == null) {
                    continue;
                }

                phi_degrees[i] = particles[i].phi_radians * 180.0 / Math.PI;

                if (particles[i].exclude) {
                    diam2D[i] = Double.NaN;
                    diam3D[i] = Double.NaN;
                    dz_top[i] = Double.NaN;
                    zintersections[i][0] = Double.NaN;
                    zintersections[i][1] = Double.NaN;
                    continue;
                }

                diam2D[i] = particles[i].r2D * 2.0;
                diam3D[i] = particles[i].r3D * 2.0;
                dz_top[i] = particles[i].dz_top;
                z_min = Math.min(z_min, particles[i].z);
                z_max = Math.max(z_max, particles[i].z);

                zintersections[i][0] = dz_top[i];
                zintersections[i][1] = particles[i].areaXY_intersect_norm();

                //if (print_diameters) {
                //    Master.log(formatter.format(diam2D[i_particle]));
                //}

                if (particles[i].above) {
                    particle_stats[isection].above++;
                } else if (particles[i].within) {
                    particle_stats[isection].within++;
                } else if (particles[i].below) {
                    particle_stats[isection].below++;
                }
                
                particle_stats[isection].AF += Math.PI * particles[i].r2D * particles[i].r2D;

                i_total++;

            }

            zeta_min_max = Math.abs(z_max - z_min);
            den2D = i_total / areaxy;

            particle_stats[isection].N3D = (int) i_total;
            particle_stats[isection].N2D = (int) i_total;
            particle_stats[isection].den3D = den3D_true;
            particle_stats[isection].den2D = den2D;
            particle_stats[isection].exclude_dmin *= 100 / i_total;
            particle_stats[isection].exclude_imax *= 100 / i_total;
            particle_stats[isection].above *= 100 / i_total;
            particle_stats[isection].within *= 100 / i_total;
            particle_stats[isection].below *= 100 / i_total;
            particle_stats[isection].d_3D_mean = Utility.average(diam3D);
            particle_stats[isection].d_3D_stdv = Utility.stdev(diam3D);
            particle_stats[isection].d_3D_min = Utility.getMin(diam3D);
            particle_stats[isection].d_3D_max = Utility.getMax(diam3D);
            particle_stats[isection].d_2D_mean = Utility.average(diam2D);
            particle_stats[isection].d_2D_stdv = Utility.stdev(diam2D);
            particle_stats[isection].d_2D_min = Utility.getMin(diam2D);
            particle_stats[isection].d_2D_max = Utility.getMax(diam2D);
            particle_stats[isection].AF /= areaxy;
            particle_stats[isection].phi_mean = Utility.average(phi_degrees);
            particle_stats[isection].phi_stdv = Utility.stdev(phi_degrees);
            particle_stats[isection].phi_min = Utility.getMin(phi_degrees);
            particle_stats[isection].phi_max = Utility.getMax(phi_degrees);
            particle_stats[isection].phi_eq2 = phi_eq2(den2D, section_thickness);
            particle_stats[isection].zeta_min_max = zeta_min_max;

            if (Fd != null) {

                histo = Utility.histogram(diam3D, histo_d_min, histo_d_max, histo_d_bwidth, histo_d_norm);

                for (int ih = 0; ih < histo.length; ih++) {
                    Fd[0][ih] = histo[ih][0];
                    Fd[isection + 1][ih] = histo[ih][1];
                }

            }

            if (Gd != null) {

                histo = Utility.histogram(diam2D, histo_d_min, histo_d_max, histo_d_bwidth, histo_d_norm);

                for (int ih = 0; ih < histo.length; ih++) {
                    Gd[0][ih] = histo[ih][0];
                    Gd[isection + 1][ih] = histo[ih][1];
                }

            }

            if (Gd_transform != null) {

                histo = Utility.histogram(diam2D, histo_d_min, histo_d_max, histo_d_bwidth, 0);

                switch (histo_d_norm) {
                    default:
                        norm_scale = 1;
                        break;
                    case 1: // probability
                        norm_scale = 1 / (i_total * histo_d_bwidth);
                        break;
                    case 2: // frequency
                        norm_scale = 1 / i_total;
                        break;
                }

                // SHOULD norm_scale = 1?
                // Goldsmith transform on histo count?
                for (int i = 0; i < histo.length; i++) {
                    histo[i][1] *= norm_scale;
                }

                histo_transform = Goldsmith_G_to_F_Transform(histo, section_thickness, D_3D_mean_true);

                for (int ih = 0; ih < histo_transform.length; ih++) {
                    Gd_transform[0][ih] = histo_transform[ih][0];
                    Gd_transform[isection + 1][ih] = histo_transform[ih][1];
                }

            }

            if (Pa != null) {

                histo = Utility.histogram(phi_degrees, histo_a_min, histo_a_max, histo_a_bwidth, histo_a_norm);

                for (int ih = 0; ih < histo.length; ih++) {
                    Pa[0][ih] = histo[ih][0];
                    Pa[isection + 1][ih] = histo[ih][1];
                }

            }

            if (Cz != null) {

                histo = Utility.histogram(dz_top, histo_z_min, histo_z_max, histo_z_bwidth, 0);

                dp_xyz = areaxy * histo_z_bwidth;

                switch (histo_z_norm) {
                    case 0: // count
                        norm_scale = 1;
                        break;
                    case 1: // density, particles/um^3
                        norm_scale = 1 / dp_xyz;
                        break;
                    //case 2: // volume fraction // using d.volumeMean is only approximate
                        //norm_scale = d.volumeMean / dp_xyz;
                        //break;
                    default:
                        return true;
                }

                for (int ih = 0; ih < histo.length; ih++) {
                    Cz[0][ih] = histo[ih][0];
                    Cz[isection + 1][ih] = histo[ih][1] * norm_scale;
                }

            }

            if (Iz != null) {

                histo = Utility.histogramIJ(zintersections, histo_z_min, histo_z_max, histo_z_bwidth);

                for (int ih = 0; ih < histo.length; ih++) {
                    Iz[0][ih] = histo[ih][0];
                    Iz[isection + 1][ih] = histo[ih][1];
                }

            }

        }

        return false;

    }
    
    public Particle[] connectOverlappingParticles(Particle[] particles) {
        
        int jmax, icount_old = 0, icount_new = 0, ioverlaps;
        double zmax, area2add, area2addTotal, r_old, r_new;

        Particle[] p_old = new Particle[particles.length];
        Particle[] p_new = new Particle[particles.length];

        System.arraycopy(particles, 0, p_old, 0, p_old.length);
        
        for (int i = 0; i < p_new.length; i++) {
            if (p_old[i].exclude) {
                p_old[i] = null;
            } else {
                icount_old++;
            }
        }
        
        for (int i = 0; i < p_new.length; i++) {
            p_new[i] = null;
        }
        
        for (int i = 0; i < p_new.length; i++) {

            zmax = Double.NEGATIVE_INFINITY;
            jmax = -1;

            for (int j = 0; j < p_old.length; j++) {
                if (p_old[j] == null) {
                    continue;
                }
                if (p_old[j].z > zmax) {
                    zmax = p_old[j].z;
                    jmax = j;
                }
            }
            
            if (jmax == -1) {
                break; // no more particles
            }
            
            p_new[i] = p_old[jmax]; // particle closest to surface, add to new array
            p_old[jmax] = null; // remove from old array
            area2addTotal = 0;
            ioverlaps = 0;
            
            for (int j = 0; j < p_old.length; j++) { // find overlapping particles, which will be below this particle
                if (p_old[j] == null) {
                    continue;
                }
                area2add = p_new[i].areaXY_intersection_area2add(p_old[j]);
                if (area2add == 0) {
                    continue; // no/little overlap, go to next particle
                } else if (Double.isNaN(area2add)) {
                    area2add = 0; // complete overlap, but does not add area (large on top of small)
                }
                area2addTotal += area2add; // overlap, add extra area and remove
                ioverlaps++;
                p_old[j] = null;
            }
            
            if (area2addTotal > 0) {
                r_old = p_new[i].r2D;
                r_new = Math.sqrt((area2addTotal + Math.PI * r_old * r_old) / Math.PI);
                p_new[i].r2D = r_new; // particle grows due to overlaps
                //Master.log("old r=" + r_old + ", new r=" + r_new + ", n = " + ioverlaps);
            }
            
            icount_new++;

        }
        
        //Master.log("overlap concat: old n=" + icount_old + ", new n=" + icount_new);
        
        Particle[] p_new2 = new Particle[icount_new];
        
        System.arraycopy(p_new, 0, p_new2, 0, p_new2.length);

        return p_new2;

    }

    public Particle[] getSectionParticles(double z0) {
        // z0 is z-midpoint of section
        int pcount;
        double dz, dz_top, phiradians;
        double x1, y1, x2, y2;

        Particle[] particles1;
        Particle[] particles2;

        if ((d == null) || (d.particles == null)) {
            return null;
        }

        double zeta, volume;

        if ((d.numParticles > 200) && (den3D_true > 0)) {
            zeta = section_thickness + D_3D_mean_true; // phi = 0
            volume = geometry.xWidth * geometry.yWidth * zeta;
            pcount = (int) (den3D_true * volume);
            pcount *= 2; // extra
        } else {
            pcount = d.numParticles;
        }

        if (Double.isFinite(pcount) && (pcount > 0)) {
            particles1 = new Particle[pcount];
        } else {
            return null;
        }

        if (roi_xy_width > 0) {
            x1 = geometry.xCenter - roi_xy_width * 0.5;
            x2 = geometry.xCenter + roi_xy_width * 0.5;
            y1 = geometry.yCenter - roi_xy_width * 0.5;
            y2 = geometry.yCenter + roi_xy_width * 0.5;
        } else {
            x1 = geometry.x1;
            x2 = geometry.x2;
            y1 = geometry.y1;
            y2 = geometry.y2;
        }

        z1 = z0 - section_thickness * 0.5;
        z2 = z0 + section_thickness * 0.5;

        if (x1 < geometry.x1) {
            Master.exit("error: out of geometry bounds: x1: " + x1);
        }

        if (x2 > geometry.x2) {
            Master.exit("error: out of geometry bounds: x2: " + x2);
        }

        if (y1 < geometry.y1) {
            Master.exit("error: out of geometry bounds: y1: " + y1);
        }

        if (y2 > geometry.y2) {
            Master.exit("error: out of geometry bounds: y2: " + y2);
        }

        if (z1 < geometry.z1) {
            Master.exit("error: out of geometry bounds: z1: " + z1);
        }

        if (z2 > geometry.z2) {
            Master.exit("error: out of geometry bounds: z2: " + z2);
        }

        for (int i = 0; i < particles1.length; i++) {
            particles1[i] = null;
        }

        pcount = 0;

        for (DiffusantParticle p : d.particles) {

            if ((p == null) || !p.insideGeometry) {
                continue;
            }

            if (pcount >= particles1.length) {
                Master.exit("error: index out of bounds: vcount for dvs: " + pcount + " vs " + particles1.length);
            }

            if (roi_xy_width > 0) {

                if ((p.x < x1 - p.radius) || (p.x > x2 - p.radius)) {
                    continue; // out of bounds
                }

                if ((p.y < y1 - p.radius) || (p.y > y2 - p.radius)) {
                    continue; // out of bounds
                }

            }

            dz_top = p.z - z2;
            phiradians = getPhiRadians();

            if ((p.z > z2) && (p.z <= z2 + p.radius)) {
                dz = p.z - z2;
                particles1[pcount++] = new Particle(p, dz, dz_top, phiradians, 'a');
            } else if ((p.z >= z1) && (p.z <= z2)) {
                dz = 0;
                particles1[pcount++] = new Particle(p, dz, dz_top, phiradians, 'w');
            } else if ((p.z >= z1 - p.radius) && (p.z < z1)) {
                dz = z1 - p.z;
                particles1[pcount++] = new Particle(p, dz, dz_top, phiradians, 'b');
            }

        }

        if (pcount == 0) {
            return null;
        }

        for (Particle p1 : particles1) { // sum zintersections

            if (p1 == null) {
                continue;
            }

            p1.areaXY_intersect_sum = 0;
            p1.areaXY_intersect_count = 0;

            for (Particle p2 : particles1) {
                if ((p2 != null) && (p2 != p1)) {
                    p1.areaXY_intersection(p2, true);
                }
            }
            
            //Master.log("overlap count =" + p1.areaXY_intersect_count);

        }

        particles2 = new Particle[pcount];

        pcount = 0;

        for (Particle p : particles1) {
            if (p == null) {
                continue;
            }
            if (!Double.isNaN(p.r2D)) {
                particles2[pcount++] = p;
            }
        }

        return particles2;

    }

    public boolean computeZstack(double dz_scan, double zPSF, int particles_per_zstack_limit, int diameters_limit, boolean particles_completely_within_section) {

        int iparticle, i_total_2D, i_total_3D;
        double dv_xyz, area, den2D;
        double z_plane, z_plane_1, z_plane_2, z_offset, z_step_section = 0;
        double z_midpoint = Double.NaN;
        double z_min, z_max, dz, zeta_min_max;
        double norm_scale;
        double avg_d_total, avg_zscans_per_particle, avg_extra_zscans;
        int predicted_zscans;
        
        boolean finished;
        
        boolean limit_zscans_per_particle = false;
        //boolean limit_zscans_per_particle = true;

        boolean particle_limit_on = particles_per_zstack_limit > 0;

        double[][] histo;

        Particle[] particles;

        double print_scale = spaceScale(print_space_units);

        int histo_d_nbins = histo_d_nbins();
        int histo_a_nbins = histo_a_nbins();
        int histo_z_nbins = histo_z_nbins();

        double[] diam2D, diam3D, dz_top, phi_degrees;

        if (!project.monteCarlo.PBC) {
            Master.exit("error: projections should be computed with PBC to get correct density measurements");
        }

        if (update_all_D_3D_stats() == null) {
            return true;
        }

        area = section_area();
        String txt;

        Master.log("computing 2D projections... ");
        Master.log("# sections = " + Integer.toString(num_sections));
        txt = "xy geometry = " + formatter.format(section_xwidth()) + " x " + formatter.format(section_ywidth()) + " " + project.spaceUnits;
        Master.log(txt);
        Master.log("xy area = " + formatter.format(area) + " " + project.spaceUnits + "^2");
        Master.log("T = " + formatter.format(section_thickness * print_scale) + " " + print_space_units);
        Master.log("zeta = " + formatter.format(zeta_true() * print_scale) + " " + print_space_units);

        Master.log("phi = " + formatter.format(phi_true) + " ± " + formatter.format(phi_CV * phi_true) + "°");

        if (num_sections > 1) {
            z_step_section = geometry.zWidth / (num_sections + 4);
            Master.log("z-distance between sections = " + formatter.format(z_step_section * print_scale) + " " + print_space_units);
        }

        particle_stats = new SectionStats[num_sections];

        for (int i = 0; i < particle_stats.length; i++) {
            particle_stats[i] = new SectionStats();
        }

        if (compute_Fd) {
            Fd = new double[num_sections + 1][histo_d_nbins];
        }

        if (compute_Gd) {
            Gd = new double[num_sections + 1][histo_d_nbins];
        }

        if (compute_Pa) {
            Pa = new double[num_sections + 1][histo_a_nbins];
        }

        if (compute_Cz) {
            Cz = new double[num_sections + 1][histo_z_nbins];
        }

        if (num_sections > 1) {
            //z_offset = 1 * (section_thickness + D_3D_mean);
            z_offset = 2 * z_step_section; // avoid z-borders
        } else {
            z_offset = 0;
        }

        int not_completely_within_section;
        int num_zscans = 1 + (int) (section_thickness / dz_scan);
        int zscans_per_particle_max = 1 + (int) (D_3D_mean_true / dz_scan);
        double avg_not_completely_within_section = 0;

        Master.log("# scans/section = " + num_zscans + ", dz-scan = " + dz_scan);
        
        avg_d_total = 0;

        for (int isection = 0; isection < num_sections; isection++) {

            if (num_sections > 1) {
                z_plane = geometry.z1 + z_offset + isection * z_step_section;
            } else if (Double.isFinite(z_midpoint) && (z_midpoint > geometry.z1) && (z_midpoint < geometry.z2)) {
                z_plane = z_midpoint;
            } else {
                z_plane = geometry.zCenter;
            }

            //Master.log("z-plane = " + zPlaneRef);

            particles = getSectionParticles(z_plane); // grab all particles for zeta = T + D // sets z1,z2

            if (particles == null) {
                Master.exit("found no particles");
                return true;
            }

            iparticle = 0;
            not_completely_within_section = 0;

            for (Particle p : particles) { // exclude particles

                if (particle_limit_on && (iparticle >= particles_per_zstack_limit)) {
                     p.exclude = true;
                     continue;
                }

                if (!p.within) {
                    
                    if (particles_completely_within_section) {
                        p.exclude = true; // no caps allowed
                        particle_stats[isection].exclude_dmin++;
                        not_completely_within_section++;
                        continue;
                    }

                    if (p.r2D < p.r2D_min) {
                        p.exclude = true;
                        particle_stats[isection].exclude_dmin++;
                        continue;
                    }

                } else { // within section

                    if (particles_completely_within_section) {

                        if ((Math.abs(p.z - z1) < p.dz_max) || (Math.abs(p.z - z2) < p.dz_max)) {
                            p.exclude = true;
                            particle_stats[isection].exclude_dmin++;
                            not_completely_within_section++;
                            continue;
                        }

                    }

                }

                if (!p.exclude) {
                    iparticle++;
                }

            }

            avg_not_completely_within_section += not_completely_within_section;

            //Master.log("z-scan # particles = " + iparticle);

            i_total_3D = iparticle;

            if (particle_limit_on && (iparticle < particles_per_zstack_limit)) {
                txt = "error: section #" + isection;
                txt += ": not enough particles: " + iparticle + " out of " + particles_per_zstack_limit;
                Master.log(txt);
                continue;
            }

            diam2D = new double[particles.length * zscans_per_particle_max];

            for (int i = 0; i < diam2D.length; i++) {
                diam2D[i] = Double.NaN;
            }

            // now do z-scan and collect 2D diameters

            i_total_2D = 0;
            finished = false;

            for (int i_zscan = 0; i_zscan < num_zscans; i_zscan++) {

                z_plane = z1 + i_zscan * dz_scan;
                z_plane_1 = z_plane - 0.5 * zPSF;
                z_plane_2 = z_plane + 0.5 * zPSF;

                for (Particle p : particles) {

                    if (p.exclude) {
                        continue;
                    }

                    if ((p.z > z_plane_2) && (p.z <= z_plane_2 + p.r3D)) {
                        dz = p.z - z_plane_2;
                    } else if ((p.z >= z_plane_1) && (p.z <= z_plane_2)) {
                        dz = 0;
                    } else if ((p.z >= z_plane_1 - p.r3D) && (p.z < z_plane_1)) {
                        dz = z_plane_1 - p.z;
                    } else {
                        continue; // out of z-range
                    }
                    
                    p.r2D = Math.sqrt(p.r3D * p.r3D - dz * dz);

                    if (p.r2D < p.r2D_min) {
                        continue; // lost  cap
                    }

                    if (limit_zscans_per_particle) {
                        
                        predicted_zscans = (int) Math.floor((2 * p.dz_max) / dz_scan);

                        if (p.n_diameters >= predicted_zscans) {
                            continue;
                        }
                        
                    }

                    diam2D[i_total_2D] = 2.0 * p.r2D;
                    p.n_diameters++;
                    i_total_2D++;
                    
                    if ((diameters_limit > 0) && (i_total_2D >= diameters_limit)) {
                        finished = true;
                        //Master.log("finished scanning at scan # " + i_zscan);
                        break;
                    }

                }
                
                if (finished) {
                    break;
                }

            }
            
            avg_d_total += i_total_2D;
            
            avg_zscans_per_particle = 0;
            avg_extra_zscans = 0;
            
            for (Particle p : particles) {

                if (p.exclude) {
                    continue;
                }

                avg_zscans_per_particle += p.n_diameters;
                predicted_zscans = (int) Math.floor((2 * p.dz_max)/ dz_scan);
                //Master.log("delta # scans =" + (predicted_zscans - plu.n_diameters));
                //Master.log(" " + predicted_zscans + ", " + plu.n_diameters);
                avg_extra_zscans += p.n_diameters - predicted_zscans;

            }
            
            avg_zscans_per_particle /= i_total_3D * 1.0;
            avg_extra_zscans /= i_total_3D * 1.0;
            //Master.log("avg # scans per particle = " + avg_zscans_per_particle);
            //Master.log("avg extra # scans per particle = " + avg_extra_zscans);
            
            //Master.log("last z-plane = " + zPlaneRef + ", " + z1 + " - " + z2);

            diam3D = new double[particles.length];
            dz_top = new double[particles.length];
            phi_degrees = new double[particles.length];

            iparticle = 0;
            z_min = Double.POSITIVE_INFINITY;
            z_max = Double.NEGATIVE_INFINITY;

            for (Particle p : particles) {

                phi_degrees[iparticle] = p.phi_radians * 180.0 / Math.PI;

                if (p.exclude) {
                    diam3D[iparticle] = Double.NaN;
                    dz_top[iparticle] = Double.NaN;
                    iparticle++;
                    continue;
                }

                diam3D[iparticle] = p.r3D * 2.0;
                dz_top[iparticle] = p.dz_top;
                z_min = Math.min(z_min, p.z);
                z_max = Math.max(z_max, p.z);

                if (p.above) {
                    particle_stats[isection].above++;
                } else if (p.within) {
                    particle_stats[isection].within++;
                } else if (p.below) {
                    particle_stats[isection].below++;
                }

                iparticle++;

            }

            zeta_min_max = Math.abs(z_max - z_min);
            den2D = i_total_2D / (area * num_zscans);

            particle_stats[isection].N3D = i_total_3D;
            particle_stats[isection].N2D = i_total_2D;

            //particle_stats[isection].den3D = i_total_3D / (area * zeta);
            particle_stats[isection].den3D = den3D_true;
            particle_stats[isection].den2D = den2D;

            particle_stats[isection].exclude_dmin *= 100 / (i_total_3D * 1.0);
            particle_stats[isection].exclude_imax *= 100 / (i_total_3D * 1.0);
            particle_stats[isection].above *= 100 / (i_total_3D * 1.0);
            particle_stats[isection].within *= 100 / (i_total_3D * 1.0);
            particle_stats[isection].below *= 100 / (i_total_3D * 1.0);
            particle_stats[isection].d_3D_mean = Utility.average(diam3D);
            particle_stats[isection].d_3D_stdv = Utility.stdev(diam3D);
            particle_stats[isection].d_3D_min = Utility.getMin(diam3D);
            particle_stats[isection].d_3D_max = Utility.getMax(diam3D);
            particle_stats[isection].d_2D_mean = Utility.average(diam2D);
            particle_stats[isection].d_2D_stdv = Utility.stdev(diam2D);
            particle_stats[isection].d_2D_min = Utility.getMin(diam2D);
            particle_stats[isection].d_2D_max = Utility.getMax(diam2D);
            particle_stats[isection].phi_mean = Utility.average(phi_degrees);
            particle_stats[isection].phi_stdv = Utility.stdev(phi_degrees);
            particle_stats[isection].phi_min = Utility.getMin(phi_degrees);
            particle_stats[isection].phi_max = Utility.getMax(phi_degrees);
            particle_stats[isection].phi_eq2 = phi_eq2(den2D, zPSF);
            particle_stats[isection].zeta_min_max = zeta_min_max;

            if (Fd != null) {

                histo = Utility.histogram(diam3D, histo_d_min, histo_d_max, histo_d_bwidth, histo_d_norm);

                for (int ih = 0; ih < histo.length; ih++) {
                    Fd[0][ih] = histo[ih][0];
                    Fd[isection + 1][ih] = histo[ih][1];
                }

            }

            if (Gd != null) {

                histo = Utility.histogram(diam2D, histo_d_min, histo_d_max, histo_d_bwidth, histo_d_norm);

                for (int ih = 0; ih < histo.length; ih++) {
                    Gd[0][ih] = histo[ih][0];
                    Gd[isection + 1][ih] = histo[ih][1];
                }

            }

            if (Pa != null) {

                histo = Utility.histogram(phi_degrees, histo_a_min, histo_a_max, histo_a_bwidth, histo_a_norm);

                for (int ih = 0; ih < histo.length; ih++) {
                    Pa[0][ih] = histo[ih][0];
                    Pa[isection + 1][ih] = histo[ih][1];
                }

            }

            if (Cz != null) {

                histo = Utility.histogram(dz_top, histo_z_min, histo_z_max, histo_z_bwidth, 0);

                dv_xyz = area * histo_z_bwidth;

                switch (histo_z_norm) {
                    case 0: // count
                        norm_scale = 1;
                        break;
                    case 1: // density, particles/um^3
                        norm_scale = 1 / dv_xyz;
                        break;
                    //case 2: // volume fraction // using dvs.volumeMean is only approximate
                        //norm_scale = dvs.volumeMean / dv_xyz;
                        //break;
                    default:
                        return true;
                }

                for (int ih = 0; ih < histo.length; ih++) {
                    Cz[0][ih] = histo[ih][0];
                    Cz[isection + 1][ih] = histo[ih][1] * norm_scale;
                }

            }

        }
        
        avg_d_total /= num_sections * 1.0;
        Master.log("# 2D diameters/section = " + avg_d_total);

        if (particles_completely_within_section) {
            Master.log("avg # particles not completely within section = " + (avg_not_completely_within_section / (num_sections * 1.0)));
        }

        return false;

    }

    public class Particle { // inner class

        public DiffusantParticle id = null;
        public double r3D = 0; // radius 3D // um
        public double r2D = 0; // projection radius 2D // um
        public double r2D_min = 0;

        public double x, y, z; // center location // um
        public double dz = 0; // distance from section above or below, 0 for within // um
        public double dz_top = Double.NaN; // distance from top of section
        public double dz_max = Double.NaN; // max dz above/below section
        public double phi_radians = Double.NaN; // cap angle limit // radians
        public int n_diameters = 0;

        public double areaXY_intersect_sum = Double.NaN; // um^2
        public int areaXY_intersect_count = 0;

        public boolean above = false; // above section
        public boolean within = false; // within section
        public boolean below = false; // below section
        public boolean exclude = false; // exclude from analysis

        public Particle(DiffusantParticle p, double DZ, double DZ_TOP, double PHI_RADIANS, char awb) {

            id = p;
            r3D = p.radius;
            r2D_min = r3D * Math.sin(PHI_RADIANS);
            x = p.x;
            y = p.y;
            z = p.z;
            dz = DZ;
            dz_top = DZ_TOP;
            dz_max = r3D * Math.cos(PHI_RADIANS);
            phi_radians = PHI_RADIANS;

            if (dz == 0) {
                r2D = r3D;
            } else {
                r2D = Math.sqrt(r3D * r3D - dz * dz);
            }

            switch (awb) {
                case 'a':
                    above = true;
                    break;
                case 'w':
                    within = true;
                    break;
                case 'b':
                    below = true;
                    break;
                default:
                    Master.exit("error: unknown AWB type: " + awb);
            }

        }

        public double areaXY_intersect_norm() {
            return areaXY_intersect_sum / (Math.PI * r2D * r2D);
        }

        public double areaXY_intersection(Particle p_other, boolean save) {
            double d2, r2D2, vr2D2, t1, t2, t3, intersection;

            double dx, dy, d, alpha;

            if (p_other.z <= z) {
                return 0; // other particle needs to be above this particle
            }

            dx = x - p_other.x;
            dy = y - p_other.y;
            d = Math.sqrt(dx * dx + dy * dy); // distance between particles
            
            alpha = (1 / (2 * 4 * r2D * p_other.r2D)) * (4 * d * d - 4 * r2D * r2D - 4 * p_other.r2D * p_other.r2D); // Hilliard

            if (d >= r2D + p_other.r2D) {
                return 0; // no overlap
            }
            
            if (d == 0) { // centered, unlikely
                if (r2D <= p_other.r2D) {
                    intersection = Math.PI * r2D * r2D;
                } else {
                    intersection = Math.PI * p_other.r2D * p_other.r2D;
                }
                if (save) {
                    areaXY_intersect_sum += intersection;
                    areaXY_intersect_count++;
                }
                return intersection;
            }

            if (d < Math.abs(r2D - p_other.r2D)) {
                // alpha < -1
                if (r2D < p_other.r2D) {
                    //Master.log("this cirle is within another circle: r < R " + a);
                    intersection = Math.PI * r2D * r2D;
                } else { // p_other.r2D < r2D
                    //Master.log("another cirle is within this circle: R < r " + a);
                    intersection = Math.PI * p_other.r2D * p_other.r2D;
                }
                if (save) {
                    areaXY_intersect_sum += intersection;
                    areaXY_intersect_count++;
                }
                return intersection;
            }

            // otherwise there is partial overlap...
            // https://mathworld.wolfram.com/Circle-CircleIntersection.html
            
            d2 = d * d;
            r2D2 = r2D * r2D;
            vr2D2 = p_other.r2D * p_other.r2D;

            t1 = (d2 + r2D2 - vr2D2) / (2 * d * r2D);
            t1 = r2D2 * Math.acos(t1);

            if (Math.abs(t1) > 1) {
                Master.exit("bad t1: greater than 1: " + t1);
            }

            t2 = (d2 + vr2D2 - r2D2) / (2 * d * p_other.r2D);
            t2 = vr2D2 * Math.acos(t2);

            if (Math.abs(t2) > 1) {
                Master.exit("bad t2: greater than 1: " + t2);
            }

            t3 = (-d + r2D + p_other.r2D);
            t3 *= (d + r2D - p_other.r2D);
            t3 *= (d - r2D + p_other.r2D);
            t3 *= (d + r2D + p_other.r2D);
            t3 = 0.5 * Math.sqrt(t3);

            intersection = t1 + t2 - t3;
            
            if (save) {
                areaXY_intersect_sum += intersection;
                areaXY_intersect_count++;
            }
            
            //if ((alpha > -1) && (alpha < 0)) {
            //    Master.log("overlap alpha =" + alpha);
            //}

            return intersection;

        }
        
        public double areaXY_intersection_area2add(Particle p_other) {
            double d2, r2D2, vr2D2, t1, t2, t3, intersection;
            double area_this, area_other;
            double dx, dy, d, alpha;

            if (p_other.z > z) {
                return 0; // other particle needs to be below this particle
            }
            
            area_this = Math.PI * r2D * r2D;
            area_other = Math.PI * p_other.r2D * p_other.r2D;
            dx = x - p_other.x;
            dy = y - p_other.y;
            d = Math.sqrt(dx * dx + dy * dy); // distance between particles
            
            alpha = (1 / (2 * 4 * r2D * p_other.r2D)) * (4 * d * d - 4 * r2D * r2D - 4 * p_other.r2D * p_other.r2D); // Hilliard

            if (d >= r2D + p_other.r2D) {
                return 0; // no overlap
            }

            if (d < Math.abs(r2D - p_other.r2D)) {
                // alpha < -1
                if (r2D < p_other.r2D) {
                    //Master.log("this cirle is within another circle: r < R " + alpha + ", " + (area_other - area_this));
                    return area_other - area_this;
                } else { // p_other.r2D < r2D
                    //Master.log("another cirle is within this circle: R < r " + alpha);
                    return Double.NaN; // flag to remove this particle
                }
            }

            // otherwise there is partial overlap...
            // https://mathworld.wolfram.com/Circle-CircleIntersection.html
            
            d2 = d * d;
            r2D2 = r2D * r2D;
            vr2D2 = p_other.r2D * p_other.r2D;

            t1 = (d2 + r2D2 - vr2D2) / (2 * d * r2D);
            t1 = r2D2 * Math.acos(t1);

            if (Math.abs(t1) > 1) {
                Master.exit("bad t1: greater than 1: " + t1);
            }

            t2 = (d2 + vr2D2 - r2D2) / (2 * d * p_other.r2D);
            t2 = vr2D2 * Math.acos(t2);

            if (Math.abs(t2) > 1) {
                Master.exit("bad t2: greater than 1: " + t2);
            }

            t3 = (-d + r2D + p_other.r2D);
            t3 *= (d + r2D - p_other.r2D);
            t3 *= (d - r2D + p_other.r2D);
            t3 *= (d + r2D + p_other.r2D);
            t3 = 0.5 * Math.sqrt(t3);

            intersection = t1 + t2 - t3;
            
            if ((alpha > -1) && (alpha < 0)) {
            //    Master.log("overlap alpha =" + alpha);
                return area_other - intersection;
            }

            return 0;

        }

    }

    public class SectionStats { // inner class
        double N3D = 0, N2D = 0;
        double den3D = 0, den2D = 0;
        double above = 0, within = 0, below = 0;
        double exclude_dmin = 0; // excluded from min diameter limit
        double exclude_imax = 0; // excluded from max intersect limit
        double d_3D_mean = 0, d_3D_stdv = 0;
        double d_3D_min = 0, d_3D_max = 0;
        double d_2D_mean = 0, d_2D_stdv = 0;
        double d_2D_min = 0, d_2D_max = 0;
        double AF = 0;
        double phi_mean = 0, phi_stdv = 0;
        double phi_min = 0, phi_max = 0;
        double phi_eq2 = 0;
        double zeta_min_max = 0;
    }

    public boolean save_section_stats(String filename, double diamscalefactor) {
        // diamscalefactor: 1 for um, 1000 for nm

        String str = "", fullName = "";

        //boolean columnNames = false;
        boolean columnNames = true; // for "all"

        try {

            fullName = directory + filename + ".dat";

            bw = new BufferedWriter(new FileWriter(fullName));

            if (columnNames) {

                str = "N3D" + "\t";
                str += "N2D" + "\t";
                str += "den3D" + "\t";
                str += "den2D" + "\t";
                str += "above%" + "\t";
                str += "within%" + "\t";
                str += "below%" + "\t";
                str += "exclude_dmin%" + "\t";
                str += "exclude_imax%" + "\t";
                str += "d_3D_mean" + "\t";
                str += "d_3D_stdv" + "\t";
                str += "d_3D_min" + "\t";
                str += "d_3D_max" + "\t";
                str += "d_2D_mean" + "\t";
                str += "d_2D_stdv" + "\t";
                str += "d_2D_min" + "\t";
                str += "d_2D_max" + "\t";
                str += "AF" + "\t";

                if (sections_stats_save_phi) {
                    str += "phi_mean" + "\t";
                    str += "phi_stdv" + "\t";
                    str += "phi_min" + "\t";
                    str += "phi_max" + "\t";
                    str += "phi_eq2" + "\t";
                }
                
                str += "zeta_min_max";

                bw.write(str);
                bw.newLine();

            }

            for (SectionStats s : particle_stats) {

                str = s.N3D + "\t";
                str += s.N2D + "\t";
                str += s.den3D + "\t";
                str += s.den2D + "\t";
                str += s.above + "\t";
                str += s.within + "\t";
                str += s.below + "\t";
                str += s.exclude_dmin + "\t";
                str += s.exclude_imax + "\t";
                str += s.d_3D_mean * diamscalefactor + "\t";
                str += s.d_3D_stdv * diamscalefactor + "\t";
                str += s.d_3D_min * diamscalefactor + "\t";
                str += s.d_3D_max * diamscalefactor + "\t";
                str += s.d_2D_mean * diamscalefactor + "\t";
                str += s.d_2D_stdv * diamscalefactor + "\t";
                str += s.d_2D_min * diamscalefactor + "\t";
                str += s.d_2D_max * diamscalefactor + "\t";
                str += s.AF + "\t";

                if (sections_stats_save_phi) {
                    str += s.phi_mean + "\t";
                    str += s.phi_stdv + "\t";
                    str += s.phi_min + "\t";
                    str += s.phi_max + "\t";
                    str += s.phi_eq2 + "\t";
                }
                
                str += s.zeta_min_max;

                bw.write(str);
                bw.newLine();

            }

            bw.close();
            bw = null;

            Master.log("created output file: " + fullName);

            return false;

        } catch (IOException e) {

            Master.log("cannot write to file: " + fullName);

            return true;

        }

    }

    public SectionStats[] section_stats_avg() {
        
        int i;

        if ((particle_stats == null) || (particle_stats.length == 0)) {
            return null;
        }
        
        int numsections = particle_stats.length;
        
        double[] temp = new double[numsections];

        SectionStats[] avg = new SectionStats[2];
        avg[0] = new SectionStats();
        avg[1] = new SectionStats();
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.N3D;
        }
        
        avg[0].N3D = Utility.average(temp);
        avg[1].N3D = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.N2D;
        }
        
        avg[0].N2D = Utility.average(temp);
        avg[1].N2D = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.den3D;
        }
        
        avg[0].den3D = Utility.average(temp);
        avg[1].den3D = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.den2D;
        }
        
        avg[0].den2D = Utility.average(temp);
        avg[1].den2D = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.above;
        }
        
        avg[0].above = Utility.average(temp);
        avg[1].above = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.within;
        }
        
        avg[0].within = Utility.average(temp);
        avg[1].within = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.below;
        }
        
        avg[0].below = Utility.average(temp);
        avg[1].below = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.exclude_dmin;
        }
        
        avg[0].exclude_dmin = Utility.average(temp);
        avg[1].exclude_dmin = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.exclude_imax;
        }
        
        avg[0].exclude_imax = Utility.average(temp);
        avg[1].exclude_imax = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_3D_mean;
        }
        
        avg[0].d_3D_mean = Utility.average(temp);
        avg[1].d_3D_mean = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_3D_stdv;
        }
        
        avg[0].d_3D_stdv = Utility.average(temp);
        avg[1].d_3D_stdv = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_3D_min;
        }
        
        avg[0].d_3D_min = Utility.average(temp);
        avg[1].d_3D_min = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_3D_max;
        }
        
        avg[0].d_3D_max = Utility.average(temp);
        avg[1].d_3D_max = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_2D_mean;
        }
        
        avg[0].d_2D_mean = Utility.average(temp);
        avg[1].d_2D_mean = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_2D_stdv;
        }
        
        avg[0].d_2D_stdv = Utility.average(temp);
        avg[1].d_2D_stdv = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_2D_min;
        }
        
        avg[0].d_2D_min = Utility.average(temp);
        avg[1].d_2D_min = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.d_2D_max;
        }
        
        avg[0].d_2D_max = Utility.average(temp);
        avg[1].d_2D_max = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.AF;
        }
        
        avg[0].AF = Utility.average(temp);
        avg[1].AF = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.phi_mean;
        }
        
        avg[0].phi_mean = Utility.average(temp);
        avg[1].phi_mean = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.phi_stdv;
        }
        
        avg[0].phi_stdv = Utility.average(temp);
        avg[1].phi_stdv = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.phi_min;
        }
        
        avg[0].phi_min = Utility.average(temp);
        avg[1].phi_min = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.phi_max;
        }
        
        avg[0].phi_max = Utility.average(temp);
        avg[1].phi_max = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.phi_eq2;
        }
        
        avg[0].phi_eq2 = Utility.average(temp);
        avg[1].phi_eq2 = Utility.stdev(temp);
        
        i = 0;
        
        for (SectionStats p : particle_stats) {
            temp[i++] = p.zeta_min_max;
        }
        
        avg[0].zeta_min_max = Utility.average(temp);
        avg[1].zeta_min_max = Utility.stdev(temp);

        return avg;

    }

    public void print_section_stats_avg() {

        double dtemp, dtemp2, xwidth, ywidth, areaxy;
        double zcap, zeta_true, dv_xyz;
        String str, tab = "   ";

        SectionStats[] avg = section_stats_avg();

        if (avg == null) {
            return;
        }

        double scale = spaceScale(print_space_units);

        double phi_radians = phi_true * Math.PI / 180.0;
        double d_phi = D_3D_mean_true * Math.sin(phi_radians);

        //zcap = all_D_3D_mean * Math.cos(phi_radians);
        zcap = D_3D_mean_true * Math.cos(phi_radians);
        zeta_true = section_thickness + zcap;

        xwidth = section_xwidth();
        ywidth = section_ywidth();
        areaxy = section_area();

        dv_xyz = xwidth * ywidth * zeta_true;

        Master.log("avg 2D diameter stats (n = " + num_sections + " sections)[expected][error%]...");

        dtemp = den3D_true * dv_xyz;
        dtemp2 = 100 * (avg[0].N3D - dtemp) / dtemp;
        str = tab + "particles/section = ";
        str += formatter.format(avg[0].N3D);
        str += " [" + formatter.format(dtemp) + "]";
        str += " [" + formatter.format(dtemp2) + "%]";
        Master.log(str);

        dtemp = den3D_true * zeta_true;
        dtemp2 = 100 * (avg[0].den2D - dtemp) / dtemp;
        str = tab + "particles/um^2 = ";
        str += formatter.format(avg[0].den2D);
        str += " [" + formatter.format(dtemp) + "]";
        str += " [" + formatter.format(dtemp2) + "%]";
        Master.log(str);

        dtemp = avg[0].N3D / dv_xyz;
        dtemp2 = 100 * (dtemp - den3D_true) / den3D_true;
        str = tab + "particles/um^3 = ";
        str += formatter.format(dtemp);
        str += " [" + formatter.format(den3D_true) + "]";
        str += " [" + formatter.format(dtemp2) + "%]";
        Master.log(str);

        Master.log(tab + "% excluded (Imax) = " + formatter.format(avg[0].exclude_imax));
        Master.log(tab + "% excluded (Dmin) = " + formatter.format(avg[0].exclude_dmin));

        str = tab + "% above = " + formatter.format(avg[0].above);
        dtemp = 100 * 0.5 * zcap / (section_thickness + zcap);
        dtemp2 = 100 * (avg[0].above - dtemp) / dtemp;
        str += " [" + formatter.format(dtemp) + "]";
        str += " [" + formatter.format(dtemp2) + "%]";
        Master.log(str);

        str = tab + "% within = " + formatter.format(avg[0].within);
        dtemp = 100 * section_thickness / (section_thickness + zcap);
        dtemp2 = 100 * (avg[0].within - dtemp) / dtemp;
        str += " [" + formatter.format(dtemp) + "]";
        str += " [" + formatter.format(dtemp2) + "%]";
        Master.log(str);

        str = tab + "% below = " + formatter.format(avg[0].below);
        dtemp = 100 * 0.5 * zcap / (section_thickness + zcap);
        dtemp2 = 100 * (avg[0].below - dtemp) / dtemp;
        str += " [" + formatter.format(dtemp) + "]";
        str += " [" + formatter.format(dtemp2) + "%]";
        Master.log(str);
        
        str = tab + "min 2D diameter (" + print_space_units + ") = ";
        str += formatter.format(avg[0].d_2D_min * scale);
        str += " [d-phi=" + formatter.format(d_phi * scale) + "]";
        Master.log(str);

        str = tab + "2D diameter (" + print_space_units + ") = ";
        str += formatter.format(avg[0].d_2D_mean * scale);
        str += " ± " + formatter.format(avg[0].d_2D_stdv * scale);

        if (phi_radians == 0) {
            dtemp = Goldsmith_d_2D_mean(section_thickness, D_3D_mean_true, D_3D_stdv_true);
            str += " [" + formatter.format(dtemp * scale);
            dtemp2 = Goldsmith_d_2D_stdv(section_thickness, D_3D_mean_true, D_3D_stdv_true);
            str +=  " ± " + formatter.format(dtemp2 * scale) + "]";
        }

        Master.log(str);

        str = tab + "3D diameter (" + print_space_units + ") = ";
        str += formatter.format(avg[0].d_3D_mean * scale);
        str += " ± " + formatter.format(avg[0].d_3D_stdv * scale);
        str += " [" + formatter.format(D_3D_mean_true * scale);
        str += " ± " + formatter.format(D_3D_stdv_true * scale) + "]";
        dtemp = 100 * (avg[0].d_3D_mean - D_3D_mean_true) / D_3D_mean_true;
        dtemp2 = 100 * (avg[0].d_3D_stdv - D_3D_stdv_true) / D_3D_stdv_true;
        str += " [" + formatter.format(dtemp);
        str += " ± " + formatter.format(dtemp2) + "]";

        Master.log(str);
        
        //str = tab + "phi Eq2 = ";
        //str += formatter.format(avg[0].phi_eq2);
        //str += " ± " + formatter.format(avg[1].phi_eq2);
        
        //Master.log(str);
        
        double kv = kv_Weibel(section_thickness, D_3D_mean_true, D_3D_stdv_true, phi_true);

        str = tab + "AF = ";
        str += formatter.format(avg[0].AF);
        str += " ± " + formatter.format(avg[1].AF);
        str += " [" + formatter.format(d.setVolumeFraction/kv) + "]";
        
        Master.log(str);

    }

    public double[][] Fd_total() { // j-bins, j=2

        double[] D_3D = d.getParameter("diameter");

        if (D_3D == null) {
            return null;
        }

        double[][] histo = Utility.histogram(D_3D, histo_d_min, histo_d_max, histo_d_bwidth, histo_d_norm);

        return histo;

    }

    public double[][] histo_stats(double[][] histo, boolean SEM) { // j-bins, j=3
        // histo: j-sections, j-bins
        double dtemp;

        if (histo == null) {
            return null;
        }

        double nsections = histo.length - 1; // first column is x-bins
        int nbins = histo[0].length;

        double[][] histo_avg = new double[nbins][3];
            // [0] xscale
            // [1] mean
            // [2] stdv

        for (int j = 0; j < histo[0].length; j++) {

            histo_avg[j][0] = histo[0][j]; // x-print_scale

            for (int i = 1; i < histo.length; i++) {
                dtemp = histo[i][j];
                histo_avg[j][1] += dtemp; // sum
                histo_avg[j][2] += dtemp * dtemp; // sum of squares
            }

        }

        for (double[] h : histo_avg) {
            if (nsections > 1) {
                h[2] = Math.sqrt((h[2] - h[1] * h[1] / nsections) / (nsections - 1.0)); // stdv
                if (SEM) {
                    h[2] /= Math.sqrt(nsections);
                }
            } else {
                h[2] = 0;
            }
            h[1] /= nsections; // mean
        }

        return histo_avg;

    }

    public double[][] Fd_stats(boolean SEM) { // j-bins, j=3
        return histo_stats(Fd, SEM);
    }

    public double[][] Gd_stats(boolean SEM) { // j-bins, j=3
        return histo_stats(Gd, SEM);
    }

    public double[][] Pa_stats(boolean SEM) { // j-bins, j=3
        return histo_stats(Pa, SEM);
    }

    public double[][] Cz_stats(boolean SEM) { // j-bins, j=3
        return histo_stats(Cz, SEM);
    }

    public double Cz_zeta(double[][] stats) {
        // stats: j-bins, j=3
        double zeta, area_sum = 0;
        String str;

        if (stats[0].length != 3) {
            Master.log("length of stats j-column is not 3");
            return Double.NaN;
        }

        double scale = spaceScale(print_space_units);

        double bwidth = stats[1][0] - stats[0][0];

        for (double[] h : stats) {
            area_sum += bwidth * h[1];
        }

        switch (histo_z_norm) {
            case 1:
                zeta = area_sum / den3D_true;
                break;
            case 2:
                zeta = area_sum / (den3D_true * d.volumeMean);
                break;
            //case 3: using dvs.volumeMean is only approximate
                //zeta = area_sum * dvs.volumeMean / (den3D_true * dvs.volumeMean);
                //break;
            default:
                return Double.NaN;
        }

        str = "equivalent zeta (" + print_space_units + ") = ";
        str += formatter.format(zeta * scale) + " ";
        str += " [" + formatter.format(zeta_true() * scale) + "]";
        Master.log(str);

        return zeta;

    }

    public double[][] Iz_stats(boolean SEM) { // j-bins, j=3
        return histo_stats(Iz, SEM);
    }

    public boolean save_histo(String filename, double[][] histo, String histo_type, double xscalefactor) {
        // histo: // j-bins, j=2
        //      j=0: xscale
        //      j=1: ydata
        // histo_type: "d"-diameter, "z"-zscale
        double xvalue, yvalue;
        String str, fullName = "";

        boolean d_probability = histo_type.equalsIgnoreCase("d") && (histo_d_norm == 1);

        if (histo == null) {
            return true;
        }

        if (histo[0].length != 2) {
            Master.log("length of stats j-column is not 2");
            return true;
        }

        if (!Double.isFinite(xscalefactor)) {
            xscalefactor = 1;
        }

        try {

            fullName = directory + filename + ".dat";

            bw = new BufferedWriter(new FileWriter(fullName));

            for (double[] h : histo) {

                xvalue = h[0];
                yvalue = h[1];

                if (xscalefactor != 1) {
                    xvalue *= xscalefactor;
                    if (d_probability) {
                        yvalue /= xscalefactor; // probability is a function of xbin
                    }
                }

                if (histo_save_xbins && (xscalefactor != 1)) {
                    xvalue *= xscalefactor;
                }

                if (histo_save_xbins) {
                    str = formatter.format(xvalue) + "\t";
                } else {
                    str = "";
                }

                str += formatter.format(yvalue);

                //Master.log(str);
                bw.write(str);
                bw.newLine();

            }

            bw.close();
            bw = null;

            Master.log("created output file: " + fullName);

            return false;

        } catch (IOException e) {

            Master.log("cannot write to file: " + fullName);

            return true;

        }

    }

    public boolean save_Fd_total(String filename, double diamscalefactor) {
        return save_histo(filename, Fd_total(), "d", diamscalefactor);
    }

    public boolean save_section_histos(String filename, double[][] histo, String histo_type, double xscalefactor) {
        // histo: j-sections, j-bins
        // histo_type: "d"-diameter, "z"-zscale
        int ibgn;
        double value;
        String str, fullName = "";

        if (!Double.isFinite(xscalefactor)) {
            return true;
        }

        if (histo == null) {
            return true;
        }

        int nsections = histo.length; // j
        int nbins = histo[0].length; // j

        boolean d_probability = histo_type.equalsIgnoreCase("d") && (histo_d_norm == 1);

        Master.log("saving " + histo_type + "-histos to disk (i=" + (nsections - 1) + " sections, j=" + nbins + " histo bins):");

        if (histo_save_xbins) {
            ibgn = 0;
        } else {
            ibgn = 1;
        }

        try {

            fullName = directory + filename + ".dat";

            bw = new BufferedWriter(new FileWriter(fullName));

            for (int j = 0; j < nbins; j++) {

                str = "";

                for (int i = ibgn; i < nsections; i++) {

                    value = histo[i][j];

                    if (xscalefactor != 1) {
                        if (i == 0) { // xbin
                            value *= xscalefactor;
                        } else if (d_probability) {
                            value /= xscalefactor; // probability is a function of xbin
                        }
                    }

                    str += formatter.format(value);

                    if (i < histo.length - 1) {
                        str += "\t";
                    }

                }

                //Master.log(str);
                bw.write(str);
                bw.newLine();

            }

            bw.close();
            bw = null;

            Master.log("created " + histo_type + "-histo output file: " + fullName);

            return false;

        } catch (IOException e) {

            Master.log("cannot write to file: " + fullName);

            return true;

        }

    }

    public boolean save_Fd(String filename, double diamscalefactor) {
        return save_section_histos(filename, Fd, "d", diamscalefactor);
    }

    public boolean save_Gd(String filename, double diamscalefactor) {
        return save_section_histos(filename, Gd, "d", diamscalefactor);
    }

    public boolean save_Pa(String filename) {
        return save_section_histos(filename, Pa, "a", 1);
    }

    public boolean save_Cz(String filename) {
        return save_section_histos(filename, Cz, "z", 1);
    }

    public boolean save_Iz(String filename) {
        return save_section_histos(filename, Iz, "z", 1);
    }

    public boolean save_stats(String filename, double[][] stats, String histo_type, double xscalefactor, double yscalefactor) {
        // stats: // j-bins, j=3
        //      j=0: xscale
        //      j=1: mean
        //      j=2: stdv
        // histo_type: "d"-diameter, "z"-zscale
        double xvalue, mn, sd;
        String str, fullName = "";

        boolean d_probability = histo_type.equalsIgnoreCase("d") && (histo_d_norm == 1);

        if (stats == null) {
            return true;
        }

        if (stats.length == 0) {
            return true;
        }

        if (stats[0].length != 3) {
            Master.log("length of stats j-column is not 3");
            return true;
        }

        if (!Double.isFinite(xscalefactor)) {
            xscalefactor = 1;
        }

        if (!Double.isFinite(yscalefactor)) {
            yscalefactor = 1;
        }

        try {

            fullName = directory + filename + ".dat";

            bw = new BufferedWriter(new FileWriter(fullName));

            for (double[] h : stats) {

                xvalue = h[0];
                mn = h[1];
                sd = h[2];

                if (xscalefactor != 1) {
                    xvalue *= xscalefactor;
                    if (d_probability) {
                        mn /= xscalefactor; // probability is a function of xbin
                        sd /= xscalefactor; // probability is a function of xbin
                    }
                }

                if (histo_save_xbins && (xscalefactor != 1)) {
                    xvalue *= xscalefactor;
                }

                if (yscalefactor != 1) {
                    mn *= yscalefactor;
                    sd *= yscalefactor;
                }

                if (histo_save_xbins) {
                    str = formatter.format(xvalue) + "\t";
                } else {
                    str = "";
                }

                str += formatter.format(mn) + "\t";
                str += formatter.format(sd);

                //Master.log(str);
                bw.write(str);
                bw.newLine();

            }

            bw.close();
            bw = null;

            Master.log("created output file: " + fullName);

            return false;

        } catch (IOException e) {

            Master.log("cannot write to file: " + fullName);

            return true;

        }

    }

    public boolean save_Fd_stats(String filename, double diamscalefactor, boolean SEM) {
        return save_stats(filename, Fd_stats(SEM), "d", diamscalefactor, 1);
    }

    public boolean save_Gd_stats(String filename, double diamscalefactor, boolean SEM) {
        return save_stats(filename, Gd_stats(SEM), "d", diamscalefactor, 1);
    }

    public boolean save_Pa_stats(String filename, boolean SEM) {
        return save_stats(filename, Pa_stats(SEM), "a", 1, 1);
    }

    public boolean save_Cz_stats(String filename, boolean SEM) {
        return save_stats(filename, Cz_stats(SEM), "z", 1, 1);
    }

    public boolean save_Iz_stats(String filename, boolean SEM) {
        return save_stats(filename, Iz_stats(SEM), "z", 1, 1);
    }

    public static double Goldsmith_d_2D_mean(double section_T, double diam_3D_mean, double diam_3D_stdv) {

        double m1 = diam_3D_mean;
        double m2 = diam_3D_mean * diam_3D_mean + diam_3D_stdv * diam_3D_stdv;

        double o1 = (section_T * m1 + 0.25 * Math.PI * m2) / (section_T + m1);

        return o1;

    }

    public static double Goldsmith_d_2D_stdv(double section_T, double diam_3D_mean, double diam_3D_stdv) {

        double m1 = diam_3D_mean;
        double m2 = diam_3D_mean * diam_3D_mean + diam_3D_stdv * diam_3D_stdv;
        double m3 = diam_3D_mean * diam_3D_mean * diam_3D_mean + 3 * diam_3D_mean * diam_3D_stdv * diam_3D_stdv;

        double o1 = (section_T * m1 + 0.25 * Math.PI * m2) / (section_T + m1);
        double o2 = (section_T * m2 + 2 * m3 / 3) / (section_T + m1);

        double stdv = Math.sqrt(o2 - o1 * o1);

        return stdv;

    }

    public static double[][] Goldsmith_G_to_F_Transform(double[][] g, double section_T, double diam_mean) {

        int nbins = g.length;

        double[][] h_array = Goldsmith_H(g, section_T);
        double[][] f = new double[nbins][2];

        double sum_g = 0, sum_h = 0;

        double binsize = (g[1][0] - g[0][0]) / 1000.0; // bin size

        double tmh = (section_T + diam_mean) / binsize;

        Master.log("T = " + section_T + ", D = " + diam_mean + ", binsize = " + binsize);

        for (int i = 0; i < nbins; i++) {
            sum_g += g[i][1];
            sum_h += h_array[i][1];
        }

        Master.log("tmh = " + tmh + ", sum_g/sum_h = " + (sum_g / sum_h));

        for (int i = 0; i < nbins; i++) {
            //f[j][1] = (sum_g/sum_h) * h_array[j][1];
            f[i][1] = tmh * h_array[i][1];
            f[i][0] = g[i][0];
        }

        return f;

    }

    private static double[][] Goldsmith_H(double[][] g, double section_T) {

        int i_index, j_index;
        int nbins = g.length;
        double sum;

        double[][] h_array = new double[nbins][2];

        double h = g[1][0] - g[0][0] / 1000.0; // bin size

        for (int i = nbins; i >= 1; i--) {
            i_index = i - 1;
            sum = g[i_index][1];
            for (int j = i + 1; j <= nbins; j++) {
                j_index = j - 1;
                sum -= Goldsmith_Aij(i, j, section_T, h) * h_array[j_index][1];
            }
            h_array[i_index][1] = sum / Goldsmith_Aij(i, i, section_T, h);
            h_array[i_index][0] = g[i_index][0];
            //Master.log("j=" + j + ", h=" + h_array[ih][1]);
        }

        return h_array;

    }

    private static double Goldsmith_Aij(int i, int j, double section_T, double h) {

        double ij0, ij1;

        if (j == i) {
            return section_T / h + Math.sqrt(i - 3.0 / 4.0);
        } else {
            ij0 = Math.sqrt((j - 0.5) * (j - 0.5) - (i - 1.0) * (i - 1.0));
            ij1 = Math.sqrt((j - 0.5) * (j - 0.5) - i * i);
            return ij0 - ij1;
        }

    }

    public void saveFileName() {

        if (save != null) {
            save.fileName(name, "");
        }

    }

    public void saveDimensions() {

        if ((save != null) && save.autoDimensions) {

            //save.xdim = project.timeUnits;
            if ((name == null) || (name.length() == 0)) {
                //save.ydim = project.concUnits;
            } else {
                //save.ydim = name + " (" + project.concUnits + ")";
            }

        }

    }

    public boolean saveInit() {

        int dataPoints = 1;

        if (save == null) {
            return false;
        }

        return save.init(name, null, -1, dataPoints);

    }

    public boolean save() { // called from simulation "run" function

        if (save == null) {
            return false;
        }

        return save.saveData(Double.NaN);

    }

    public boolean saveFinish() {

        if (save == null) {
            return false;
        }

        return save.finish(name, null, -1);

    }

    @Override
    public boolean addUser(ParamVector pv) {

        super.addUser(pv);

        if (save != null) {
            save.addUser(pv);
        }

        return true;

    }

    @Override
    public boolean createVector(boolean close) {

        if (!super.createVector(false)) {
            return false;
        }

        if (save != null) {
            addBlankParam();
            save.createVector(true);
            addVector(save.getVector());
            save.addUser(this);
        }

        if (close) {
            closeVector();
        }

        return true;

    }

    @Override
    public void updateVector(ParamObject[] v) {

        super.updateVector(v);

        if (save != null) {
            save.updateVector(v);
        }

    }

    // do not use this function directly, use set() instead
    @Override
    public boolean setMyParams(ParamObject o, double v) {

        if (o == null) {
            return false;
        }

        if (!(o.paramVector instanceof DiffusantParticlesProjectionKeiding)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("C0")) {
            if (v < 0) {
                return false;
            }
            //C0 = plu;
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

        if (!(o.paramVector instanceof DiffusantParticlesProjectionKeiding)) {
            return false;
        }

        String n = o.getName();

        if (n.equalsIgnoreCase("name")) {
            name = s;
            return true;
        }
        //if (n.equalsIgnoreCase("displaySelect")) {
        //displaySelect = s;
        //return true;
        //}
        return super.setMyParams(o, s);
    }

}
