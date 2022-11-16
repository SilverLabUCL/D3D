package ucl.silver.d3d.core;

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
public class DiffusantCalretininCaFaas extends DiffusantCalmodulinCaFaas {
    // see DiffusantCalmodulinCaFaas
    //
    // "Resolving the Fast Kinetics of Cooperative Binding: Ca2+ Buffering by Calretinin"
    // PLoS Biol. 2007 Nov;5(11):e311.
    // Faas GC, Schwaller B, Vergara JL, Mody I.
    // Room Temperature ~25C
    //
    // Calretinin is a calcium-binding protein belonging to the calmodulin superfamily
    //
    public DiffusantCalretininCaFaas(Project p, String NAME, double CTOTAL, double diffusionCoefficient, CoordinatesVoxels c,
            int ICa, int ICaTR) {

        super(p, NAME, CTOTAL, diffusionCoefficient, c, ICa, ICaTR, "");

        // Ctotal = 0.1; // mM
        
        // D = 0.02 um^2/ms
        // Faas et al... "a diffusion coefficient of approximately 20 um^2/s
        // assuming a similar mobility as the closely related CaBP CB [29,35]"
        
        iCa = ICa;
        iCaTR = ICaTR;

        reaction = true;
        
        name = "Calretinin-Ca";
        K1on = 3.6; // on reaction rate (1/(mM*ms))
        K1off = 0.053; // off reaction rate (1/ms)
        K2on = 310.0; // on reaction rate (1/(mM*ms))
        K2off = 0.040; // off reaction rate (1/ms)
        Kd1 = K1off / K1on;
        Kd2 = K2off / K2on;

        temp = 25; // Room Temperature ~25C

        createVector(true);

    }

}
