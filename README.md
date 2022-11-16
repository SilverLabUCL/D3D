# D3D
D3D is a 3D reaction-diffusion simulator created by The Silver Lab at University College London (http://silverlab.org/) written in Java.

To simulate diffusion, D3D uses either an explicit finite-difference method similar to that described by Crank (1975) or a Monte Carlo algorithm for non-overlapping hard spheres that explicitly accounts for excluded volume effects (Cichocki and Hinsen, 1990).

For more information contact j.rothman@ucl.ac.uk

# D3D Features Include:
1. An easy-to-use graphical interface for viewing simulation parameters, geometries, point-spread functions (PSFs), diffusants, sources and detectors.

2. The ability to preview simulations in a 2D viewer.

3. Functions for specifying arbitrary 3D geometries.

4. Functions for specifying illumination and detection PSFs for photolysis reactions (see PSF.java) including uncaging and fluorescence recovery after photobleaching (FRAP).

5. Automatic generation of log files that include nearly all simulation parameters.

6. Master class functions for creating user-defined simulations (see Master.java).

7. Demo initialisation classes that include several simulation examples (D3D > Project > Init).

8. Functions for creating 2D projections of spherical particles for Monte Carlo simulations, including demo functions for recreating the 2D projections used in the preprint 'Validation of a stereological method for estimating particle size and density from 2D projections with high accuracy' https://doi.org/10.1101/2022.10.21.513285 (D3D > Project > Init > InitMC_Projection_Demo).

# Published D3D Studies Include:
1. Glutamate release in the synaptic cleft between cerebellar mossy fiber terminals (MFTs) and granule cell (GC) dendrites (Nielsen et al., 2004).

2. Glutamate uncaging in the synaptic cleft of the cerebellar MFT, to investigate the desensitization properties of post-synaptic AMPA receptors of GCs (DiGregorio et al., 2007).

3. The spatiotemporal distribution of Ca2+ in the vicinity of voltage-gate Ca2+ channel (VGCC) clusters at the calyx of Held terminal, including Ca2+ diffusion and binding with buffers and fluorescence dye, and the detection of fluorescence using a Gaussian confocal PSF (Nakamura et al., 2015; Nakamura et al., 2018).

4. Brownian motion of synaptic vesicles in cerebellar MFTs, including steric interactions (Rothman et al., 2016).

5. Synaptic vesicle dynamics near active zones, including connectors, tethers, docking and release (Rothman et al., 2016).

6. 2D projections of spherical particles in planar and thick sections (unpublished preprint Rothman et al., 2022).

# References
B Cichocki, K Hinsen (1990) Dynamic computer simulation of concentrated hard sphere suspensions: I. Simulation technique and mean square displacement data. Physica A 166:473–491. https://doi.org/10.1016/0378-4371(90)90068-4

J Crank (1975) The Mathematics of Diffusion. Oxford UK: Clarendon Press.

DA DiGregorio, JS Rothman, TA Nielsen, RA Silver (2007) Desensitization properties of AMPA receptors at the cerebellar mossy fiber granule cell synapse. Journal of Neuroscience 27:8344–8357. https://doi.org/10.1523/JNEUROSCI.2399-07.2007

Y Nakamura, H Harada, N Kamasawa, K Matsui, JS Rothman, R Shigemoto, RA Silver, DA DiGregorio, T Takahashi (2015) Nanoscale distribution of presynaptic Ca(2+) channels and its impact on vesicular release during development. Neuron 85:145–158. https://doi.org/10.1016/j.neuron.2014.11.019

Nakamura Y, Reva M, DiGregorio DA (2018) Variations in Ca2+ Influx Can Alter Chelator-Based Estimates of Ca2+ Channel-Synaptic Vesicle Coupling Distance. J Neurosci. 38(16):3971-3987. https://doi.org/10.1523/JNEUROSCI.2061-17.2018.

TA Nielsen, DA DiGregorio, RA Silver (2004) Modulation of glutamate mobility reveals the mechanism underlying slow-rising AMPAR EPSCs and the diffusion coefficient in the synaptic cleft. Neuron 42:757–771. https://doi.org/10.1016/j.neuron.2004.04.003

JS Rothman, L Kocsis, E Herzog, Z Nusser, RA Silver (2016) Physical determinants of vesicle mobility and supply at a central synapse. Elife. 2016 Aug 19;5. pii: e15133. https://doi.org/10.7554/eLife.15133

JS Rothman, C Borges-Merjane, N Holderith, P Jonas, RA Silver (2022) Validation of a stereological method for estimating particle size and density from 2D projections with high accuracy. bioRxiv. 2022 Oct 25. https://doi.org/10.1101/2022.10.21.513285
# D3D GUI
![alt tag](http://ucl.ac.uk/silverlab/images/D3D_MonteCarlo_FRAP.png)

Screnshot of a D3D Monte Carlo FRAP simulation showing mitochondria (dark grey), synaptic vesicles (green and light gray) and the illumination point-spread function (blue). See Rothman et al., 2016.
