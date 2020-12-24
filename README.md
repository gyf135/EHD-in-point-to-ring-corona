# EHD-in-point-to-ring-corona
The UDF and FLUENT case to simulate EHD flow in the point-to-ring corona configuration

This code accompanies the paper:
Experimental and numerical investigation of electrohydrodynamic flow in a point-to-ring corona discharge
By: Yifei Guan, Ravi Sankar Vaddi, Alberto Aliseda, and Igor Novosselov
URL: https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.3.043701

To implement the UDF properly:
1. Open the FLUENT case file
2. Compile UDF
3. Setup everything properly as in the example. (The corresponding names of the UDS is listed in the UDF file)
Here something to notice is to unselect the equations for other UDS except for potential and charge density.
The diffusivity of charge density should be specified.
The time-dependent term for potential should be unselected.
The uds1_flux function is actually for uds2 which is the charge density.
4. Un-load the libudf
5. Initialize the flow field (zero velocities should be fine)
6. Run for 1 iteration (No need to run the whole time step)
7. Load the libudf
8. Unselect the charge equation in solution controls.
9. Run it until the potential is steady-state
10. Activate the charge equation solver in solution controls
11. Run it at time step size 1e-7, then incrementally increase the time step size as the charge density is becoming steady state.
