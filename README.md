# P53 / Mdm2 Modeling using an EKF

- This repository contains the code for a project for the university course "Think Mathematically, Act Algorithmically: Optimization Techniques (Complementary Studies)" in the winter semester 24/25.

The project uses an Extended Kalman Filter to attempt to recover theoretical parameter values of a system of DDEs modeling p53 and Mdm2 interaction.

## Structure of the repo 
### Main files 
- Simulink_generate_measurement_data.m (run to generate the measurement data needed for the EKF files, need to adjust the loaded Simulink model in the "load" function, depending on if delay parameter tau is assumed constant or not - if constant load "EKF_constantTau_simulation", if non-constant load "EKF_nonconstant_pi")
- EKF_constantTau_dx_ks_tau.m (run if tau is assumed constant and after generating the measurement data)
- EKF_dx_ks_tau_nonconstantTau.m (run if tau is assumed non-constant and after generating measurement data)

### Additional files (**Do not need to be run separately**)
- params.m (Contains theoretical parameter values)
- Files in "data" folder (contain measurement data for constant and non-constant tau cases)
