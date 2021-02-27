# H-k-stacking

This repository contains software to perform H-k stacking of receiver functions (RF) to estimate crustal thickness and Vp/Vs ratio as in Zhu and Kanamori (2000).

-time_to_hk.c converts a single time-domain RF to an H-k-domain RF.

-hk_stacking.c performs stacking and normalization of multiple H-k-domain RF, and outputs a file containing H, Vp/Vs, and Amplitude which can be easily plotted.

-bootstrap.c computes standard error and correlation of the H-k stack via the bootstrap method.



REQUIREMENTS

    -SAC library. Check the SAC manual (https://ds.iris.edu/files/sac-manual/) before compiling.

ARGUMENTS

    time_to_hk [RF file] [output file] [minimum H] [maximum H] [minimum Vp/Vs] [maximum Vp/Vs] [P wave velocity] [weight1] [weight2] [weight3]

    hk_stacking [output file] [H-k RF 1] [H-k RF 2] ... [H-k RF n]

    bootstrap [output file] [H-k RF 1] [H-k RF 2] ... [H-k RF n]

  
  ******* Additional parameters to be set in  hkstacking.c as constants *******
  
    MAX (maximum length of data)
    
    HSTEP (resolution of H)
    
    KSTEP (resolution of Vp/Vs)
    
    RAYP (SAC header variable containing the ray parameter)
   
