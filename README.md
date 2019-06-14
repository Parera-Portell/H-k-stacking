# H-k-stacking

This repository contains C programs to perform H-k stacking of receiver functions (RF) to estimate crustal thickness and Vp/Vs ratio as in Zhu and Kanamori (2000), doi:10.1029/1999jb900322.
The program hkstacking.c converts a single time-domain RF to a H-k-domain RF.
The program hk_sum.c performs stacking and normalization of multiple H-k-domain RF, and outputs a file containing H, Vp/Vs, and Amplitude which can be easily plotted.

Hkstacking.c relies on the SAC library to read sac files and header variables. Be sure to check the SAC manual (https://ds.iris.edu/files/sac-manual/) before trying to compile. Ray parameter must be already set in the sac header.

Arguments for hkstacking:
  program [RF file] [output file] [minimum H] [maximum H] [minimum kappa] [maximum kappa] [RF sample rate] [P wave velocity] [weight1] [weight2] [weight3]

** Additional parameters (MAX length of data array, SEC seconds before time 0 in RF, STEP number of iterations, and RAYP header variable containing ray parameter) should be defined in hkstacking.c as constants. **

Arguments for hk_sum:
  program [output file] [H-k RF 1] [H-k RF 2] ... [H-k RF N]
  
Both programs are short and simple, easily controlled through a bash script to do batch processing and plotting.
 
