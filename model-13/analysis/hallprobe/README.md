SI-DIPOLE-BC ANALYSIS

Murilo B. Alves - LNLS-FACS - June 28, 2019

In the beginning the analysis was performed with x=0 in the initial RK trajectory. The average BC model deflection angles was included in the Matlab script 'sirius_nominal_dipole_trajectory'. The matching x0 point (obtained from extrapolation of straight lines) obtained was about 79um different from th x0 aligned in the real machine, which is 7.7030mm.

Based on this difference, the analysis was done with initial Rung-Kutta trajectory at x=+79um and the results are in folder x0-0p079mm. With these results another average model was produced and its deflection angles were included in the Matlab script 'sirius_nominal_dipole_trajectory'.

The x0 obtained was (7.703 +- 0.009) mm, which is different from the installation x0 only by 92nm. With the Matlab script it was produced the reference trajectories 'trajectory-bc-pos.txt' and 'trajectory-bc-neg.txt'. We introduced the reference trajectory in the inputs and performed the analysis over the reference trajectory instead of Runge-Kutta. The results are in x0-0p079mm-reftraj folder.
