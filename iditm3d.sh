#!/bin/bash
mpiexec -n $1 ./iditm3d -ts_type arkimex -ksp_type fgmres -pc_type jacobi -ts_rtol 1.0e-6 -ts_atol 1.0e-6 \
-snes_monitor -snes_rtol 1.0e-6 -snes_atol 1.0e-6 #-ksp_converged_reason -ksp_monitor_true_residual -ts_max_snes_failures -1
