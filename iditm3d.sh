#!/bin/bash
mpiexec -n $1 ./iditm3d -ts_type arkimex -ksp_type gcr -pc_type asm -ts_rtol 1.0e-8 -ts_atol 1.0e-8 \
-snes_monitor -ksp_converged_reason -ksp_monitor_true_residual -ts_max_snes_failures -1
