#!/bin/bash
mpiexec -n $1 ./iditm3d -ts_type arkimex -snes_tyep ngmres -ksp_type cgs -pc_type asm \
-ts_rtol 1.0e-4 -ts_atol 1.0e-4 -snes_monitor -snes_rtol 1.0e-4 -snes_atol 1.0e-4 \
-snes_converged_reason #-snes_test_jacobian \
#-snes_test jacobian_view #-ksp_converged_reason -ksp_monitor_true_residual -ts_max_snes_failures -1
