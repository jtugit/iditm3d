#!/bin/bash
mpiexec -n $1 ./iditm3d -ts_type cn -ksp_type fgmres -pc_type jacobi -ts_rtol 1.0e-8 -ts_atol 1.0e-8 \
-snes_monitor -ksp_monitor -ts_max_snes_failures -1
