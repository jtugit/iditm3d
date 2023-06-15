#!/bin/bash
mpiexec -n $1 valgrind \
 --tool=memcheck --trace-children=yes --leak-check=full --show-reachable=yes \
 --undef-value-errors=yes --show-leak-kinds=all ./iditm3d -ksp_type fgmres -pc_type bjacobi  \
-ksp_rtol 1.0e-8 -ksp_atol 1.0e-8 -snes_monitor -snes_rtol 1.0e-8 -snes_atol 1.0e-8 -snes_stol 1e-8 \
-ksp_converged_reason -snes_converged_reason
