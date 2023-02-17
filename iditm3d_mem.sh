#!/bin/bash
mpiexec -n $1 valgrind \
 --tool=memcheck --trace-children=yes --leak-check=full --show-reachable=no \
--track-origins=yes  --show-leak-kinds=all ./iditm3d -ts_type arkimex \
-snes_tyep ngmres -ksp_type fgmres -pc_type jacobi -ts_rtol 1.0e-4 -ts_atol \
1.0e-4 -snes_monitor -snes_rtol 1.0e-4 -snes_atol 1.0e-4 -ksp_converged_reason
