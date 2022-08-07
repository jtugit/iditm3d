#!/bin/bash
mpiexec -n $1 valgrind \
 --tool=memcheck --trace-children=yes --leak-check=full --show-reachable=no \
--track-origins=yes  --show-leak-kinds=all ./iditm3d
