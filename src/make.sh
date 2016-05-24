#!/bin/sh

g++ SCMatrix.cpp -c
g++ SCVector.cpp -c
g++ utils.cpp -c

mpicxx mm_single_proc.cpp utils.o -o mm_single_proc
mpicxx mm_MPI.cpp utils.o -o mm_MPI
mpicxx mm_MPI_update.cpp utils.o -o mm_MPI_update
mpicxx mm_MPI_SCMatrix.cpp SCMatrix.o -o mm_MPI_SCMatrix

echo "Compilation succeed."
