#!/bin/sh
#BSUB -J test
#BSUB -u y.ma@mfm.tu-darmstadt.de
#BSUB -B -N
# Bitte achten Sie auf vollständige Pfad-Angaben:

#
#BSUB -n 1       # Anzahl der Rechenkerne
#BSUB -M 3000    # Hauptspeicher in MByte für alle Threads
#BSUB -W 0:10     # in Stunden und Minuten, oder '#BSUB -W 10' - nur Minuten
#
#BSUB -a openmp   # Spezifiziert spezielle Unterstützung für OpenMP

# -------------------------------
# Anschließend schreiben Sie Ihre eigenen Befehle, wie z.B.
module unload
module load gcc
module load gnuplot

# Spezifikation von OMP_NUM_THREADS kann entfallen wenn
#   obige Unterstützung "-a openmp" verwendet wurde

gfortran -O3 -ffree-form dom1.f -o d1.out
./d1.out
gfortran -O3 -ffree-form dom2.f -o d2.out
./d2.out
gfortran -O3 -ffree-form dom3.f -o d3.out
./d3.out
gfortran -O3 -ffree-form dom4.f -o d4.out
./d4.out


