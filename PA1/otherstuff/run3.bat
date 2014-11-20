#
# Batch submission script for Prog. Assignment 1
# If this script is placed in file RunPA1.bat,
# to submit batch job, use "qsub RunPA1.bat"
# Status of batch job can be checked using "qstat -u osu...." 
#
#
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=12
#PBS -N PA1
#PBS -S /bin/ksh
#PBS -j oe
cd $TMPDIR
#
# Assumes the source program is in file pa2-p1-sol.c in a
# subdirectory called 5441/Au14/PA1 in the OSC home directory
# Change it appropriately to match your directory and file name
#
cp $HOME/5441/Au14/PA1/pa3.c .

echo " "; echo "--------------------------------------------------------------"
echo "Compiling Problem 3"
echo "--------------------------------------------------------------"; echo " "
gcc -O3 -fopenmp -o pa1-p3-gcc pa3.c
icc -fast -openmp -o pa1-p3-icc pa3.c
echo " "; echo "ICC: Run 1"; echo " "
./pa1-p3-icc
echo " "; echo "GCC: Run 1"; echo " "
./pa1-p3-gcc
echo " "; echo "ICC: Run 2"; echo " "
./pa1-p3-icc
echo " "; echo "GCC: Run 2"; echo " "
./pa1-p3-gcc
echo " "; echo "ICC: Run 3"; echo " "
./pa1-p3-icc
echo " "; echo "GCC: Run 3"; echo " "
./pa1-p3-gcc

