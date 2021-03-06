# to submit batch job, use "qsub RunPA2.bat"
# Status of batch job can be checked using "qstat -u osu...." 

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
cp /nfs/13/osu8844/PA2/2/* .

echo " "; echo "--------------------------------------------------------------"
echo "Compiling Programming Assignment 2"
echo "--------------------------------------------------------------"; echo " "
gcc -O3 -fopenmp -o pa2_a pa2_pt_a.c  -lm
gcc -O3 -fopenmp -o pa2_b pa2_pt_b.c  -lm
gcc -O3 -fopenmp -o pa2_c pa2_pt_c.c  -lm
gcc -O3 -fopenmp -o pa2_d pa2_pt_d.c  -lm
gcc -O3 -fopenmp -o pa2_e pa2_pt_e.c  -lm

echo " "; echo "GCC: Run a"; echo " "
./pa2_a
echo " "; echo "GCC: Run b"; echo " "
./pa2_b
echo " "; echo "GCC: Run c"; echo " "
./pa2_c
echo " "; echo "GCC: Run d"; echo " "
./pa2_d
echo " "; echo "GCC: Run e"; echo " "
./pa2_e






