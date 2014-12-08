#PBS -N mpi_ping_pong
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=12
#PBS -j oe
export OMP_NUM_THREADS=12
export MV2_ENABLE_AFFINITY=0
cd $TMPDIR

#
# Copy file
#

cp $HOME/5441/Au14/PA4/ping_pong.c .

#
# Compile
#
module load mvapich2
mpicc -o ping_pong.out -O3 ping_pong.c

#
# Run
#

mpiexec -pernode ping_pong.out
