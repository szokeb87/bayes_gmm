#! /bin/sh

# This shell script works for a quad board machine with quad core CPUs
# running OpenMPI Version 2.1.

export PATH="$PATH:/home/balint/openmpi/bin/"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/balint/openmpi/lib:/usr/lib/x86_64-linux-gnu"
export C_INCLUDE_PATH="$C_INCLUDE_PATH:/home/balint/openmpi/include"

echo "localhost cpu=16" > OpenMPIhosts

test -f bayes_gmm_mpi.err  && mv -f bayes_gmm_mpi.err  bayes_gmm_mpi.err.bak
test -f bayes_gmm_mpi.out  && mv -f bayes_gmm_mpi.out  bayes_gmm_mpi.out.bak

rm -f core core.*

make -f makefile.mpi >bayes_gmm_mpi.out 2>&1 && \
  $HOME/openmpi/bin/mpirun --hostfile OpenMPIhosts ${PWD}/bayes_gmm_mpi >>bayes_gmm_mpi.out 2>bayes_gmm_mpi.err

RC=$?

case $RC in
  0) exit 0 ;;
  esac
exit 1;
~
