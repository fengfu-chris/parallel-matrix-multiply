#!/bin/sh
m=2000
p=1000
n=3000
np=5

echo -e "\n====== Run tests ======\n"
echo -e "\nTest matrix multiply with MPI"
echo "np=$np, m=$m, p=$p, n=$n"
mpirun -np $np mm_MPI $m $p $n

echo -e "\nTest matrix multiply with MPI-updated"
echo "np=$np, m=$m, p=$p, n=$n"
mpirun -np $np mm_MPI_update $m $p $n

echo -e "\nTest matrix multiply with MPI + SCMatrix"
echo "np=$np, m=$m, p=$p, n=$n"
mpirun -np $np mm_MPI_SCMatrix $m $p $n

echo -e "\n===== Test complete =====\n"
