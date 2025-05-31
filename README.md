CPU code is first compiled with g++ -std=c++17 -03 gotohCPU.cpp -pthread -o gotoh
then it is run with ./gotoh seqA.fasta seqB.fasta n

here for n u can put any number u want i tested with 4


GPU compilation and running:

run first:
/usr/local/cuda/bin/nvcc gotohCUDA.cu -o gotohCUDA -arch=sm_75
then: ./gotohCUDA 
