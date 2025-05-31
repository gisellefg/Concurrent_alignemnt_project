CPU code:

first compiled with g++ -std=c++17 -03 gotohCPU.cpp -pthread -o gotoh
then it is run with ./gotoh seqA.fasta seqB.fasta n

here for n u can put any number u want i tested with 4


GPU compilation and running: (remember to uncomment main function)

/usr/local/cuda/bin/nvcc gotohCUDA.cu -o gotohCUDA -arch=sm_75
./gotohCUDA 


Running main.cpp (runs GPU and CPU at the same time and compares)

 nvcc -std=c++17 -O3 -I$CUDA_HOME/include -L$CUDA_HOME/lib64 gotohCUDA.cu main.cpp gotohCPU.cpp -o benchmark
 ./benchmark

 Example output:
Scores (CPU/GPU): 208  208 Times taken: 14.2469 ms (CPU) 8.09587 ms (GPU). Speedup: 1.75978