CPU code:

first compiled with g++ -std=c++17 -03 gotohCPU.cpp -pthread -o gotoh
then it is run with ./gotoh seqA.fasta seqB.fasta n

here for n u can put any number u want i tested with 4


GPU compilation and running: (remember to uncomment main function and fastaReader)

/usr/local/cuda/bin/nvcc gotohCUDA.cu -o gotohCUDA -arch=sm_75
./gotohCUDA 

For linear space CPU:
/usr/local/cuda/bin/nvcc gotoh_linearspace_CPU.cu fasta_reader.cpp -o gotoh_linearspace_CPU -arch=sm_75


For linear space GPU:
/usr/local/cuda/bin/nvcc gotoh_linearspace_GPU.cu fasta_reader.cpp -o gotoh_linearspace_GPU -arch=sm_75



Running main.cpp (runs GPU and CPU at the same time and compares)


(Making sure CUDA_HOME is set might be needed: 
 export CUDA_HOME=/usr/local/cuda
 export PATH=$CUDA_HOME/bin:$PATH
 export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH)

Compile the benchmark with:
 nvcc -std=c++17 -O3 -I$CUDA_HOME/include -L$CUDA_HOME/lib64 gotohCUDA.cu main.cpp gotohCPU.cpp -o benchmark
Run the program with
 ./benchmark


Example output:
Scores (CPU/GPU): 208  208 Times taken: 14.2469 ms (CPU) 8.09587 ms (GPU). Speedup: 1.75978