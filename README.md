This project implements the Gotoh algorithm for pairwise DNA sequence alignment using affine gap penalties. It includes multithreaded CPU, CUDA-based GPU, and linear space variants that allow alignment of longer sequences.

To compile and run the CPU version, use:

    g++ -std=c++17 -03 gotohCPU.cpp -pthread -o gotoh

Then run the program with:

    ./gotoh seqA.fasta seqB.fasta n

Here, `seqA.fasta` and `seqB.fasta` are input sequences in FASTA format. Replace `n` with the number of threads you want to use (e.g., 4).


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


For detailed explanation of the algorithm, implementation, and results, see the accompanying report: DNA_sequence_alignment.pdf.
