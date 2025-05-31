#include <iostream>
#include <fstream> //we need this for fasta input and output
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
using namespace std; //i used this bcs i keep forgetting to add std::

struct Alignment { 
    string alignedA, alignedB; 
    int score; 
}; //so we can have aligned strings and their score grouped

enum Matrix { matM, matI, matD }; //i use this later to mark which DP part is being done in the backtracking
/*matM is matrix for (mis)match, matI is inserttion (so gap in B) and matD is deletion which marks gap in A*/

//this function is a function to read fasta file, i used GPT to write it
string fastaReader(const string &path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open FASTA file: " + path);
    string line, seq;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '>') continue;
        seq += line;
    }
    return seq;
}

class Threads {
    vector<thread> workers;
    queue<function<void()>> tasks; //work left to be done
    mutex lock;
    condition_variable cvQueue;
    //threads will sleep when there are no tasks and wake up when theres a new task or there is shutdown of the pool of threads
    bool shutDown = false;

public:
    Threads(size_t n) {
        for (size_t i = 0; i < n; i++)
            workers.emplace_back([this] { workerLoop(); });
    }//n threads that run worker loop

    ~Threads() {
        {
            lock_guard<mutex> lk(lock);
            shutDown = true;
        }
        cvQueue.notify_all();
        for (auto &w : workers) w.join();
    }//signals shutdown, notifies and joins threads

    
    void submit(function<void()> f) {
        {
            lock_guard<mutex> lk(lock);
            tasks.push(move(f));
        }
        cvQueue.notify_one();
    }

private:
    void workerLoop() {
        while (true) {
            function<void()> task;
            {
                unique_lock<mutex> lk(lock);
                cvQueue.wait(lk, [&] { return shutDown || !tasks.empty(); });
                if (shutDown && tasks.empty()) {
                    return;
                }
                task = move(tasks.front());
                tasks.pop();
            }
            task();
        }
    }
};

class CountdownLatch {
    mutex lock;
    condition_variable cv;
    int count;
public:
    explicit CountdownLatch(int n){
        count = n; //initializing count as number of threads
    }
    void wait() {
        unique_lock<mutex> lk(lock);
        cv.wait(lk, [&]{ return count == 0; });//we block until all parts are done (count is then 0)
    }
    void decrement() {
        lock_guard<mutex> lk(lock);
        if (--count == 0) {
            cv.notify_all();
        }
    }//after finishing their part, workers decrement
};

Alignment gotoh_align(
    const string &A,
    const string &B,
    int openGap,
    int extendGap,
    const vector<vector<int>> &submat
) {
    const int m = A.size(), n = B.size();//seq lengths
    const int NEG_INF = numeric_limits<int>::min()/2;//i googled how to put a neg inf to mark impossible state


//now I start with allocating and initializinf DP matrices

    vector<vector<int>> M(m+1, vector<int>(n+1, NEG_INF));
    vector<vector<int>> I(m+1, vector<int>(n+1, NEG_INF));
    vector<vector<int>> D(m+1, vector<int>(n+1, NEG_INF));

    //traceback matrices record what was previously chosen at a cell
    //i initialize them to invalid values
    vector<vector<char>> traceM(m+1, vector<char>(n+1, -1));
    vector<vector<char>> traceI(m+1, vector<char>(n+1, -1));
    vector<vector<char>> traceD(m+1, vector<char>(n+1, -1));

    //base cases 
    M[0][0] = 0;//zero cost of aligning two emtpty prefixes
    for (int i = 1; i <= m; i++){
        I[i][0] = -openGap - (i-1)*extendGap;//cost of aligning A[...i] to empty B wuth one open gap and extensions
        traceI[i][0] = 1;//bcs all gaps in I column 0 are extensions
    }
    for (int j = 1; j <= n; j++){
        D[0][j] = -openGap - (j-1)*extendGap;//analog to I
        traceD[0][j] = 1;
    }

    //here is first put a number for numthreads but after googling there is a hardware concurrency default so Klara feel free to override this to test diff thread counts
    const int num_threads = thread::hardware_concurrency();
    Threads pool(num_threads);

    //here i deal wiht the anti diagonals (so cells where i+j is diagonal)
    for (int diag = 1; diag <= m + n; diag++) {
        int i_min = max(1, diag - n);//these are just to compute the range that is valid
        int i_max = min(m, diag - 1);
        int total = i_max - i_min + 1;//nb of cells on that diagonal
        if (total <= 0) continue;

        CountdownLatch latch(num_threads);
        //now i divide total into numthread chunks plus some will get an extra cell if rem is not 0
        int base = total / num_threads;
        int rem = total % num_threads;
        int start = i_min;


        
        for (int t = 0; t < num_threads; t++) {
            int chunk = base ;
            if (t < rem) {
                chunk++;
            } 
            int cs = start;
            int ce = min(cs + chunk - 1, i_max);  // to prevent overflow
            start += chunk;
            //here i want to submit chunk tasks
            pool.submit([&, cs, ce, diag, chunk] {
                if (chunk > 0) {
                    for (int i = cs; i <= ce; i++) {
                        int j = diag - i;
                        if (j < 1 || j > n) continue;//so we dont go out of bound            
                        //writing matrix I
                        int openI = M[i-1][j] - (openGap + extendGap);
                        int extI = I[i-1][j] - extendGap;
                        I[i][j] = max(openI, extI);
                        if (openI > extI) {
                            traceI[i][j] = 0;
                        }
                        else{
                            traceI[i][j] = 1;
                        } 

                        //matrix D
                        int openD = M[i][j-1] - (openGap + extendGap);
                        int extD = D[i][j-1] - extendGap;
                        D[i][j] = max(openD, extD);
                        if (openD > extD) {
                            traceD[i][j] = 0;
                        }
                        else{
                            traceD[i][j] = 1;
                        } 

                        //matrix M
                        int sub = submat[A[i-1]][B[j-1]];
                        int diagS = M[i-1][j-1] + sub;
                        M[i][j] = max({diagS, I[i][j], D[i][j]});
                        
                        if (M[i][j] == diagS) traceM[i][j] = 0;
                        else if (M[i][j] == I[i][j]) traceM[i][j] = 1;
                        else traceM[i][j] = 2;
                    
                }

                }
                latch.decrement();
            });
        }
        latch.wait(); //here we have to wait for the current diagonal to complete
    }

    //this part is backtracking to find the alignment
    //first i want to find a correct starting point
    int final_score = max({M[m][n], I[m][n], D[m][n]});
    Matrix current;
    if (final_score == M[m][n]) {
        current = matM;
    } else if (final_score == I[m][n]) {
        current = matI;
    } else {
        current = matD;
    }
    string Aaligned, Baligned;
    int i = m, j = n; //we start here because the final score is in M[m][n]

    while (i > 0 || j > 0) {//we loop until they are both 0
        if (current == matM) {
            //egde cases
            if (i == 0) {
                current = matD;
                continue;
            }
            if (j == 0) {
                current = matI;
                continue;
            }
            switch (traceM[i][j]) {
                case 0: // if we came from M[i-1][j-1]
                    Aaligned += A[i-1];
                    Baligned += B[j-1];
                    i--;
                    j--;
                    break;
                case 1: // if we came from I[i][j] we switch to I without taking characters
                    current = matI;
                    break;
                case 2: // similar to prev case
                    current = matD;
                    break;
            }
        }
        else if (current == matI) {
            //first forcing gap extension in this case
            if (i == 0) {
                current = matD;
                continue;
            }
        //take one character from A and gap in B and then see traceI to check if we stay in I (which is extension) or if we return to M (which is open)
            Aaligned += A[i-1];
            Baligned += '-';
            if (traceI[i][j] == 0) {
                current = matM;
            }
            else{
                current = matI;
            }
            i = max(0, i-1);
        }
        else { // case similar to prev
            if (j == 0) {
                current = matI;
                continue;
            }
            Baligned += B[j-1];
            Aaligned += '-';
            if (traceD[i][j] == 0) {
                current = matM;
            }
            else{
                current = matD;
            }
            j = max(0, j-1);
        }
    }

    //we reverse so that the strings we get are from left to right 
    reverse(Aaligned.begin(), Aaligned.end());
    reverse(Baligned.begin(), Baligned.end());

    return {Aaligned, Baligned, final_score};
}

//i used GPT to give me simple main function 
int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " seqA.fasta seqB.fasta\n";
        return 1;
    }

    string A = fastaReader(argv[1]);
    string B = fastaReader(argv[2]);

    // Simple match/mismatch matrix for nucleotides
    vector<vector<int>> submat(128, vector<int>(128, -1));
    for (char c : {'A','C','G','T'}) submat[c][c] = 1;

    Alignment result = gotoh_align(A, B, 10, 1, submat);

    cout << "Score: " << result.score << "\n";
    cout << "A: " << result.alignedA << "\n";
    cout << "B: " << result.alignedB << "\n";

    return 0;
}