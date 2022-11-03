#include <iostream>
#include <iomanip>
#include <fstream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <unordered_map>
#include <vector>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>

#define SAVE_RESULTS 0
#define MAIN_PE 0

using namespace std;

typedef vector<int> EdgeList;
typedef map<int, EdgeList> AdjList;
typedef unordered_map<int, int> NodeTranslator;

void printGraph(AdjList &graph, int N) {
    cout << "N = " << N << endl;
    for (auto v : graph) {
        cout << v.first << ": ";
        for (auto edge : v.second)
            cout << edge << " ";
        cout << endl;
    }
}

template <typename T>
void printArray(vector<T> a, int N) {
    for (int i=0; i<N; i++)
        cout << a[i] << " ";
    cout << endl;
}

AdjList loadGraph(string filename, NodeTranslator &nodeList, int &N, int &E) {
    string line;
    ifstream infile(filename);

    AdjList adjList;
    nodeList.clear();
    int count = 0;

    E = 0;
    // Get each edge
    while(getline(infile, line)) {
        // Skip commented lines
        if (ispunct(line.front()))
            continue;

        stringstream ss(line);
        int from, dest;
        ss >> from >> dest;
        
        // Add to node translator if not present
        if (nodeList.find(from) == nodeList.end()) {
            nodeList[from] = count++;
            EdgeList edges;
            adjList.insert({nodeList[from], edges});
        }
        if (nodeList.find(dest) == nodeList.end()) {
            nodeList[dest] = count++;
            EdgeList edges;
            adjList.insert({nodeList[dest], edges});
        }
        adjList[nodeList[from]].push_back(nodeList[dest]);
        E++;

        //cout << from << " - " << nodeList[from] << " | " << dest << " - " << nodeList[dest] << endl;
    }

    N = adjList.size();
    //cout << "count: " << count << " | N = " << N << endl;
    return adjList;
}

vector<double> determineBC(AdjList const &graph, int const N, int const data_start, int const data_end) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    vector<double> bc(N, 0);
    vector<double> bc_local(N,0);

    int j, v, w, src;
    double pairwise, min_bc = 1, max_bc = 0;
    list<int>::iterator it;

    for (src=data_start; src<=data_end; src++) {
        //printf("Rank %d / %d - Thread %d of %d - src = %d\n", world_rank, world_size, omp_get_thread_num( ), omp_get_max_threads( ), src);

        // Setup variables for loop
        stack<int> visited;
        vector<list<int>> pred(N); // Array of predecessor lists
        vector<int> sigma(N, 0);
        vector<int> dist(N, -1);
        queue<int> frontier;
        vector<double> delta(N,0);

        sigma[src] = 1;
        dist[src] = 0;
        frontier.push(src);

        while (!frontier.empty()) {
            v = frontier.front();
            frontier.pop();
            visited.push(v);
            //cout << "v = " << v << endl;

            // Get all neighbors of v
            for (j=0; j<graph.at(v).size(); j++) {
                w = graph.at(v)[j];
                // has w been found before?
                if (dist[w] < 0) {
                    frontier.push(w);
                    dist[w] = dist[v] + 1;
                }
                // shortest path to w via v
                if (dist[w] == (dist[v] + 1)) {
                    sigma[w] += sigma[v];
                    pred[w].push_back(v);
                }
            }
        }

        // Stack contains vertices in order of non-increasing distance form s
        while (!visited.empty()) {
            w = visited.top();
            visited.pop();
            for (it = pred[w].begin(); it != pred[w].end(); it++) {
                v = *it;
                pairwise = ((double) sigma[v] / (double) sigma[w]) * (1.0 + delta[w]);
                delta[v] += pairwise;
            }

            if (w != src) {
                bc_local[w] += delta[w];
            }
            
        }
    }

    // implicit barrier at end
    MPI_Reduce(bc_local.data(), bc.data(), N, MPI_DOUBLE, MPI_SUM, MAIN_PE, MPI_COMM_WORLD);

    // Only the main PE will normalize the data
    if (world_rank == MAIN_PE) {
        //printf("Rank %d / %d: Rescale Begun\n", world_rank, world_size);

        // Rescale the BC by dividing by (N-1)(N-2)/2
        for (int i = 0; i < N; i++) {
            bc[i] /= ((N-1)*(N-2)/2);
            if (bc[i] < min_bc) min_bc = bc[i];
            if (bc[i] > max_bc) max_bc = bc[i];
        }

        //printf("Rank %d / %d: Normalization Begun\n", world_rank, world_size);
        // Normalize BC between 0 and 1 (Check for divide by zero)
        if (max_bc != min_bc) {
            for (auto &g : bc)
                g = (g - min_bc) / (max_bc - min_bc);
        }
    }

    if (world_rank == MAIN_PE)
        printf("Min = %.5g, Max = %.5g, ", min_bc, max_bc);

    //printf("Rank %d / %d: Exiting BC\n", world_rank, world_size);

    return bc;
}

int main(int argc, char* argv[]) {
    // ==================== //
    //  MPI Initialization  //
    // ==================== //
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // ==================== //
    //     Cmd line args    //
    // ==================== //
	// #Trials, *graphname*
	// ./bc_hybrid 5 small_example
    int trials = atoi(argv[1]);
    string graphName(argv[2]);
    string ifname = "./graphs/" + graphName + ".txt";

    // ==================== //
    //     Graph Loading    //
    // ==================== //
    if (world_rank == MAIN_PE) {
        printf("----- Loading graph -----\n");
        printf("File: %s\n", ifname.c_str());
    }
    NodeTranslator nodeList;
    double load_start, load_end;
    int N = 0, E = 0;
    load_start = MPI_Wtime();
    AdjList graph = loadGraph(ifname, nodeList, N, E);
    load_end = MPI_Wtime();
    if (world_rank == MAIN_PE) {
        printf("N = %d, E=%d\n", N, E);
        printf("Graph Load time: %f\n", load_end - load_start);
    }
    //printGraph(graph, N);

    // ==================== //
    //   Load Distribution  //
    // ==================== //
    int data_start, data_end; // FIRST and LAST node to calculate SP weights for
    // Use world_rank to determine role in computations
    data_start = (N / world_size) * world_rank;
    data_end = (N / world_size) * (world_rank + 1) - 1;
    // Last node computes a the remainder in case M_SIZE not evenly divisible by WORLD_SIZE
    if (world_rank == world_size - 1)
        data_end = N - 1;
    // Print a message to signal node is running
    //printf("Processor %s, rank %d / %d: data_start = %d, data_end = %d\n", processor_name, world_rank, world_size, data_start, data_end);

    // ==================== //
    //      BC Testing      //
    // ==================== //
    if (world_rank == MAIN_PE) {
        printf("\n----- BC Results -----\n");
        printf("Nodes = %d\n", world_size);
    }
    double start, end;
    vector<double> bc;
    for (int i = 0; i < trials; i++) {
        MPI_Barrier(MPI_COMM_WORLD); // Ensure all nodes start at the same time
        if (world_rank == MAIN_PE) {
            printf("Trial %d, ", i);
            start = MPI_Wtime();
        }
        bc = determineBC(graph, N, data_start, data_end);
        if (world_rank == MAIN_PE) {
            end = MPI_Wtime();
            printf("Time: %f\n", end - start);
        }
    }
    //printArray(bc, N);

    // ==================== //
    //     Save Results     //
    // ==================== //
    #if SAVE_RESULTS
        if (world_rank == MAIN_PE) {
            double save_start, save_end;
            save_start = MPI_Wtime();
            string ofname = "./outputs/" + graphName + "-mpi-" + to_string(world_size) + ".txt";
            ofstream outfile;
            outfile.open(ofname, ios::out | ios::trunc);

            vector<int> nodesOrdered(N,0);
            for (auto i : nodeList)
                nodesOrdered[i.second] = i.first;

            printf("\n----- Saving Data -----\n");
            printf("BC Data: %s\n", ofname.c_str());
            outfile << "# of Nodes: " << world_size << endl;
            outfile << "Betweenness runtime: " << (end - start) << endl;
            outfile << "BC results:" << endl;
            outfile << "Node # \t| Betweenness" << endl;
            for (int i = 0; i < N; i++) {
                outfile << nodesOrdered[i] << " \t| " << bc[i] << endl;
            }
            outfile.close();
            save_end = MPI_Wtime();
            printf("Save time: %f\n", save_end - save_start);
        }
    #endif

    MPI_Finalize();
    return 0;
}