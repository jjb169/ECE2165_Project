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

#include "wrapper.cpp"

#define MAIN_PE 0

using namespace std;

template <typename T>
void printArray(vector<T> a, int N) {
    for (int i=0; i<N; i++)
        cout << a[i] << " ";
    cout << endl;
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

    wrapper::wrapper_init(0.01, world_rank, world_size);

    if (world_rank == MAIN_PE) {
        printf("\n----- Toy Results -----\n");
        printf("Nodes = %d\n", world_size);
    }

    int N = 5;
    vector<double> vec = {0, 0.5, 1, 1.5, 2};
    if (world_rank == MAIN_PE) {
        printArray(vec, N);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    vec[world_rank] += 1;
    /*
    if (world_rank == MAIN_PE) {
        for (int src=1; src<world_size; src++) {
            double buf[5] = {0};
            wrapper::FI_MPI_Recv(buf, N, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i=0; i<N; i++)
                vec[i] += buf[i];
        }
    } else {
        wrapper::FI_MPI_Send(vec.data(), N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    */

    wrapper::FI_MPI_Reduce(vec.data(), vec.data(), N, MPI_DOUBLE, MPI_SUM, MAIN_PE, MPI_COMM_WORLD);

    if (world_rank == MAIN_PE) {
        printArray(vec, N);
    }

    MPI_Finalize();
    return 0;
}