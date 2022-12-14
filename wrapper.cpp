#include <cstring>
#include <random>
#include <vector>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>

#include "secded.h"

#define MAIN_PE 0
#define PRINT_OUT 0

using namespace std;

namespace wrapper {

    typedef mt19937                        ENG;     // Mersenne Twister
    typedef uniform_int_distribution<>     iDIST;   // Uniform Integer Distribution
    typedef uniform_real_distribution<>  dDIST;     // Uniform Double Distribution

    class RNG {
        private:
            static ENG e;
            iDIST idist;
            dDIST ddist;
        public:
            static void seed(int s) {e.seed(s);}
            RNG() : idist(0,1), ddist(0,1) {};
            void set_int(int imin, int imax) { idist = iDIST(imin, imax); }
            void set_double(double dmin, double dmax) { ddist = dDIST(dmin, dmax); }
            int gen_int() { return idist(e); }
            double gen_double() { return ddist(e); }

        /* 
        USAGE:
        RNG::seed(ANY_SEED); // or wrapper_init()
        RNG myrng();
        myrng.set_int(0,10);
        int val = myrng.gen_int();
        myrng.set_double(0,10);
        double val = myrng.gen_double();
        */
    };
    ENG RNG::e;

    double lambda = 0;
    static int world_rank = -1;
    static int world_size = -1;

    static int total_faults = 0;
    static int total_sec = 0;
    static int total_ded = 0;
    static int total_retrans = 0;

    /*
    initialize the random elements and set the bit error chance
    */
    void wrapper_init(double bit_error_rate, int rank, int size) {
        random_device rd;
        RNG::seed(rd());

        lambda = bit_error_rate;
        world_rank = rank;
        world_size = size;

        if (world_rank == MAIN_PE) printf("Wrapper initialized with lambda = %f\n", lambda);
    }

    /*
    Reset the error counts, retransmission, etc
    */
    void wrapper_reset() {
        total_faults = 0;
        total_retrans = 0;
        total_sec = 0;
        total_ded = 0;
    }

    // Combines the total_fault count from every node
    void combine_faults() {
        int local_faults = total_faults;
        //MPI_Reduce(&local_faults, &local_faults, 1, MPI_INT, MPI_SUM, MAIN_PE, MPI_COMM_WORLD);
        total_faults = local_faults;
    }

    // Provides the Fault Injector runtime information
    void wrapper_output(int &faults, int &sec, int &ded, int &retrans) {
        faults = total_faults;
        sec = total_sec;
        ded = total_ded;
        retrans = total_retrans;
    }

    /*
    Check if a bit error occurred
    */
    bool check_fault() {
        RNG myrng;
        double val = myrng.gen_double();

        return (val < lambda);
    }

    void FI_MPI_Send(const void *sendbuf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {

        int type_size, total_size; // in bytes
        MPI_Type_size(datatype, &type_size);
        total_size = type_size * count;
        // Allocate new send buffer that has an added byte for each double's check byte
        int new_size = count * (type_size + 1); // in bytes
        unsigned char *sendbuf_ft = (unsigned char *) malloc(new_size);
        const double *buf = static_cast<const double *>(sendbuf);

        int call_num = 0;
        int passed = 0;
        do {
            // Loop through 8 bytes at a time
            unsigned char *it = sendbuf_ft;
            for (int i=0; i<count; i++) {
                // Copy data to new buffer
                memcpy(it, &buf[i], type_size);
                // send to get check byte added
                addCheckByte(buf[i], it);

                // Loop through each byte for a fault/flip
                for (int b=0; b < (type_size + 1); b++) {
                    // 8 bits per byte
                    for (int bit_idx=0; bit_idx < 8; bit_idx++) {
                        if (check_fault()) {
                            it[b] ^= (0x1<<bit_idx);
                            total_faults++;
                            if (PRINT_OUT == 1) printf("FAULT! Node %d, Call %d, Chunk %d, Byte %d, Bit %d\n", world_rank, call_num, i, b, bit_idx);
                        }
                    }
                }

                // Increment byte array to next entry
                it += (type_size + 1);
            }

            // Call the api with the new array and size
            MPI_Send(sendbuf_ft, new_size, MPI_UNSIGNED_CHAR, dest, tag, comm);

            // Check if successful or retransmission necessary
            MPI_Recv(&passed, 1, MPI_INT, dest, tag, comm, MPI_STATUS_IGNORE);
            call_num++;
        } while (passed == 0);

        free(sendbuf_ft);
        return;
    }

    void FI_MPI_Recv(void *recvbuf, int count, MPI_Datatype datatype, int source, int tag, 
    MPI_Comm comm, MPI_Status *status) {

        // Create new buffer large enough for the data
        int type_size, total_size; // in bytes
        MPI_Type_size(datatype, &type_size);
        total_size = type_size * count;
        // Allocate new recv buffer that has an added byte for each double's check byte
        int new_size = count * (type_size + 1); // in bytes
        unsigned char *recvbuf_ft = (unsigned char *) malloc(new_size);
        double *buf = static_cast<double *>(recvbuf);

        int call_num = 0;
        vector<int> valid(count, 0);
        int passed = 0; // 0 means request retransmission
        do {
            if (PRINT_OUT == 1) printf("MPI Call: From %d to %d Call %d\n", source, world_rank, call_num);
            call_num++;

            // Set passed to 1 then check for errors
            passed = 1;

            // Call the api with the new array and size
            MPI_Recv(recvbuf_ft, new_size, MPI_UNSIGNED_CHAR, source, tag, comm, status);

            // Loop through each chunk and check validity
            unsigned char *it = recvbuf_ft;
            for (int i=0; i < count; i++) {
                // check if valid data already saved
                if (valid[i] == 0) {
                    int decode = checkReceivedData(it);
                    if (decode == 0) {
                        // Data is good, set in true output
                        buf[i] = charToDouble(it);
                        valid[i] = 1;
                    } else if (decode == 1) {
                        // SEC occurred, still good
                        total_sec++;
                        buf[i] = charToDouble(it);
                        valid[i] = 1;
                    } else {
                        // DED occurred, not good
                        total_ded++;
                        if (PRINT_OUT == 1) printf("Uncorrectable error from Node %d, Chunk %d\n", source, i);
                        buf[i] = 0;
                        passed = 0; // retransmission needed
                    }
                }
                
                it += (type_size + 1);
            }

            // Return retransmission request to sender
            MPI_Send(&passed, 1, MPI_INT, source, tag, comm);

        } while (passed == 0);

        total_retrans += (call_num - 1);

        free(recvbuf_ft);
        return;
    }

    void FI_MPI_Reduce(void *sendbuf, void *recvbuf, int count, 
    MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
        int type_size; // in bytes
        MPI_Type_size(datatype, &type_size);

        if (world_rank == root) {
            // Create temp buffer to hold results
            double *tempbuf = (double *) malloc(type_size * count);
            // Recast real result buffer to doubles
            double *out_buf = static_cast<double *>(recvbuf);
            double *in_buf = static_cast<double *>(sendbuf);

            for (int i=0; i<world_size; i++) {
                if (i == world_rank) {
                    for (int i=0; i<count; i++)
                        out_buf[i] += in_buf[i];
                    continue;
                }

                // Recieve data from node i
                FI_MPI_Recv(tempbuf, count, datatype, i, 0, comm, MPI_STATUS_IGNORE);
                // Add data to real buffer
                for (int i=0; i<count; i++)
                    out_buf[i] += tempbuf[i];
            }
            free(tempbuf);
        } else {
            FI_MPI_Send(sendbuf, count, datatype, MAIN_PE, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }
}