#include <random>
#include <vector>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>

#include "secded.h"

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

    /*
    initialize the random elements and set the bit error chance
    */
    void wrapper_init(double bit_error_rate) {
        random_device rd;
        RNG::seed(rd());

        lambda = bit_error_rate;
    }

    /*
    Check if a bit error occurred
    */
    bool check_fault() {
        RNG myrng;
        double val = myrng.gen_double();

        return (val < lambda);
    }

    /*
    obsolete
    */
    vector<int> get_flipped_bits(int length, int type_size) {
        vector<int> flipped(fault_flip_count, 0);
        RNG myrng;
        myrng.set_int(0, (length * type_size) - 1);
        for (int &i : flipped) {
            i = myrng.gen_double(); 
        }

        return flipped;
    }

    void FI_MPI_Send(const void *sendbuf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {

        int type_size, total_size; // in bytes
        MPI_Type_size(datatype, &type_size);
        total_size = type_size * count;
        // Allocate new send buffer that has an added byte for each double's check byte
        int new_size = count * (type_size + 1); // in bytes
        unsigned char *sendbuf_ft = malloc(new_size);
        double *buf = static_cast<double *>(sendbuf);

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
                        printf("FAULT! Flipping Byte %d, Bit %d\n", b, bit_idx);
                    }
                }
            }

            // Increment byte array to next entry
            it += (type_size + 1);
        }

        // Call the api with the new array and size
        MPI_Send(sendbuf_ft, new_size, MPI_UNSIGNED_CHAR, dest, tag, comm);

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
        unsigned char *recvbuf_ft = malloc(new_size);
        double *buf = static_cast<double *>(recvbuf);

        // Call the api with the new array and size
        MPI_Recv(recvbuf_ft, new_size, MPI_UNSIGNED_CHAR, source, tag, comm, status);

        // Loop through each chunk and check validity
        unsigned char *it = recvbuf_ft;
        for (int i=0; i < count; i++) {
            if (checkReceivedData(it)) {
                // Data is good, set in true output
                buf[i] = charToDouble(it);
            } else {
                printf("Uncorrectable error in Byte %d", i);
                //TODO what happens now?
            }
            it += (type_size + 1);
        }

        free(recvbuf_ft);
        return;
    }

    void FI_MPI_Reduce(void *sendbuf, void *recvbuf, int count, 
    MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
        
        

        // Call the actual communication now
        MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);

        return;
    }
}