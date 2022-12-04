#include <random>
#include <vector>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>

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

    double fault_chance = 0;
    int fault_flip_count = 1;

    void wrapper_init(double fc, int ffc) {
        random_device rd;
        RNG::seed(rd());

        fault_chance = fc;
        fault_flip_count = ffc;
    }

    bool check_fault() {
        RNG myrng;
        myrng.set_double(0,1);
        double val = myrng.gen_double();

        return (fault_chance > val);
    }

    vector<int> get_flipped_bits(int length, int type_size) {
        vector<int> flipped(fault_flip_count, 0);
        RNG myrng;
        myrng.set_int(0, (length * type_size) - 1);
        for (int &i : flipped) {
            i = myrng.gen_double(); 
        }

        return flipped;
    }

    void FI_MPI_Send() {

        return;
    }

    void FI_MPI_Recv() {

        return;
    }

    void FI_MPI_Reduce(void *sendbuf, void *recvbuf, int count, 
        MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
        
        // Inject faults into the sendbuf before communicating
        if (check_fault) {
            // fault occurred, get the index of bits flipped
            int type_size;
            MPI_Type_size(datatype, &type_size);
            vector<int> flipped = get_flipped_bits(count, type_size);

            uint8_t *buf = static_cast<uint8_t *>(sendbuf);
            for (int idx : flipped) {
                int byte_idx = idx / 8;
                int bit_idx = idx % 8;
                buf[byte_idx] ^= (1<<bit_idx);
                printf("Flipping bit %d => (%d,%d)\n", idx, byte_idx, bit_idx);
            }
        }

        // Call the actual communication now
        MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);

        return;
    }
}