#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int MPI_FlattreeColective(const void *sendbuf, void *recvbuf, int count,
                          MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int error = MPI_SUCCESS;
    void *sum;

    error = MPI_Send(sendbuf, count, datatype, root, world_rank, comm);

    if (root == world_rank) {
        for (int i = 0; i < world_size; ++i) {
            MPI_Recv(sum, count, datatype, MPI_ANY_SOURCE, MPI_ANY_TAG,
                     comm, MPI_STATUS_IGNORE);
            *(int *) recvbuf += *(int *) sum;
        }
    }
    return error;
}

void MPI_BinomialBcast(void *buffer, int count, MPI_Datatype datatype,
                       int root, MPI_Comm comm) {

    int size, rank, src, dst;

    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size == 1) return;

    int relative_rank = (rank >= root) ? rank - root : rank - root + size;

    int mask = 0x1;
//    printf("mask %d, size %d, rank %d \n", mask, size, rank);
    while (mask < size) {
        if (relative_rank & mask) {
            src = rank - mask;

            if (src < 0)
                src += size;

            MPI_Recv(buffer, count, datatype, src, 0, comm, &status);

            break;
        }
        mask <<= 1;
//        printf(" updated mask %d, rank %d \n", mask, rank);
    }


    mask >>= 1;

    while (mask > 0) {
        if (relative_rank + mask < size) {
            dst = rank + mask;

            if (dst >= size) {
                dst -= size;
            }

            MPI_Send(buffer, count, datatype, dst, 0, comm);
        }
        mask >>= 1;
    }
}

int main(int argc, char *argv[]) {

    int i, done = 0, n = 0, count;
    double PI25DT = 3.141592653589793238462643;
    double pi, x, y, z;
    // Unique rank is assigned to each process in a communicator
    int rank;

    // Total number of ranks
    int size;

    // Initializes the MPI execution environment
    MPI_Init(&argc, &argv);

    // Get this process' rank (process within a communicator)
    // MPI_COMM_WORLD is the default communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the total number ranks in this communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // If not rank zero, send your message to be printed
    while (!done) {
        if (rank == 0) {
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d", &n);
        }

        MPI_BinomialBcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (n == 0) break;

        count = 0;

        for (i = (rank + 1); i <= n; i += size) {
            // Get the random numbers between 0 and 1
            x = ((double) rand()) / ((double) RAND_MAX);
            y = ((double) rand()) / ((double) RAND_MAX);

            // Calculate the square root of the squares
            z = sqrt((x * x) + (y * y));

            // Check whether z is within the circle
            if (z <= 1.0)
                count++;

        }
        int local_count = 0;

        MPI_FlattreeColective(&count, &local_count, 1, MPI_INT, MPI_SUM, 0,
                              MPI_COMM_WORLD);

        if (rank == 0) {
            pi = ((double) local_count / (double) n) * 4.0;
            printf("pi is approximately %.16f, Error is %.16f\n",
                   pi, fabs(pi - PI25DT));
        }

    }

    MPI_Finalize();

}


