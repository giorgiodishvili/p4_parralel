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
    int *sum;

    error = MPI_Send(sendbuf, count, datatype, root, world_rank, comm);

    if (root == world_rank) {
        for (int i = 0; i < world_size; ++i) {
            MPI_Recv(sum, count, datatype, MPI_ANY_SOURCE, MPI_ANY_TAG,
                     comm, MPI_STATUS_IGNORE);

            *(int *) recvbuf += *(int*) sum;
        }
    }
    return error;
}

int main(int argc, char *argv[]) {

    int i, done = 0, n, count;
    double PI25DT = 3.141592653589793238462643;
    double myPi, pi, x, y, z;
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

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
        myPi = ((double) count / (double) n) * 4.0;

        MPI_FlattreeColective(&myPi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0,
                              MPI_COMM_WORLD);

        if (rank == 0)
            printf("pi is approximately %.16f, Error is %.16f\n",
                   pi, fabs(pi - PI25DT));

    }

    MPI_Finalize();

}


