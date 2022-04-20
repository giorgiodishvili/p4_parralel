#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
//#include <opencl-c-base.h>

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

    int size, rank, dst;

    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size == 1) return;

    int index = 1;

    if (rank != 0) {
//        printf("Rank %d, step %d \n", rank, step);
        MPI_Recv(buffer, count, datatype, MPI_ANY_SOURCE, 0, comm, &status); // receive n
        MPI_Recv(&index, count, datatype, MPI_ANY_SOURCE, 0, comm, &status); // receive index
//        printf("Rank %d, index %d \n", rank, index);
    }
//    int step = pow(2, index - 1) - rank;
//    int step = 1;
//    while (step < (size-index)){
//        step<<=1;
//    }
//    printf("after while Rank %d, step %d \n", rank, step);


    int step = size - rank - 1;
//    int step = size - (size / index) + 1 ;
//    int step = pow(2, index - 1);
    // 0 -> 1
    // 1 -> 3. not 0 -> 2
    // 3 -> 7, 2 -> 6
//    if (rank == 7) {
//        printf("im 7777 startingg with step %d \n", step);
//    }
    int iterations = index - 1;

    while (step > 0 && step <= size) {

        long long power = pow(2, iterations); //
        iterations++;
        dst = rank + power;
        if (dst > 0 && dst < size && dst > rank ) {
//            printf("Before Send Rank %d,step %d, dst %d, sent %d \n", rank, step, dst, *(int *) buffer);
            MPI_Send(buffer, count, datatype, dst, 0, comm);
            index++;
            MPI_Send(&index, count, datatype, dst, 0, comm);
            printf("After Send Rank %d, step %d, dst %d, index %d \n", rank, step, dst, index);
        }
        step >>= 1;

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
        printf("Recevied n %d, rank %d \n", n, rank);
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


