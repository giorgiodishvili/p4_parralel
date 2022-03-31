#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[]) {

    int i, done = 0, n, count;
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

    // Pack these values together into a string
    int buffer_len = 150;
    char buffer[buffer_len];
    // If not rank zero, send your message to be printed
    while (!done) {
        if (rank == 0) {
            // Always print from rank 0
            printf("Enter the number of points: (0 quits) \n");
            scanf("%d", &n);

            for (int i = 1; i < size; i++) {
                // Takes buffer, size, type, source, tag, communicator, and status
                MPI_Send(&n, buffer_len, MPI_INT, i, i, MPI_COMM_WORLD);
            }
            // Takes buffer, size, type, source, tag, communicator, and status
            MPI_Recv(buffer, buffer_len, MPI_CHAR, 0, MPI_ANY_TAG,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Print our received message
            printf("%s\n", buffer);
        } else {
            printf("here %d", rank);

            MPI_Recv(&n, buffer_len, MPI_INT, rank, rank,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("here");
            if (n == 0) break;

            count = 0;

            for (i = 1; i <= n; i++) {
                // Get the random numbers between 0 and 1
                x = ((double) rand()) / ((double) RAND_MAX);
                y = ((double) rand()) / ((double) RAND_MAX);

                // Calculate the square root of the squares
                z = sqrt((x * x) + (y * y));

                // Check whether z is within the circle
                if (z <= 1.0)
                    count++;
            }
            pi = ((double) count / (double) n) * 4.0;
            sprintf(buffer, "pi is approx. %.16f, Error is %.16f\n", pi, fabs(pi - PI25DT));
            MPI_Send(buffer, buffer_len, MPI_CHAR, 0, rank, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();

}

//
//#include <stdio.h>
//#include <mpi.h>
//int main(int argc, char **argv) {
//    int ierr, num_procs, my_id;
//
//    /* find out MY process ID, and how many processes were started. */
//
//    ierr = MPI_Init(&argc, &argv);
//
//    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
//    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
//
//    printf("Hello world! I'm process %i out of %i processes\n",
//           my_id, num_procs);
//
//    ierr = MPI_Finalize();
//}