#include "prog.h"

// Create a test array of the given size with a fixed pattern
void create_test_array(double* test_array, int array_length) {
    for (int i = 0; i < array_length; i++) {
        if (i % 7 == 0) {
            test_array[i] = 0.4;
        }
        else if (i % 5 == 0) {
            test_array[i] = 0.8;
        }
        else if (i % 3 == 0) {
            test_array[i] = 0.5;
        }
        else if (i % 2 == 0) {
            test_array[i] = 0.2;
        }
    }  
}

// Prints the provided array of doubles
void print_double_array(double* array, int size) {
    printf("[");
    for (int row = 0; row < size; row ++) {
        for (int column = 0; column < size; column++) {
            printf("%f, ", array[row * size + column]);
        }
        printf("\n");
    }
    printf("]");
}

// Get the start and end row indices for this process
void get_start_and_end_row_indices(int* start_row, int* end_row, int rank, int process_count, int size) {
    // Find the excess number of rows
    int excess_rows = size % process_count;

    // Find the number of rows to calculate per process
    int rows_per_process = (size - excess_rows) / process_count;

    /*
        Calculate the start and end row indices for this process, based
        on the rank.
        The excess rows are divided up as equally as
        possible over the processes.
    */
    if (rank < excess_rows) {
        *start_row = rank * rows_per_process + rank;
        *end_row = *start_row + rows_per_process + 1;
    }
    else {
        *start_row = rank * rows_per_process + excess_rows;
        *end_row = *start_row + rows_per_process;
    }
}

// Send the designated rows to a process
void send_rows(double* test_array, int start_row, int end_row) {
    
}

// Average the four values surrounding the input index
double average(double* test_array, int size, int index) {
    double sum = test_array[index + 1] + test_array[index - 1] + test_array[index + size] + test_array[index - size];
    return sum / 4.0;
}

int main(int argc, char** argv)
{
    // Initialise the MPI environment
    MPI_Init(&argc, &argv);

    // The size of the array is passed in as a command line argument
    int size = atoi(argv[1]);
    int array_length = size * size;

    // Get the unique rank of this process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Process 0 sets up the test array and distributes the workload
    if (rank == 0) {
        // Create the test array
        double* test_array = malloc(array_length * sizeof(double));
        create_test_array(test_array, array_length);

        // get the number of processes in the environment
        int process_count;
        MPI_Comm_size(MPI_COMM_WORLD, &process_count);

        // Each process works on set rows based on its rank
        for (int process = 0; process < process_count; process++) {
            int start_row;
            int end_row;
            get_start_and_end_row_indices(&start_row, &end_row, process, process_count, size);
            
        }
    }

    // Close the MPI environment
    MPI_Finalize();
    return 0;
}
