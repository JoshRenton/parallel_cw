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

// Get the start and end indices for this process
void get_start_and_end_indices(int* start_index, int* end_index, int array_length) {
    // Get the number of processes in the cluster
    int size_of_cluster;
    MPI_Comm_size(MPI_COMM_WORLD, &size_of_cluster);

    // Get the unique rank of this process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Find the excess number of values
    int excess_values = array_length % size_of_cluster;

    // Find the number of values to calculate per process
    int values_per_process = (array_length - excess_values) / size_of_cluster;

    /*
        Calculate the start and end index for this process, based
        on the rank.
        The excess values are divided up equally over the processes.
    */
    if (rank < excess_values) {
        *start_index = rank * values_per_process + rank;
        *end_index = *start_index + values_per_process + 1;
    }
    else {
        *start_index = rank * values_per_process + excess_values;
        *end_index = *start_index + values_per_process;
    }
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

    /*
        Create the test array of doubles.
        Each process creates the test array individually.
        The test array will still be identical for all processes.
    */
    double* test_array = malloc(array_length * sizeof(double));
    create_test_array(test_array, array_length);

    /*
        Each process works on a set amount of values, based
        on its rank.
    */
    int start_index;
    int end_index;
    get_start_and_end_indices(&start_index, &end_index, array_length);

    printf("\n%d, %d\n", start_index, end_index);

    // Close the MPI environment
    MPI_Finalize();
    return 0;
}
