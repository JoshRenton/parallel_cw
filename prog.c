#include "prog.h"

// Currently broken when using 2 processes as the process should not get either the row above or the row below

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
void get_start_and_end_indices(int* start_index, int* end_index, int rank, int process_count, int size) {
    // Decrement rank and process_count by 1 as process 0 does not compute values
    process_count--;
    rank--;

    // Find the excess number of rows
    int excess_rows = size % process_count;

    // Find the number of rows to calculate per process
    int rows_per_process = (size - excess_rows) / process_count;

    int start_row;
    int end_row;

    /*
        Calculate the start and end row indices for this process, based
        on the rank.
        The excess rows are divided up as equally as
        possible over the processes.
    */
    if (rank < excess_rows) {
        start_row = rank * rows_per_process + rank;
        end_row = start_row + rows_per_process + 1;
    }
    else {
        start_row = rank * rows_per_process + excess_rows;
        end_row = start_row + rows_per_process;
    }

    // Convert row indices to element indices
    *start_index = start_row * size;
    *end_index = end_row * size;
}

// Send the designated rows to a process
void send_rows(double* test_array, int start_index, int end_index, int rank) {
    int sent_array_length = end_index - start_index;
    double* sent_array = malloc(sent_array_length * sizeof(double));

    // Copy the needed row values to sent_array
    for (int i = 0; i < sent_array_length; i++) {
        sent_array[i] = test_array[start_index + i];
    }

    // Send the data values to the corresponding process
    printf("\nSending rows to process %d", rank);
    if (MPI_Send(sent_array, sent_array_length, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD) != 0) {
        printf("\nFailed to send rows to process %d\n", rank);
    }
    free(sent_array);
}

// Average the four values surrounding the input index
double average(double* test_array, int size, int index) {
    double sum = test_array[index + 1] + test_array[index - 1] + test_array[index + size] + test_array[index - size];
    return sum / 4.0;
}

/* 
    Send a row of values to one process and recieve a row of values
    from another process. The recieved row of values replaces
    the sent row of values.
*/
void get_adjacent_row(int send_to_rank, int recieve_from_rank, double* buffer, int buffer_size) {
    MPI_Sendrecv_replace(buffer, buffer_size, MPI_DOUBLE, send_to_rank, 1, 
    recieve_from_rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

int main(int argc, char** argv)
{
    // Initialise the MPI environment
    if (MPI_Init(&argc, &argv) != 0) {
        exit(1);
    }

    // The size of the array is passed in as a command line argument
    int size = atoi(argv[1]);
    int array_length = size * size;

    // Get the unique rank of this process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the number of processes in the environment
    int process_count;
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);

    // Process 0 sets up the test array and distributes the workload
    if (rank == 0) {
        // Create the test array
        double* test_array = malloc(array_length * sizeof(double));
        create_test_array(test_array, array_length);

        // Each process works on set rows based on its rank
        for (int process = 1; process < process_count; process++) {
            int start_index;
            int end_index;
            // Get the start and end indices of the values for this process
            get_start_and_end_indices(&start_index, &end_index, process, process_count, size);
            // Send the initial rows of values to this process
            send_rows(test_array, start_index, end_index, process);
        }

        free(test_array);
    }
    // All processes other than process 0
    else {
        // Check for an incoming message from process 0
        MPI_Status status;
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // Get the size of the message
        int buffer_size;
        MPI_Get_count(&status, MPI_DOUBLE, &buffer_size);
        printf("\n%d\n", buffer_size);

        double* values = malloc(buffer_size * sizeof(double));

        // Recieve the message containing the rows from process 0
        MPI_Recv(values, buffer_size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double* row_above;
        double* row_below;

        // Process 1 does not have a row above
        if (rank == 1) {
            row_below = malloc(size * sizeof(double));

            // Copy the top row of values into row_below
            for (int i = 0; i < size; i++) {
                row_below[i] = values[i];
            }
            
            // Only communicate with process 1 rank higher
            get_adjacent_row(rank + 1, rank + 1, row_below, size);
            printf("\nProcess %d recieved values from process %d\n", rank, rank + 1);
            free(row_below);
        }
        // Process with highest rank has no row below
        else if (rank == process_count - 1) {
            row_above = malloc(size * sizeof(double));

            // Copy the bottom row of values into row_above
            for (int i = buffer_size - size; i < buffer_size; i++) {
                row_above[i - (buffer_size - size)] = values[i];
            }

            // Only communicate with process 1 rank lower
            get_adjacent_row(rank - 1, rank - 1, row_above, size);
            printf("\nProcess %d recieved values from process %d\n", rank, rank - 1);
            free(row_above);
        }
        else {
            row_above = malloc(size * sizeof(double));
            row_below = malloc(size * sizeof(double));

            // Copy the top row of values into row_above
            for (int i = buffer_size - size; i < buffer_size; i++) {
                row_above[i - (buffer_size - size)] = values[i];
            }
            /*
                Send bottom row to process 1 rank higher.
                Recieve row above from process 1 rank lower.
            */
            get_adjacent_row(rank + 1, rank - 1, row_above, size);
            printf("\nProcess %d recieved values from process %d\n", rank, rank - 1);

            // Copy the bottom row of values into row_below
            for (int i = 0; i < size; i++) {
                row_below[i] = values[i];
            }
            /*
                Send top row to process 1 rank lower.
                Recieve row below from process 1 rank higher.
            */
            get_adjacent_row(rank - 1, rank + 1, row_below, size);
            printf("\nProcess %d recieved values from process %d\n", rank, rank + 1);
            free(row_above);
            free(row_below);
        }

        free(values);
    }

    MPI_Finalize();
    return 0;
}
