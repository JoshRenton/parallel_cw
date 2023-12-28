#include "prog.h"

// Returns the absolute value of a double
double absolute(double value) {
    if (value < 0.0) {
        return 0.0 - value;
    }
    else {
        return value;
    }
}

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
        else {
            test_array[i] = 0.0;
        }
    }  
}

// Prints the provided array of doubles
void print_double_array(double* array, int width, int height) {
    printf("\n[");
    for (int row = 0; row < height; row ++) {
        for (int column = 0; column < width; column++) {
            printf("%f, ", array[row * width + column]);
        }
        printf("\n");
    }
    printf("]\n");
}

// Get the start and end indices for this process
void get_start_and_end_indices(int* start_index, int* end_index, int rank, int process_count, int size) {
    // Decrement rank and process_count by 1 as process 0 does not compute values
    process_count--;
    rank--;

    // Find the excess number of rows
    int excess_rows = (size - 2) % process_count;

    // Find the number of rows to calculate per process
    int rows_per_process = ((size - 2) - excess_rows) / process_count;

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

    // Shift the start and end rows up by one to ignore the first and last rows.
    start_row++;
    end_row++;

    // Convert row indices to element indices
    *start_index = start_row * size;
    *end_index = end_row * size;
}

// Send the given rows to the given process
void send_rows(double* sent_array, int sent_array_length, int rank) {
    // Send the data values to the corresponding process
    printf("\nSending rows to process %d", rank);
    if (MPI_Send(sent_array, sent_array_length, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD) != 0) {
        printf("\nFailed to send rows to process %d\n", rank);
    }
    free(sent_array);
}

// Average the given values and return if within the given precision
int average_values(double* values, double* row_above, double* row_below, int size, int buffer_size, double* precision) {
    /*
        Stores the previous values of the row above
        so that the averaged values do not affect the
        the calculation of neighbouring values.
    */
    double* prev_row_state = malloc((size - 2) * sizeof(double));

    int rows = buffer_size / size;

    int within_precision = 1;

    for (int row = 0; row < rows; row++) {
        for (int column = 1; column < size - 1; column++) {
            int index = row * size + column;
            double sum;
            // The case of the process only having one row to work on
            if (rows == 1) {
                if (column == 1) {
                    sum = values[index - 1] + values[index + 1] + row_above[column] + row_below[column];
                }
                else {
                    sum = prev_row_state[column - 2] + values[index + 1] + row_above[column] + row_below[column];
                }
            }
            else if (row == 0) {
                if (column == 1) {
                    sum = values[index - 1] + values[index + 1] + row_above[column] + values[index + size];
                }
                else {
                    sum = prev_row_state[column - 2] + values[index + 1] + row_above[column] + values[index + size];
                }
            }
            else if (row == rows - 1) {
                if (column == 1) {
                    sum = values[index - 1] + values[index + 1] + prev_row_state[column - 1] + row_below[column];
                }
                else {
                    sum = prev_row_state[column - 2] + values[index + 1] + prev_row_state[column - 1] + row_below[column];
                }
            }
            else {
                if (column == 1) {
                    sum = values[index - 1] + values[index + 1] + prev_row_state[column - 1] + values[index + size];
                }
                else {
                    sum = prev_row_state[column - 2] + values[index + 1] + prev_row_state[column - 1] + values[index + size];
                }
            }
            // Store current value before updating
            prev_row_state[column - 1] = values[index];

            // Update the old value with the new averaged value
            values[index] = sum / 4.0;

            // Check if within the required precision
            if (within_precision != 0) {
                if (absolute(values[index] - prev_row_state[column - 1]) > *precision) {
                    within_precision = 0;
                }
            }
        }
    }
    free(prev_row_state);
    return within_precision;
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

// Get the top row from process 0
void get_top_row(double* row_above, int size) {
    MPI_Recv(row_above, size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// Get the bottom row from process 0
void get_bottom_row(double* row_below, int size) {
    MPI_Recv(row_below, size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

/*
    For process 0.
    Recieves the final values from a process and
    writes them to the test array.
*/
void recieve_and_write(double* test_array, int process, int* current_index) {
    // Check for an incoming message from process
    MPI_Status status;
    MPI_Probe(process, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    // Get the size of the message
    int buffer_size;
    MPI_Get_count(&status, MPI_DOUBLE, &buffer_size);

    double* values = malloc(buffer_size * sizeof(double));

    // Recieve the message containing the rows from process
    MPI_Recv(values, buffer_size, MPI_DOUBLE, process, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Write averaged values to test array
    for (int i = *current_index; i < *current_index + buffer_size; i++) {
        test_array[i] = values[i - *current_index];
    }

    *current_index += buffer_size;

    printf("\nRecieved values from process %d up to index %d\n", process, *current_index);

    free(values);
}

// Check if the given process is within the required precision, and update within_precision if necessary
void check_precision(int rank, int* within_precision) {
    int* buffer = malloc(sizeof(int));
    MPI_Recv(buffer, 1, MPI_INT, rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // If within_precision is not 0, set it to 0 to indicate that another iteration is needed
    if (*buffer == 0 && *within_precision != 0) {
        *within_precision = 0;
    }
    free(buffer);
}

int main(int argc, char** argv)
{
    // Initialise the MPI environment
    if (MPI_Init(&argc, &argv) != 0) {
        exit(1);
    }

    // The size of the array is passed in as the first command line argument
    int size = atoi(argv[1]);
    int array_length = size * size;

    // The precision is passed in as the second command line argument
    double precision = strtod(argv[2], NULL);

    // Get the unique rank of this process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the number of processes in the environment
    int process_count;
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);

    // This will be turned to 0 to indicate that a process is not below the required precision
    int within_precision;

    // Process 0 sets up the test array and distributes the workload among the other processes
    if (rank == 0) {
        // Create the test array
        double* test_array = malloc(array_length * sizeof(double));
        create_test_array(test_array, array_length);

        // Each process works on an (almost) equal number of consecutive rows based on its rank
        for (int process = 1; process < process_count; process++) {
            int start_index;
            int end_index;
            // Get the start and end indices of the values for this process
            get_start_and_end_indices(&start_index, &end_index, process, process_count, size);
            
            int sent_array_length = end_index - start_index;
            double* sent_array = malloc(sent_array_length * sizeof(double));

            // Copy the needed row values to sent_array
            for (int i = 0; i < sent_array_length; i++) {
                sent_array[i] = test_array[start_index + i];
            }
            // Send the initial rows of values to this process
            send_rows(sent_array, sent_array_length, process);
        }

        // Send top row of array to process 1
        double* top_row = malloc(size * sizeof(double));
        for (int i = 0; i < size; i++) {
            top_row[i] = test_array[i];
        }
        MPI_Send(top_row, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        free(top_row);

        // Send bottom row of array to highest ranking process
        double* bottom_row = malloc(size * sizeof(double));
        for (int i = array_length - size; i < array_length; i++) {
            bottom_row[i - (array_length - size)] = test_array[i];
        }
        MPI_Send(bottom_row, size, MPI_DOUBLE, process_count - 1, 0, MPI_COMM_WORLD);
        free(bottom_row);

        do {
            within_precision = 1;
            // Check if any process is not within the required precision
            for (int process = 1; process < process_count; process++) {
                check_precision(process, &within_precision);
            }
            // Communicate to all other processes if they should continue or not
            for (int process = 1; process < process_count; process++) {
                MPI_Send(&within_precision, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            }
        }
        while (within_precision != 1);

        int current_index = size;
        // Recieve all averaged values from each process
        for (int process = 1; process < process_count; process++) {
            recieve_and_write(test_array, process, &current_index);
        }

        print_double_array(test_array, size, size);

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

        double* values = malloc(buffer_size * sizeof(double));

        // Recieve the message containing the rows to work on from process 0
        MPI_Recv(values, buffer_size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double* row_above;
        double* row_below;

        // Process 1 gets the top row from process 0
        if (rank == 1) {
            // Get the top row of array from process 0, this is only done once as it will not change
            row_above = malloc(size * sizeof(double));
            get_top_row(row_above, size);
            // printf("\nProcess %d recieved top row\n", rank);

            row_below = malloc(size * sizeof(double));
            // If there are only 2 process, get the bottom row from process 0
            if (process_count == 2) {
                get_bottom_row(row_below, size);
            }

            int stop;

            do {
                // Copy the bottom row of values into row_below
                if (process_count != 2) {
                    for (int i = buffer_size - size; i < buffer_size; i++) {
                        row_below[i - (buffer_size - size)] = values[i];
                    }  
                    // Exchange bottom row with top row of process 1 rank higher
                    get_adjacent_row(rank + 1, rank + 1, row_below, size);
                }
                // printf("\nProcess %d recieved values from process %d\n", rank, rank + 1);

                stop = average_values(values, row_above, row_below, size, buffer_size, &precision);

                // Send value of stop to process 0 and check if this should stop
                MPI_Sendrecv_replace(&stop, 1, MPI_INT, 0, 0, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            while (stop != 1);

            // Send resulting values back to process 0
            send_rows(values, buffer_size, 0);

            // printf("\n%d\n", rank);
            // print_double_array(values, size, buffer_size / size);
        }
        // The highest ranking process gets the bottom row from process 0
        else if (rank == process_count - 1) {
            // Get the bottom row of array from process 0, this is only done once as it will not change
            row_below = malloc(size * sizeof(double));
            get_bottom_row(row_below, size);
            // printf("\nProcess %d recieved bottom row\n", rank);

            row_above = malloc(size * sizeof(double));

            int stop;

            do {
                // Copy the top row of values into row_above
                for (int i = 0; i < size; i++) {
                    row_above[i] = values[i];
                }

                // Exchange top row with bottom row of process 1 rank lower
                get_adjacent_row(rank - 1, rank - 1, row_above, size);
                // printf("\nProcess %d recieved values from process %d\n", rank, rank - 1);

                stop = average_values(values, row_above, row_below, size, buffer_size, &precision);

                // Send value of stop to process 0 and check if this should stop
                MPI_Sendrecv_replace(&stop, 1, MPI_INT, 0, 0, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            while (stop != 1);

            // Send resulting values back to process 0
            send_rows(values, buffer_size, 0);

            // printf("\n%d\n", rank);
            // print_double_array(values, size, buffer_size / size);
        }
        else {
            row_above = malloc(size * sizeof(double));
            row_below = malloc(size * sizeof(double));

            int stop;

            do {
                // Copy the bottom row of values into row_above
                for (int i = buffer_size - size; i < buffer_size; i++) {
                    row_above[i - (buffer_size - size)] = values[i];
                }
                /*
                    Send bottom row to process 1 rank higher.
                    Recieve row above from process 1 rank lower.
                */
                get_adjacent_row(rank + 1, rank - 1, row_above, size);
                // printf("\nProcess %d recieved values from process %d\n", rank, rank - 1);

                // Copy the top row of values into row_below
                for (int i = 0; i < size; i++) {
                    row_below[i] = values[i];
                }
                /*
                    Send top row to process 1 rank lower.
                    Recieve row below from process 1 rank higher.
                */
                get_adjacent_row(rank - 1, rank + 1, row_below, size);
                // printf("\nProcess %d recieved values from process %d\n", rank, rank + 1);

                stop = average_values(values, row_above, row_below, size, buffer_size, &precision);

                // Send value of stop to process 0 and check if this should stop
                MPI_Sendrecv_replace(&stop, 1, MPI_INT, 0, 0, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            while (stop != 1);

            // Send resulting values back to process 0
            send_rows(values, buffer_size, 0);

            // printf("\n%d\n", rank);
            // print_double_array(values, size, buffer_size / size);
        }
        // Free allocated memory
        free(row_below);
        free(row_above);
    }

    // Close MPI environment
    MPI_Finalize();
    return 0;
}
