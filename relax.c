#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int upper_bound(int rank, int process_count, int length) {
    if((length - 2) % process_count == 0) {
        return (rank + 1) * ((length - 2) / process_count);
    }
	else if(process_count > (length - 2) / 2) {
		if(rank < (length - 2) % process_count) {
			return 2 * (rank + 1);
		}
		else {
			return 2 * ((length - 2) % process_count) + (1 + rank - ((length - 2) % process_count));
		}
	}
    else if(((length - 2) % process_count) - rank > 0) {
        return (rank + 1) * (((length - 2) / process_count) + 1);
    }
    else {
        return ((rank + 1) * ((length - 2) / process_count) + ((length - 2) % process_count));
    }
}

void init_complex(double *line, int length, int r) {
	for(int c = 0; c < length; ++c) {
		if(r%2==0){
			if(c%2==0){
				line[c] = 1;
			}
			else {
				line[c] = 0;
			}
		}
		else {
			if(c%2==0){
				line[c] = 0;
			}
			else {
				line[c] = 1;
			}
		}
	}
}

void init_simple(double *line, int length, int r) {
	for(int c = 0; c < length; ++c) {
		if(r==0 || c == 0){
			line[c] = 1;
		}
        else {
            line[c] = 0;
        }
	}
}

void print_sub_matrix(double **sub_matrix, int length, int width) {

	int r, c;

	for(r = 0; r < width; ++r){
		for(c = 0; c < length; ++c) {
			printf("%f ", sub_matrix[r][c]);
		}
		printf("\r\n");
	}
}

void print_line(double *line, int length) {
	for(int c = 0; c < length; ++c) {
		printf("%f ", line[c]);
	}
	printf("\r\n");
}

void copy_buffer(double *line, double *buffer, int length) {
    for(int i = 1; i < length - 1; ++i) line[i] = buffer[i-1];
}

int main(int argc, char **argv) {
    int rc = MPI_Init(&argc, &argv);

    // Calculating
    double start, finish;
    start = MPI_Wtime();


    // Edit length to change the dimensions of the array. Dimensions are length x length.
    int length = (44 * 4) + 2;
    char is_simple = 1;
    // Allows you to choose the matrix to work with. If 1 then the example matrix within the spec will be used, otherwise the following will be used.
    /*
    0 1 0 1 0 ...
    1 0 1 0 1 ...
    0 1 0 1 0 ...
    1 0 1 0 1 ...
    0 1 0 1 0 ...
    ...
    */
    // The former tends to be better for testing larger matrices, the latter is better for testing matrices which require large amounts of iterations
    
    // Check if initializing MPI was succesful
    if(rc != MPI_SUCCESS) {
        // Unsucessful, abort.
        printf("Failed to start MPI program\r\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    

    /* Information about the number of threads, the rank number for this thread, and the name of the node.
       rank and process_count are used in calculating the section sizes, and the latter is used in determining if the matrix is in an end state.
       A sum of the number of threads ready to finish is taken, if that sum is equal to the process_count then iterations stop.
    */
    int rank, process_count, name_length;
    char name[MPI_MAX_PROCESSOR_NAME];
    name_length = MPI_MAX_PROCESSOR_NAME;

    MPI_Get_processor_name(name, &name_length);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    
    if((length - 2) / process_count == 0) {
        MPI_Finalize();
        printf("The amount of threads is greater than available indexes\n");
        return 1;
    }

    /* Creates indexes in the effective matrix (the matrix being processd without the top or bottom row)
       These will be used for iterating through the sub matrix. The width is essentially the length for the
       first index.
    */
    int lower, upper, width;
    lower = upper_bound(rank - 1, process_count, length);
    upper = upper_bound(rank, process_count, length);
    width = upper - lower;

    
    double **sub_matrix, **outer_values;
    sub_matrix = malloc(width * sizeof(*sub_matrix));
    outer_values = malloc(2 * sizeof(*outer_values));
    outer_values[0] = malloc(length * sizeof(outer_values[0]));
    outer_values[1] = malloc(length * sizeof(outer_values[1]));

    // Starts initialising values in the sub matrix and the outer values. This is done in the order where lines appear in the matrix.
    for(int r = 0; r < width + 2; ++r) {
        // Is r within the submatrix's index range?
        if(r < width ) {
            // Yes, allocate rows for the matrix.
            sub_matrix[r] = malloc(length * sizeof(*sub_matrix[r]));
        }
        for(int c = 0; c < length; ++c) {
            // Is the matrix simple?
            if(is_simple == 1) {
                // Yes. Initialize values.
                // Is it the first line?
                if(r == 0)
                {
                    // Yes, handle the outer values.
                    init_simple(outer_values[0], length, r + lower);
                } 
                // Is it the last line?
                else if(r == width + 1) {
                    // Yes, handle the outer values
                    init_simple(outer_values[1], length, r + lower);
                }
                else {
                    // It's all the values inbetween, initialize them.
                    init_simple(sub_matrix[r-1], length, lower + r);
                }
            }
            else {
                // No. It's complex.
                // Is it the first line?
                if(r == 0)
                {
                    // Yes, initialize the outer values
                    init_complex(outer_values[0], length, r + lower);
                } 
                else if(r == width + 1) {
                    // Yes, initialize the outer values
                    init_complex(outer_values[1], length, r + lower);
                }
                else {
                    // It's all the values inbetween, initialize them.
                    init_complex(sub_matrix[r-1], length, lower + r);
                }
            }
            
        }
    }

    // Buffer used to store new values temporarily before they're added into their relavent location in the sub matrix.
    double **buffer;
    buffer = malloc(2 * sizeof(*buffer));
    buffer[0] = malloc((length - 2) * sizeof(buffer[0]));
    buffer[1] = malloc((length - 2) * sizeof(buffer[1]));

    // Count stores the sum of all stop values across threads. It's value is available in all threads. Once count is the same as process_count relaxation finishes.
    int count = 0;
    // Used to indicate if all values within the sub matrix are under the maximum difference.
    int stop = 0;
    // The difference worked towards.
    double max_diff = 0.0001;

    // Declaring that the thread is ready to start relaxation, stating where it's working on alongside its rank and node name.
    printf("Node: %s\nRank: %d\n\tStarting at %d\n\tEnding at %d\n\n======================\n", name, rank, lower, upper);
    MPI_Status stat;
    while(count < 44) {
        stop = 1;
        // Is the submatrix a line?
        if(width != 1){
            // No, treat as a submatrix of (x > 2) x length.
            for(int r = 0; r < width; ++r) {
                // Is the buffer ready to be copied to the submatrix?
                if(r > 1) {
                    // Yes. Copy the buffer to the corresponding row.
                    copy_buffer(sub_matrix[r-2], buffer[r%2], length);
                }
                for(int c = 1; c < length - 1; ++c) {
                    // Does the value involve the last outer value row?
                    if(r == width - 1){
                        // Yes, use it in the calculation.
                        buffer[r%2][c-1] = (outer_values[1][c] + sub_matrix[r - 1][c] + sub_matrix[r][c + 1] + sub_matrix[r][c - 1]) / 4;
                    // Does the value involve the first outer value row?
                    } else if(r == 0) {
                        // Yes, use it in the calculation.
                        buffer[r%2][c-1] = (outer_values[0][c] + sub_matrix[r + 1][c] + sub_matrix[r][c + 1] + sub_matrix[r][c - 1]) / 4;
                    } else {
                        // No outer values used, calculate normally.
                        buffer[r%2][c-1] = (sub_matrix[r + 1][c] + sub_matrix[r - 1][c] + sub_matrix[r][c + 1] + sub_matrix[r][c - 1]) / 4;
                    }
                    // Is the difference too big between iterations to finish?
                    if(fabs(buffer[r%2][c-1] - sub_matrix[r][c]) > max_diff) {
                        // Yes, set stop to 0.
                        stop = 0;
                    }
                }
            }
            // Copy the remaining rows to the sub matrix.
            copy_buffer(sub_matrix[width-2], buffer[(width-2)%2], length);
            copy_buffer(sub_matrix[width-1], buffer[(width-1)%2], length);
        }
        else {
            // It's a single row sub matrix.
            for(int c = 1; c < length - 1; ++c) {
                // Calculate the new value
                buffer[0][c] = (outer_values[1][c] + outer_values[0][c] + sub_matrix[0][c + 1] + sub_matrix[0][c - 1]) / 4;
                // Is the difference too big between iterations to finish?
                if(fabs(sub_matrix[0][c] - buffer[0][c]) > max_diff) {
                    // Yes, set stop to 0
                    stop = 0;
                }
            }
            // Copy the buffer to the row. Yes I know this is lazy.
            copy_buffer(sub_matrix[0], buffer[0], length);
        }
        
        // Start communicating with other ranks. Get relevant rows to replace outer values.

        // Is this the first rank?
        if(rank == 0) {
            // Yes. Just send the last row, and receive their first row.
            MPI_Sendrecv(sub_matrix[width - 1], length, MPI_DOUBLE, rank + 1, rank + 1, outer_values[1], length, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD, &stat);
        }
        // Is this the last rank?
        else if(rank == process_count - 1) {
            // Yes. Just send the first row, and receive their last row.
            MPI_Sendrecv(sub_matrix[0], length, MPI_DOUBLE, rank - 1, rank - 1, outer_values[0], length, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD, &stat);
        }
        else {
            // It's all the ranks in between.
            // Is the rank even?
            if(rank % 2 == 0) {
                // Send last row and replace the last outer values
                MPI_Sendrecv(sub_matrix[width - 1], length, MPI_DOUBLE, rank + 1, rank + 1, outer_values[1], length, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD, &stat);
                // Send first row and replace the first outer values
                MPI_Sendrecv(sub_matrix[0], length, MPI_DOUBLE, rank - 1, rank - 1, outer_values[0], length, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD, &stat);
            }
            else {
                // Send first row and replace the first outer values
                MPI_Sendrecv(sub_matrix[0], length, MPI_DOUBLE, rank - 1, rank - 1, outer_values[0], length, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD, &stat);
                // Send last row and replace the last outer values
                MPI_Sendrecv(sub_matrix[width - 1], length, MPI_DOUBLE, rank + 1, rank + 1, outer_values[1], length, MPI_DOUBLE, rank + 1, rank, MPI_COMM_WORLD, &stat);
            }
            
        }
        // Get all threads to wait before performing the coming calculation.
        MPI_Barrier(MPI_COMM_WORLD);
        // Get the sum of stop values and compare to the loop condition. If all are ready to stop exit.
        MPI_Allreduce(&stop, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    // Start cleaning up.
    // Release buffer memory, no longer needed.
    free(buffer[0]);
    free(buffer[1]);
    free(buffer);

    // Remove comments below to print values.

    /*
    if(rank == 0) {
        for(int i = 0; i < process_count; ++i) {
            if(i == 0){
                print_line(outer_values[0], length);
                print_sub_matrix(sub_matrix, length, width);
            }
            else if(i == process_count - 1) {
                MPI_Recv(&(sub_matrix[0][0]), length * width, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
                MPI_Recv(&(outer_values[0][0]), length * 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
                print_sub_matrix(sub_matrix, length, width);
                print_line(outer_values[1], length);
            }
            else {
                MPI_Recv(&(sub_matrix[0][0]), length * width, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
                print_sub_matrix(sub_matrix, length, width);
            }
        }
    }
    else if(rank == process_count - 1) {
        MPI_Send(&(sub_matrix[0][0]), length * width, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&(outer_values[0][0]), length * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Send(&(sub_matrix[0][0]), length * width, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    */

    // Release everything else, no longer needed.
    for(int r = 0; r < width; ++r) {
        free(sub_matrix[r]);
    }
    free(sub_matrix);
    free(outer_values[0]);
    free(outer_values[1]);
    free(outer_values);

    finish = MPI_Wtime();
    
    MPI_Finalize();

    // Get the first rank to print the time taken.
    if(rank == 0) {
        double final_time = finish - start;
        printf("Node: %s\nRank: %d\n\tTime taken %lf\n\n================================\n", name, rank, final_time);
    }
    return 0;
}