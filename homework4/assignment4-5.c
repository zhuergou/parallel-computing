 /***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Yuchen Liang && Diya LI       *******************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
//#include</home/diya/codes/assignment4/clcg4.h>
#include</gpfs/u/home/PCP7/PCP7lngc/barn/homework4/clcg4.h>
#include<pthread.h>
#include<mpi.h>


/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define BGQ 1 // when running BG/Q, comment out when running on kratos
#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime
#endif


#define ALIVE 1
#define DEAD  0
#define UNIVERSE 16384
#define UNIVERSE_test 10
#define CHILDREN 63
#define CHILDREN_test 255
#define TICK 20
#define TICK 128
#define THRESHOLD 0.25

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

// You define these

int semaphore1 = 0;
int semaphore2 = 0;
int semaphore3 = 0;

typedef struct package{
	int start;
	int row;
	int column;
	int ** cell;
	int ** cell2;
	int ** count;
	int * ghostup;
	int * ghostdown;
	int mpi_myrank;
	int row_per_rank;

} variable;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these

int update_counter(int **cell, int j, int m, int row, int column, int *ghostup, int *ghostdown);

void * whattodo(void *arg);

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[]){

	int i;
	int j;
	int m;
	int mpi_myrank;
	int mpi_commsize;
	int rc;
	int zero = 0;
	int one = 0;
	int p;
	int q;
	pthread_t tid[CHILDREN];
	
// Example MPI startup and using CLCG4 RNG
	MPI_Init( &argc, &argv);
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

	int row_per_rank = UNIVERSE/mpi_commsize;


    double time_in_secs = 0, time_in_secs2 = 0;
    double processor_frequency = 1600000000.0;
    unsigned long long start_cycles=0, start_cycles2 = 0;
    unsigned long long end_cycles=0, end_cycles2 = 0;

    start_cycles = GetTimeBase();

	// board for each rank
	int **cell = malloc(row_per_rank *sizeof(int *));
	for (i = 0; i < row_per_rank; i++){
		cell[i] = malloc(UNIVERSE *sizeof(int));
	}

	// counting borad for each rank
	int **count =  malloc(row_per_rank *sizeof(int *));
	for (i = 0; i < row_per_rank; i++){
		count[i] = malloc(UNIVERSE *sizeof(int));
	}

	// another board
	int **cell2 =  malloc(row_per_rank *sizeof(int *));
	for (i = 0; i < row_per_rank; i++){
		cell2[i] = malloc(UNIVERSE *sizeof(int));
	}

	// ghost row for each rank
	int ghostup[UNIVERSE];
	int ghostdown[UNIVERSE];

	MPI_Request myRequest;

// Init 16,384 RNG streams - each rank has an independent stream


	InitDefault();


// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.


	printf("Rank %d of %d has been started and a first Random Value of %lf\n", mpi_myrank, mpi_commsize, GenVal(mpi_myrank));

// Insert your code

//	initial value
	for (i = 0; i < row_per_rank; i++){
		for(j = 0; j < UNIVERSE; j++){
			if (GenVal(mpi_myrank * row_per_rank + i) > 0.5){
				cell[i][j] = 1;
			}else{
				cell[i][j] = 0;
			}
		}
	}

#ifdef DEBUG_MODE
	for (i = 0; i < row_per_rank; i++){
		for(j = 0; j < UNIVERSE; j++){
			cell[i][j] = rand()%2;
		}
	}
#endif

	if (mpi_myrank == 0){
		for (i = 0; i < row_per_rank; i++){
			for (j = 0; j < UNIVERSE; j++){
				if (cell[i][j] == 0){
					zero += 1;
				}else{
					one += 1;
				}
			}
		}
		//printf("%d\n", zero);
		//printf("%d\n", one);
		//printf("##################initial##################\n");
	}


	MPI_Barrier(MPI_COMM_WORLD);

	// number of row for each thread
	int row_per_thread = row_per_rank/(CHILDREN + 1);

	// if it has children, then put it into the new thread
	if (CHILDREN != 0){
		for(j = 0; j < CHILDREN; j++){
			variable *t;
			t = (variable *)malloc(sizeof(variable));
			t -> start = row_per_thread * (j + 1);
			t -> row = row_per_thread;
			t -> column = UNIVERSE;
			t -> cell = cell;
			t -> count = count;
			t -> cell2 = cell2;
			t -> ghostup = ghostup;
			t -> ghostdown = ghostdown;
			t -> mpi_myrank = mpi_myrank;
			t -> row_per_rank = row_per_rank;

			rc = pthread_create(&tid[j], NULL, whattodo, t);
			if (rc != 0){
				fprintf(stderr, "MIAN: Could not create thread (%d)\n", rc);
				return EXIT_FAILURE;
			}
		}
	}

	for (i = 0; i < TICK; i++){
		//receive and send ghost row from main thread
		MPI_Irecv(ghostup, UNIVERSE, MPI_INT, (mpi_myrank - 1)%mpi_commsize, 0, MPI_COMM_WORLD, &myRequest);
		MPI_Isend(cell[row_per_rank - 1], UNIVERSE, MPI_INT, (mpi_myrank+1)%mpi_commsize, 0, MPI_COMM_WORLD, &myRequest);

		MPI_Irecv(ghostdown, UNIVERSE, MPI_INT, (mpi_myrank + 1)%mpi_commsize, 1, MPI_COMM_WORLD, &myRequest);
		MPI_Isend(cell[0], UNIVERSE, MPI_INT, (mpi_myrank - 1)%mpi_commsize, 1, MPI_COMM_WORLD, &myRequest);

		MPI_Barrier(MPI_COMM_WORLD);
		pthread_mutex_lock(&mutex);
		semaphore1 += 1;
		pthread_mutex_unlock(&mutex);
		while (1){
			if (semaphore1 %(CHILDREN + 1) == 0){
				break;
			}
		}

		for (j = 0; j < row_per_thread; j++){
			for (m = 0; m < UNIVERSE; m++){
				// if smaller than threshold, randomly assign value, if it greater than threshold, then apply the rule.
				if (GenVal(mpi_myrank * row_per_rank + j) < THRESHOLD){
					if (GenVal(mpi_myrank * row_per_rank + j) < 0.5){
						cell2[j][m] = 0;
					}else{
						cell2[j][m] = 1;
					}

				}else{
					count[j][m] = update_counter(cell, j, m, row_per_rank, UNIVERSE, ghostup, ghostdown);
					if(count[j][m] < 2){
						cell2[j][m] = 0;
					}
					else if(count[j][m] == 2){
						cell2[j][m] = cell[j][m];
					}
					else if(count[j][m] == 3){
						cell2[j][m] = 1;
					}
					else{
						cell2[j][m] = 0;
					}
				}
			}
		}

		// here we wait for other threads before update the cell board
		pthread_mutex_lock(&mutex);
		semaphore2 += 1;
		pthread_mutex_unlock(&mutex);
		while (1){
			if (semaphore2 %(CHILDREN + 1) == 0){
				break;
			}
		}

		// after all the thread finish updating the counting board, we update the counting board.
		for (j = 0; j < row_per_thread; j++){
			for (m = 0; m < UNIVERSE; m++){
				cell[j][m] = cell2[j][m];
			}
		}

		// we make sure all the thread finish updating the cell board in this tick
		pthread_mutex_lock(&mutex);
		semaphore3 += 1;
		pthread_mutex_unlock(&mutex);
		while (1){
			if (semaphore3 %(CHILDREN + 1) == 0){
				break;
			}
		}

		// result after each tick

#ifdef DEBUG_MODE
		if (mpi_myrank == 0){
			zero = 0;
			one = 0;
			for (p = 0; p < row_per_rank; p++){
				for (q = 0; q < UNIVERSE; q++){
					if (cell[p][q] == 0){
						zero += 1;
					}else{
						one += 1;
					}
				}
			}
			printf("%d\n", zero);
			printf("%d\n", one);
			printf("######### %d  ############\n", i);
		}
		if (i == TICK - 1){
            int number = 0;
            int ii,jj;
			for (ii = 0; ii < row_per_rank; ii++){
				for (jj = 0; jj < UNIVERSE; jj++){
					if (cell[ii][jj] != 0){
						number += 1;
					}
				}
			}
			printf("#########myRank: %d \t number: %d ############\n", mpi_myrank, number);

		}
#endif	

	}


	for (i = 0; i < CHILDREN; i++){
		rc = pthread_join(tid[i], NULL);
		if(rc != 0){
			fprintf(stderr, "MAIN: Could not join thread (%d)\n", rc);
		}
	}

	MPI_Barrier( MPI_COMM_WORLD );

    end_cycles = GetTimeBase();
    time_in_secs = ((double)(end_cycles - start_cycles)) / processor_frequency;


   
    printf("execution time: %f\n", time_in_secs);
	fflush(stdout);
    

    //Perform an output of Cell Universe using MPI_file_write_at

#ifdef DEBUG_MODE
            int number = 0;
            int ii,jj;
			for (ii = 0; ii < row_per_rank; ii++){
				for (jj = 0; jj < UNIVERSE; jj++){
					if (cell[ii][jj] != 0){
						number += 1;
					}
				}
			}
			printf("#########myRank: %d \t number: %d ############\n", mpi_myrank, number);

#endif

    start_cycles2 = GetTimeBase();

#ifdef DEBUG_MODE
    char *filenamefix;
    itoa(mpi_commsize,filenamefix);
    printf("%s\n",filenamefix);
#endif

    char *filename = "data.txt";

// convert 123 to string [buf]
	MPI_File fh;
	MPI_Status status;
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);

	 /* Allocate and initialize a buffer (buff1) containing NUM_INTS integers,
     where the integer in location i is set to i. */
    // if(mpi_myrank == 0){
        int ii,jj;
        char *buff1;
        int result = 0;
        buff1 = (char *)malloc(UNIVERSE*sizeof(char));

        for (ii = 0; ii < row_per_rank; ii++){
				for (jj = 0; jj < UNIVERSE; jj++){
                    if(cell[ii][jj]!=0)
                    buff1[jj] = '1';
                    else
                    buff1[jj] = '0';
            }
        //printf("%s\n", buff1);

        MPI_Offset offset = (UNIVERSE * row_per_rank) * mpi_myrank + UNIVERSE * ii;
        //printf("\n offset: %d\n", offset);
        result = MPI_File_write_at(fh, offset, buff1, UNIVERSE, MPI_CHAR, &status);
        if(result != MPI_SUCCESS)
            sample_error(result, "MPI_File_write_at");
    }


#ifdef DEBUG_MODE
	char buf[256];
	int count = 100;
	each process read a block at different offset
	MPI_File_read_at(fh, offset, buf, count, MPI_CHAR, &status);
	printf("rankID = %d, offset = %lld, \n buf = %s\n", rankID, offset, buf);

	int write_buf[128];
	sprintf(write_buf, "rankID = %d, totalTaskNum = %d\n", mpi_myrank, mpi_commsize);
	each process write a block at different offset
	MPI_File_write_at(fh, offset, write_buf, sizeof(write_buf), MPI_INT, &status);
    	MPI_File_write_at(fh, offset, **cell, sizeof(cell), MPI_INT, &status);
#endif

	MPI_File_close(&fh);//after open, fh has the communicator info
	free(buff1);


    //Collect I/O performance stats using MPI_Wtime.

    //Print Compute (and I/O if needed) Performance Stats!!

    end_cycles2 = GetTimeBase();
    time_in_secs2 = ((double)(end_cycles2 - start_cycles2)) / processor_frequency;



    if(mpi_myrank == 0){
    printf("Parallel I/O time: %f\n", time_in_secs2);
    }

	//free all the dynamic allocated memory
	for (i = 0; i < row_per_rank; i++){
		free(cell[i]);
	}
	free(cell);

	for (i = 0; i < row_per_rank; i++){
		free(count[i]);
	}
	free(count);

	for (i = 0; i < row_per_rank; i++){
		free(cell2[i]);
	}
	free(cell2);

// END -Perform a barrier and then leave MPI
	MPI_Finalize();
	return 0;
}
/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

int update_counter(int **cell, int j, int m, int row, int column, int *ghostup, int *ghostdown){
	int rst = 0;
	if(j == 0){
		rst += cell[j][(m - 1)%column];
		rst += cell[j][(m + 1)%column];
		rst += cell[j + 1][(m - 1)%column];
		rst += cell[j + 1][(m)%column];
		rst += cell[j + 1][(m + 1)%column];
		rst += ghostup[(m - 1)%column];
		rst += ghostup[(m)%column];
		rst += ghostup[(m + 1)%column];
	}
	else if(j == (row - 1)){
		rst += cell[j - 1][(m - 1)%column];
		rst += cell[j - 1][(m)%column];
		rst += cell[j - 1][(m + 1)%column];
		rst += cell[j][(m - 1)%column];
		rst += cell[j][(m + 1)%column];
		rst += ghostdown[(m - 1)%column];
		rst += ghostdown[(m)%column];
		rst += ghostdown[(m + 1)%column];
	}
	else{

		rst += cell[j - 1][(m - 1)%column];
		rst += cell[j - 1][(m)%column];
		rst += cell[j - 1][(m + 1)%column];
		rst += cell[j][(m - 1)%column];
		rst += cell[j][(m + 1)%column];
		rst += cell[j + 1][(m - 1)%column];
		rst += cell[j + 1][(m)%column];
		rst += cell[j + 1][(m + 1)%column];
	}
	return rst;

}


void * whattodo(void *arg){
	variable * arg_temp;
	arg_temp = (variable *)arg;
	int ** cell = arg_temp -> cell;
	int ** count  = arg_temp -> count;
	int ** cell2 = arg_temp -> cell2;
	int start = arg_temp -> start;
	int row = arg_temp -> row;
	int mpi_myrank = arg_temp -> mpi_myrank;
	int row_per_rank = arg_temp -> row_per_rank;
	int column = arg_temp -> column;
	int * ghostup = arg_temp -> ghostup;
	int * ghostdown = arg_temp -> ghostdown;

	free(arg);
	int i;
	int j;
	int k;
	int m;

	for (k = 0; k < TICK; k++){
		// make sure execute the children after the main thread has distributed the new ghost row
		pthread_mutex_lock(&mutex);
		semaphore1 += 1;
		pthread_mutex_unlock(&mutex);
		while (1){
			if (semaphore1 %(CHILDREN + 1) == 0){
				break;
			}
		}

		// the same logic as the main thread
		for (i = start; i < (start + row); i++){
			for (j = 0; j < column; j++){
				if (GenVal(mpi_myrank * row_per_rank + i) < THRESHOLD){
					if (GenVal(mpi_myrank * row_per_rank + i) < 0.5){
						cell2[i][j] = 0;
					}else{
						cell2[i][j] = 1;
					}

				}else{
					count[i][j] = update_counter(cell, i, j, row_per_rank, UNIVERSE, ghostup, ghostdown);
					if(count[i][j] < 2){
						cell2[i][j] = 0;
					}
					else if(count[i][j] == 2){
						cell2[i][j] = cell[i][j];
					}
					else if(count[i][j] == 3){
						cell2[i][j] = 1;
					}
					else{
						cell2[i][j] = 0;
					}
				}
			}
		}

		// here we wait for all other threads
		pthread_mutex_lock(&mutex);
		semaphore2 += 1;
		pthread_mutex_unlock(&mutex);
		while (1){
			if (semaphore2 %(CHILDREN + 1) == 0){
				break;
			}
		}

		// after all threads finish the update counting board, we update the cell board.
		for (j = start; j < (start + row); j++){
			for (m = 0; m < UNIVERSE; m++){
				cell[j][m] = cell2[j][m];
			}
		}

		// make sure all thread finish updating the cell board
		pthread_mutex_lock(&mutex);
		semaphore3 += 1;
		pthread_mutex_unlock(&mutex);
		while (1){
			if (semaphore3 %(CHILDREN + 1) == 0){
				break;
			}
		}


	}
	return NULL;
}
int sample_error(int error, char *string)
{
  fprintf(stderr, "Error %d in %s\n", error, string);
  MPI_Finalize();
  exit(-1);
}
