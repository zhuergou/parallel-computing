#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

#define BGQ 1

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime
#endif


long long global_array[1073741824];

// handwritten MPI_P2P_Reduce
void MPI_P2P_Reduce(long long *send_data, long long *recv_data, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
	int comm_rank, comm_size;
	MPI_Comm_rank(comm, &comm_rank);
	MPI_Comm_size(comm, &comm_size);
	MPI_Request myRequest;
	MPI_Status myStatus;

	int i;

	//use this to control stride
	int index = 2;
	long long tmp = 0;

	while (index <= comm_size){
		for(i = 0; i < comm_size; i  = i + index){
			if (comm_rank == i){
				MPI_Irecv(&tmp, 1, MPI_LONG_LONG, i + (index/2), 0, comm, &myRequest);
				MPI_Wait(&myRequest, &myStatus);
				*send_data = *send_data + tmp;
				tmp = 0;		
			}
			else if (comm_rank == i + (index / 2)){
				MPI_Isend(send_data, 1, MPI_LONG_LONG, i, 0, comm, &myRequest);
				MPI_Wait(&myRequest, &myStatus);
			}
		}
		index = index * 2;
	}
	if (comm_rank == root){
		*recv_data = *send_data;
	}
}



int main(int argc, char *argv[]){
	int world_size;
	int world_rank;
	long long i;
	long long local_sum = 0;
	long long local_sum_2 = 0;
	long long global_sum_1 = 0;
	long long global_sum_2 = 0;
	double time_in_secs_1 = 0;
	double time_in_secs_2 = 0;
	double processor_frequency = 1600000000.0;
	unsigned long long start_1 = 0;
	unsigned long long end_1 = 0;	
	unsigned long long start_2 = 0;
	unsigned long long end_2 = 0;
	unsigned long long start_3 = 0;
	unsigned long long end_3 = 0;

	for (i = 0; i < 1073741824; i++){
		global_array[i] = i;
	}

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Get the number of process
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	long long elements_per_proc = 1073741824/world_size;
	

	start_1 = GetTimeBase();
	for (i = world_rank *elements_per_proc; i < (world_rank + 1) *elements_per_proc; i++){
		local_sum += global_array[i];
	}
	end_1 = GetTimeBase();
	
	local_sum_2 = local_sum;
	
	start_2 = GetTimeBase();
	MPI_P2P_Reduce(&local_sum, &global_sum_1, 1, MPI_LONG_LONG,  0, MPI_COMM_WORLD);
	end_2 = GetTimeBase();
	
	start_3 = GetTimeBase();
	MPI_Reduce(&local_sum_2, &global_sum_2, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	end_3 = GetTimeBase();

	time_in_secs_1 = ((double)(end_1 + end_2 - start_1 - start_2))/ processor_frequency;
	time_in_secs_2 = ((double)(end_1 + end_3 - start_1 - start_3))/ processor_frequency;


	if(world_rank == 0){
		printf("%lld %f\n", global_sum_1, time_in_secs_1);
		printf("%lld %f\n", global_sum_2, time_in_secs_2);
	}

	// Finalize mpi
	MPI_Finalize();
	return 0;
}
