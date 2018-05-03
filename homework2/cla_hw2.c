#include<stdio.h>
#include<string.h>
#include<mpi.h>
#include<stdlib.h>

int g[1048576];
int p[1048576];
int gg[1048576];
int gp[1048576];
int sg[1048576];
int sp[1048576];
int ssg[1048576];
int ssp[1048576];
int sssg[1048576];
int sssp[1048576];
int c[1048577];
char hex_input_a[262145] = {0};
char hex_input_b[262145] = {0};
char rst[262150];

// reverse the array
void reverse(int *s, long len){
	long i;
	long j;
	int tmp;
	for (i = 0, j = len - 1; i < j; i++, j--){
		tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
	} 
}

// main cla function
void cla(int *partial1, int *partial2, long elements_per_proc, long elements_per_proc2, int *result, int world_rank, int world_size){

	long m;
	long i;
	long j;
	long k;
	int tmp1;
	long w;
	MPI_Request myRequest;
	MPI_Status myStatus;

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//calculate g_i and p_i
	for (i = 0; i < elements_per_proc; i++){
		g[i] = partial1[i] & partial2[i];
		p[i] = partial1[i] | partial2[i];
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif


	// calculate gg_j and gp_j
	i = 0;
	for (j = 0; j < elements_per_proc/16; j++){
		gg[j] = 0;
		gp[j] = 1;
		for (k = 0; k < 16; k++){
			tmp1 = 1;
			for (m = k + 1; m < 16; m++){
				tmp1 = tmp1 & p[i + m];
			}
			tmp1 = tmp1 & g[i + k];
			gg[j] = gg[j] | tmp1;
		}

		for (k = 0; k < 16; k++){
			gp[j] = gp[j] & p[i + k];
		}
		
		i = i + 16;
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// calculate sg_j and sp_j
	i = 0;
	for (j = 0; j < elements_per_proc/256; j++){
		sg[j] = 0;
		sp[j] = 1;
		for (k = 0; k < 16; k++){
			tmp1 = 1;
			for (m = k + 1; m < 16; m++){
				tmp1 = tmp1 & gp[i + m];
			}
			tmp1 = tmp1 & gg[i + k];
			sg[j] = sg[j] | tmp1;
		}

		for (k = 0; k < 16; k++){
			sp[j] = sp[j] & gp[i + k];
		}
		
		i = i + 16;
	}


#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// calculate ssg_j and ssp_j
	i = 0;
	for (j = 0; j < elements_per_proc/4096; j++){
		ssg[j] = 0;
		ssp[j] = 1;
		for (k = 0; k < 16; k++){
			tmp1 = 1;
			for (m = k + 1; m < 16; m++){
				tmp1 = tmp1 & sp[i + m];
			}
			tmp1 = tmp1 & sg[i + k];
			ssg[j] = ssg[j] | tmp1;
		}

		for (k = 0; k < 16; k++){
			ssp[j] = ssp[j] & sp[i + k];
		}
		
		i = i + 16;
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// calculate sssg_j and sssp_j
	i = 0;
	for (j = 0; j < elements_per_proc/65536; j++){
		sssg[j] = 0;
		sssp[j] = 1;
		for (k = 0; k < 16; k++){
			tmp1 = 1;
			for (m = k + 1; m < 16; m++){
				tmp1 = tmp1 & ssp[i + m];
			}
			tmp1 = tmp1 & ssg[i + k];
			sssg[j] = sssg[j] | tmp1;
		}

		for (k = 0; k < 16; k++){
			sssp[j] = sssp[j] & ssp[i + k];
		}
		
		i = i + 16;
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// calculate sssc_j
	
	int token = 0;

	if(world_rank != 0){

		MPI_Irecv(&token, 1, MPI_INT, world_rank - 1, 0, MPI_COMM_WORLD, &myRequest);
		MPI_Wait(&myRequest, &myStatus);
	}

	c[0] = token;
	for(i = 1; i <= elements_per_proc/65536; i++){
		c[i * 65536] = sssg[i - 1] | (sssp[i - 1] & c[(i - 1)* 65536]);
	}

	int sendValue = c[elements_per_proc];

	if(world_rank != (world_size - 1)){
	
		MPI_Isend(&sendValue, 1, MPI_INT, world_rank + 1, 0, MPI_COMM_WORLD, &myRequest);
		MPI_Wait(&myRequest, &myStatus);
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//calculate ssc_j
	for(i = 1; i <= 16; i++){
		for(j = 0; j < elements_per_proc/65536; j++){
			c[i * 4096 + j * 65536] = ssg[(i - 1) + j * 16] | (ssp[(i - 1) + j * 16] & c[(i - 1) * 4096 + j * 65536]);

		}
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//calculate sc_j
	for(k = 1; k <= 16; k++){
		for(i = 0; i < 16; i++){
			for(j = 0; j < elements_per_proc/65536; j++){
				c[k * 256 + i * 4096 + j * 65536] = sg[(k - 1) + i * 16 + j * 256] | (sp[(k - 1) + i * 16 + j * 256] & c[(k - 1) * 256 + i * 4096 + j * 65536]);
			}
		}
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//calculate gc_i
	for(m = 1; m <= 16; m++){
		for(k = 0; k < 16; k++){
			for(i = 0; i < 16; i++){
				for(j = 0; j < elements_per_proc/65536; j++){
					c[m * 16 + k * 256 + i * 4096 + j * 65536] = gg[(m - 1) + k * 16 + i * 256 + j * 4096] | (gp[(m - 1) + k * 16 + i * 256 + j * 4096] & c[(m - 1) * 16 + k * 256 + i * 4096 + j * 65536]);
				}
			}
		}
	}


#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//calculate c_i
	for (w = 1; w <= 16; w++){
		for(m = 0; m < 16; m++){
			for(k = 0; k < 16; k++){
				for(i = 0; i < 16; i++){
					for(j = 0; j < elements_per_proc/65536; j++){
						c[w + m * 16 + k * 256 + i * 4096 + j * 65536] = g[(w - 1) + m * 16 + k * 256 + i * 4096 + j * 65536] | (p[(w - 1) + m * 16 + k * 256 + i * 4096 + j * 65536] & c[(w - 1) + m  * 16 + k * 256 + i * 4096 + j * 65536]);
					}
				}
			}
		}
	
	}

#ifdef Barrier
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	//calculate sum_i
	for(i = 0; i < elements_per_proc; i++){
		result[i] = partial1[i]^partial2[i]^c[i]; 
	}

	/*
	if (c[elements_per_proc] == 1){
		result[elements_per_proc] = 1;
	}
	*/


}

//convert binary result back to hex output
void convertHex(int * result, char *rst){

	long i;
	int temp;
	long j = 0;
	for(i = 0; i < 1048576; i = i + 4){
		temp = 0;
		temp = temp + result[i]* 8 + result[i + 1]*4 + result[i + 2]*2 + result[i + 3];
		if (temp < 10){
			rst[j] = temp + '0';
		}else{
			rst[j] = temp - 10 + 'A';
		}
		j++;
	}
	rst[j] = '\0';
}
 

// read the data, and put it in to the array
int dataRead(int * hexNumber, char * hexString){
    long count = 0L;
	long i = 0L;

	while(hexString[count] != '\0'){
		switch(hexString[count]){
			case '0': hexNumber[i] = 0; hexNumber[i + 1] = 0; 
					  hexNumber[i + 2] = 0; hexNumber[i + 3] = 0; break;
			case '1': hexNumber[i] = 0; hexNumber[i + 1] = 0; 
                      hexNumber[i + 2] = 0; hexNumber[i + 3] = 1; break;
			case '2': hexNumber[i] = 0; hexNumber[i + 1] = 0; 
                      hexNumber[i + 2] = 1; hexNumber[i + 3] = 0; break;
			case '3': hexNumber[i] = 0; hexNumber[i + 1] = 0;
                      hexNumber[i + 2] = 1; hexNumber[i + 3] = 1; break;					case '4': hexNumber[i] = 0; hexNumber[i + 1] = 1;
                      hexNumber[i + 2] = 0; hexNumber[i + 3] = 0; break;
			case '5': hexNumber[i] = 0; hexNumber[i + 1] = 1;
                      hexNumber[i + 2] = 0; hexNumber[i + 3] = 1; break;
			case '6': hexNumber[i] = 0; hexNumber[i + 1] = 1;
                      hexNumber[i + 2] = 1; hexNumber[i + 3] = 0; break;	
			case '7': hexNumber[i] = 0; hexNumber[i + 1] = 1;
                      hexNumber[i + 2] = 1; hexNumber[i + 3] = 1; break;
			case '8': hexNumber[i] = 1; hexNumber[i + 1] = 0;
                      hexNumber[i + 2] = 0; hexNumber[i + 3] = 0; break;
			case '9': hexNumber[i] = 1; hexNumber[i + 1] = 0;
                      hexNumber[i + 2] = 0; hexNumber[i + 3] = 1; break;
			case 'a':
			case 'A': hexNumber[i] = 1; hexNumber[i + 1] = 0;
                      hexNumber[i + 2] = 1; hexNumber[i + 3] = 0; break;	
			case 'b':
	        case 'B': hexNumber[i] = 1; hexNumber[i + 1] = 0;
		              hexNumber[i + 2] = 1; hexNumber[i + 3] = 1; break;
			case 'c':
			case 'C': hexNumber[i] = 1; hexNumber[i + 1] = 1;
		              hexNumber[i + 2] = 0; hexNumber[i + 3] = 0; break;
			case 'd':
		    case 'D': hexNumber[i] = 1; hexNumber[i + 1] = 1;
		              hexNumber[i + 2] = 0; hexNumber[i + 3] = 1; break;
			case 'e':
		    case 'E': hexNumber[i] = 1; hexNumber[i + 1] = 1;
		              hexNumber[i + 2] = 1; hexNumber[i + 3] = 0; break;
			case 'f':
		    case 'F': hexNumber[i] = 1; hexNumber[i + 1] = 1;
		              hexNumber[i + 2] = 1; hexNumber[i + 3] = 1; break;
		}
		i = i + 4;
		count = count + 1;
	}
	reverse(hexNumber, i);

}

int main(int argc, char *argv[]){

	int *result_all = NULL;
	int *num1 = NULL;
	int *num2 = NULL;
	int world_size;
	int world_rank;
	double start;
	double end;

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Get the number of process
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	start = MPI_Wtime();

	long elements_per_proc = 1048576L/world_size;
	int *result = NULL;
	int *partial1 = NULL;
	int *partial2 = NULL;
	result = malloc(elements_per_proc * sizeof(int));
	partial1 = malloc(elements_per_proc * sizeof(int));
	partial2 = malloc(elements_per_proc * sizeof(int));

	if (world_rank == 0){
	        if ((result_all = malloc(1048576L * sizeof(int))) == NULL){
			perror("allocate failed");
			return EXIT_FAILURE;
		}
	        if((num1 = (int *)malloc(1048576L * sizeof(int))) == NULL){
			perror("allocate failed");
			return EXIT_FAILURE;
		}
	        if((num2 = (int *)malloc(1048576L * sizeof(int))) == NULL){
			perror("allocate failed");
			return EXIT_FAILURE;
		}

		FILE *my_input_file = NULL;
		if ((my_input_file = fopen(argv[1], "r")) == NULL){
			printf("Failed to open input data file: %s\n", argv[1]);
		
		}
		fscanf(my_input_file, "%s %s", hex_input_a, hex_input_b);
		fclose(my_input_file);
		dataRead(num1, hex_input_a);
		dataRead(num2, hex_input_b);
	}
	


	MPI_Scatter(num1, elements_per_proc, MPI_INT, partial1, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Scatter(num2, elements_per_proc, MPI_INT, partial2, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

	cla(partial1, partial2, elements_per_proc, elements_per_proc, result, world_rank, world_size);


/*  this is for ripple adder

	long i;
	int g2;
	int p2;
	c[0] = 0;
	for(i = 0; i < 1048576; i++){
		result_all[i] = num1[i]^num2[i]^c[i];
		g2 = num1[i] & num2[i];
		p2 = num1[i] | num2[i];
		c[i + 1] = g2 | (p2 & c[i]);
	}
*/


	MPI_Gather(result, elements_per_proc, MPI_INT, result_all, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

	if (world_rank == 0){
		FILE *my_output_file = NULL;
		reverse(result_all, 1048576);
		convertHex(result_all, rst);
		if ((my_output_file = fopen(argv[2], "w")) == NULL){
			printf("Failed to open input data file: %s \n", argv[2]);
		}
		fprintf(my_output_file, "%s\n", rst);
		fclose(my_output_file);
		free(result_all);
		result_all = NULL;
		free(num1);
		num1 = NULL;
		free(num2);
		num2 = NULL;
		end = MPI_Wtime();
//		printf("the running time is %f\n", end - start);
	}
	free(result);
	free(partial1);
	free(partial2);


	MPI_Finalize();
	return 0;
}
