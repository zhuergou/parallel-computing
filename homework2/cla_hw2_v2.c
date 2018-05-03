#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<mpi.h>

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
int dataRead(int * hexNumber, char * path){
	FILE *fp = NULL;
        long count = 0L;
	long i = 0L;
	fp = fopen(path,"r");
	if (NULL == fp){
		perror("ERROR: <can't open>");
		return EXIT_FAILURE;
	}
	char hexString[262150L];
//	scanf("%s",hexString);
	while((hexString[count] = fgetc(fp)) != EOF){
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

	fclose(fp);
	return 0;
	}

int main(int argc, char *argv[]){

	int *result_all = NULL;
	int *num1 = NULL;
	int *num2 = NULL;
	int *c = NULL;
	int g;
	int p;
	int i;
	double start;
	double end;

	MPI_Init(&argc, &argv);
	
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	start = MPI_Wtime();

	if ((c = malloc(1048577L * sizeof(int))) == NULL){
		perror("allocate failed");
		return EXIT_FAILURE;
	}

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

	dataRead(num1, argv[1]);
	dataRead(num2, argv[2]);
	
	c[0] = 0;
	for (i = 0; i < 1048576; i++){
		result_all[i] = num1[i]^num2[i]^c[i];
		g = num1[i] & num2[i];
		p = num1[i] | num2[i];
		c[i + 1] = g | (p & c[i]);	
	}
		
	char rst[262150];
	reverse(result_all, 1048576);
	convertHex(result_all, rst);
	
	if (world_rank == 0){
		printf("%s\n",rst);
		fflush(stdout);
	}
	free(result_all);
	result_all = NULL;
	free(num1);
	num1 = NULL;
	free(num2);
	num2 = NULL;
	free(c);
	c = NULL;

	end = MPI_Wtime();

	printf("the running time is %f\n", end - start);
	
	MPI_Finalize();

	return 0;
}
