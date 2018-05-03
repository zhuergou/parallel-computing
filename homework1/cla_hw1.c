#include<stdio.h>
#include<string.h>

// reverse the array
void reverse(int *s, int len){
	int i;
	int j;
	int tmp;
	for (i = 0, j = len - 1; i < j; i++, j--){
		tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
	} 
}

// main cla function
void cla(int * num1, int *num2, int end1, int end2, int *result){
	int g[256];
	int p[256];
	int gg[64];
	int gp[64];
	int sg[16];
	int sp[16];
	int ssg[4];
	int ssp[4];
	// initial c with the length of 257, so as to use c[0] to indicate c_{i-1}
	int c[257];
	int m;
	int i;
	int j;
	int k;
	//calculate g_i and p_i
	for (i = 0; i < end1; i++){
		g[i] = num1[i] & num2[i];
		p[i] = num1[i] | num2[i];
	}

	// calculate gg_j and gp_j
	i = 0;
	for (j = 0; j < 64; j++){
		gg[j] = g[i + 3] | (p[i + 3] & g[i + 2]) | (p[i + 3] & p[i + 2] & g[i + 1]) | (p[i + 3] & p[i + 2] & p[i + 1] & g[i]);
		gp[j] = p[i + 3] & p[i + 2] & p[i + 1] & p[i];
		i = i + 4;
	}

	//calculate sg_k and sp_k
	i = 0;
    for (j = 0; j < 16; j++){
		sg[j] = gg[i + 3] | (gp[i + 3] & gg[i + 2]) | (gp[i + 3] & gp[i + 2] & gg[i + 1]) | (gp[i + 3] & gp[i + 2] & gp[i + 1] & gg[i]);
	    sp[j] = gp[i + 3] & gp[i + 2] & gp[i + 1] & gp[i];
		i = i + 4;
	}

	// calculate ssg_l and ssp_l
	i = 0;
	for (j = 0; j < 4; j++){
		ssg[j] = sg[i + 3] | (sp[i + 3] & sg[i + 2]) | (sp[i + 3] & sp[i + 2] & sg[i + 1]) | (sp[i + 3] & sp[i + 2] & sp[i + 1] & sg[i]);
	    ssp[j] = sp[i + 3] & sp[i + 2] & sp[i + 1] & sp[i];
		i = i + 4;
	}

	// calculate ssc_l
	c[0] = 0;
	for(i = 1; i <= 4; i++){
		c[i * 64] = ssg[i - 1] | (ssp[i - 1] & c[(i - 1)* 64]);
	}

	//calculate sc_k
	for(i = 1; i <= 4; i++){
		for(j = 0; j < 4; j++){
			c[i * 16 + j * 64] = sg[(i - 1) + j * 4] | (sp[(i - 1) + j * 4] & c[(i - 1) * 16 + j * 64]);	
		}
	}

	//calculate gc_j
	for(k = 1; k <= 4; k++){
		for(i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				c[k * 4 + i * 16 + j * 64] = gg[(k - 1) + i * 4 + j * 16] | (gp[(k - 1) + i * 4 + j * 16] & c[(k - 1) * 4 + i * 16 + j * 64]);			
			}
		}
	}

	//calculate c_i
	for(m = 1; m <= 4; m++){
		for(k = 0; k < 4; k++){
			for(i = 0; i < 4; i++){
				for(j = 0; j < 4; j++){
					c[m + k * 4 + i * 16 + j * 64] = g[(m - 1) + k * 4 + i * 16 + j * 64] | (p[(m - 1) + k * 4 + i * 16 + j * 64] & c[(m - 1) + k * 4 + i * 16 + j * 64]);
				}
			}
		}
	}
	
	//calculate sum_i
	for(i = 0; i < 256; i++){
		result[i] = num1[i]^num2[i]^c[i]; 
	}

	/*
	if (c[256] == 1){
		result[256] = 1;
	}
	*/

	reverse(result, 256);

}

//convert binary result back to hex output
void convertHex(int * result, char *rst){

	int i;
	int temp;
	int j = 0;
	for(i = 0; i < 256; i = i + 4){
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
int dataRead(int * hexNumber){
    int count = 0;
	int i = 0;
	char hexString[70];
	scanf("%s",hexString);
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
	return i;
}

int main(){
	int num1[260];	
	int end1 = dataRead(num1);
	int num2[260];
	int end2 = dataRead(num2);
	int result[256];
	char rst[130];
	cla(num1,num2,end1,end2,result);
	convertHex(result, rst);

	printf("%s\n",rst);
}
