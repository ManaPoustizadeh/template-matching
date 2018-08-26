#include <stdio.h>
#include <math.h>
#include "bmpread.h"
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#pragma optimize( "2", on )


bool piece_check(unsigned char *I, unsigned char *T, int i, int j, int T_width, int T_height, int I_width){
	long double dif = 0.0;
	double dif_avg = 0;
//#pragma omp parallel for reduction(+:dif)
	for (int m = i; m <= (i + T_height - 1); m++){
		for (int n = j; n <= j + 3 * T_width - 1; n++){
			dif += abs(I[3 * m*I_width + n] - T[3 * (m - i)*T_width + n - j]);
		}
		dif_avg = dif/(T_width*(m-i+1));
		if (dif_avg > 20)
			return false;
	}
	dif /= T_width*T_height * 3;
	if (dif < 8)
		return true;
	return false;
}

int tm(unsigned char *I, unsigned char *T, int I_width, int I_height, int T_width, int T_height){
	int cnt = 0;
#pragma omp parallel for
	for (int i = 0; i <= (I_height - T_height); i++){
#pragma omp parallel for reduction(+:cnt)
		for (int j = 0; j <= 3 * (I_width - T_width); j += 3){
			if (piece_check(I, T, i, j, T_width, T_height, I_width))
				cnt++;
		}
	}
	//printf("%d\n", cnt);
	return cnt;
}

unsigned char *rotate(unsigned char* T, int T_width, int T_height){
	unsigned char *rotationR, *rotationG, *rotationB, *r,*g, *b, *rgb2;
	int T_size = T_height*T_width * 3;
	rotationR = (unsigned char*)malloc(sizeof(unsigned char) * T_size/3);
	rotationG = (unsigned char*)malloc(sizeof(unsigned char) * T_size / 3);
	rotationB = (unsigned char*)malloc(sizeof(unsigned char) * T_size / 3);
	r = (unsigned char*)malloc(sizeof(unsigned char) * T_size/3);
	g = (unsigned char*)malloc(sizeof(unsigned char) * T_size / 3);
	b = (unsigned char*)malloc(sizeof(unsigned char) * T_size / 3);
	rgb2 = (unsigned char*)malloc(sizeof(unsigned char) * T_size);
	for (int i = 0; i < T_size; i+=3){
		r[i / 3] = T[i];
		g[i / 3] = T[i + 1];
		b[i / 3] = T[i + 2];
	}
	for (int j = 0; j < T_height; j++){
		for (int i = 0; i < T_width; i++){
			rotationR[i * T_height + T_height - j - 1] = r[j * T_width + i];
			rotationG[i * T_height + T_height - j - 1] = g[j * T_width + i];
			rotationB[i * T_height + T_height - j - 1] = b[j * T_width + i];
		}
	}
	for (int i = 0; i < T_size; i+=3){
		rgb2[i] = rotationR[i/3];
		rgb2[i + 1] = rotationG[i/3];
		rgb2[i + 2] = rotationB[i/3];
	}
	return rgb2;
}

void omp_check() {
	printf("------------ Info -------------\n");
#ifdef _DEBUG
	printf("[!] Configuration: Debug.\n");
#pragma message ("Change configuration to Release for a fast execution.")
#else
	printf("[-] Configuration: Release.\n");
#endif // _DEBUG
#ifdef _M_X64
	printf("[-] Platform: x64\n");
#elif _M_IX86 
	printf("[-] Platform: x86\n");
#pragma message ("Change platform to x64 for more memory.")
#endif // _M_IX86 
#ifdef _OPENMP
	printf("[-] OpenMP is on.\n");
	printf("[-] OpenMP version: %d\n", _OPENMP);
#else
	printf("[!] OpenMP is off.\n");
	printf("[#] Enable OpenMP.\n");
#endif // _OPENMP
	printf("[-] Maximum threads: %d\n", omp_get_max_threads());
	printf("[-] Nested Parallelism: %s\n", omp_get_nested() ? "On" : "Off");
#pragma message("Enable nested parallelism if you wish to have parallel region within parallel region.")
	printf("===============================\n");
}


int main(int argc, char *argv[]) {
	if (argc % 2 == 0){
		printf("Invalid number of inputs!\n");
		exit(0);
	}
	omp_set_nested(1);
	//omp_check();
	//int I_width, I_height, T_width, T_height;
	unsigned char *I, *T, *X1;
	bmpread_t img, temp;
	double start, end;
	int cnt;
	start = omp_get_wtime();

	for (int i = 1; i < argc; i+=2){
		if (!bmpread(argv[i], BMPREAD_BYTE_ALIGN | BMPREAD_ANY_SIZE, &img)){
			fprintf(stderr, "%s: error loading bitmap file\n", argv[i]);
			//return EXIT_FAILURE;
		}
		if (!bmpread(argv[i+1], BMPREAD_BYTE_ALIGN | BMPREAD_ANY_SIZE, &temp)){
			fprintf(stderr, "%s: error loading bitmap file\n", argv[i+1]);
			//return EXIT_FAILURE;
		}
		//printf("height and width are: %d %d \n", img.width, img.height);
		//printf("height and width are: %d %d \n", temp.width, temp.height);

		cnt = tm(img.data, temp.data, img.width, img.height, temp.width, temp.height);
		
		X1 = rotate(temp.data, temp.width, temp.height);
		cnt += tm(img.data, X1, img.width, img.height, temp.height, temp.width);
		
		printf("%d\n", cnt);
	}
	
	end = omp_get_wtime();
	//printf("Elapsed Time: %f secs\n", end-start);
	/*I = img.data;
	T = temp.data;
	I_width = img.width;
	I_height = img.height;
	T_width = temp.width;
	T_height = temp.height;
	for (int i = 0; i < 3 * T_width* T_height; i += 3)
	{
		printf("(%d, ", T[i]);
		printf("%d, ", T[i + 1]);
		printf("%d)\n", T[i + 2]);
	}
	unsigned char *X = rotate(T, T_width, T_height);*/
	/*for (int i = 0; i < 3 * I_width* I_height; i += 3)
	{
	printf("(%d, ", I[i]);
	printf("%d, ", I[i + 1]);
	printf("%d)\n", I[i + 2]);
	}
	printf("\n T: \n");*/
	/*for (int i = 0; i < 3* T_width* T_height; i+=3)
	{
	printf("%d, ", X[i]);
	printf("%d, ", X[i+1]);
	printf("%d)\n", X[i+2]);
	}*/
}
