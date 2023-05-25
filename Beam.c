#define  _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>
#include <sys/types.h>
#include <string.h> 
#include <math.h>

/*My method uses a banded matrix of size (N-1) to compute values of w(x) for x != 0, using ghost points at each boundary.*/


long read_input(double *L, long *N, double *om, char *fname) {
	
	FILE* fptr=fopen(fname,"r");

	if (fptr==NULL) return 1;

	if (3!=fscanf(fptr,"%lf %ld %lf", L, N, om)) {
  		return 1;
	}

	fclose(fptr);
	return 0;
}

long read_coefficients(char *fname, double mu[], double K[], double q[]) {

	
	FILE * fp;
    	char * line = NULL;
	int i = 0, j = 0;
	size_t len = 0;
    	ssize_t read;

    	fp = fopen(fname, "r");
   	
	if (fp == NULL)
		return 1;

	read = getline(&line, &len, fp);
	char * token = strtok(line, " ");
        token = strtok(NULL, " ");
	K[0] = atof(token);
	j = 0;

  	while ((read = getline(&line, &len, fp)) != -1) {
		token = strtok(line, " ");
   		while( token != NULL ) {

             		if (j==0) {
				mu[i] = atof(token);
			}

             		if (j==1) {
				K[i+1] = atof(token);
			}

             		if (j==2) {
				q[i] = atof(token);
			}

                   	token = strtok(NULL, " ");
			j += 1;
              	}
		j = 0;
		i += 1;
    	}
    	fclose(fp);
    	if (line)
        	free(line);
	return 0;
}

void write_output(double *w, double dx, long N) {
	FILE *fptr;
	fptr = fopen("output.txt", "w");
	if(fptr == NULL) {
      		printf("File write error");   
      		exit(1);             
   	}
	
	fprintf(fptr, "0.0 0.0 \n");
	
	double x = dx, dxfour = pow(dx, 4);

	for (int i = 0; i < N - 1; i++) {
		fprintf(fptr, "%lf %lf \n", x, w[i] * -dxfour);
		x += dx;
	}

	fclose(fptr);
}

struct band_mat{
	
 	long ncol;       
  	long nbrows;    
  	long nbands_up;
  	long nbands_low; 
  	double *array;  
  	long nbrows_inv;
  	double *array_inv;
  	int *ipiv;       

};

typedef struct band_mat band_mat;

int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
 	
	bmat->nbrows = nbands_lower + nbands_upper + 1;
 	bmat->ncol   = n_columns;
 	bmat->nbands_up = nbands_upper;
 	bmat->nbands_low= nbands_lower;
 	bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
 	bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
 	bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
 	bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
 	
	if (bmat->array==NULL||bmat->array_inv==NULL) {
 		return 0;
	}  
  	
	long i;
  	
	for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    		bmat->array[i] = 0.0;
  	}
  	
	return 1;

}; 

double *getp(band_mat *bmat, long row, long column) {
  	int bandno = bmat->nbands_up + row - column;
  	if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    		printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    		exit(1);
  	} 	
	return &bmat->array[bmat->nbrows*column + bandno];
}

double getv(band_mat *bmat, long row, long column) {
  	return *getp(bmat,row,column);
}

double setv(band_mat *bmat, long row, long column, double val) {
  	*getp(bmat,row,column) = val;
  	return val;
}

int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  	int i,bandno;
  	for(i=0;i<bmat->ncol;i++) { 
    		for (bandno=0;bandno<bmat->nbrows;bandno++) {
      			bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
  		}
    		x[i] = b[i];
  	}

  	long nrhs = 1;
  	long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
  	int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
  	return info;
}

void finalise_band_mat(band_mat *bmat) {
  	free(bmat->array);
  	free(bmat->array_inv);
  	free(bmat->ipiv);
}

void print_arr(double * x, long size) {
	long i = 0;
	for(; i < size; i++) {
		printf("%lf\n", x[i]); 
	}
}

int main(void) {

	double L, om;
	long N;
	char *fname = "input.txt";
	if (read_input(&L, &N, &om, fname)) {
    		printf("File read error\n");
   		return 1;
 	}
	
	double *mu, *K, *q;
		
	mu = malloc(sizeof(double)*N-1);
	K = malloc(sizeof(double)*N);
	q = malloc(sizeof(double)*N-1);

	read_coefficients("coefficients.txt", mu, K, q);
	
	double dx = L / (double) (N-1) ;
	
	double a = pow(dx, 4) * pow(om, 2);
	
	long ncols = N - 1, nbands_low = 2, nbands_up = 2;
	
	band_mat bmat;
	
	init_band_mat(&bmat, nbands_low, nbands_up, ncols);
	
	long i = 1;

	/* Set matrix values, looping from i = 1 (k = 2) to i = N - 2 (k = N - 1), to exclude w0 (already known) and  */
	/* to avoid evaluating K[-1], K[N] */
	for (;i < ncols-1; i++) {
		{setv(&bmat, i, i-1, 2*(K[i-1] + K[i]));};
		if (i > 1) {setv(&bmat, i, i-2, -K[i-1]);};
		setv(&bmat, i, i, -K[i-1] - 4*K[i] - K[i+1] - a*mu[i]);
		if (i < ncols-2) {setv(&bmat, i, i+2, -K[i+1]);};
		setv(&bmat, i, i+1, 2*(K[i] + K[i+1]));
	}

	/* Assign values for row 1 */	
	
	setv(&bmat, 0, 0, -2*K[0] -4*K[1] -K[2] -a*mu[0]);
	setv(&bmat, 0, 1, 2*(K[1] + K[2]));
	setv(&bmat, 0, 2, -K[2]);
		
	/* Correct value for entry [N-3][N-3] */
	
	setv(&bmat, N-3, N-3, -K[N-3] -4*K[N-2] -a*mu[N-3]);
	
	/* Assign values for row N-2 */
	
	setv(&bmat, N-2, N-4, 2*K[N-1]);
	setv(&bmat, N-2, N-3, -4*K[N-1]);
	setv(&bmat, N-2, N-2, 2*K[N-1] + -a*mu[N-2]);

	double * w = malloc(sizeof(double)*ncols);

	solve_Ax_eq_b(&bmat, w, q);

	write_output(w, dx, N);
	
	finalise_band_mat(&bmat);
	free(mu);
	free(K);
	free(q);	
	free(w);
	return(0);
}	
