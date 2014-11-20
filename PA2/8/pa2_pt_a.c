/*
 * worksharing across the rows
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

#define N 1024
#define diag 5.0
#define recipdiag 0.2
#define odiag -1.0
#define eps  1.0E-6
#define maxiter 10000
#define nprfreq 10

#define NUM_THREADS 8
#define CHUNK_SIZE 32

double clkbegin, clkend;
double t;
double rtclock();

void init(double*,double*,double*);
double rhocalc(double*);
void update(double*,double*,double*,double*); 
void copy(double*,double*);
int main (int argc, char * argv[])
{
  double * b = malloc(sizeof(double)*(N+2)*(N+2));
  double * xold = malloc(sizeof(double)*(N+2)*(N+2));
  double * xnew = malloc(sizeof(double)*(N+2)*(N+2));
  double * resid = malloc(sizeof(double)*(N+2)*(N+2));
  double rhoinit,rhonew; 
  int i,j,iter,u;

  init(xold,xnew,b);
  rhoinit = rhocalc(b);

  clkbegin = rtclock();
  for(iter=0;iter<maxiter;iter++){
    // This can of course not be parallelized. The number of iterations will be the same no matter how we divide the work up amongst the processors.
    update(xold,xnew,resid,b);
    rhonew = rhocalc(resid);
    if(rhonew<eps){
      clkend = rtclock();
      t = clkend-clkbegin;
      printf("Solution converged in %d iterations\n",iter);
      printf("Final residual norm = %f\n",rhonew);
      printf("Solution at center and four corners of interior N/2 by N/2 grid : \n");
      i=(N+2)/4; j=(N+2)/4; printf("xnew[%d][%d]=%f\n",i,j,xnew[i*(N+2)+j]);
      i=(N+2)/4; j=3*(N+2)/4; printf("xnew[%d][%d]=%f\n",i,j,xnew[i*(N+2)+j]);
      i=(N+1)/2; j=(N+1)/2; printf("xnew[%d][%d]=%f\n",i,j,xnew[i*(N+2)+j]);
      i=3*(N+2)/4; j=(N+2)/4; printf("xnew[%d][%d]=%f\n",i,j,xnew[i*(N+2)+j]);
      i=3*(N+2)/4; j=3*(N+2)/4; printf("xnew[%d][%d]=%f\n",i,j,xnew[i*(N+2)+j]);
      printf("Sequential Jacobi: Matrix Size = %d; %.1f GFLOPS; Time = %.3f sec; \n",
	     N,13.0*1e-9*N*N*(iter+1)/t,t); 
      break;
    } 
    copy(xold,xnew);
    if((iter%nprfreq)==0)
      printf("Iter = %d Resid Norm = %f\n",iter,rhonew);
  }
} 

void init(double * xold, double * xnew, double * b)
{ 
  //Do not parallelize this.
  int i,j;
  for(i=0;i<N+2;i++){
    for(j=0;j<N+2;j++){
      xold[i*(N+2)+j]=0.0;
      xnew[i*(N+2)+j]=0.0;
      b[i*(N+2)+j]=i+j; 
    }
  }
}

double rhocalc(double * A)
{
	//Parallelize this
	double tmp, temp;
	int i,j;
	tmp = 0.0;
	temp = 0.0;
	for(i=1;i<N+1;i++)
	{
		#pragma omp parallel for shared(A) num_threads(NUM_THREADS)\
			schedule(dynamic, CHUNK_SIZE)
		for(j=1;j<N+1;j++)
			#pragma omp atomic
			tmp+=A[i*(N+2)+j]*A[i*(N+2)+j];
			// This is a terrible solution
	}
	return(sqrt(tmp));
}

void update(double * xold,double * xnew,double * resid, double * b)
{
  //Parallelize this
  int i,j;
  for(i=1;i<N+1;i++)
	#pragma omp parallel for shared(xnew, xold) num_threads(NUM_THREADS)\
		schedule(dynamic, CHUNK_SIZE)
    for(j=1;j<N+1;j++){
      xnew[i*(N+2)+j]=b[i*(N+2)+j]-odiag*(xold[i*(N+2)+j-1]+xold[i*(N+2)+j+1]+\
	    xold[(i+1)*(N+2)+j]+xold[(i-1)*(N+2)+j]);
      xnew[i*(N+2)+j]*=recipdiag;
    }
  for(i=1;i<N+1;i++)
  	#pragma omp parallel for shared(xnew, resid) num_threads(NUM_THREADS)\
		schedule(dynamic, CHUNK_SIZE)
    for(j=1;j<N+1;j++)
    {
      resid[i*(N+2)+j]=b[i*(N+2)+j]-diag*xnew[i*(N+2)+j]-odiag*(xnew[i*(N+2)+j+1]+\
           xnew[i*(N+2)+j-1]+xnew[(i-1)*(N+2)+j]+xnew[(i+1)*(N+2)+j]);
    } 
} 
  
void copy(double * xold, double * xnew)
{
	//Parallelize this
	int i,j;
	for(i=1;i<N+1;i++)
		#pragma omp parallel for shared(xnew, xold) num_threads(NUM_THREADS)\
			schedule(dynamic, CHUNK_SIZE)
		for(j=1;j<N+1;j++)
			xold[i*(N+2)+j]=xnew[i*(N+2)+j];
}

double rtclock()
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

