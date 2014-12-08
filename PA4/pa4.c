#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include "mpi.h"

#define N 1024
#define diag 5.0
#define recipdiag 0.2
#define odiag -1.0
#define eps  1.0E-6
#define maxiter 10000
#define nprfreq 10

#define MASTER 0
#define NONE 0
#define BEGIN 1
#define LTAG 2
#define RTAG 3
#define DONE 4

double clkbegin, clkend;
double t;
double rtclock();

void init(double*,double*,double*);
double rhocalc(double*);
void update(double*,double*,double*,double*); 
void copy(double*,double*);

int main (int argc, char * argv[])
{

	/* MPI variables */
	MPI_Status status
	int 	taskid,
			numworkers,
			numtasks,
			averow, rows, offset, extra
			dest, source, 
			left, right, 
			msgtype,
			k;
			
	/* MPI setup*/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	
	/************************* MASTER CODE *******************************/
	if (taskid == MASTER)
	{
		/* Initialize method variables */
		double * b = malloc(sizeof(double)*(N+2)*(N+2));
		double * xold = malloc(sizeof(double)*(N+2)*(N+2));
		double * xnew = malloc(sizeof(double)*(N+2)*(N+2));
		double * resid = malloc(sizeof(double)*(N+2)*(N+2));
		double rhoinit,rhonew; 
		int i,j,iter,u;
		init(xold,xnew,b);
		rhoinit = rhocalc(b);
		
		/* Information about the distribution of work. */
		averow = N/numtasks; 
		extra = N%numtasks;
		offset = 0;
		for (k = 1; k <= numtasks; k++)
		{
			rows = (i <= extra) ? averow + 1 : averow;
			if(i == 1)
				left = NONE;
			else 
				left = i-1;
			if(i == numworkers)
				right == NONE;
			else
				right = i+1;
			dest = i; /* Who will we send the current information to?*/
			MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&xold, 1, MPI_DOUBLE, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&b, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&resid, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);			
		}
	}
	/************************* WORKER CODE *******************************/
	if (taskid != MASTER)
	{
	}
	MPI_Barrier
	clkbegin = rtclock();
	for(iter=0;iter<maxiter;iter++){
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
 double tmp;
 int i,j;

 tmp = 0.0;
 for(i=1;i<N+1;i++)
   for(j=1;j<N+1;j++)
     tmp+=A[i*(N+2)+j]*A[i*(N+2)+j];
 return(sqrt(tmp));
}

void update(double * xold, double * xnew,double * resid, double * b)
{
 int i,j;
 for(i=1;i<N+1;i++)
   for(j=1;j<N+1;j++){
     xnew[i*(N+2)+j]=b[i*(N+2)+j]-odiag*(xold[i*(N+2)+j-1]+xold[i*(N+2)+j+1]+xold[(i+1)*(N+2)+j]+xold[(i-1)*(N+2)+j]);
     xnew[i*(N+2)+j]*=recipdiag;
   }
 for(i=1;i<N+1;i++)
   for(j=1;j<N+1;j++){
     resid[i*(N+2)+j]=b[i*(N+2)+j]-diag*xnew[i*(N+2)+j]-odiag*(xnew[i*(N+2)+j+1]+xnew[i*(N+2)+j-1]+xnew[(i-1)*(N+2)+j]+xnew[(i+1)*(N+2)+j]);
   } 
} 
  
void copy(double * xold, double * xnew)
{ 
 int i,j;
 for(i=1;i<N+1;i++)
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


