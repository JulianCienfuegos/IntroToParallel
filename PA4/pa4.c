/*
 * Everyone gets the initial data copy. 
 * Then everyone does the update computation for this portion. 
 * The updated data is sent to master.
 * MASTER computes the residual.
 * MASTER checks the residual.
 * If no good, receive left, send right, receive right, send left.
 * Then new is copied to old.
 * Repeat
 */
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include "mpi.h"

#define N 256
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
void update(double*,double*,double*,double*, int, int); 
void copy(double*,double*, int, int);

int main (int argc, char * argv[])
{
		
/* Get ready to use MPI */
MPI_Status status;
MPI_Request request;
int 	taskid,
		numworkers,
		numtasks,
		averow, rows, master_rows, offset, extra,
		dest, source, 
		left, right, 
		msgtype,
		k;
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
numworkers = numtasks-1;

/* 
 * Initialize variables for the PDE solver.
 * Everyone has a local copy of all of these variables.
 * Each process will update only a portion, however, and the MASTER
 * nope is responsible for fitting the pieces into their appropriate 
 * places in the MASTER copy.
 */
 
double * b = malloc(sizeof(double)*(N+2)*(N+2));
double * xold = malloc(sizeof(double)*(N+2)*(N+2));
double * xnew = malloc(sizeof(double)*(N+2)*(N+2));
double * resid = malloc(sizeof(double)*(N+2)*(N+2));
double rhoinit,rhonew; 
int i,j,iter,u, go;
init(xold,xnew,b);
rhoinit = rhocalc(b);

/************************* MASTER CODE *******************************/
if (taskid == MASTER)
{
	/*
	 * Information about the distribution of work. 
	 * The MASTER has to do a little bit of extra work
	 */
	averow = (N+2)/numtasks; 
	extra = (N+2)%numtasks;
	offset = 0;
	for (k = 0; k < numtasks; k++)
	{
		rows = (k <= extra-1) ? averow + 1 : averow;
		if(k == 0)
		{
			left = NONE;
			master_rows = rows;
		}
		else 
			left = k-1;
		if(k == numworkers)
			right == NONE;
		else
			right = k+1;
		if(k!=0)
		{
			dest = k; 
			MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
		}	
		offset+=rows;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	clkbegin = rtclock();
	for (iter = 0; iter< maxiter; iter++)
	{
		/* Perform a calculation*/
		int start = 1;
		update(xold, xnew, resid, b, start, master_rows);
		
		/*Gather results*/
		for(k = 1; k <= numworkers; k++)
		{
			msgtype = DONE;
			MPI_Recv(&offset, 1, MPI_INT, k, msgtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows, 1, MPI_INT, k, msgtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&resid[offset*N], rows*(N+2), MPI_DOUBLE, k, msgtype, MPI_COMM_WORLD, &status);
		}
		/* Check error. We could potentially use MPI_Reduce here, but there is no time. */
		rhonew = rhocalc(resid);
		if((iter%nprfreq)==0)
		printf("Iter = %d Resid Norm = %f\n",iter,rhonew);
		if(rhonew<eps)
		{
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
			printf("MPI Jacobi: Matrix Size = %d; %.1f GFLOPS; Time = %.3f sec; \n", N,13.0*1e-9*N*N*(iter+1)/t,t); 
			break;
		} 
		
		/* Send/Recv with neighbor */
		if(numworkers > 0)
		{
			MPI_Send(&xnew[(master_rows-1)*(N+2)], N+2, MPI_DOUBLE, right, BEGIN, MPI_COMM_WORLD);
			MPI_Recv(&xnew[(master_rows+1)*(N+2)], N+2, MPI_DOUBLE, right, BEGIN, MPI_COMM_WORLD, &status);
		}		
		/*Copy*/
		copy(xold, xnew, start, master_rows);
	}
	
	MPI_Abort(MPI_COMM_WORLD,0); // Instead of finalize, because we want all of the processes to terminate once an answer is found.
/*********************** END OF MASTER CODE **************************/	
}
/************************* WORKER CODE *******************************/

if (taskid != MASTER)
{
/*
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Recv(&offset, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
	MPI_Recv(&rows, 1, MPI_INT, MSTER, BEGIN, MPI_COMM_WORLD, &status);
	MPI_Recv(&left, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
	MPI_Recv(&right, 1, MPI_INT, MASTER, BEGIN, MPI_COMM_WORLD, &status);
	for(iter=0;iter<maxiter;iter++)
	{
		// Perform a calculation
		update(xold,xnew,resid,b);
		//Send Results
		
		//Send/recv with neighbor
		
		//Copy		
		copy(xold,xnew);
	}
	
	MPI_Finalize();
*/
	
}

/*********************** END OF WORKER CODE **************************/

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

void update(double * xold, double * xnew,double * resid, double * b, int start, int rows)
{
 int i,j;
 for(i=start;i<start+rows;i++)
   for(j=1;j<N+1;j++){
     xnew[i*(N+2)+j]=b[i*(N+2)+j]-odiag*(xold[i*(N+2)+j-1]+xold[i*(N+2)+j+1]+xold[(i+1)*(N+2)+j]+xold[(i-1)*(N+2)+j]);
     xnew[i*(N+2)+j]*=recipdiag;
   }
 for(i=start;i<start+rows;i++)
   for(j=1;j<N+1;j++){
     resid[i*(N+2)+j]=b[i*(N+2)+j]-diag*xnew[i*(N+2)+j]-odiag*(xnew[i*(N+2)+j+1]+xnew[i*(N+2)+j-1]+xnew[(i-1)*(N+2)+j]+xnew[(i+1)*(N+2)+j]);
   } 
} 
  
void copy(double * xold, double * xnew, int start, int rows)
{ 
 int i,j;
 for(i=start;i<start+rows;i++)
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


