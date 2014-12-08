#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MSGSIZE 1
#define SERVER 0 // This is the processor that serves the ball                                                                                                                                                \
																																														   
#define PLAYER2 1
#define NITER 10
#define BEGIN 1
#define COUNT 1

double clkbegin, clkend;
double rtclock();

int main(int argc, char * argv[])
{
int taskid,
numworkers,
numtasks,
i; // Loop Variable                                                                                                                                                                                        

double bandwidth,
time,
time2;


double msg1[MSGSIZE]={0};
double msg2[MSGSIZE]={0}; // Messages to pass.                                                                                                                                                              \
																																														   

MPI_Status status;
MPI_Request request;
MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
numworkers = numtasks-1;
if(taskid == SERVER)
{
	/************************* Master Code ********************************/
	printf("**************A FRIENDLY GAME OF PING PONG**************\n");
	printf("The message size is %u\n", MSGSIZE);
	MPI_Barrier(MPI_COMM_WORLD);
	clkbegin = rtclock();
	for(i = 0; i < NITER; i++)
	{
	MPI_Isend(&msg1, MSGSIZE, MPI_DOUBLE, PLAYER2, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Irecv(&msg2, MSGSIZE, MPI_DOUBLE, PLAYER2, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Isend(&msg2, MSGSIZE, MPI_DOUBLE, PLAYER2, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Irecv(&msg1, MSGSIZE, MPI_DOUBLE, PLAYER2, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	}
	clkend = rtclock();
	time = (clkend - clkbegin)/(4*NITER);
	bandwidth = 8*MSGSIZE/time;
	printf("I am the SERVER because my taskid is %u\n", taskid);
	
																																					  
	printf("The total time is %f seconds\n", time);
	printf("Bandwidth is %f\n", bandwidth);
	
	MPI_Irecv(&time2, 1, MPI_DOUBLE, PLAYER2, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	bandwidth = 8*MSGSIZE/time2;
	printf("\nThe total time for PLAYER2 is %f\n", time2);
	printf("Bandwidth is %f\n", bandwidth);
	printf("***********************GAME OVER************************\n");
	MPI_Finalize();
} 
if(taskid != SERVER) {
	/************************* Worker Code ********************************/
	MPI_Barrier(MPI_COMM_WORLD);
	clkbegin = rtclock();
	for(i = 0; i < NITER; i++)
	{
	MPI_Irecv(&msg1, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Isend(&msg1, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Irecv(&msg2, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Isend(&msg2, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	}
	clkend = rtclock();
	time2 = (clkend - clkbegin)/(4*NITER);
	MPI_Isend(&time2, 1, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	
	MPI_Finalize();
	
	} if(taskid != SERVER) {
	/************************* Worker Code ********************************/
	MPI_Barrier(MPI_COMM_WORLD);
	clkbegin = rtclock();
	for(i = 0; i < NITER; i++)
	{
	MPI_Irecv(&msg1, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Isend(&msg1, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Irecv(&msg2, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	MPI_Isend(&msg2, MSGSIZE, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Wait(&request, &status);
	}
	clkend = rtclock();
	time2 = (clkend - clkbegin)/(4*NITER);
	MPI_Isend(&time2, 1, MPI_DOUBLE, SERVER, BEGIN, MPI_COMM_WORLD, &request);
	MPI_Finalize();
}
}
/************************* Timer Code *********************************/
double rtclock()
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

