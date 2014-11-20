#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#define N 4096
#define Niter 3
#define threshold 0.0000001
#define blocksize 512

double A[N][N], AA[N][N], B[N][N], BB[N][N];
int main(){
printf("blocksize = %d\n", blocksize);
double rtclock();
void pa1p2(int n, double a[][n], double b[][n]);
void pa1p2opt(int n, double a[][n], double b[][n], int is_sym);
void compare(int n, double wref[][n], double w[][n]);

double clkbegin, clkend;
double t;
double rtclock();

int i,j,it;

  for(i=0;i<N;i++)
   for(j=0;j<N;j++) A[i][j] = (i+2*j)/(2*N);

  clkbegin = rtclock();
  for(it=0;it<Niter;it++) { pa1p2(N,A,B); pa1p2(N,B,A);}
  clkend = rtclock();
  t = clkend-clkbegin;
  if (B[N/2][N/2]*B[N/2][N/2] < -100.0) printf("%f\n",B[N/2][N/2]);
  printf("Problem 1 Reference Version: Matrix Size = %d; %.2f GFLOPS; Time = %.3f sec; \n",
          N,4.0*1e-9*N*N*Niter/t,t);

  for(i=0;i<N;i++)
   for(j=0;j<N;j++) AA[i][j] = (i+2*j)/(2*N);

  clkbegin = rtclock();
  for(it=0;it<Niter;it++) {pa1p2opt(N,AA,BB, 0);pa1p2opt(N,BB,AA, 1);}
  clkend = rtclock();
  t = clkend-clkbegin;
  if (BB[N/2][N/2]*BB[N/2][N/2] < -100.0) printf("%f\n",BB[N/2][N/2]);
  printf("Problem 1 Optimized Version: Matrix Size = %d; %.2f GFLOPS; Time = %.3f sec; \n",
          N,4.0*1e-9*N*N*Niter/t,t);
  compare(N,B,BB);
  compare(N,A,AA);
}

void pa1p2(int n, double a[][n], double b[][n])
{ int i,j;
  for(i=0;i<n;i++)
   for(j=0;j<n;j++)
    b[i][j] = 0.5*(a[i][j] + a[j][i]);
}

void pa1p2opt(int n, double a[][n], double b[][n], int is_sym)
{
/*
 * For this assignment I try a combination of outerloop unrolling, which 
 * get more rows in the l1 cache at once, and tiling, to make sure that
 * whole tiles can be accessed efficiently. 
 * For this function, a blocksize of 32 works best because all 3 
 * 2d blocks of 8 byte words need 24576 bytes of cache memory, and the L1
 * cache size is 32Kb.
 * The second part of this function uses an additional argument "is_sym"
 * to decide whether the matrix we are reding from is symmetric. If it 
 * is, the we don't have to access a_ij and a_ji, because they are the 
 * same. We dont need tiling if we are just copying one array to another. 
 */	
	 int i,j, it, jt;
	if(is_sym == 0)
	{
	for(it = 0; it < n; it += blocksize)
		for(jt = 0; jt < n; jt += blocksize)
			for(i=it;i<it+blocksize;i+=4)
				for(j=jt;j<jt+blocksize;j++)
				{
					b[i][j] = 0.5*(a[i][j] + a[j][i]);
					b[i+1][j] = 0.5*(a[i+1][j] + a[j][i+1]);
					b[i+2][j] = 0.5*(a[i+2][j] + a[j][i+2]);
					b[i+3][j] = 0.5*(a[i+3][j] + a[j][i+3]);
				}
	}
	else
	{
	for(it = 0; i < n; i ++)
		for(j = 0; j < n; j ++)
		{	
			b[i][j] = a[i][j];
			b[i][j+1] = a[i][j+1];
			b[i][j+2] = a[i][j+2];
			b[i][j+3] = a[i][j+3];
		}
	}
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

void compare(int n, double wref[][n], double w[][n])
{
double maxdiff,this_diff;
int numdiffs;
int i,j;
  numdiffs = 0;
  maxdiff = 0;
  for (i=0;i<n;i++)
   for (j=0;j<n;j++)
    {
     this_diff = wref[i][j]-w[i][j];
     if (this_diff < 0) this_diff = -1.0*this_diff;
     if (this_diff>threshold)
      { numdiffs++;
        if (this_diff > maxdiff) maxdiff=this_diff;
      }
    }
   if (numdiffs > 0)
      printf("%d Diffs found over threshold %f; Max Diff = %f\n",
               numdiffs,threshold,maxdiff);
   else
      printf("No differences found between base and test versions\n");
}
