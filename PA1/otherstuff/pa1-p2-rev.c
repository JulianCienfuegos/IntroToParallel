#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#define N 4096
#define Niter 3
#define threshold 0.0000001

double A[N][N], B[N][N], BB[N][N];
int main(){
double rtclock();
void pa1p2(int n, double a[][n], double b[][n]);
void pa1p2opt(int n, double a[][n], double b[][n]);
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
   for(j=0;j<N;j++) A[i][j] = (i+2*j)/(2*N);

  clkbegin = rtclock();
  for(it=0;it<Niter;it++) {pa1p2opt(N,A,BB);pa1p2opt(N,BB,A);}
  clkend = rtclock();
  t = clkend-clkbegin;
  if (BB[N/2][N/2]*BB[N/2][N/2] < -100.0) printf("%f\n",BB[N/2][N/2]);
  printf("Problem 1 Optimized Version: Matrix Size = %d; %.2f GFLOPS; Time = %.3f sec; \n",
          N,4.0*1e-9*N*N*Niter/t,t);
  compare(N,B,BB);

}

void pa1p2(int n, double a[][n], double b[][n])
{ int i,j;
  for(i=0;i<n;i++)
   for(j=0;j<n;j++)
    b[i][j] = 0.5*(a[i][j] + a[j][i]);
}

void pa1p2opt(int n, double a[][n], double b[][n])
// Initially identical to reference; make your changes to optimize this code
{ int i,j;
  for(i=0;i<n;i++)
   for(j=0;j<n;j++)
    b[i][j] = 0.5*(a[i][j] + a[j][i]);
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
