#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#define N 4096
#define Niter 10
#define threshold 0.0000001
#define blocksize 32

double A[N][N], x[N],y[N],z[N],yy[N],zz[N];
int main(){
double rtclock();
void pa1p1(int n, double a[][n], double p[n], double q[n], double r[n]);
void pa1p1opt(int n, double a[][n], double p[n], double q[n], double r[n]);
void compare(int n, double wref[n], double w[n]);

double clkbegin, clkend;
double t;
double rtclock();

int i,j,it;

  for(i=0;i<N;i++)
   { 
     x[i] = i; 
     y[i]= 0; z[i] = 1.0;
     yy[i]= 0; zz[i] = 1.0;
     for(j=0;j<N;j++) A[i][j] = (i+2*j)/(2*N);
   }

  clkbegin = rtclock();
  for(it=0;it<Niter;it++) pa1p1(N,A,x,y,z);
  clkend = rtclock();
  t = clkend-clkbegin;
  if (y[N/2]*y[N/2] < -100.0) printf("%f\n",y[N/2]);
  printf("Problem 1 Reference Version: Matrix Size = %d; %.1f GFLOPS; Time = %.3f sec; \n",
          N,4.0*1e-9*N*N*Niter/t,t);

  clkbegin = rtclock();
  for(it=0;it<Niter;it++) pa1p1opt(N,A,x,yy,zz);
  clkend = rtclock();
  t = clkend-clkbegin;
  if (yy[N/2]*yy[N/2] < -100.0) printf("%f\n",yy[N/2]);
  printf("Problem 1 Optimized Version: Matrix Size = %d; %.1f GFLOPS; Time = %.3f sec; \n",
          N,4.0*1e-9*N*N*Niter/t,t);
  compare(N,y,yy);
  compare(N,z,zz);

}

void pa1p1(int n, double a[][n], double p[n], double q[n], double r[n])
{  int i, j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
      q[i] = q[i] + a[i][j]*p[j];
      r[i] = r[i] + a[j][i]*p[j];
    }
}

void pa1p1opt(int n, double a[][n], double p[n], double q[n], double r[n])
// Initially identical to reference; make your changes to optimize this code
{ int i,j, it, jt;
	
	for(it = 0; it < n; it+=blocksize)
		for(jt = 0; jt < n; jt += blocksize)	
		{
			for(i=it;i<it+blocksize;i+=4)
				for(j=jt;j<jt+blocksize;j++)
				{
				  q[i] += a[i][j]*p[j];
				  q[i+1] += a[i+1][j]*p[j];		  
				  q[i+2] += a[i+2][j]*p[j];
				  q[i+3] += a[i+3][j]*p[j];

				}
			for(i=it;i<it+blocksize;i+=4)
			    for(j=jt;j<jt+blocksize;j++)
				{	  
				  r[j] += a[i][j]*p[i]; 
  				  r[j] += a[i+1][j]*p[i+1]; 
				  r[j] += a[i+2][j]*p[i+2]; 
				  r[j] += a[i+3][j]*p[i+3]; 			  
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

void compare(int n, double wref[n], double w[n])
{
double maxdiff,this_diff;
int numdiffs;
int i;
  numdiffs = 0;
  maxdiff = 0;
  for (i=0;i<n;i++)
    {
     this_diff = wref[i]-w[i];
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
