extern "C" __global__ void
mmkernel( float* a, float* b, float* c,
  int pitch_a, int pitch_b, int pitch_c,
  int n, int m, int p )
{
    int tx = threadIdx.x;
    int i = blockIdx.x*64 + tx;
    int j = blockIdx.y*2;
    __shared__ float cb0[32], cb1[32];

    float sum0 = 0.0, sum1 = 0., sum2 = 0.0, sum3 = 0.0;
    for( int ks = 0; ks < p; ks += 32 ){
      cb0[tx] = c[ks+tx+pitch_c*j];
      cb1[tx] = c[ks+tx+pitch_c*(j+1)];
      __syncthreads();
      for( int k = ks; k < ks+32; ++k ){
        float rb0 = b[i+pitch_b*k];
        float rb1 = b[i+32+pitch_b*k];

        sum0 += rb0 * cb0[k - ks];
        sum1 += rb0 * cb1[k - ks];

        sum2 += rb1*cb0[k - ks];
        sum3 += rb1*cb1[k - ks];

      }
      __syncthreads();
    }
    a[i+pitch_a*j] = sum0;
    a[i+pitch_a*(j+1)] = sum1;
    a[i+32+pitch_a*j] = sum2;
    a[i+32+pitch_a*(j+1)] = sum3;
}

