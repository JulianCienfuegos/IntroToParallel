extern "C" __global__ void
mmkernel( float* a, float* b, float* c,
  int pitch_a, int pitch_b, int pitch_c,
  int n, int m, int p )
{
    int tx = threadIdx.x;
    int i = blockIdx.x*32 + tx;
    int j = blockIdx.y*2;
    __shared__ float cb0[32], cb1[32];

    float sum0 = 0.0, sum1 = 0.0;
    for( int ks = 0; ks < p; ks += 32 ){
      cb0[tx] = c[ks+tx+pitch_c*j];
      cb1[tx] = c[ks+tx+pitch_c*(j+1)];
      __syncthreads();
      for( int k = ks; k < ks+32; ++k ){
        float rb = b[i+pitch_b*k];
        sum0 += rb * cb0[k-ks];
        sum1 += rb * cb1[k-ks];
      }
      __syncthreads();
    }
    a[i+pitch_a*j] = sum0;
    a[i+pitch_a*(j+1)] = sum1;
}
