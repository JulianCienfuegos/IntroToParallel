extern "C" __global__ void
mmkernel( float* a, float* b, float* c,
  int pitch_a, int pitch_b, int pitch_c,
  int n, int m, int p )
{
    int i = blockIdx.x*32 + threadIdx.x;
    int j = blockIdx.y;

    float sum = 0.0;
    for( int k = 0; k < p; ++k ) // p iterations
      sum += b[i+pitch_b*k] * c[k+pitch_c*j]; // 2 flops. 2p mem accesses per thread.
    a[i+pitch_a*j] = sum; //1 more global memory access
}
