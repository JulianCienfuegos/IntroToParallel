extern "C" __global__ void
mmkernel( float* a, float* b, float* c,
  int pitch_a, int pitch_b, int pitch_c,
  int n, int m, int p )
{
	// This is a shared memory version of k1.
    int tx = threadIdx.x;
    int i = blockIdx.x*32 + tx;
    int j = blockIdx.y;
    __shared__ float cb[32];

    float sum = 0.0;
    for( int ks = 0; ks < p; ks += 32 ){
      cb[tx] = c[ks+tx+pitch_c*j];
      for( int k = ks; k < ks+32; ++k )
        sum += b[i+pitch_b*k] * cb[k-ks];
    }
    a[i+pitch_a*j] = sum;
}
