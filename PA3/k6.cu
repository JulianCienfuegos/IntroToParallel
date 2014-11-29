/*
We want to use shared memory, and access coalesced global memory.
*/

extern "C" __global__ void
mmkernel( float* a, float* b, float* c,
        int pitch_a, int pitch_b, int pitch_c,
        int n, int m, int p )
{
        int i = blockIdx.x*64 + threadIdx.x;
        int j = blockIdx.y*2;
        float sum0 = 0.0, sum1 = 0.0, sum2=0.0, sum3=0.0;

        for(int k = 0; k < p; ++k)
        {
                float c_tmp = c[k*pitch_c+i];
                float c_tmp2 = c[k*pitch_c+i+32];
                float b_tmp = b[j+pitch_b*k];
                float b_tmp2 = b[j+1+pitch_b*k];
                sum0 += b_tmp*c_tmp;
                sum1 += b_tmp2*c_tmp;
                sum2 += b_tmp*c_tmp2;
                sum3 += b_tmp2*c_tmp2;
        }

        a[j+pitch_a*i] = sum0;
        a[j+pitch_a*(i+32)] = sum1;
        a[j+1+pitch_a*i]=sum2;
        a[j+1+pitch_a*(i+32)] =sum3;
}


