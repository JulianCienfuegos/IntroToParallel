extern "C" __global__ void
mmkernel( float* a, float* b, float* c,
        int pitch_a, int pitch_b, int pitch_c,
        int n, int m, int p )
{
        int i = blockIdx.x*128 + threadIdx.x;
        int j = blockIdx.y*2;
        float sum0 = 0.0, sum1 =0.0, sum2 = 0.0, sum3 = 0.0;
        for( int k = 0; k < p; ++k )
        {
                float b_tmp0 = b[j+pitch_b*k];
                float b_tmp1 = b[j+1+pitch_b*k];
                sum0 += b_tmp0*c[k*pitch_c+i];
                sum1 += b_tmp0*c[k*pitch_c + i + 32];
                sum2 += b_tmp0*c[k*pitch_c + i + 64];
                sum3 += b_tmp0*c[k*pitch_c + i + 96];
                
                sum4 += b_tmp1*c[k*pitch_c+i];
				sum5 += b_tmp1*c[k*pitch_c + i + 32];
                sum6 += b_tmp1*c[k*pitch_c + i + 64];
                sum7 += b_tmp1*c[k*pitch_c + i + 96];
        }
        a[i+pitch_a*j] = sum0;
        a[i+32+pitch_a*j] = sum1;
        a[i+64+pitch_a*j] = sum2;
        a[i+96+pitch_a*j] = sum3;
        
        a[i+pitch_a*(j+1)] = sum4;
        a[i+32+pitch_a*(j+1)] = sum5;
        a[i+64+pitch_a*(j+1)] = sum6;
        a[i+96+pitch_a*(j+1)] = sum7;
}


