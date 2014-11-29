/*
 *       Copyright (C) 2008, STMicroelectronics, Incorporated.
 *       All rights reserved.
 *
 *         THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT
 *  WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT
 *  NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR
 *  FITNESS FOR A PARTICULAR PURPOSE. 
 *
 * Driver for the CUDA matmul kernels
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

extern void exit(int);

float *a, *b, *c;
float *agold;

void
mmgold( float* a, int an, float* b, int bn, float* c, int cn, int n1, int n2, int n3 )
{
    int i, j, k;
    for( j = 0; j < n2; ++j )
	for( k = 0; k < n3; ++k )
	    for( i = 0; i < n1; ++i )
		a[j+i*an] += b[i+k*bn]*c[j+k*cn];
}

void
docheck( float* a, float* agold, int an, int n1, int n2 )
{
    int i, j, err;
    err = 0;
    for( i = 0; i < n1; ++i ){
	for( j = 0; j < n2; ++j ){
	    if( a[i+j*an] != agold[i+j*an] ){
		++err;
		if( err < 10 ){
		    fprintf( stderr, "a[%d,%d] = %f, should be %f\n",
			i, j, a[i+j*an], agold[i+j*an] );
		}
	    }
	}
    }
    if( err ){
	fprintf( stderr, "%d errors found\n", err );
    }else{
	fprintf( stderr, "no errors found\n" );
    }
}

void
help()
{
    fprintf( stderr, "%s -bin binfile -size arraysize -mat matrixsize -block gangsize1 gangsize2 -thread threadsize1 threadsize2 [-debug -check -print]\n", "mm" );
}

typedef int CUresult;
typedef unsigned long CUdeviceptr;
typedef int CUdevice;
typedef void* CUmodule;
typedef void* CUcontext;
typedef void* CUfunction;

extern CUresult cuInit( unsigned int flags );
extern CUresult cuDeviceGet( CUdevice* dev, int index );
extern CUresult cuCtxCreate( CUcontext* ctx, unsigned int Flags, CUdevice dev );
extern CUresult cuModuleLoad( CUmodule* mod, const char* filename );
extern CUresult cuModuleGetFunction( CUfunction* func, CUmodule mod, char* funcname );
extern CUresult cuMemAlloc( CUdeviceptr* devptr, unsigned int bytes );
extern CUresult cuMemcpyHtoD( CUdeviceptr dstptr, void* srcptr, unsigned int bytes );
extern CUresult cuMemcpyDtoH( void* dstptr, CUdeviceptr srcptr, unsigned int bytes );
extern CUresult cuCtxSynchronize( void );

extern CUresult cuParamSetv( CUfunction func, unsigned int offset, void* ptr, unsigned int bytes );
extern CUresult cuLaunchGrid( CUfunction func, int width, int height );
extern CUresult cuMemFree( CUdeviceptr devptr );
extern CUresult cuModuleUnload( CUmodule mod );
extern CUresult cuCtxDestroy( CUcontext ctx );
int
main( int argc, char* argv[] )
{
    CUdevice device;
    CUcontext context;
    CUmodule module;
    CUfunction func;

    char* binfile;
    int nsize, matsize, memsize, debug, check, print, gs[2], ts[2];
    int i, j, ar;
    int iter, niter = 1;
    CUdeviceptr ap, bp, cp, ip;
    union{
	int i[16];
	long l[8];
    }args;
    struct timeval t1, t2, ta, tb;
    long msec1, msec2, mseca, msecb;
    float flop, mflop, gflop, mflopa, gflopa;
    float gflops[10];

    if( argc < 7 ){
	help();
	return 1;
    }

    ar = 0;
    matsize = 0;
    nsize = 0;
    debug = 0;
    check = 0;
    print = 0;
    binfile = NULL;
    gs[0] = gs[1] = 1;
    ts[0] = ts[1] = 1;
    while( ++ar < argc ){
	if( strncmp( argv[ar], "-bin", 4 ) == 0 ){
	    binfile = argv[++ar];
	}else if( strncmp( argv[ar], "-block", 6 ) == 0 ||
		  strncmp( argv[ar], "-gang", 5 ) == 0 ){
	    int ng = 0;
	    gs[0] = gs[1] = 1;
	    while( ar < argc-1 && argv[ar+1][0] != '-' && ng < 2 )
		gs[ng++] = atoi( argv[++ar] );
	}else if( strncmp( argv[ar], "-thread", 7 ) == 0 ){
	    int nt = 0;
	    ts[0] = ts[1] = 1;
	    while( ar < argc-1 && argv[ar+1][0] != '-' && nt < 2 )
		ts[nt++] = atoi( argv[++ar] );
	}else if( strncmp( argv[ar], "-size", 5 ) == 0 ){
	    nsize = atoi( argv[++ar] );
	}else if( strncmp( argv[ar], "-mat", 4 ) == 0 ){
	    matsize = atoi( argv[++ar] );
	}else if( strncmp( argv[ar], "-debug", 6 ) == 0 ){
	    debug = 1;
	}else if( strncmp( argv[ar], "-check", 6 ) == 0 ){
	    check = 1;
	}else if( strncmp( argv[ar], "-print", 6 ) == 0 ){
	    print = 1;
	}else if( strncmp( argv[ar], "-iter", 5 ) == 0 ){
	    niter = atoi( argv[++ar] );
	    if( niter < 1 ) niter = 1;
	    if( niter > 10 ) niter = 10;	/* sanity? */
	}
    }
    if( niter < 1 ) niter = 1;

    if( nsize == 0 ) nsize = matsize;
    if( matsize == 0 ) matsize = nsize;
    if( nsize <= 0 && matsize <= 0 ){
	help();
	fprintf( stderr, "no size given\n" );
	return 1;
    }

    if( binfile == NULL ){
	help();
	fprintf( stderr, "no binfile given\n" );
	return 1;
    }

    if( ts[0] <= 1 ){
	help();
	fprintf( stderr, "thread size must be > 1\n" );
    }

    printf( "binfile=%s  array=%dx%d  matrix=%dx%d  block=<%dx%d>  thread=<%dx%d>\n",
	binfile, nsize, nsize, matsize, matsize, gs[0], gs[1], ts[0], ts[1] );


    memsize = nsize*nsize*sizeof(float);
    a = (float*) malloc( memsize );
    b = (float*) malloc( memsize );
    c = (float*) malloc( memsize );

    for( j = 0; j < matsize; ++j ){
	for( i = 0; i < matsize; ++i ){
	    a[i+nsize*j] = 0.0;	/* a(i,j) = 0 */
	    b[i+nsize*j] = i*2+j*2;	/* b(i,j) = i+j*2 */
	    c[i+nsize*j] = 0.0;	/* c(i,j) = 0 */
	}
	/* c(j,j+1) = c(j,j-1) = 2 */
	if( j > 0 )
	    c[j-1+nsize*j] = -1.0;
	if( j < matsize-1 )
	    c[j+1+nsize*j] = 2.0;
	
    }

    /* initialize runtime */
    cuInit(0);
    cuDeviceGet( &device, 0 );
    cuCtxCreate( &context, 0, device );

    flop = matsize;
    flop = flop*flop*flop*2.0;

    printf( 
	    "matrix = %dx%d\n"
	    " array = %dx%d\n"
	    "  grid = %dx%d\n"
	    " block = %dx%dx1\n"
	    " flops = %.0f\n",
	matsize, matsize, nsize, nsize, gs[0], gs[1], ts[0], ts[1], flop );
    for( iter = 0; iter < niter; ++iter ){
	gettimeofday( &t1, NULL );
	msec1 = t1.tv_sec * 1000000 + t1.tv_usec;

	/* load module */
	cuModuleLoad( &module, binfile );

	cuModuleGetFunction( &func, module, "mmkernel" );

	/* allocate memory */
	bp = 0; ap = 0; cp = 0;
	cuMemAlloc( &bp, memsize );
	cuMemAlloc( &ap, memsize );
	cuMemAlloc( &cp, memsize );

	/* upload memory */
	cuMemcpyHtoD( bp, b, memsize );
	cuMemcpyHtoD( cp, c, memsize );
	cuMemcpyHtoD( ap, a, memsize );
	cuCtxSynchronize();

	gettimeofday( &ta, NULL );
	mseca = ta.tv_sec * 1000000 + ta.tv_usec;

	/* set function arguments */
	args.l[0] = ap;		/*      i[0-1] */
	args.l[1] = bp;		/*      i[2-3] */
	args.l[2] = cp;		/*      i[4-5] */
	args.i[6] = nsize;		/* l[3] */
	args.i[7] = nsize;		/* l[3] */
	args.i[8] = nsize;		/* l[4] */
	args.i[9] = matsize;	/* l[4] */
	args.i[10] = matsize;	/* l[5] */
	args.i[11] = matsize;	/* l[5] */

	cuParamSetv( func, 0, &args, 48 );
	cuParamSetSize( func, 48 );

	cuFuncSetBlockShape( func, ts[0], ts[1], 1 );
	cuLaunchGrid( func, gs[0], gs[1] );
	cuCtxSynchronize();

	gettimeofday( &tb, NULL );
	msecb = tb.tv_sec * 1000000 + tb.tv_usec;
	msecb -= mseca;

	/* download results */
	cuMemcpyDtoH( a, ap, memsize );
	if( print ){
	    char* format = "  a[%3d,%3d] = %9.3f";
	    for( j = 0; j < 4; ++j ){
		for( i = 0; i < 2; ++i ){
		    printf( format, i, j, a[i+j*nsize] );
		}
		printf( "  " );
		for( i = matsize-2; i < matsize; ++i ){
		    printf( format, i, j, a[i+j*nsize] );
		}
		printf( "\n" );
	    }
	    for( j = matsize-4; j < matsize; ++j ){
		for( i = 0; i < 2; ++i ){
		    printf( format, i, j, a[i+j*nsize] );
		}
		printf( "  " );
		for( i = matsize-2; i < matsize; ++i ){
		    printf( format, i, j, a[i+j*nsize] );
		}
		printf( "\n" );
	    }
	}

	/* free memory */
	cuMemFree( ap );
	cuMemFree( bp );
	cuMemFree( cp );

	/* unload module and quit */
	cuModuleUnload( module );

	gettimeofday( &t2, NULL );
	msec2 = t2.tv_sec * 1000000 + t2.tv_usec;

	msec2 -= msec1;
	mflop = flop / (float)msec2;
	gflop = mflop / 1000.;
	mflopa = flop / (float)msecb;
	gflopa = mflopa / 1000.;
	printf( 
		"  msec = %10ld   GFLOPS = %7.2f, %7.2f (kernel)\n",
	    msec2, gflop, gflopa );
	gflops[iter] = gflopa;
    }

    cuCtxDestroy( context );
    if( check ){
	agold = (float*) malloc( memsize );
	for( j = 0; j < matsize; ++j )
	    for( i = 0; i < matsize; ++i )
		agold[i+nsize*j] = 0.0;
	mmgold( agold, nsize, b, nsize, c, nsize, matsize, matsize, matsize );
	docheck( a, agold, nsize, matsize, matsize );
    }
    exit(0);
}
