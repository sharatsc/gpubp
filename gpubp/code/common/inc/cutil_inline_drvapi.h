#ifndef _CUTIL_INLINE_FUNCTIONS_DRVAPI_H_
#define _CUTIL_INLINE_FUNCTIONS_DRVAPI_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// We define these calls here, so the user doesn't need to include __FILE__ and __LINE__
// The advantage is the developers gets to use the inline function so they can debug
#define cutilDrvSafeCallNoSync(err)     __cuSafeCallNoSync  (err, __FILE__, __LINE__)
#define cutilDrvSafeCall(err)           __cuSafeCall        (err, __FILE__, __LINE__)
#define cutilDrvCtxSync()               __cuCtxSync         (__FILE__, __LINE__)
#define cutilDrvCheckMsg(msg)           __cuCheckMsg        (msg, __FILE__, __LINE__)
#define cutilDrvAlignOffset(off, dptr)  __cuAlignMemOffset  (offset, dptr)

inline int __cuAlignMemOffset(int offset, void *dptr)
{
   // adjust offset to meet device memory alignment requirements
   return (offset + __alignof(dptr) - 1) & ~(__alignof(dptr) - 1);
}

// These are the inline versions for all of the CUTIL functions
inline void __cuSafeCallNoSync( CUresult err, const char *file, const int line )
{
    if( CUDA_SUCCESS != err) {
        fprintf(stderr, "cuSafeCallNoSync() Driver API error = %04d from file <%s>, line %i.\n",
                err, file, line );
        exit(-1);
    }
}
inline void __cuSafeCall( CUresult err, const char *file, const int line )
{
    __cuSafeCallNoSync( err, file, line );
}

inline void __cuCtxSync(const char *file, const int line )
{
    CUresult err = cuCtxSynchronize();
    if( CUDA_SUCCESS != err ) {
        fprintf(stderr, "cuCtxSynchronize() API error = %04d in file <%s>, line %i.\n",
                err, file, line );
        exit(-1);
    }
}

inline void __cuCheckMsg( const char * msg, const char *file, const int line )
{
    CUresult err = cuCtxSynchronize();
    if( CUDA_SUCCESS != err) {
		fprintf(stderr, "cutilDrvCheckMsg -> %s", msg);
        fprintf(stderr, "cutilDrvCheckMsg -> cuCtxSynchronize API error = %04d in file <%s>, line %i.\n",
                err, file, line );
        exit(-1);
    }
}


#if __DEVICE_EMULATION__
    inline void cutilDeviceInitDrv(int cuDevice, int ARGC, char **ARGV) { } 
#else
    inline void cutilDeviceInitDrv(int cuDevice, int ARGC, char ** ARGV) 
    {
        cuDevice = 0;
        int deviceCount = 0;
        CUresult err = cuInit(0);
        if (CUDA_SUCCESS == err)
            cutilDrvSafeCallNoSync(cuDeviceGetCount(&deviceCount));
        if (deviceCount == 0) {
            fprintf(stderr, "CUTIL DeviceInitDrv error: no devices supporting CUDA\n");
            exit(-1);
        }
        int dev = 0;
        cutGetCmdLineArgumenti(ARGC, (const char **) ARGV, "device", &dev);
	    if (dev < 0) dev = 0;
        if (dev > deviceCount-1) dev = deviceCount - 1;
        cutilDrvSafeCallNoSync(cuDeviceGet(&cuDevice, dev));
        char name[100];
        cuDeviceGetName(name, 100, cuDevice);
        if (cutCheckCmdLineFlag(ARGC, (const char **) ARGV, "quiet") == CUTFalse)
            fprintf(stderr, "Using device %d: %s\n", dev, name);
    }
#endif

//! Check for CUDA context lost
inline void cutilDrvCudaCheckCtxLost(const char *errorMessage, const char *file, const int line ) 
{
    CUresult err = cuCtxSynchronize();
    if( CUDA_ERROR_INVALID_CONTEXT != err) {
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i\n",
                errorMessage, file, line );
        exit(-1);
    }
    err = cuCtxSynchronize();
    if( CUDA_SUCCESS != err) {
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i\n",
                errorMessage, file, line );
        exit(-1);
    } 
}

// General check for CUDA GPU SM Capabilities
inline bool cutilDrvCudaCapabilities(int major_version, int minor_version)
{
    int major, minor;
    int dev;
    char device_name[256];

#ifdef __DEVICE_EMULATION__
    printf("> Compute Device Emulation Mode \n");
#endif

    cutilDrvSafeCallNoSync( cuDeviceGet(&dev, 0) );
    cutilDrvSafeCallNoSync( cuDeviceComputeCapability(&major, &minor, dev));
    cutilDrvSafeCallNoSync( cuDeviceGetName(device_name, 256, dev) ); 

    if(major >= major_version && minor >= minor_version)
    {
        printf("> Compute SM %d.%d Device Detected\n", major, minor);
        printf("> Device %d: <%s>\n", dev, device_name);
        return true;
    }
    else
    {
        printf("There is no device supporting CUDA compute capability %d.%d.\n", major_version, minor_version);
        printf("TEST PASSED\n");
        return false;
    }
}

#endif // _CUTIL_INLINE_FUNCTIONS_DRVAPI_H_
