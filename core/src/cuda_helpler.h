#pragma once

#ifdef __CUDACC__
    #include <cuda_runtime.h>
    #define CUDA_AVALABLE 1
#else
    #define CUDA_AVALABLE 0
#endif

inline bool isCudaAvailable()
{
#if CUDA_AVAILABLE
    return true;
#else
    return false;
#endif

}
