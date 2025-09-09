#ifndef GOMAT_H
#define GOMAT_H

#include "matrix.h"
#include "vector.h"
#include "linear_solver.h"
#include "matrix_view.h"
#include "util.h"
#include "matrix_stream.h"
#include "vector_stream.h"
#include "mem_alloc.h"

// for macOS 
#ifdef __APPLE__
    #define METAL_SUPPORT 0    
    #define ACCELERATE_SUPPORT 1  
#endif

namespace gomat {
    const char* version() { 
        return "0.0.1"; 
    }
}


#endif

