#pragma once
// 跨平台动态库导出宏定义
#if defined(_WIN32)
    // Windows 平台
    #ifdef GOMAT_EXPORT
        #define GOMAT_API __declspec(dllexport)  // 编译动态库时导出符号
    #else
        #define GOMAT_API __declspec(dllimport)  // 使用动态库时导入符号
    #endif
#else
    // Unix-like 平台 (macOS/Linux)
    #ifdef GOMAT_EXPORT
        #define GOMAT_API __attribute__((visibility("default")))  // 导出符号
    #else
        #define GOMAT_API  // 导入时不需要特殊属性
    #endif
#endif

#ifdef __APPLE__
    #define GOMAT_ALIGNED __attribute__((aligned(16)))  // 内存对齐优化
    #define GOMAT_HOT __attribute__((hot))             // 标记热点函数
    #define GOMAT_COLD __attribute__((cold))           // 标记冷门函数
#else
    #define GOMAT_ALIGNED
    #define GOMAT_HOT
    #define GOMAT_COLD
#endif

// 头文件包含
#include "matrix.h"
#include "vector.h"
#include "linear_solver.h"
#include "matrix_view.h"
#include "util.h"
#include "matrix_stream.h"
#include "vector_stream.h"

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