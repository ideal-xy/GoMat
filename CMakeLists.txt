cmake_minimum_required(VERSION 3.15)
set(GOMAT_VERSION "0.0.1")
project(GoMat VERSION ${GOMAT_VERSION} LANGUAGES CXX)

# 全局设置
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)


add_subdirectory(core)
add_subdirectory(series)

# 可选的测试与调试目录（避免缺失依赖导致配置失败）
if(EXISTS ${CMAKE_SOURCE_DIR}/tests)
  add_subdirectory(tests)
endif()
if(EXISTS ${CMAKE_SOURCE_DIR}/debug)
  add_subdirectory(debug)
endif()

