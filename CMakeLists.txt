cmake_minimum_required(VERSION 3.15)
set(GOMAT_VERSION "0.0.1")
project(GoMat VERSION ${GOMAT_VERSION} LANGUAGES CXX)

# 全局设置
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)


add_subdirectory(core)
add_subdirectory(dataframe)
add_subdirectory(series)


add_subdirectory(tests)
add_subdirectory(debug)

