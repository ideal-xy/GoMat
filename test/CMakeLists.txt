
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
# 本项目的实现并没有利用第三方库，上面的两个库都是为了性能测试对比使用的（这俩是结合使用的）

add_executable(tes main.cpp)

#link
target_include_directories(tes PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(tes PRIVATE gomat Eigen3::Eigen OpenMP::OpenMP_CXX)

#自动运行测试
add_custom_command(TARGET tes POST_BUILD
    COMMAND $<TARGET_FILE:tes>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "正在运行 GoMat 测试..."
)
