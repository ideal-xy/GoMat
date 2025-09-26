
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
# 本项目的实现并没有利用第三方库，上面的两个库都是为了性能测试对比使用的（这俩是结合使用的）

add_executable(tes main.cpp)


# target_include_directories(tes PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(tes
    PRIVATE
    gomat_core          # 链接核心库
    gomat_dataframe     # 链接数据分析库
    Eigen3::Eigen 
    OpenMP::OpenMP_CXX   
)