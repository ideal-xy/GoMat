# in dataframe/CMakeLists.txt


add_library(gomat_dataframe STATIC
    src/dataframe.cpp
    src/io.cpp # 例如，CSV读写功能也在这里
)

# 指定头文件目录
target_include_directories(gomat_dataframe PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

#！
# 因为 DataFrame 是由 Series 组成的，所以必须包含 Series 的头文件并链接它。
# 同样使用 PUBLIC，让最终用户也能访问到构成 DataFrame 的 Series 类型。
target_link_libraries(gomat_dataframe PUBLIC gomat_series)

