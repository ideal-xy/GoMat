# 创建测试可执行文件
add_executable(tes main.cpp)

# 关键链接配置
target_link_libraries(tes PRIVATE gomat)

# 确保包含路径正确，这一行确保我们可以找到gomat.h 
target_include_directories(tes PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
)


# 添加编译完成后的自动运行
add_custom_command(TARGET tes POST_BUILD
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/tes
    COMMENT "Running tests..."
)