# --- 创建测试可执行文件 ---
add_executable(tes main.cpp)

# --- 链接库 ---
# 将我们的 gomat 库链接到 tes 可执行文件。
# PRIVATE 意味着链接关系只对 tes 本身有效。
target_link_libraries(tes PRIVATE gomat)

# --- 构建后自动运行测试 ---
# 这个命令会在每次成功编译 tes 目标后，自动执行它。
# 非常适合快速验证代码。
add_custom_command(TARGET tes POST_BUILD
    COMMAND $<TARGET_FILE:tes>
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "正在运行 GoMat 测试..."
)
