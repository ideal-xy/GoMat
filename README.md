# Introoduction

This is a lightweight library for linear algebra, with C++ language ,constructed with CMake.

But it's immature, with a lot of bugs and disadvantages.

My next goal is to fix there bugs and gradually complete it.（2025.08.25）

# Update Log

1.把 `m_data` 的存储形式改为列优先，这是为了方便矩阵乘法

2.手写 SIMD 指令集实现了矩阵乘法

3.修复了 `loadFromCsv` 成员函数的 bug

4.使用SIMD指令集实现大矩阵加法

5.修复了特殊的矩阵填充函数



