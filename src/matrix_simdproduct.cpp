#ifdef __x86_64__
#include "matrix.h"
#include <thread>
#include <immintrin.h>


namespace gomat{

inline void transpose_4x4(__m256d& in0, __m256d& in1, __m256d& in2, __m256d& in3,
                          __m256d& out0, __m256d& out1, __m256d& out2, __m256d& out3) {
    __m256d tmp0 = _mm256_unpacklo_pd(in0, in1);
    __m256d tmp1 = _mm256_unpackhi_pd(in0, in1);
    __m256d tmp2 = _mm256_unpacklo_pd(in2, in3);
    __m256d tmp3 = _mm256_unpackhi_pd(in2, in3);

    out0 = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
    out1 = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
    out2 = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
    out3 = _mm256_permute2f128_pd(tmp1, tmp3, 0x31);
}

inline void kernel_4x4(size_t start_row, size_t start_col, size_t k_start, size_t k_end, 
    Matrix& result, const double* packed_left, const double* packed_right, const Matrix& right)
{
    __m256d res_col_0 = _mm256_load_pd(&result(start_row, start_col));
    __m256d res_col_1 = _mm256_load_pd(&result(start_row, start_col + 1));
    __m256d res_col_2 = _mm256_load_pd(&result(start_row, start_col + 2));
    __m256d res_col_3 = _mm256_load_pd(&result(start_row, start_col + 3));

    // 转置到行主序计算
    __m256d row0, row1, row2, row3;
    // transpose_4x4(res_col_0, res_col_1, res_col_2, res_col_3, row0, row1, row2, row3);

    __m256d tmp0 = _mm256_unpacklo_pd(res_col_0, res_col_1);
    __m256d tmp1 = _mm256_unpackhi_pd(res_col_0, res_col_1);
    __m256d tmp2 = _mm256_unpacklo_pd(res_col_2, res_col_3);
    __m256d tmp3 = _mm256_unpackhi_pd(res_col_2, res_col_3);

    row0 = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
    row1 = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
    row2 = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
    row3 = _mm256_permute2f128_pd(tmp1, tmp3, 0x31);

    size_t k_len = k_end - k_start;

    size_t i = 0;
    for (; i+3 < k_len; i+=4) 
    {
        __m256d left_col_i = _mm256_load_pd(&packed_left[i * 4]);
        __m256d right_vec = _mm256_load_pd(&packed_right[i * 4]);
        __m256d l0 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 0);
        __m256d l1 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 1);
        __m256d l2 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 2);
        __m256d l3 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 3);
        row0 = _mm256_fmadd_pd(l0, right_vec, row0);
        row1 = _mm256_fmadd_pd(l1, right_vec, row1);
        row2 = _mm256_fmadd_pd(l2, right_vec, row2);
        row3 = _mm256_fmadd_pd(l3, right_vec, row3);

        
        left_col_i = _mm256_load_pd(&packed_left[(i+1) * 4]);
        right_vec = _mm256_load_pd(&packed_right[(i+1) * 4]);
        l0 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 0);
        l1 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 1);
        l2 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 2);
        l3 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 3);
        row0 = _mm256_fmadd_pd(l0, right_vec, row0);
        row1 = _mm256_fmadd_pd(l1, right_vec, row1);
        row2 = _mm256_fmadd_pd(l2, right_vec, row2);
        row3 = _mm256_fmadd_pd(l3, right_vec, row3);

        left_col_i = _mm256_load_pd(&packed_left[(i+2) * 4]);
        right_vec = _mm256_load_pd(&packed_right[(i+2) * 4]);
        l0 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 0);
        l1 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 1);
        l2 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 2);
        l3 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 3);
        row0 = _mm256_fmadd_pd(l0, right_vec, row0);
        row1 = _mm256_fmadd_pd(l1, right_vec, row1);
        row2 = _mm256_fmadd_pd(l2, right_vec, row2);
        row3 = _mm256_fmadd_pd(l3, right_vec, row3);

        left_col_i = _mm256_load_pd(&packed_left[(i+3) * 4]);
        right_vec = _mm256_load_pd(&packed_right[(i+3) * 4]);
        l0 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 0);
        l1 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 1);
        l2 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 2);
        l3 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 3);
        row0 = _mm256_fmadd_pd(l0, right_vec, row0);
        row1 = _mm256_fmadd_pd(l1, right_vec, row1);
        row2 = _mm256_fmadd_pd(l2, right_vec, row2);
        row3 = _mm256_fmadd_pd(l3, right_vec, row3);

    }

    for (; i < k_len; i++) 
    {
        __m256d left_col_i = _mm256_load_pd(&packed_left[i * 4]);
        __m256d right_vec = _mm256_load_pd(&packed_right[i * 4]);
        __m256d l0 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 0);
        __m256d l1 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 1);
        __m256d l2 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 2);
        __m256d l3 = _mm256_broadcast_sd(reinterpret_cast<const double*>(&left_col_i) + 3);
        row0 = _mm256_fmadd_pd(l0, right_vec, row0);
        row1 = _mm256_fmadd_pd(l1, right_vec, row1);
        row2 = _mm256_fmadd_pd(l2, right_vec, row2);
        row3 = _mm256_fmadd_pd(l3, right_vec, row3);
    }

    // 转置回列主序
    tmp0 = _mm256_unpacklo_pd(row0, row1);
    tmp1 = _mm256_unpackhi_pd(row0, row1);
    tmp2 = _mm256_unpacklo_pd(row2, row3);
    tmp3 = _mm256_unpackhi_pd(row2, row3);

    res_col_0 = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
    res_col_1 = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
    res_col_2 = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);
    res_col_3 = _mm256_permute2f128_pd(tmp1, tmp3, 0x31);

    _mm256_store_pd(&result(start_row, start_col), res_col_0);
    _mm256_store_pd(&result(start_row, start_col + 1), res_col_1);
    _mm256_store_pd(&result(start_row, start_col + 2), res_col_2);
    _mm256_store_pd(&result(start_row, start_col + 3), res_col_3);
}

inline void kernel_4x4(size_t start_row, size_t start_col, size_t k_start, size_t k_end, 
                       Matrix& result, const double* packed_left, const double* packed_right)                
{
    __m256d res_col_0 = _mm256_load_pd(&result(start_row, start_col));
    __m256d res_col_1 = _mm256_load_pd(&result(start_row, start_col + 1));
    __m256d res_col_2 = _mm256_load_pd(&result(start_row, start_col + 2));
    __m256d res_col_3 = _mm256_load_pd(&result(start_row, start_col + 3));

    for (size_t i = 0; i < k_end - k_start; i++) 
    {
        _mm_prefetch(&packed_left[(i + 2) * 4], _MM_HINT_T0);
        _mm_prefetch(&packed_right[(i + 2) * 4], _MM_HINT_T0);
        __m256d left_col = _mm256_load_pd(&packed_left[i * 4]);
        __m256d right_vec = _mm256_load_pd(&packed_right[i * 4]);
        res_col_0 = _mm256_fmadd_pd(left_col, _mm256_broadcast_sd(reinterpret_cast<const double*>(&right_vec)), res_col_0);
        res_col_1 = _mm256_fmadd_pd(left_col, _mm256_broadcast_sd(reinterpret_cast<const double*>(&right_vec)+1), res_col_1);
        res_col_2 = _mm256_fmadd_pd(left_col, _mm256_broadcast_sd(reinterpret_cast<const double*>(&right_vec)+2), res_col_2);
        res_col_3 = _mm256_fmadd_pd(left_col, _mm256_broadcast_sd(reinterpret_cast<const double*>(&right_vec)+3), res_col_3);
    }

    _mm256_store_pd(&result(start_row, start_col), res_col_0);
    _mm256_store_pd(&result(start_row, start_col + 1), res_col_1);
    _mm256_store_pd(&result(start_row, start_col + 2), res_col_2);
    _mm256_store_pd(&result(start_row, start_col + 3), res_col_3);
}

inline void kernel_small(size_t start_row, size_t start_col, size_t k_start, size_t k_end, size_t r, size_t c, Matrix& result, 
    const Matrix& left, const Matrix& right,const double* packed_right)
{
    const size_t rows = left.getRows();
    __m256d res_cols[4] = { _mm256_setzero_pd(), _mm256_setzero_pd(), _mm256_setzero_pd(), _mm256_setzero_pd() };

    // 加载结果矩阵
    for (size_t j = 0; j < c; ++j)
    {
        res_cols[j] = _mm256_load_pd(&result(start_row, start_col + j));
    }

    for (size_t k = k_start; k < k_end; k++)
    {
        __m256d left_vec = _mm256_set_pd(
        (r > 3 && start_row + 3 < rows) ? left(start_row + 3, k) : 0.0,
        (r > 2 && start_row + 2 < rows) ? left(start_row + 2, k) : 0.0,
        (r > 1 && start_row + 1 < rows) ? left(start_row + 1, k) : 0.0,
        (r > 0 && start_row + 0 < rows) ? left(start_row + 0, k) : 0.0
      );

        __m256d right_vec = _mm256_load_pd(&packed_right[(k - k_start) * 4]);

        for (size_t j = 0; j < c; ++j) 
        {
            __m256d rb = _mm256_broadcast_sd(reinterpret_cast<const double*>(&right_vec) + j);
            res_cols[j] = _mm256_fmadd_pd(left_vec, rb, res_cols[j]);
        }
    }

    // 存储结果
    for (size_t j = 0; j < c; ++j)
    {
        _mm256_storeu_pd(&result(start_row, start_col + j), res_cols[j]);
    }
}


Matrix Matrix::multiplyWithMulThread(const Matrix& other) const
{
    const size_t rows = getRows();
    const size_t cols = other.getCols();
    const size_t invariant = getCols();
    const size_t KERNEL_BLOCK_ROWS = 4;
    const size_t KERNEL_BLOCK_COLS = 4;
    const size_t KERNEL_BLOCK_K = 64;


    Matrix result(rows, cols, true); 

    const unsigned int nums = 8;
    std::vector<std::thread> threads;
    threads.reserve(nums);
    size_t rows_per_thread = (rows + nums - 1) / nums;

    for (unsigned int t = 0; t < nums; ++t) 
    {
        const size_t start_row_for_thread = t * rows_per_thread;
        const size_t end_row_for_thread = std::min(start_row_for_thread + rows_per_thread, rows);

        if (start_row_for_thread >= end_row_for_thread) 
        {
            continue; 
        }
        threads.emplace_back([=, &result, &other, this] 
        {
            std::vector<double,SpecialAllocator<double,32>> packed_left(KERNEL_BLOCK_ROWS * KERNEL_BLOCK_K);
            std::vector<double,SpecialAllocator<double,32>> packed_right(KERNEL_BLOCK_K * KERNEL_BLOCK_COLS);
            for (size_t k = 0;k < invariant;k+= KERNEL_BLOCK_K)
            {
                size_t k_end = std::min(k + KERNEL_BLOCK_K,invariant);
                size_t bk_actual = k_end  - k;
                for (size_t i = start_row_for_thread;i < end_row_for_thread;i+= KERNEL_BLOCK_ROWS)
                {
                    const size_t i_end = std::min(i + KERNEL_BLOCK_ROWS, end_row_for_thread);
                    const size_t num_rows_left = i_end - i;
  
                    if (num_rows_left == KERNEL_BLOCK_ROWS)
                    {
                        for (size_t kk = 0; kk < bk_actual; ++kk)
                        {
                            __m256d vec = _mm256_load_pd(&(*this).at(i, k + kk));
                            _mm256_store_pd(&packed_left[kk * 4], vec);
                        }
                    }

                    for (size_t j = 0; j < cols; j += KERNEL_BLOCK_COLS)
                    {
                        const size_t j_end = std::min(j + KERNEL_BLOCK_COLS, cols);
                        const size_t num_cols = j_end - j;

                        if (num_rows_left == KERNEL_BLOCK_ROWS && num_cols == KERNEL_BLOCK_COLS)
                        {
                            for (size_t kk = 0; kk < bk_actual; ++kk)
                            {
                                packed_right[kk * 4 + 0] = other.at(k + kk, j + 0);
                                packed_right[kk * 4 + 1] = other.at(k + kk, j + 1);
                                packed_right[kk * 4 + 2] = other.at(k + kk, j + 2);
                                packed_right[kk * 4 + 3] = other.at(k + kk, j + 3);
                            }

                            kernel_4x4(i, j, k, k_end, result, packed_left.data(), packed_right.data()); 
                        }
                        else
                        {
                            for (size_t kk = 0; kk < bk_actual; ++kk) 
                            {
                                packed_right[kk * 4 + 0] = (j + 0 < cols) ? other.at(k + kk, j + 0) : 0.0;
                                packed_right[kk * 4 + 1] = (j + 1 < cols) ? other.at(k + kk, j + 1) : 0.0;
                                packed_right[kk * 4 + 2] = (j + 2 < cols) ? other.at(k + kk, j + 2) : 0.0;
                                packed_right[kk * 4 + 3] = (j + 3 < cols) ? other.at(k + kk, j + 3) : 0.0;
                            }
                            kernel_small(i, j, k, k_end, num_rows_left, num_cols, result, *this,other, packed_right.data()); 
                        }     
                    }
                 
                }
            }
        });

    }
    // 3. 等待所有线程完成
    for (auto& thread : threads) 
    { 
        thread.join();
    }

    return result;
}


} // namespace 


#elif defined(__aarch64__)
#include "matrix.h"
#include <arm_neon.h>
#include <thread>
#include <algorithm>
#include <vector>

namespace gomat{

inline void kernel_4x4(size_t start_row, size_t start_col, size_t k_start, size_t k_end,
    Matrix& result, const double* packed_left, const double* packed_right)
{
    // 寄存器分别存储结果矩阵 C 的 4x4 块，每列使用两个 float64x2_t 向量
    // res_col_J_01 表示第 J 列的 0,1行; res_col_J_23 表示第 J 列的 2,3行
    float64x2_t res_col_0_01 = vld1q_f64(&result(start_row + 0, start_col + 0));
    float64x2_t res_col_0_23 = vld1q_f64(&result(start_row + 2, start_col + 0));
    float64x2_t res_col_1_01 = vld1q_f64(&result(start_row + 0, start_col + 1));
    float64x2_t res_col_1_23 = vld1q_f64(&result(start_row + 2, start_col + 1));
    float64x2_t res_col_2_01 = vld1q_f64(&result(start_row + 0, start_col + 2));
    float64x2_t res_col_2_23 = vld1q_f64(&result(start_row + 2, start_col + 2));
    float64x2_t res_col_3_01 = vld1q_f64(&result(start_row + 0, start_col + 3));
    float64x2_t res_col_3_23 = vld1q_f64(&result(start_row + 2, start_col + 3));

    for (size_t i = 0; i < k_end - k_start; i++)
    {
        // 预取数据到缓存
        __builtin_prefetch(&packed_left[(i + 4) * 4], 0, 3);
        __builtin_prefetch(&packed_right[(i + 4) * 4], 0, 3);

        // 加载左矩阵 A 的一列 (4个元素)
        float64x2_t left_col_01 = vld1q_f64(&packed_left[i * 4 + 0]);
        float64x2_t left_col_23 = vld1q_f64(&packed_left[i * 4 + 2]);

        // 加载右矩阵 B 的一行 (4个元素) 的指针
        const double* right_ptr = &packed_right[i * 4];

        // 将 B 矩阵行中的元素广播到单独的向量中
        float64x2_t right_bcast_0 = vdupq_n_f64(right_ptr[0]);
        float64x2_t right_bcast_1 = vdupq_n_f64(right_ptr[1]);
        float64x2_t right_bcast_2 = vdupq_n_f64(right_ptr[2]);
        float64x2_t right_bcast_3 = vdupq_n_f64(right_ptr[3]);

        // 执行外积并累加: C_col += A_col * B_scalar
        // 更新 C 矩阵的第0列
        res_col_0_01 = vfmaq_f64(res_col_0_01, left_col_01, right_bcast_0);
        res_col_0_23 = vfmaq_f64(res_col_0_23, left_col_23, right_bcast_0);

        // 更新 C 矩阵的第1列
        res_col_1_01 = vfmaq_f64(res_col_1_01, left_col_01, right_bcast_1);
        res_col_1_23 = vfmaq_f64(res_col_1_23, left_col_23, right_bcast_1);

        // 更新 C 矩阵的第2列
        res_col_2_01 = vfmaq_f64(res_col_2_01, left_col_01, right_bcast_2);
        res_col_2_23 = vfmaq_f64(res_col_2_23, left_col_23, right_bcast_2);

        // 更新 C 矩阵的第3列
        res_col_3_01 = vfmaq_f64(res_col_3_01, left_col_01, right_bcast_3);
        res_col_3_23 = vfmaq_f64(res_col_3_23, left_col_23, right_bcast_3);
    }

    // 将结果写回内存
    vst1q_f64(&result(start_row + 0, start_col + 0), res_col_0_01);
    vst1q_f64(&result(start_row + 2, start_col + 0), res_col_0_23);
    vst1q_f64(&result(start_row + 0, start_col + 1), res_col_1_01);
    vst1q_f64(&result(start_row + 2, start_col + 1), res_col_1_23);
    vst1q_f64(&result(start_row + 0, start_col + 2), res_col_2_01);
    vst1q_f64(&result(start_row + 2, start_col + 2), res_col_2_23);
    vst1q_f64(&result(start_row + 0, start_col + 3), res_col_3_01);
    vst1q_f64(&result(start_row + 2, start_col + 3), res_col_3_23);
}

// 针对边缘小块的 NEON 内核函数
inline void kernel_small(size_t start_row, size_t start_col, size_t k_start, size_t k_end, size_t r, size_t c, Matrix& result,
    const Matrix& left, const double* packed_right)
{
    const size_t rows = left.getRows();
    float64x2_t res_cols_01[4] = { vdupq_n_f64(0.0), vdupq_n_f64(0.0), vdupq_n_f64(0.0), vdupq_n_f64(0.0) };
    float64x2_t res_cols_23[4] = { vdupq_n_f64(0.0), vdupq_n_f64(0.0), vdupq_n_f64(0.0), vdupq_n_f64(0.0) };

    // 加载结果矩阵的现有值
    // 注意：这里即使 c 或 r 不足4，也加载完整的向量。
    // 这要求 result 矩阵的内存分配能保证边界安全。
    for (size_t j = 0; j < c; ++j)
    {
        res_cols_01[j] = vld1q_f64(&result(start_row, start_col + j));
        res_cols_23[j] = vld1q_f64(&result(start_row + 2, start_col + j));
    }

    for (size_t k = k_start; k < k_end; k++)
    {
        // 手动构建左矩阵 A 的列向量，并处理边界情况
        double left_vals[4];
        left_vals[0] = (r > 0 && start_row + 0 < rows) ? left(start_row + 0, k) : 0.0;
        left_vals[1] = (r > 1 && start_row + 1 < rows) ? left(start_row + 1, k) : 0.0;
        left_vals[2] = (r > 2 && start_row + 2 < rows) ? left(start_row + 2, k) : 0.0;
        left_vals[3] = (r > 3 && start_row + 3 < rows) ? left(start_row + 3, k) : 0.0;
        
        float64x2_t left_vec_01 = vld1q_f64(&left_vals[0]);
        float64x2_t left_vec_23 = vld1q_f64(&left_vals[2]);
        
        const double* right_ptr = &packed_right[(k - k_start) * 4];

        for (size_t j = 0; j < c; ++j)
        {
            float64x2_t rb = vdupq_n_f64(right_ptr[j]);
            res_cols_01[j] = vfmaq_f64(res_cols_01[j], left_vec_01, rb);
            res_cols_23[j] = vfmaq_f64(res_cols_23[j], left_vec_23, rb);
        }
    }

    // 存储结果
    for (size_t j = 0; j < c; ++j)
    {
        vst1q_f64(&result(start_row, start_col + j), res_cols_01[j]);
        vst1q_f64(&result(start_row + 2, start_col + j), res_cols_23[j]);
    }
}


// 主乘法函数
Matrix Matrix::multiplyWithMulThread(const Matrix& other) const
{
    const size_t rows = getRows();
    const size_t cols = other.getCols();
    const size_t invariant = getCols();
    const size_t KERNEL_BLOCK_ROWS = 4;
    const size_t KERNEL_BLOCK_COLS = 4;
    const size_t KERNEL_BLOCK_K = 64;

    Matrix result(rows, cols, true);

    const unsigned int nums = 8;
    std::vector<std::thread> threads;
    threads.reserve(nums);
    size_t rows_per_thread = (rows + nums - 1) / nums;

    for (unsigned int t = 0; t < nums; ++t)
    {
        const size_t start_row_for_thread = t * rows_per_thread;
        const size_t end_row_for_thread = std::min(start_row_for_thread + rows_per_thread, rows);

        if (start_row_for_thread >= end_row_for_thread)
        {
            continue;
        }
        threads.emplace_back([=, &result, &other, this]
        {
            std::vector<double, SpecialAllocator<double, 32>> packed_left(KERNEL_BLOCK_ROWS * KERNEL_BLOCK_K);
            std::vector<double, SpecialAllocator<double, 32>> packed_right(KERNEL_BLOCK_K * KERNEL_BLOCK_COLS);
            for (size_t k = 0; k < invariant; k += KERNEL_BLOCK_K)
            {
                size_t k_end = std::min(k + KERNEL_BLOCK_K, invariant);
                size_t bk_actual = k_end - k;
                for (size_t i = start_row_for_thread; i < end_row_for_thread; i += KERNEL_BLOCK_ROWS)
                {
                    const size_t i_end = std::min(i + KERNEL_BLOCK_ROWS, end_row_for_thread);
                    const size_t num_rows_left = i_end - i;
  
                    if (num_rows_left == KERNEL_BLOCK_ROWS)
                    {
                        for (size_t kk = 0; kk < bk_actual; ++kk)
                        {
                            // 将左矩阵的 4x1 列向量打包
                            packed_left[kk * 4 + 0] = (*this).at(i + 0, k + kk);
                            packed_left[kk * 4 + 1] = (*this).at(i + 1, k + kk);
                            packed_left[kk * 4 + 2] = (*this).at(i + 2, k + kk);
                            packed_left[kk * 4 + 3] = (*this).at(i + 3, k + kk);
                        }
                    }

                    for (size_t j = 0; j < cols; j += KERNEL_BLOCK_COLS)
                    {
                        const size_t j_end = std::min(j + KERNEL_BLOCK_COLS, cols);
                        const size_t num_cols = j_end - j;

                        if (num_rows_left == KERNEL_BLOCK_ROWS && num_cols == KERNEL_BLOCK_COLS)
                        {
                            for (size_t kk = 0; kk < bk_actual; ++kk)
                            {
                                // 打包右矩阵的 1x4 行向量
                                packed_right[kk * 4 + 0] = other.at(k + kk, j + 0);
                                packed_right[kk * 4 + 1] = other.at(k + kk, j + 1);
                                packed_right[kk * 4 + 2] = other.at(k + kk, j + 2);
                                packed_right[kk * 4 + 3] = other.at(k + kk, j + 3);
                            }

                            kernel_4x4(i, j, k, k_end, result, packed_left.data(), packed_right.data()); 
                        }
                        else
                        {
                            for (size_t kk = 0; kk < bk_actual; ++kk) 
                            {
                                packed_right[kk * 4 + 0] = (j + 0 < cols) ? other.at(k + kk, j + 0) : 0.0;
                                packed_right[kk * 4 + 1] = (j + 1 < cols) ? other.at(k + kk, j + 1) : 0.0;
                                packed_right[kk * 4 + 2] = (j + 2 < cols) ? other.at(k + kk, j + 2) : 0.0;
                                packed_right[kk * 4 + 3] = (j + 3 < cols) ? other.at(k + kk, j + 3) : 0.0;
                            }
                            kernel_small(i, j, k, k_end, num_rows_left, num_cols, result, *this, packed_right.data()); 
                        }     
                    }
                }
            }
        });
    }

    for (auto& thread : threads) 
    { 
        thread.join();
    }

    return result;
}




}




#endif