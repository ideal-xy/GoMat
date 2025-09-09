
#ifndef MEM_ALLOC_H
#define MEM_ALLOC_H
// 我们希望我们给std::vector分配的内存都是对齐的，因此需要实现一个分配器
#if defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h> 
#elif defined(__aarch64__)
    #include <arm_neon.h>   
#endif

#include <memory>
namespace gomat{

template <typename T, size_t Alignment>
class SpecialAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    // 传播特性
    using propagate_on_container_copy_assignment = std::true_type;
    using propagate_on_container_move_assignment = std::true_type;
    using propagate_on_container_swap = std::true_type;
    using is_always_equal = std::true_type;

    SpecialAllocator() noexcept = default;

    template <typename U>
    SpecialAllocator(const SpecialAllocator<U, Alignment>&) noexcept {}

   T* allocate(size_t n) 
   {
        const size_t bytes = n * sizeof(T);
        if (bytes == 0) 
        {
            return nullptr;
        }

        void* ptr = nullptr;

    #if defined(__x86_64__) || defined(_M_X64)
        ptr = _mm_malloc(bytes, Alignment);

    #elif defined(_WIN32)
        ptr = _aligned_malloc(bytes, Alignment);

    #else 
        if (posix_memalign(&ptr, Alignment, bytes) != 0) 
        {
            ptr = nullptr; 
        }
    #endif

        if (!ptr) {
            throw std::bad_alloc();
        }
        return static_cast<T*>(ptr);
    }

    void deallocate(T* ptr, size_t /*n*/) noexcept {
        if (!ptr) {
            return;
        }

    #if defined(__x86_64__) || defined(_M_X64)
        _mm_free(ptr);
    #elif defined(_WIN32)
        _aligned_free(ptr);
    #else
        free(ptr);
    #endif
    }


    template <typename U>
    struct rebind {
        using other = SpecialAllocator<U, Alignment>;
    };

    bool operator==(const SpecialAllocator&) const noexcept { return true; }
    bool operator!=(const SpecialAllocator&) const noexcept { return false; }
};

}


#endif