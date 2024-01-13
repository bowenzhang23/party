#pragma once

#include <cinttypes>
#include <vector>

/**
 * @brief Print memory layout
 *
 * @param mem Pointer to the begin
 * @param rows how many rows to show
 */
void print_mem(uint8_t* mem, unsigned rows);

/**
 * @brief Translate vector of bytes to vector of given data type
 *
 * @tparam T Desired data type
 * @param byte_vec The original vector of bytes
 * @return std::vector<T> The results vector of T
 */
template <typename T>
inline std::vector<T> from_bytevec(const std::vector<uint8_t>& byte_vec)
{
    std::vector<T> array;

    constexpr std::size_t sz = sizeof(T);
    union floating_point_u {
        T      f_val;
        int8_t bytes[sz];
    } fu;
    for (std::size_t i = 0; i < byte_vec.size(); i += sz) {
        for (std::size_t j = 0; j < sz; ++j) {
            fu.bytes[sz - j - 1] = byte_vec[i + j];
        }
        array.push_back(fu.f_val);
    }
    return array;
}