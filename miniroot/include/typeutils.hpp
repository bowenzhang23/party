#include <algorithm>
#include <cinttypes>

enum class Type {
    Float,
    Double,
    Byte,
    Short,
    Integer,
    Long,
};

void print_mem(uint8_t* mem, unsigned short rows);

template <typename T>
std::vector<T> from_bytevec(const std::vector<uint8_t>& byte_vec)
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
