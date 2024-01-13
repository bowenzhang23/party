#pragma once

#include "zlib.h"

#include <cinttypes>
#include <iostream>

/**
 * @brief Information about how to retrive the bytes.
 */
struct ZipInfo_t {
    int32_t nbytes;
    int32_t key_len;
    int32_t obj_len;
    int     seekpos;
    bool    compressed;

    constexpr ZipInfo_t(int32_t nb, int32_t nk, int32_t no, int sp,
                        bool c) noexcept
        : nbytes(nb), key_len(nk), obj_len(no), seekpos(sp), compressed(c)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const ZipInfo_t& zi)
{
    os << std::hex << "nbytes = " << zi.nbytes << '\n'
       << "key_len = " << zi.key_len << '\n'
       << "obj_len = " << zi.obj_len << '\n'
       << "seekpos = " << zi.seekpos << '\n'
       << "compressed = " << zi.compressed << '\n';

    return os;
}

/**
 * @brief Compress the input data using zlib
 *
 * @tparam T Type of output data
 * @param data Input data
 * @param data_size Size of input data
 * @param comp Output data
 * @param comp_size Size of output data
 * @return int Status
 */
template <typename T>
[[maybe_unused]] bool zlib_compress(const T* data, int32_t data_size, T*& comp,
                                    int32_t& comp_size)
{
    z_stream stream;
    stream.zalloc = Z_NULL;
    stream.zfree  = Z_NULL;
    stream.opaque = Z_NULL;

    if (deflateInit(&stream, Z_DEFAULT_COMPRESSION) != Z_OK) {
        std::cout << "Failed to initialize zlib for compression.\n";
        return false;
    }

    stream.next_in  = reinterpret_cast<Bytef*>(const_cast<T*>(data));
    stream.avail_in = data_size;

    do {
        stream.avail_out = comp_size - stream.total_out;
        stream.next_out  = reinterpret_cast<Bytef*>(comp + stream.total_out);

        if (deflate(&stream, Z_FINISH) == Z_STREAM_ERROR) {
            std::cout << "Compression error with zlib.\n";
            deflateEnd(&stream);
            return false;
        }

        // Resize the buffer if needed
        if (stream.avail_out == 0) {
            comp_size *= 2;
            T* temp = new T[comp_size];
            std::copy(comp, comp + stream.total_out, temp);
            delete[] comp;
            comp = temp;
        }
    } while (stream.avail_out == 0);

    deflateEnd(&stream);

    if (comp_size != stream.total_out) {
        T* comp_final = new T[stream.total_out];
        std::copy(comp, comp + stream.total_out, comp_final);
        delete[] comp;
        comp      = comp_final;
        comp_size = stream.total_out;
    }
    return true;
}

/**
 * @brief Uncompress the input data using zlib
 *
 * @tparam T Type of output data
 * @param data Input data
 * @param data_size Size of input data
 * @param decomp Output data
 * @param decomp_size Size of output data
 * @return int Status
 */
template <typename T>
[[maybe_unused]] bool zlib_uncompress(const T* data, int32_t data_size,
                                      T*& decomp, int32_t& decomp_size)
{
    z_stream stream;
    stream.zalloc = Z_NULL;
    stream.zfree  = Z_NULL;
    stream.opaque = Z_NULL;

    if (inflateInit(&stream) != Z_OK) {
        std::cout << "Failed to initialize zlib for decompression.\n";
        return false;
    }

    stream.next_in  = reinterpret_cast<Bytef*>(const_cast<T*>(data));
    stream.avail_in = data_size;

    do {
        stream.avail_out = decomp_size - stream.total_out;
        stream.next_out  = reinterpret_cast<Bytef*>(decomp + stream.total_out);

        if (inflate(&stream, Z_FINISH) == Z_STREAM_ERROR) {
            std::cout << "Decompression error with zlib.\n";
            inflateEnd(&stream);
            return false;
        }

        // resize
        if (stream.avail_out == 0) {
            decomp_size *= 2;
            T* temp = new T[decomp_size];
            std::copy(decomp, decomp + stream.total_out, temp);
            delete[] decomp;
            decomp = temp;
        }
    } while (stream.avail_out == 0);

    inflateEnd(&stream);

    if (decomp_size != stream.total_out) {
        T* decomp_final = new T[stream.total_out];
        std::copy(decomp, decomp + stream.total_out, decomp_final);
        delete[] decomp;
        decomp      = decomp_final;
        decomp_size = stream.total_out;
    }

    return true;
}