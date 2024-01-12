#include <cinttypes>
#include <iostream>

#include "zlib.h"

struct zip_info {
    int32_t nbytes;
    int32_t key_len;
    int32_t obj_len;
    int     seekpos;
    bool    compressed;

    constexpr zip_info(int32_t nb, int32_t nk, int32_t no, int sp,
                       bool c) noexcept
        : nbytes(nb), key_len(nk), obj_len(no), seekpos(sp), compressed(c)
    {
    }
};

std::ostream& operator<<(std::ostream& os, const zip_info& zi)
{
    os << std::hex << "nbytes = " << zi.nbytes << '\n'
       << "key_len = " << zi.key_len << '\n'
       << "obj_len = " << zi.obj_len << '\n'
       << "seekpos = " << zi.seekpos << '\n'
       << "compressed = " << zi.compressed << '\n';

    return os;
}

template <typename T>
int zlib_unzip(const T* data, int32_t data_size, T*& decomp,
               int32_t& decomp_size)
{
    z_stream stream;
    stream.zalloc = Z_NULL;
    stream.zfree  = Z_NULL;
    stream.opaque = Z_NULL;

    if (inflateInit(&stream) != Z_OK) {
        std::cout << "Failed to initialize zlib for decompression."
                  << std::endl;
        return -1;
    }

    stream.next_in  = reinterpret_cast<Bytef*>(const_cast<T*>(data));
    stream.avail_in = data_size;

    do {
        stream.avail_out = decomp_size - stream.total_out;
        stream.next_out  = reinterpret_cast<Bytef*>(decomp + stream.total_out);

        if (inflate(&stream, Z_FINISH) == Z_STREAM_ERROR) {
            std::cerr << "Decompression error with zlib." << std::endl;
            inflateEnd(&stream);
            return -1;
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

    // resize to final out size
    T* finaldecomp = new T[stream.total_out];
    std::copy(decomp, decomp + stream.total_out, finaldecomp);
    delete[] decomp;
    decomp      = finaldecomp;
    decomp_size = stream.total_out;

    return 0;
}