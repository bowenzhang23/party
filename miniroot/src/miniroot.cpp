#include "miniroot.hpp"
#include <cassert>
#include <cstring>
#include <iostream>

Miniroot::Miniroot(const std::string& filename)
    : m_iname(0), m_fs(filename, std::ios::binary | std::ios::in)
{
    std::cout << "Read " << filename << '\n';
    ResetStringBuffer();
    m_fh.read(m_fs);
}

void Miniroot::Read()
{
    auto    cur       = (int) m_fs.beg; // 0
    int32_t increment = m_fh.f_begin;   // 100

    while (!m_fs.eof()) {
        cur += increment;
        if (cur == m_fh.f_end) break;
        m_fs.seekg(cur);
        m_lrh.read(m_fs);
        increment = m_lrh.nbytes;

        // class name
        if (!ReadByIname() || m_sbuf[0] != 'T') continue;
        std::string query_name(m_sbuf);
        // object name
        if (!ReadByIname()) continue;
        query_name += "_";
        query_name += m_sbuf;
        // object title
        if (!ReadByIname()) continue;
        query_name += "_";
        query_name += m_sbuf;

        if (query_name.find("TBasket") != std::string::npos) {
            std::cout << "Reading .. " << query_name << '\n';
            m_basket_map[query_name].emplace_back(
                m_lrh.nbytes, m_lrh.key_len, m_lrh.obj_len,
                cur + m_lrh.key_len + UNZIP_HEADER_LEN,
                m_lrh.obj_len > m_lrh.nbytes - m_lrh.key_len);
        }
    }
}

uint8_t* Miniroot::GetUncompressedBytes(const zip_info& info)
{
    std::cout << info;
    m_fs.seekg(info.seekpos);
    int32_t data_size   = info.nbytes - info.key_len - UNZIP_HEADER_LEN;
    int32_t decomp_size = info.obj_len;

    uint8_t* data = (uint8_t*) malloc(sizeof(uint8_t) * data_size);
    m_fs.read((char*) data, sizeof(uint8_t) * data_size);
    if (!info.compressed) return data;

    uint8_t* decomp = (uint8_t*) malloc(sizeof(uint8_t) * decomp_size);
    zlib_unzip(data, data_size, decomp, decomp_size);
    assert(decomp_size == info.obj_len);
    std::free(data);

    return decomp;
}

bool Miniroot::ReadByIname()
{
    ResetStringBuffer();
    m_fs.read((char*) &m_iname, sizeof(int8_t));
    if (m_iname <= 0) return false;
    m_fs.read(m_sbuf, sizeof(int8_t) * m_iname);
    return true;
}

void Miniroot::ResetStringBuffer()
{
    std::memset(m_sbuf, 0, sizeof(char) * STR_BUF_SIZE);
}
