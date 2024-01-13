#include "miniroot.hpp"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <memory>

Miniroot::Miniroot(const std::string& filename)
    : m_iname(0), m_fs(filename, std::ios::binary | std::ios::in)
{
    std::cout << "Read " << filename << '\n';
    ResetStringBuffer();
    Read();
    ValidateBasketMap();
    SetBranches();
}

std::vector<uint8_t> Miniroot::GetBytes(const std::string& branch_name)
{
    const auto&          info_vec = m_basket_map.at(branch_name);
    std::vector<uint8_t> decomp;
    for (const auto& info : info_vec) {
        auto di = std::unique_ptr<uint8_t>(GetUncompressedBytes(info));
        decomp.insert(decomp.end(), di.get(), di.get() + info.obj_len);
    }
    return decomp;
}

void Miniroot::Read()
{
    m_fh.read(m_fs);
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
        std::string class_name(m_sbuf);
        // object name
        if (!ReadByIname()) continue;
        std::string object_name(m_sbuf);
        // object title
        if (!ReadByIname()) continue;
        std::string object_title(m_sbuf);
        if (class_name == TBASKET) {
            std::string path(object_title);
            path.push_back('/');
            path.append(object_name);

            m_basket_map[path].emplace_back(
                m_lrh.nbytes,  // nbytes including key
                m_lrh.key_len, // nbytes key only
                m_lrh.obj_len, // uncompressed object length
                cur + m_lrh.key_len + UNZIP_HEADER_LEN,      // seekpos
                m_lrh.obj_len > m_lrh.nbytes - m_lrh.key_len // compressed
            );
        }
    }
}

uint8_t* Miniroot::GetUncompressedBytes(const zip_info& info)
{
    m_fs.seekg(info.seekpos);
    int32_t  header_len  = info.compressed ? UNZIP_HEADER_LEN : 0;
    int32_t  data_size   = info.nbytes - info.key_len - header_len;
    int32_t  decomp_size = info.obj_len;
    uint8_t* data        = new uint8_t[data_size];
    m_fs.read((char*) data, sizeof(uint8_t) * data_size);
    if (!info.compressed) return data;

    uint8_t* decomp = new uint8_t[decomp_size];
    zlib_unzip(data, data_size, decomp, decomp_size);
    assert(decomp_size == info.obj_len);

    delete[] data;
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

void Miniroot::ValidateBasketMap()
{
    int nbaskets = 0;
    for (const auto& p : m_basket_map) {
        if (nbaskets == 0)
            nbaskets = p.second.size();
        else
            assert(nbaskets == p.second.size());
    }
}

void Miniroot::SetBranches()
{
    m_branches.clear();
    for (const auto& p : m_basket_map) { m_branches.push_back(p.first); }
    std::sort(m_branches.begin(), m_branches.end());
}
