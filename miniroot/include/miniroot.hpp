#include "fileheader.hpp"
#include "logicalrecordheader.hpp"
#include "unzip.hpp"

#include <fstream>
#include <map>
#include <string>
#include <vector>

class Miniroot
{
public:
    using BasketMap = std::map<std::string, std::vector<zip_info>>;
    explicit Miniroot(const std::string& filename);
    void                    Read();
    uint8_t*                GetUncompressedBytes(const zip_info& info);
    inline const BasketMap& GetData() const { return m_basket_map; }

private:
    bool ReadByIname();
    void ResetStringBuffer();

private:
    static constexpr std::size_t STR_BUF_SIZE     = 100;
    static constexpr int         UNZIP_HEADER_LEN = 9;
    int8_t                       m_iname;
    char                         m_sbuf[STR_BUF_SIZE];
    fileheader_t                 m_fh;
    logicalrecordheader_t        m_lrh;
    std::ifstream                m_fs;
    BasketMap                    m_basket_map;
};