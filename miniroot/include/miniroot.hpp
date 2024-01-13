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
    using StrVec    = std::vector<std::string>;
    explicit Miniroot(const std::string& filename);
    inline const StrVec& GetBranches() const { return m_branches; }
    std::vector<uint8_t> GetBytes(const std::string& branch_name);
    template <typename T>
    std::vector<T> Get(const std::string& branch_name);

private:
    void                    Read();
    uint8_t*                GetUncompressedBytes(const zip_info& info);
    inline const BasketMap& GetData() const { return m_basket_map; }
    bool                    ReadByIname();
    void                    ResetStringBuffer();
    void                    ValidateBasketMap();
    void                    SetBranches();

private:
    static constexpr std::size_t      STR_BUF_SIZE     = 100;
    static constexpr int              UNZIP_HEADER_LEN = 9;
    static constexpr std::string_view TBASKET          = "TBasket";
    int8_t                            m_iname;
    char                              m_sbuf[STR_BUF_SIZE];
    fileheader_t                      m_fh;
    logicalrecordheader_t             m_lrh;
    std::ifstream                     m_fs;
    BasketMap                         m_basket_map;
    StrVec                            m_branches;
};

template <typename T>
inline std::vector<T> Miniroot::Get(const std::string& branch_name)
{
    return from_bytevec<T>(GetBytes(branch_name));
}
