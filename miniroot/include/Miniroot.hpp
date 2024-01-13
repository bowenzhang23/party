#pragma once

#include "FileHeader.hpp"
#include "LogicalRecordHeader.hpp"
#include "Unzip.hpp"
#include "TypeUtils.hpp"

#include <fstream>
#include <map>
#include <string>
#include <vector>

/**
 * @brief A Class that read simple root files
 *
 * Read root file storing TTree (single but can be easily extended to multiple),
 * assuming the types of all branches are trivial and known.
 *
 */
class Miniroot
{
public:
    using BasketMap = std::map<std::string, std::vector<ZipInfo_t>>;
    using StrVec    = std::vector<std::string>;

    /**
     * @brief Construct a new Miniroot object
     *
     * @param filename Path to root file
     */
    explicit Miniroot(const std::string& filename);

    /**
     * @brief Get the list of branch names
     */
    inline const StrVec& GetBranches() const { return m_branches; }

    /**
     * @brief Get the data as a continuous storage (stl vector) of bytes,
     * given the name of a specific branch
     *
     * @param branch_name Name of branch
     * @return std::vector<uint8_t> Bytes
     */
    std::vector<uint8_t> GetBytes(const std::string& branch_name);

    /**
     * @brief Get the data as a continuous storage (stl vector) of given data
     * type
     *
     * @tparam T Data type
     * @param branch_name Name of branch
     * @return std::vector<T> An array of type T
     */
    template <typename T>
    std::vector<T> Get(const std::string& branch_name);

private:
    /**
     * @brief Main method to read the file
     */
    void Read();

    /**
     * @brief Get the uncompressed data
     *
     * @param info Information about how to retrive the bytes. See ZipInfo_t
     * @return uint8_t*
     */
    uint8_t* GetUncompressedBytes(const ZipInfo_t& info);

    /**
     * @brief Get the basket map
     */
    inline const BasketMap& GetData() const { return m_basket_map; }

    /**
     * @brief Read block of memory by current iname (m_iname)
     *
     * @return true Where iname > 0, can read
     * @return false Where iname <= 0, cannot read
     */
    bool ReadByIname();

    /**
     * @brief Reset memory of m_sbuf
     */
    void ResetStringBuffer();

    /**
     * @brief Validate if all branches has the same number of baskets
     */
    void ValidateBasketMap();

    /**
     * @brief Set the list of branches (alphabet sorted)
     */
    void SetBranches();

private:
    static constexpr std::size_t      STR_BUF_SIZE     = 100;
    static constexpr int              UNZIP_HEADER_LEN = 9;
    static constexpr std::string_view TBASKET          = "TBasket";
    int8_t                            m_iname;
    char                              m_sbuf[STR_BUF_SIZE];
    FileHeader_t                      m_fh;
    LogicalRecordHeader_t             m_lrh;
    std::ifstream                     m_fs;
    BasketMap                         m_basket_map;
    StrVec                            m_branches;
};

template <typename T>
inline std::vector<T> Miniroot::Get(const std::string& branch_name)
{
    return from_bytevec<T>(GetBytes(branch_name));
}
