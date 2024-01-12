#include <arpa/inet.h> // for ntohl()
#include <array>
#include <cinttypes>
#include <fstream>
#include <iomanip>
#include <iostream>

struct logicalrecordheader_t {
    int32_t nbytes;       // Length of compressed object (in bytes)
    int16_t version;      // TKey version identifier
    int32_t obj_len;      // Length of uncompressed object
    int32_t datime;       // Date and time when object was written to file
    int16_t key_len;      // Length of the key structure (in bytes)
    int16_t cycle;        // Cycle of key
    int32_t seek_key;     // Pointer to record itself (consistency check)
    int32_t seek_key_64;  // Pointer to record itself (consistency check)
    int32_t seek_pdir;    // Pointer to directory header
    int32_t seek_pdir_64; // Pointer to directory header

    void read(std::ifstream& fs);
};

std::ostream& operator<<(std::ostream& os, const logicalrecordheader_t& lrh);
