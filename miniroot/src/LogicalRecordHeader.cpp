#include "LogicalRecordHeader.hpp"

#include <arpa/inet.h> // for ntohl()
#include <iostream>

void LogicalRecordHeader_t::read(std::ifstream& fs)
{
    fs.read((char*) &nbytes, sizeof(int32_t));
    fs.read((char*) &version, sizeof(int16_t));
    fs.read((char*) &obj_len, sizeof(int32_t));
    fs.read((char*) &datime, sizeof(int32_t));
    fs.read((char*) &key_len, sizeof(int16_t));
    fs.read((char*) &cycle, sizeof(int16_t));

    nbytes  = std::abs((int) ntohl(nbytes)); // negative value might exist
    version = ntohs(version);
    obj_len = ntohl(obj_len);
    datime  = ntohl(datime);
    key_len = ntohs(key_len);
    cycle   = ntohs(cycle);

    // see https://root.cern.ch/doc/master/classTKey.html
    fs.read((char*) &seek_key, sizeof(int32_t));
    if (version > 1000)
        fs.read((char*) &seek_key_64, sizeof(int32_t));
    else
        seek_key_64 = 0;
    fs.read((char*) &seek_pdir, sizeof(int32_t));
    if (version > 1000)
        fs.read((char*) &seek_pdir_64, sizeof(int32_t));
    else
        seek_pdir_64 = 0;

    seek_key     = ntohl(seek_key);
    seek_key_64  = ntohl(seek_key_64);
    seek_pdir    = ntohl(seek_pdir);
    seek_pdir_64 = ntohl(seek_pdir_64);
}

std::ostream& operator<<(std::ostream& os, const LogicalRecordHeader_t& lrh)
{
    os << " --- Logical Record Header --- \n"
       << "nbytes = " << lrh.nbytes << '\n'
       << "version = " << lrh.version << '\n'
       << "obj_len = " << lrh.obj_len << '\n'
       << "datime = " << lrh.datime << '\n'
       << "key_len = " << lrh.key_len << '\n'
       << "cycle = " << lrh.cycle << '\n'
       << "seek_key = " << lrh.seek_key << '\n'
       << "seek_key_64 = " << lrh.seek_key_64 << '\n'
       << "seek_pdir = " << lrh.seek_pdir << '\n'
       << "seek_pdir_64 = " << lrh.seek_pdir_64 << '\n';
    os << std::dec;
    return os;
}
