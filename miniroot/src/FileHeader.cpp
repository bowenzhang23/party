#include "FileHeader.hpp"

#include <arpa/inet.h> // for ntohl()
#include <iostream>

void FileHeader_t::read(std::ifstream& fs)
{
    fs.read((char*) &identifier, sizeof(int32_t));
    fs.read((char*) &f_version, sizeof(int32_t));
    fs.read((char*) &f_begin, sizeof(int32_t));
    fs.read((char*) &f_end, sizeof(int32_t));
    fs.read((char*) &f_seekfree, sizeof(int32_t));
    fs.read((char*) &f_nbytesfree, sizeof(int32_t));
    fs.read((char*) &nfree, sizeof(int32_t));
    fs.read((char*) &f_nbytesname, sizeof(int32_t));
    fs.read((char*) &f_units, sizeof(int8_t));
    fs.read((char*) &f_compress, sizeof(int32_t));
    fs.read((char*) &f_seekinfo, sizeof(int32_t));
    fs.read((char*) &f_nbytesinfo, sizeof(int32_t));
    fs.read((char*) &f_uuid, sizeof(int8_t) * 18);

    identifier   = ntohl(identifier);
    f_version    = ntohl(f_version);
    f_begin      = ntohl(f_begin);
    f_end        = ntohl(f_end);
    f_seekfree   = ntohl(f_seekfree);
    f_nbytesfree = ntohl(f_nbytesfree);
    nfree        = ntohl(nfree);
    f_nbytesname = ntohl(f_nbytesname);
    f_compress   = ntohl(f_compress);
    f_seekinfo   = ntohl(f_seekinfo);
    f_nbytesinfo = ntohl(f_nbytesinfo);
}

std::ostream& operator<<(std::ostream& os, const FileHeader_t& fh)
{
    os << " --- File Header --- \n"
       << "identifier = " << fh.identifier << '\n'
       << "f_version = " << fh.f_version << '\n'
       << "f_begin = " << fh.f_begin << '\n'
       << "f_end = " << fh.f_end << '\n'
       << "f_seekfree = " << fh.f_seekfree << '\n'
       << "f_nbytesfree = " << fh.f_nbytesfree << '\n'
       << "nfree = " << fh.nfree << '\n'
       << "f_nbytesname = " << fh.f_nbytesname << '\n'
       << "f_units = " << (int) fh.f_units << '\n'
       << "f_compress = " << fh.f_compress << '\n'
       << "f_seekinfo = " << fh.f_seekinfo << '\n'
       << "f_nbytesinfo = " << fh.f_nbytesinfo << '\n'
       << "f_uuid = ";

    for (char i = 0; i < 17; ++i) { os << (int) fh.f_uuid[i] << ':'; }
    os << (int) fh.f_uuid[18] << '\n';
    os << std::dec;
    return os;
}