#include "ttreemetadata.hpp"

void ttreemetadata_t::read(std::stringstream& fs)
{
    fs.read((char*) &entries, sizeof(int64_t));
    fs.read((char*) &tot_bytes, sizeof(int64_t));
    fs.read((char*) &zip_bytes, sizeof(int64_t));
    fs.read((char*) &saved_bytes, sizeof(int64_t));
    fs.read((char*) &flushed_bytes, sizeof(int64_t));
    fs.read((char*) &weight, sizeof(int64_t));
    fs.read((char*) &timer_interval, sizeof(int32_t));
    fs.read((char*) &scan_field, sizeof(int32_t));
    fs.read((char*) &update, sizeof(int32_t));
    fs.read((char*) &default_entry_offset_len, sizeof(int32_t));
    fs.read((char*) &n_cluster_range, sizeof(int32_t));
    fs.read((char*) &max_entries, sizeof(int64_t));
    fs.read((char*) &max_entry_loop, sizeof(int64_t));
    fs.read((char*) &max_virtual_size, sizeof(int64_t));
    fs.read((char*) &auto_save, sizeof(int64_t));
    fs.read((char*) &auto_flush, sizeof(int64_t));
    fs.read((char*) &estimate, sizeof(int64_t));

    entries                  = be64toh(entries);
    tot_bytes                = be64toh(tot_bytes);
    zip_bytes                = be64toh(zip_bytes);
    saved_bytes              = be64toh(saved_bytes);
    flushed_bytes            = be64toh(flushed_bytes);
    weight                   = be64toh(weight);
    timer_interval           = ntohl(timer_interval);
    scan_field               = ntohl(scan_field);
    update                   = ntohl(update);
    default_entry_offset_len = ntohl(default_entry_offset_len);
    n_cluster_range          = ntohl(n_cluster_range);
    max_entries              = be64toh(max_entries);
    max_entry_loop           = be64toh(max_entry_loop);
    max_virtual_size         = be64toh(max_virtual_size);
    auto_save                = be64toh(auto_save);
    auto_flush               = be64toh(auto_flush);
    estimate                 = be64toh(estimate);
}

std::ostream& operator<<(std::ostream& os, const ttreemetadata_t& tmd)
{
    os << " --- TTree Metadata --- \n"
       << std::hex << "entries = " << tmd.entries << '\n'
       << "tot_bytes = " << tmd.tot_bytes << '\n'
       << "zip_bytes = " << tmd.zip_bytes << '\n'
       << "saved_bytes = " << tmd.saved_bytes << '\n'
       << "flushed_bytes = " << tmd.flushed_bytes << '\n'
       << "weight = " << tmd.weight << '\n'
       << "timer_interval = " << tmd.timer_interval << '\n'
       << "scan_field = " << tmd.scan_field << '\n'
       << "update = " << tmd.update << '\n'
       << "default_entry_offset_len = " << tmd.default_entry_offset_len << '\n'
       << "n_cluster_range = " << tmd.n_cluster_range << '\n'
       << "max_entries = " << tmd.max_entries << '\n'
       << "max_entry_loop = " << tmd.max_entry_loop << '\n'
       << "max_virtual_size = " << tmd.max_virtual_size << '\n'
       << "auto_save = " << tmd.auto_save << '\n'
       << "auto_flush = " << tmd.auto_flush << '\n'
       << "estimate = " << tmd.estimate << '\n';

    os << std::dec;
    return os;
}