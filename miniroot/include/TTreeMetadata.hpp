#pragma once

#include <cinttypes>
#include <sstream>

/**
 * @brief TTree Metadata of TFile
 *
 * No direct reference shows the exact layout, 
 * This is not used at the moment,
 * However the information can still be extracted.
 */
struct TTreeMetadata_t {
    int64_t entries;
    int64_t tot_bytes;
    int64_t zip_bytes;
    int64_t saved_bytes;
    int64_t flushed_bytes;
    int64_t weight; // double
    int32_t timer_interval;
    int32_t scan_field;
    int32_t update;
    int32_t default_entry_offset_len;
    int32_t n_cluster_range;
    int64_t max_entries;
    int64_t max_entry_loop;
    int64_t max_virtual_size;
    int64_t auto_save;
    int64_t auto_flush;
    int64_t estimate;

    /**
     * @brief Read field values from string stream
     *
     * @param ss string stream
     */
    void read(std::stringstream& ss);
};

std::ostream& operator<<(std::ostream& os, const TTreeMetadata_t& lrh);
