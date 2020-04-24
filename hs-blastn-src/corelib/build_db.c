#include "build_db.h"

#include <ctype.h>
#include <string.h>
#include <sys/stat.h>   
#include <unistd.h>

#include "cstr_util.h"
#include "db_format.h"
#include "ksort.h"
#include "line_reader.h"
#include "seqdb.h"
#include "fasta.h"

static void
pack_one_seq(int* id_in_seqdb,
    size_t* seq_offset_in_seqdb,
    size_t* hdr_offset_in_seqdb,
    kstring_t* seq_name,
    kstring_t* seq_data,
    FILE* hdr_file,
    FILE* seq_info_file,
    FILE* packed_seq_file,
    size_t* ambig_offset_in_seqdb,
    FILE* ambig_subseq_file)
{
    CSeqInfo seq_info;
    seq_info.seq_offset = *seq_offset_in_seqdb;
    seq_info.seq_size = ks_size(*seq_data);
    seq_info.hdr_offset = *hdr_offset_in_seqdb;
    seq_info.ambig_offset = *ambig_offset_in_seqdb;
    seq_info.ambig_size = 0;

    /// dump header
    if (id_in_seqdb) {  // rename sequence
        char buf[64];
        u64_to_fixed_width_string_r(*id_in_seqdb, buf, HBN_DIGIT_WIDTH);
        ++(*id_in_seqdb);
        buf[HBN_DIGIT_WIDTH] = '\0';
        hbn_fwrite(buf, 1, HBN_DIGIT_WIDTH+1, hdr_file);
        seq_info.hdr_size = HBN_DIGIT_WIDTH;
        *hdr_offset_in_seqdb += (HBN_DIGIT_WIDTH + 1);
    } else { // use the orignal header
        size_t n = 0;
        while (n < ks_size(*seq_name)) {
            int c = ks_A(*seq_name, n);
            if (isspace(c)) break;
            ++n;
        }
        hbn_assert(n > 0);
        seq_info.hdr_size = n;
        if (n == ks_size(*seq_name)) {
            kputc('\0', seq_name);
            hbn_fwrite(ks_s(*seq_name), 1, n+1, hdr_file);
            ks_pop_back(*seq_name);
        } else {
            int c = ks_A(*seq_name, n);
            ks_A(*seq_name, n) = '\0';
            hbn_fwrite(ks_s(*seq_name), 1, n+1, hdr_file);
            ks_A(*seq_name, n) = c;
        }
        *hdr_offset_in_seqdb += (n + 1);
    }

    /// dump packed sequence
    const size_t ubytes = (seq_info.seq_size + 3) >> 2;
    u8* es = (u8*)calloc(ubytes, 1);
    size_t i = 0;
    int c, c1;
    u8 ec, ec1;
    while (i < seq_info.seq_size) {
        c = ks_A(*seq_data, i);
        ec = nst_nt16_table[c];
        if (ec > 3) {
            CAmbigSubseq ambig = { i, c, 1 };
            ec1 = 0;
            _set_pac(es, i, ec1);
            ++i;
            while (i < seq_info.seq_size) {
                c1 = ks_A(*seq_data, i);
                ec1 = nst_nt16_table[c1];
                if (ec != ec1) break;
                ec1 = 0;
                _set_pac(es, i, ec1);
                ++ambig.count;
                ++i;
            }
            ++seq_info.ambig_size;
            hbn_fwrite(&ambig, sizeof(CAmbigSubseq), 1, ambig_subseq_file);
            ++(*ambig_offset_in_seqdb);
        } else {
            _set_pac(es, i, ec);
            ++i;
        }
    }
    hbn_fwrite(es, 1, ubytes, packed_seq_file);

    /// dump sequence info
    hbn_fwrite(&seq_info, sizeof(CSeqInfo), 1, seq_info_file);

    *seq_offset_in_seqdb += (ubytes << 2);
    free(es);
}

static void
pack_one_file(const char* file_path,
    int* id_in_seqdb,
    size_t* seq_offset_in_seqdb,
    size_t* hdr_offset_in_seqdb,
    FILE* hdr_file,
    FILE* seq_info_file,
    FILE* packed_seq_file,
    size_t* ambig_offset_in_seqdb,
    FILE* ambig_subseq_file,
    const int min_seq_size,
    size_t* total_seq_count,
    size_t* total_res_count)
{
    HBN_LOG("pack %s", file_path);
    size_t file_seq_count = 0;
    size_t file_res_count = 0;
    HbnFastaReader* reader = HbnFastaReaderNew(file_path);
    HbnFastaReaderSkipErrorFormatedSequences(reader);
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        if (!HbnFastaReaderReadOneSeq(reader)) continue;
        if (ks_size(reader->sequence) < min_seq_size) continue;
        pack_one_seq(id_in_seqdb,
            seq_offset_in_seqdb,
            hdr_offset_in_seqdb,
            &reader->name,
            &reader->sequence,
            hdr_file,
            seq_info_file,
            packed_seq_file,
            ambig_offset_in_seqdb,
            ambig_subseq_file);
        ++file_seq_count;
        file_res_count += ks_size(reader->sequence);
    }
    HbnFastaReaderFree(reader);

    char buf1[64], buf2[64];
    HBN_LOG("pack %s sequences (%s)", u64_to_string_comma(file_seq_count, buf1), u64_to_string_datasize(file_res_count, buf2));
    *total_seq_count += file_seq_count;
    *total_res_count += file_res_count;
}

void
build_volume_info(const char* seqdb_dir, 
    const char* seqdb_title, 
    CSeqDBInfo dbinfo,
    const size_t volume_size,
    const int volume_seq_cnt)
{
    char path[HBN_MAX_PATH_LEN];
    make_ascii_volume_info_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(ascii_out, path, "w");
    make_bin_volume_info_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(bin_out, path, "wb");

    dump_seqdb_info(fprintf, ascii_out, dbinfo);
    hbn_fwrite(&dbinfo, sizeof(CSeqDBInfo), 1, bin_out);

    size_t i = 0;
    size_t last_i = 0;
    size_t s = 0;
    int curr_seq_cnt = 0;
    const int seq_count = dbinfo.num_seqs;
    CSeqInfo* seq_info_list = load_seq_infos(seqdb_dir, seqdb_title, 0, seq_count);

    for (; i < seq_count; ++i) {
        s += seq_info_list[i].seq_size;
        ++curr_seq_cnt;
        if (s >= volume_size || curr_seq_cnt >= volume_seq_cnt) {
            dbinfo.db_size = s;
            dbinfo.num_seqs = i + 1 - last_i;
            dbinfo.seq_start_id = last_i;
            dbinfo.seq_offset_from = seq_info_list[last_i].seq_offset;
            dbinfo.seq_offset_to = seq_info_list[i].seq_offset + seq_info_list[i].seq_size;
            dbinfo.hdr_offset_from = seq_info_list[last_i].hdr_offset;
            dbinfo.hdr_offset_to = seq_info_list[i].hdr_offset + seq_info_list[i].hdr_size + 1;
            dbinfo.ambig_offset_from = seq_info_list[last_i].ambig_offset;
            dbinfo.ambig_offset_to = seq_info_list[i].ambig_offset + seq_info_list[i].ambig_size;
            dump_seqdb_info(fprintf, ascii_out, dbinfo);
            hbn_fwrite(&dbinfo, sizeof(CSeqDBInfo), 1, bin_out);
            last_i = i + 1;
            s = 0;
            curr_seq_cnt = 0;
        }
    }
    
    if (s > 0) {
        --i;
        dbinfo.db_size = s;
        dbinfo.num_seqs = i + 1 - last_i;
        dbinfo.seq_start_id = last_i;
        dbinfo.seq_offset_from = seq_info_list[last_i].seq_offset;
        dbinfo.seq_offset_to = seq_info_list[i].seq_offset + seq_info_list[i].seq_size;
        dbinfo.hdr_offset_from = seq_info_list[last_i].hdr_offset;
        dbinfo.hdr_offset_to = seq_info_list[i].hdr_offset + seq_info_list[i].hdr_size + 1;
        dbinfo.ambig_offset_from = seq_info_list[last_i].ambig_offset;
        dbinfo.ambig_offset_to = seq_info_list[i].ambig_offset + seq_info_list[i].ambig_size;
        dump_seqdb_info(fprintf, ascii_out, dbinfo);
        hbn_fwrite(&dbinfo, sizeof(CSeqDBInfo), 1, bin_out);
    }

    free(seq_info_list);
    hbn_fclose(ascii_out);
    hbn_fclose(bin_out);
}

int
hbndb_is_built(const char* input,
    const char* seqdb_dir,
    const char* seqdb_title,
    const int min_seq_size,
    const int max_file_seqs,
    const size_t max_file_res,
    const int rename_seq)
{
    char path[HBN_MAX_PATH_LEN];
    path[0] = '\0';
    if (seqdb_dir) sprintf(path, "%s/", seqdb_dir);
    if (seqdb_title) strcat(path, seqdb_title);
    if (path[0] != '\0') strcat(path, ".hbndb_is_built");
    else sprintf(path, "hbndb_is_built");
    if (access(path, F_OK) != 0) return 0;

    hbn_dfopen(old_in, path, "rb");
    int old_input_len;
    hbn_fread(&old_input_len, sizeof(int), 1, old_in);
    char* old_input = (char*)malloc(old_input_len+1);
    hbn_fread(old_input, 1, old_input_len, old_in);
    old_input[old_input_len] = '\0';
    int old_min_seq_size;
    hbn_fread(&old_min_seq_size, sizeof(int), 1, old_in);
    int old_max_file_seqs;
    hbn_fread(&old_max_file_seqs, sizeof(int), 1, old_in);
    size_t old_max_file_res;
    hbn_fread(&old_max_file_res, sizeof(size_t), 1, old_in);
    int old_rename_seq;
    hbn_fread(&old_rename_seq, sizeof(int), 1, old_in);
    size_t old_input_bytes;
    hbn_fread(&old_input_bytes, sizeof(size_t), 1, old_in);
    hbn_fclose(old_in);
    size_t new_input_bytes = hbn_file_size(input);
    int r = 1;
    if (strcmp(input, old_input)) {
        HBN_LOG("rebuild hbndb");
        r = 0;
    } else if (min_seq_size != old_min_seq_size) {
        HBN_LOG("rebuild hbndb");
        r = 0;
    } else if (max_file_seqs != old_max_file_seqs) {
        HBN_LOG("rebuild hbndb");
        r = 0;
    } else if (max_file_res != old_max_file_res) {
        HBN_LOG("rebuild hbndb");
        r = 0;
    } else if (rename_seq != old_rename_seq) {
        HBN_LOG("rebuild hbndb");
        r = 0;
    } else if (new_input_bytes != old_input_bytes) {
        HBN_LOG("rebuild hbndb, old bytes: %d, new bytes: %d, input: %s", old_input_bytes, new_input_bytes, input);
        r = 0;
    }

    free(old_input);
    return r;
}

void
hbndb_make_built(const char* input,
    const char* seqdb_dir,
    const char* seqdb_title,
    const int min_seq_size,
    const int max_file_seqs,
    const size_t max_file_res,
    const int rename_seq)
{
    char path[HBN_MAX_PATH_LEN];
    path[0] = '\0';
    if (seqdb_dir) sprintf(path, "%s/", seqdb_dir);
    if (seqdb_title) strcat(path, seqdb_title);
    if (path[0] != '\0') strcat(path, ".hbndb_is_built");
    else sprintf(path, "hbndb_is_built");
    hbn_dfopen(out, path, "wb");

    int input_len = strlen(input);
    hbn_fwrite(&input_len, sizeof(int), 1, out);
    hbn_fwrite(input, 1, input_len, out);
    hbn_fwrite(&min_seq_size, sizeof(int), 1, out);
    hbn_fwrite(&max_file_seqs, sizeof(int), 1, out);
    hbn_fwrite(&max_file_res, sizeof(size_t), 1, out);
    hbn_fwrite(&rename_seq, sizeof(int), 1, out);
    size_t input_file_size = hbn_file_size(input);
    hbn_fwrite(&input_file_size, sizeof(size_t), 1, out);

    hbn_fclose(out);
}

void
build_db(const char* input,
    const char* seqdb_dir,
    const char* seqdb_title,
    const int min_seq_size,
    const int max_file_seqs,
    const size_t max_file_res,
    const int rename_seq)
{
    EDbFormat fmt = hbn_guess_db_format(input);
    if (fmt == eDbFormatEmptyFile) {
        HBN_LOG("file %s is empty", input);
        return;
    }

    char path[HBN_MAX_PATH_LEN];
    make_seq_info_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(seq_info_file, path, "wb");
    make_packed_seq_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(packed_seq_file, path, "wb");
    make_header_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(hdr_file, path, "wb");
    make_ambig_subseq_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(ambig_subseq_file, path, "wb");
    int _read_id = 0;
    int* read_id = rename_seq ? (&_read_id) : NULL;
    size_t seq_offset_in_seqdb = 0;
    size_t hdr_offset_in_seqdb = 0;
    size_t ambig_offset_in_seqdb = 0;
    size_t seq_count = 0;
    size_t res_count = 0;

    if (fmt == eDbFormatUnknown) {
        HbnLineReader* line_reader = HbnLineReaderNew(input);
        while (!HbnLineReaderAtEof(line_reader)) {
            HbnLineReaderReadOneLine(line_reader);
            kstring_t* line = &line_reader->line;
            if (ks_empty(*line)) continue;
            kputc('\0', line);
            if (truncate_both_end_spaces(ks_s(*line)) == 0) continue;
            pack_one_file(ks_s(*line), 
                read_id, 
                &seq_offset_in_seqdb, 
                &hdr_offset_in_seqdb, 
                hdr_file,
                seq_info_file,
                packed_seq_file,
                &ambig_offset_in_seqdb,
                ambig_subseq_file,
                min_seq_size,
                &seq_count,
                &res_count);
        }
        HbnLineReaderFree(line_reader);
    } else {
        pack_one_file(input, 
            read_id, 
            &seq_offset_in_seqdb, 
            &hdr_offset_in_seqdb, 
            hdr_file,
            seq_info_file,
            packed_seq_file,
            &ambig_offset_in_seqdb,
            ambig_subseq_file,
            min_seq_size,
            &seq_count,
            &res_count);
    }    

    hbn_fclose(seq_info_file);
    hbn_fclose(packed_seq_file);
    hbn_fclose(hdr_file);
    hbn_fclose(ambig_subseq_file);

    CSeqDBInfo volinfo;
    volinfo.seq_start_id = 0;
    volinfo.num_seqs = seq_count;
    volinfo.db_size = res_count;
    volinfo.seq_offset_from = 0;
    volinfo.seq_offset_to = seq_offset_in_seqdb;
    volinfo.hdr_offset_from = 0;
    volinfo.hdr_offset_to = hdr_offset_in_seqdb;
    volinfo.ambig_offset_from = 0;
    volinfo.ambig_offset_to = ambig_offset_in_seqdb;
    build_volume_info(seqdb_dir, seqdb_title, volinfo, max_file_res, max_file_seqs);
}