#include "seqdb.h"

#include "../ncbi_blast/setup/blast_encoding.h"

void
make_ambig_subseq_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "ambig");
}

CAmbigSubseq*
load_ambig_subseqs(const char* data_dir, const char* db_name, const size_t from, const size_t to)
{
    hbn_assert(from >= 0);
    hbn_assert(from <= to);
    const size_t bytes_from = sizeof(CAmbigSubseq) * from;
    const size_t bytes_to = sizeof(CAmbigSubseq) * to;
    char path[HBN_MAX_PATH_LEN];
    make_ambig_subseq_path(data_dir, db_name, path);
    const size_t total_bytes = hbn_file_size(path);
    hbn_assert(bytes_to <= total_bytes);
    size_t n = to - from;
    CAmbigSubseq* ambig_subseq_list = (CAmbigSubseq*)malloc(bytes_to - bytes_from);
    hbn_dfopen(in, path, "rb");
    fseek(in, bytes_from, SEEK_SET);
    hbn_fread(ambig_subseq_list, sizeof(CAmbigSubseq), n, in);
    hbn_fclose(in);
    return ambig_subseq_list;
}

void
make_header_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "hdr");
}

char* 
load_seq_headers(const char* data_dir, const char* db_name, const size_t from, const size_t to)
{
    hbn_assert(from >= 0);
    hbn_assert(from <= to);
    char path[HBN_MAX_PATH_LEN];
    make_header_path(data_dir, db_name, path);
    const size_t total_bytes = hbn_file_size(path);
    hbn_assert(to <= total_bytes, "to = %zu, total_bytes = %zu", to, total_bytes);
    const size_t n = to - from;
    char* hdr = (char*)malloc(n);
    hbn_dfopen(in, path, "rb");
    fseek(in, from, SEEK_SET);
    hbn_fread(hdr, 1, n, in);
    hbn_fclose(in);
    return hdr;
}

void
make_seq_info_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "seq_info");
}

void
make_seqdb_summary_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "summary");
}

CSeqInfo*
load_seq_infos(const char* data_dir, const char* db_name, const size_t from, const size_t to)
{
    hbn_assert(from >= 0);
    hbn_assert(from <= to);
    const size_t bytes_from = sizeof(CSeqInfo) * from;
    const size_t bytes_to = sizeof(CSeqInfo) * to;
    char path[HBN_MAX_PATH_LEN];
    make_seq_info_path(data_dir, db_name, path);
    const size_t total_bytes = hbn_file_size(path);
    hbn_assert(bytes_to <= total_bytes);
    size_t n = to - from;
    CSeqInfo* seq_info_list = (CSeqInfo*)malloc(sizeof(CSeqInfo) * n);
    hbn_dfopen(in, path, "rb");
    fseek(in, bytes_from, SEEK_SET);
    hbn_fread(seq_info_list, sizeof(CSeqInfo), n, in);
    hbn_fclose(in);
    return seq_info_list;
}

u8*
seqdb_load_pac(const char* seqdb_dir, const char* seqdb_title, const size_t res_from, size_t res_to)
{
    hbn_assert((res_from&3) == 0);
    hbn_assert(res_from <= res_to);
    const size_t bytes_from = res_from >> 2;
    const size_t bytes_to = (res_to + 3) >> 2;
    const size_t bytes_cnt = bytes_to - bytes_from;
    char path[HBN_MAX_PATH_LEN];
    make_packed_seq_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(in, path, "rb");
    u8* pac = (u8*)calloc(bytes_cnt, sizeof(u8));
    fseek(in, bytes_from, SEEK_SET);
    hbn_fread(pac, sizeof(u8), bytes_cnt, in);
    hbn_fclose(in);
    return pac;
}

u8*
seqdb_load_unpac(const char* seqdb_dir, const char* seqdb_title, const size_t res_from, size_t res_to)
{
    u8* pac = seqdb_load_pac(seqdb_dir, seqdb_title, res_from, res_to);
    size_t res_cnt = res_to - res_from;
    u8* unpac = (u8*)calloc(res_cnt, sizeof(u8));
    for (size_t i = 0; i < res_cnt; ++i) unpac[i] = _get_pac(pac, i);
    free(pac);
    return unpac;
}

void
make_ascii_volume_info_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "volume_info.txt");
}

void
make_bin_volume_info_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "volume_info.bin");
}

CSeqDBInfo
seqdb_load_volume_info(const char* data_dir, const char* db_name, int volume_index)
{
    char path[HBN_MAX_PATH_LEN];
    make_bin_volume_info_path(data_dir, db_name, path);
    const size_t file_size = hbn_file_size(path);
    hbn_assert(file_size % sizeof(CSeqDBInfo) == 0);
    const int nvol = file_size / sizeof(CSeqDBInfo);
    hbn_assert(volume_index >= 0, ": %d", volume_index);
    hbn_assert(volume_index < nvol);
    CSeqDBInfo* vol_info_list = (CSeqDBInfo*)malloc(file_size);
    hbn_dfopen(in, path, "rb");
    hbn_fread(vol_info_list, sizeof(CSeqDBInfo), nvol, in);
    hbn_fclose(in);
    CSeqDBInfo vol_info = vol_info_list[volume_index];
    free(vol_info_list);
    return vol_info;
}

void
make_packed_seq_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "pac");
}

size_t seqdb_seq_offset(const CSeqDB* seqdb, const int seq_id)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs, "seq_id = %d, num_seqs = %d",
			seq_id, seqdb->dbinfo.num_seqs);
    return seqdb->seq_info_list[seq_id].seq_offset;
}

size_t seqdb_seq_size(const CSeqDB* seqdb, const int seq_id)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs, 
        "seq_id = %d, num_seqs = %d", seq_id, seqdb->dbinfo.num_seqs);
    return seqdb->seq_info_list[seq_id].seq_size;    
}

const char* seqdb_seq_name(const CSeqDB* seqdb, const int seq_id)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs);
    return seqdb->seq_header_list + seqdb->seq_info_list[seq_id].hdr_offset;
}

int seqdb_num_seqs(const CSeqDB* seqdb)
{
    return seqdb->dbinfo.num_seqs;
}

size_t seqdb_size(const CSeqDB* seqdb)
{
    return seqdb->dbinfo.db_size;
}

int seqdb_offset_to_seq_id(const CSeqDB* seqdb, const size_t offset)
{
    int left = 0, mid = 0, right = seqdb->dbinfo.num_seqs;
    while (left < right) {
        mid = (left + right) >> 1;
        if (offset >= seqdb->seq_info_list[mid].seq_offset) {
            if (mid == seqdb->dbinfo.num_seqs - 1) break;
            if (offset < seqdb->seq_info_list[mid+1].seq_offset) break;
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return mid;    
}

size_t seqdb_max_offset(const CSeqDB* seqdb)
{
    return seqdb->dbinfo.seq_offset_to - seqdb->dbinfo.seq_offset_from;
}

int seqdb_load_num_volumes(const char* seqdb_dir, const char* seqdb_title)
{
    char path[HBN_MAX_PATH_LEN];
    make_bin_volume_info_path(seqdb_dir, seqdb_title, path);
    size_t file_size = hbn_file_size(path);
    hbn_assert(file_size > 0);
    hbn_assert(file_size % sizeof(CSeqDBInfo) == 0);
    const int num_vols = file_size / sizeof(CSeqDBInfo);
    hbn_assert(num_vols > 1);
    return num_vols - 1;
}

int seqdb_load_num_reads(const char* seqdb_dir, const char* seqdb_title)
{
    CSeqDBInfo dbinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, 0);
    return dbinfo.num_seqs;
}

void
seqdb_extract_subsequence(const CSeqDB* seqdb,
    const int seq_id,
    const size_t from,
    const size_t to,
    const int strand,
    vec_u8* seq)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs);
    size_t start = seqdb_seq_offset(seqdb, seq_id);
    size_t size = seqdb_seq_size(seqdb, seq_id);
    hbn_assert(from <= to && to <= size);
    kv_clear(*seq);
    if (strand == FWD) {
        size_t sfrom = start + from;
        size_t sto = start + to;
        for (size_t i = sfrom; i < sto; ++i) {
            u8 c = _get_pac(seqdb->packed_seq, i);
            kv_push(u8, *seq, c);
        }        
    } else {
        size_t sfrom = start + (size - to);
        size_t sto = start + (size - from);
        size_t i = sto;
        while (i > sfrom) {
            --i;
            u8 c = _get_pac(seqdb->packed_seq, i);
            c = 3 - c;
            kv_push(u8, *seq, c);
        }
    }
}

void
seqdb_extract_sequence(const CSeqDB* seqdb,
    const int seq_id,
    const int strand,
    vec_u8* seq)
{
    const size_t seq_size = seqdb_seq_size(seqdb, seq_id);
    seqdb_extract_subsequence(seqdb, seq_id, 0, seq_size, strand, seq);
}

static const Uint1 BLASTNA_REVERSE_COMPLEMENT_TABLE[16] = 
		{3, 2, 1, 0, 5, 4, 7, 6, 8, 9, 13, 12, 11, 10, 14, 15};

void
seqdb_recover_subsequence_ambig_res(const CSeqDB* seqdb,
    const int seq_id,
    const size_t from,
    const size_t to,
    const int strand,
    u8* seq)
{
    CSeqInfo seqinfo = seqdb->seq_info_list[seq_id];
    CAmbigSubseq* amb_array = seqdb->ambig_subseq_list + seqinfo.ambig_offset;
    if (strand == FWD) {
        for (int i = 0; i < seqinfo.ambig_size; ++i) {
            CAmbigSubseq amb = amb_array[i];
            int res = amb.ambig_residue;
            hbn_assert((res >= 'A' && res <= 'Z') || (res >= 'a' && res <= 'z'));
            res = nst_nt16_table[res];
            hbn_assert(res >= 0 && res < 16);
            size_t pos = amb.offset;
            for (int k = 0; k < amb.count; ++k, ++pos) {
                if (pos >= from && pos < to) seq[pos - from] = res;
            }
        }
    } else {
        hbn_assert(strand == REV);
        size_t size = seqdb_seq_size(seqdb, seq_id);
        hbn_assert(to <= size);
        size_t offset_delta = size - to;
        for (int i = 0; i < seqinfo.ambig_size; ++i) {
            CAmbigSubseq amb = amb_array[i];
            int res = amb.ambig_residue;
            hbn_assert((res >= 'A' && res <= 'Z') || (res >= 'a' && res <= 'z'));
            res = nst_nt16_table[res];
            hbn_assert(res >= 0 && res < 16);
            res = BLASTNA_REVERSE_COMPLEMENT_TABLE[res];
            size_t pos = amb.offset;
            for (int k = 0; k < amb.count; ++k, ++pos) {
                if (pos >= from && pos < to) {
                    size_t rev_pos = size - 1 - pos - offset_delta;
                    seq[rev_pos] = res;
                }
            }
        }    
    }
}

void
seqdb_recover_sequence_ambig_res(const CSeqDB* seqdb,
    const int seq_id,
    const int strand,
    u8* seq)
{
    const size_t seq_size = seqdb_seq_size(seqdb, seq_id);
    seqdb_recover_subsequence_ambig_res(seqdb, seq_id, 0, seq_size, strand, seq);
}

void
seqdb_extract_subsequence_with_ambig_res(const CSeqDB* seqdb,
    const int seq_id,
    const size_t from,
    const size_t to,
    const int strand,
    vec_u8* seq_v)
{
    seqdb_extract_subsequence(seqdb, seq_id, from, to, strand, seq_v);
    seqdb_recover_subsequence_ambig_res(seqdb, seq_id, from, to, strand, kv_data(*seq_v));
}

void
seqdb_extract_sequence_with_ambig_res(const CSeqDB* seqdb,
    const int seq_id,
    const int strand,
    vec_u8* seq)
{
    const size_t seq_size = seqdb_seq_size(seqdb, seq_id);
    seqdb_extract_subsequence_with_ambig_res(seqdb, seq_id, 0, seq_size, strand, seq);
}

void
seqdb_extract_raw_subsequence(const CSeqDB* seqdb,
    const int seq_id,
    const size_t from,
    const size_t to,
    const int strand,
    vec_u8* seq)
{
    seqdb_extract_subsequence(seqdb, seq_id, from, to, strand, seq);
    size_t pos = 0;
    for (size_t i = from; i < to; ++i, ++pos) {
        int c = kv_A(*seq, pos);
        hbn_assert(c >= 0 && c < 16);
        c = BLASTNA_TO_IUPACNA[c];
        kv_A(*seq, pos) = c;
    }
}

void
seqdb_extract_raw_sequence(const CSeqDB* seqdb,
    const int seq_id,
    const int strand,
    vec_u8* seq)
{
    const size_t seq_size = seqdb_seq_size(seqdb, seq_id);
    seqdb_extract_raw_subsequence(seqdb, seq_id, 0, seq_size, strand, seq);    
}

void
seqdb_extract_raw_subsequence_with_ambig_res(const CSeqDB* seqdb,
    const int seq_id,
    const size_t from,
    const size_t to,
    const int strand,
    vec_u8* seq)
{
    seqdb_extract_subsequence_with_ambig_res(seqdb, seq_id, from, to, strand, seq);
    size_t pos = 0;
    for (size_t i = from; i < to; ++i, ++pos) {
        int c = kv_A(*seq, pos);
        hbn_assert(c >= 0 && c < 16);
        c = BLASTNA_TO_IUPACNA[c];
        kv_A(*seq, pos) = c;
    }
}

void
seqdb_extract_raw_sequence_with_ambig_res(const CSeqDB* seqdb,
    const int seq_id,
    const int strand,
    vec_u8* seq)
{
    const size_t seq_size = seqdb_seq_size(seqdb, seq_id);
    seqdb_extract_raw_subsequence_with_ambig_res(seqdb, seq_id, 0, seq_size, strand, seq);    
}

CSeqDB*
seqdb_load(const char* seqdb_dir, const char* seqdb_title, int vol_id)
{
    CSeqDB* vol = (CSeqDB*)calloc(1, sizeof(CSeqDB));
    ++vol_id;
    vol->dbinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, vol_id);
    vol->packed_seq = seqdb_load_pac(seqdb_dir, seqdb_title, vol->dbinfo.seq_offset_from, vol->dbinfo.seq_offset_to);
    vol->seq_header_list = load_seq_headers(seqdb_dir, seqdb_title, vol->dbinfo.hdr_offset_from, vol->dbinfo.hdr_offset_to);
    vol->seq_info_list = load_seq_infos(seqdb_dir, seqdb_title, vol->dbinfo.seq_start_id, vol->dbinfo.seq_start_id + vol->dbinfo.num_seqs);
    vol->ambig_subseq_list = load_ambig_subseqs(seqdb_dir, seqdb_title, vol->dbinfo.ambig_offset_from, vol->dbinfo.ambig_offset_to);

    const size_t hdr_offset_from = vol->dbinfo.hdr_offset_from;
    const size_t seq_offset_from = vol->dbinfo.seq_offset_from;
    const size_t ambig_offset_from = vol->dbinfo.ambig_offset_from;
    for (int i = 0; i < vol->dbinfo.num_seqs; ++i) {
        hbn_assert(vol->seq_info_list[i].hdr_offset >= hdr_offset_from);
        vol->seq_info_list[i].hdr_offset -= hdr_offset_from;
        hbn_assert(vol->seq_info_list[i].seq_offset >= seq_offset_from);
        vol->seq_info_list[i].seq_offset -= seq_offset_from;
        vol->seq_info_list[i].ambig_offset -= ambig_offset_from;
    }

    return vol;
}

CSeqDB*
seqdb_load_unpacked(const char* seqdb_dir, const char* seqdb_title, int vol_id)
{
    CSeqDB* vol = (CSeqDB*)calloc(1, sizeof(CSeqDB));
    ++vol_id;
    vol->dbinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, vol_id);
    vol->unpacked_seq = seqdb_load_unpac(seqdb_dir, seqdb_title, vol->dbinfo.seq_offset_from, vol->dbinfo.seq_offset_to);
    vol->seq_header_list = load_seq_headers(seqdb_dir, seqdb_title, vol->dbinfo.hdr_offset_from, vol->dbinfo.hdr_offset_to);
    vol->seq_info_list = load_seq_infos(seqdb_dir, seqdb_title, vol->dbinfo.seq_start_id, vol->dbinfo.seq_start_id + vol->dbinfo.num_seqs);
    vol->ambig_subseq_list = load_ambig_subseqs(seqdb_dir, seqdb_title, vol->dbinfo.ambig_offset_from, vol->dbinfo.ambig_offset_to);

    const size_t hdr_offset_from = vol->dbinfo.hdr_offset_from;
    const size_t seq_offset_from = vol->dbinfo.seq_offset_from;
    const size_t ambig_offset = vol->dbinfo.ambig_offset_from;
    for (int i = 0; i < vol->dbinfo.num_seqs; ++i) {
        hbn_assert(vol->seq_info_list[i].hdr_offset >= hdr_offset_from);
        vol->seq_info_list[i].hdr_offset -= hdr_offset_from;
        hbn_assert(vol->seq_info_list[i].seq_offset >= seq_offset_from);
        vol->seq_info_list[i].seq_offset -= seq_offset_from;  
        vol->seq_info_list[i].ambig_offset -= ambig_offset;  
    }
    return vol;
}

CSeqDB*
seqdb_load_unpacked_with_ambig_res(const char* seqdb_dir, const char* seqdb_title, int vol_id)
{
    CSeqDB* vol = seqdb_load_unpacked(seqdb_dir, seqdb_title, vol_id);
    for (int i = 0; i < seqdb_num_seqs(vol); ++i) {
        u8* seq = vol->unpacked_seq + seqdb_seq_offset(vol, i);
        seqdb_recover_sequence_ambig_res(vol, i, FWD, seq);
    }
    return vol;
}

CSeqDB*
CSeqDBFree(CSeqDB* vol)
{
    free(vol->packed_seq);
    free(vol->unpacked_seq);
    free(vol->seq_header_list);
    free(vol->seq_info_list);
    free(vol->ambig_subseq_list);
    free(vol);
    return NULL;
}

CSeqDB*
CSeqDBNew()
{
    CSeqDB* db = (CSeqDB*)calloc(1, sizeof(CSeqDB));
    return db;
}

int
extract_sequence_block_from_packed_seqdb(const CSeqDB* pdb,
    int* next_seq_id,
    pthread_mutex_t* seq_id_lock,
    const int seq_batch_size,
    const BOOL extract_fwd_sequence,
    const BOOL extract_rev_sequence,
    const BOOL extract_ambig_sequence,
    BLAST_SequenceBlk* seq_blk,
    BlastQueryInfo* seq_info)
{
    hbn_assert(pdb->packed_seq);
    int from = 0, to = 0;
    if (seq_id_lock) pthread_mutex_lock(seq_id_lock);
    from = *next_seq_id;
    *next_seq_id += seq_batch_size;
    if (seq_id_lock) pthread_mutex_unlock(seq_id_lock);
    if (from >= pdb->dbinfo.num_seqs) return 0;
    to = hbn_min(from + seq_batch_size, pdb->dbinfo.num_seqs);
    int num_seq = to - from;

    int length = 0;
    for (int i = from; i < to; ++i) {
        int size = seqdb_seq_size(pdb, i);
        if (extract_fwd_sequence) length += size;
        if (extract_rev_sequence) length += size;
    }

    BlastContextInfo ctx_info;
    int ctx_idx = 0;
    int seq_idx = 0;
    int max_length = 0;
    int min_length = I32_MAX;
    seq_blk->sequence = (Uint1*)realloc(seq_blk->sequence, length);
    if (extract_ambig_sequence) {
        seq_blk->sequence_nomask = (Uint1*)realloc(seq_blk->sequence_nomask, length);
    }
    seq_blk->length = length;
    kv_dinit(vec_u8, seq);

    for (int i = from; i < to; ++i) {
        ctx_info.query_length = seqdb_seq_size(pdb, i);
        ctx_info.eff_searchsp = 0;
        ctx_info.length_adjustment = 0;
        ctx_info.query_index = i;
        ctx_info.frame = 0;
        ctx_info.is_valid = TRUE;
        ctx_info.segment_flags = 0;

        max_length = hbn_max(max_length, ctx_info.query_length);
        min_length = hbn_min(min_length, ctx_info.query_length);

        if (extract_fwd_sequence) {
            ctx_info.query_offset = seq_idx;
            ctx_info.frame = FWD;
            seq_info->contexts[ctx_idx++] = ctx_info;
            seqdb_extract_sequence(pdb, i, FWD, &seq);
            hbn_assert(ctx_info.query_length == kv_size(seq));
            memcpy(seq_blk->sequence + seq_idx, kv_data(seq), kv_size(seq));
            if (extract_ambig_sequence) {
                seqdb_recover_sequence_ambig_res(pdb, i, FWD, kv_data(seq));
                hbn_assert(ctx_info.query_length == kv_size(seq));
                memcpy(seq_blk->sequence_nomask + seq_idx, kv_data(seq), kv_size(seq));
            }
            seq_idx += ctx_info.query_length;
        }

        if (extract_rev_sequence) {
            ctx_info.query_offset = seq_idx;
            ctx_info.frame = REV;
            seq_info->contexts[ctx_idx++] = ctx_info;
            seqdb_extract_sequence(pdb, i, REV, &seq);
            hbn_assert(ctx_info.query_length == kv_size(seq));
            memcpy(seq_blk->sequence + seq_idx, kv_data(seq), kv_size(seq));
            if (extract_ambig_sequence) {
                seqdb_recover_sequence_ambig_res(pdb, i, REV, kv_data(seq));
                hbn_assert(ctx_info.query_length == kv_size(seq));
                memcpy(seq_blk->sequence_nomask + seq_idx, kv_data(seq), kv_size(seq));
            }
            seq_idx += ctx_info.query_length;
        }
    }
    hbn_assert(seq_idx == length);
    kv_destroy(seq);

    seq_info->first_context = 0;
    seq_info->last_context = ctx_idx - 1;
    seq_info->num_queries = num_seq;
    seq_info->max_length = max_length;
    seq_info->min_length = min_length;  

    return num_seq;
}

static void
s_extract_subsequence_without_ambig_res_from_unpacked_seqdb(
    const text_t* pdb, const int pid, 
    const size_t from, const size_t to, vec_u8* primer)
{
    hbn_assert(pdb != NULL);
    kv_clear(*primer);
    const u8* s = pdb->unpacked_seq + seqdb_seq_offset(pdb, pid);
    const u8 code_table[16] = {0,1,2,3,0,0,0,0,0,0,0,0,0,0,0,0xf};
    for (size_t i = from; i < to; ++i) {
        u8 c = s[i];
        c = code_table[c];
        kv_push(u8, *primer, c);
    }
}

int
extract_sequence_block_from_unpacked_seqdb(const CSeqDB* updb,
    int* next_seq_id,
    pthread_mutex_t* seq_id_lock,
    const int seq_batch_size,
    const BOOL extract_fwd_sequence,
    const BOOL extract_rev_sequence,
    const BOOL extract_ambig_sequence,
    BLAST_SequenceBlk* seq_blk,
    BlastQueryInfo* seq_info)
{
    hbn_assert(updb->unpacked_seq);
    int from = 0, to = 0;
    if (seq_id_lock) pthread_mutex_lock(seq_id_lock);
    from = *next_seq_id;
    *next_seq_id += seq_batch_size;
    if (seq_id_lock) pthread_mutex_unlock(seq_id_lock);
    if (from >= updb->dbinfo.num_seqs) return 0;
    to = hbn_min(from + seq_batch_size, updb->dbinfo.num_seqs);
    int num_seq = to - from;

    int length = 0;
    for (int i = from; i < to; ++i) {
        int size = seqdb_seq_size(updb, i);
        if (extract_fwd_sequence) length += size;
        if (extract_rev_sequence) length += size;
    }

    BlastContextInfo ctx_info;
    int ctx_idx = 0;
    int seq_idx = 0;
    int max_length = 0;
    int min_length = I32_MAX;
    seq_blk->sequence = (Uint1*)realloc(seq_blk->sequence, length);
    if (extract_ambig_sequence) {
        seq_blk->sequence_nomask = (Uint1*)realloc(seq_blk->sequence_nomask, length);
    }
    seq_blk->length = length;
    kv_dinit(vec_u8, seq);

    for (int i = from; i < to; ++i) {
        ctx_info.query_length = seqdb_seq_size(updb, i);
        ctx_info.eff_searchsp = 0;
        ctx_info.length_adjustment = 0;
        ctx_info.query_index = i;
        ctx_info.frame = 0;
        ctx_info.is_valid = TRUE;
        ctx_info.segment_flags = 0;

        max_length = hbn_max(max_length, ctx_info.query_length);
        min_length = hbn_min(min_length, ctx_info.query_length);

        if (extract_fwd_sequence) {
            ctx_info.query_offset = seq_idx;
            ctx_info.frame = FWD;
            seq_info->contexts[ctx_idx++] = ctx_info;
            s_extract_subsequence_without_ambig_res_from_unpacked_seqdb(
                updb, i, 0, ctx_info.query_length, &seq);
            hbn_assert(ctx_info.query_length == kv_size(seq));
            memcpy(seq_blk->sequence + seq_idx, kv_data(seq), kv_size(seq));
            if (extract_ambig_sequence) {
                const u8* ambig_seq = updb->unpacked_seq + seqdb_seq_offset(updb, i);
                memcpy(seq_blk->sequence_nomask + seq_idx, ambig_seq, ctx_info.query_length);
            }
            seq_idx += ctx_info.query_length;
        }

        if (extract_rev_sequence) {
            ctx_info.query_offset = seq_idx;
            ctx_info.frame = REV;
            seq_info->contexts[ctx_idx++] = ctx_info;
            if (!extract_fwd_sequence) {
                s_extract_subsequence_without_ambig_res_from_unpacked_seqdb(
                    updb, i, 0, ctx_info.query_length, &seq);                
            }
            size_t p = 0, q = kv_size(seq);
            while (q) {
                --q;
                u8 c = 3 - kv_A(seq, q);
                seq_blk->sequence[seq_idx + p] = c;
                ++p;
            }
            hbn_assert(p == ctx_info.query_length);
            if (extract_ambig_sequence) {
                const u8* fwd_ambig_from = updb->unpacked_seq + seqdb_seq_offset(updb, i);
                const u8* fwd_ambig_to = fwd_ambig_from + ctx_info.query_length;
                const u8* as = fwd_ambig_to;
                p = 0;
                while (as > fwd_ambig_from) {
                    --as;
                    u8 c = *as;
                    hbn_assert(c < BLASTNA_SIZE);
                    c = BLASTNA_REVERSE_COMPLEMENT_TABLE[c];
                    seq_blk->sequence_nomask[seq_idx + p] = c;
                    ++p;
                }
                hbn_assert(p == ctx_info.query_length);
            }
            seq_idx += ctx_info.query_length;
        }
    }
    kv_destroy(seq);
    hbn_assert(seq_idx == length);

    seq_info->first_context = 0;
    seq_info->last_context = ctx_idx - 1;
    seq_info->num_queries = num_seq;
    seq_info->max_length = max_length;
    seq_info->min_length = min_length;  

    return num_seq;
}