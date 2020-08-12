#include "cmdline_args.h"
#include "primer_map_one_volume.h"
#include "../hbnmap/hbn_build_seqdb.h"
#include "../hbnmap/hbn_results.h"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/seqdb.h"

int main(int argc, char* argv[])
{
    HbnProgramOptions* opts = (HbnProgramOptions*)calloc(1, sizeof(HbnProgramOptions));
    ParseHbnProgramCmdLineArguments(argc, argv, opts);
    hbn_build_seqdb(opts, INIT_QUERY_DB_TITLE, INIT_SUBJECT_DB_TITLE);
    HbnOptionsHandle* opts_handle = HbnOptionsHandleNew(eMegablast);
    HbnOptionsHandle_Update(opts, opts_handle);
    hbn_dfopen(out, opts->output, "w");
    if (opts->outfmt == eSAM) {
        print_sam_prolog(out, 
            kSamVersion, 
            HBN_PACKAGE_VERSION, 
            opts->rg_info,
            opts->rg_sample,
            argc, 
            argv);
    }
    const int num_query_vols = seqdb_load_num_volumes(opts->db_dir, INIT_QUERY_DB_TITLE);
    const int num_subject_vols = seqdb_load_num_volumes(opts->db_dir, INIT_SUBJECT_DB_TITLE);
    const int query_vol_stride = opts->num_nodes;
    const int subject_vol_stride = 1;
    char job_name[256];

    for (int svid = 0; svid < num_subject_vols; svid += subject_vol_stride) {
        CSeqDB* svol = seqdb_load_unpacked_with_ambig_res(opts->db_dir, INIT_SUBJECT_DB_TITLE, svid);
        int qvid = 0;
        HBN_LOG("Searching against S%s", u64_to_fixed_width_string(svid, HBN_DIGIT_WIDTH));
        for (; qvid < num_query_vols; qvid += query_vol_stride) {
            CSeqDB* qvol = seqdb_load(opts->db_dir, INIT_QUERY_DB_TITLE, qvid);
            char qibuf[64], sjbuf[64];
            u64_to_fixed_width_string_r(qvid, qibuf, HBN_DIGIT_WIDTH);
            u64_to_fixed_width_string_r(svid, sjbuf, HBN_DIGIT_WIDTH);
            sprintf(job_name, "Q%s_vs_S%s", qibuf, sjbuf);
            hbn_timing_begin(job_name);
            primer_map_one_volume(qvol, svol, opts, opts_handle, out);
            hbn_timing_end(job_name);
        }
        svol = CSeqDBFree(svol);
    }
    free(opts);
    HbnOptionsHandleFree(opts_handle);
    hbn_fclose(out);
    return 0;
}