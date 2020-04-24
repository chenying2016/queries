#include "hbn_job_control.h"
#include "cmdline_args.h"
#include "hbn_build_seqdb.h"
#include "hbn_task_struct.h"
#include "hbn_find_subseq_hit.h"
#include "map_one_volume.h"
#include "hbn_results.h"
#include "../../corelib/hbn_package_version.h"

int main(int argc, char* argv[])
{
    HbnProgramOptions* opts = (HbnProgramOptions*)calloc(1, sizeof(HbnProgramOptions));
    ParseHbnProgramCmdLineArguments(argc, argv, opts);
    hbn_build_seqdb(opts, INIT_QUERY_DB_TITLE, INIT_SUBJECT_DB_TITLE);

    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/%s", opts->db_dir, kBackupResultsDir);
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }

    hbn_task_struct* task_struct = hbn_task_struct_new(opts);
    if (opts->outfmt == eSAM) {
        print_sam_prolog(task_struct->out, 
            kSamVersion, 
            HBN_PACKAGE_VERSION, 
            opts->rg_info,
            opts->rg_sample,
            argc, 
            argv);
    }
    const int num_query_vols = seqdb_load_num_volumes(opts->db_dir, task_struct->query_db_title);
    const int num_subject_vols = seqdb_load_num_volumes(opts->db_dir, task_struct->subject_db_title);
    const int query_vol_stride = opts->num_nodes;
    const int subject_vol_stride = 1;
    char job_name[256];

    for (int svid = 0; svid < num_subject_vols; svid += subject_vol_stride) {
        int qvid = (task_struct->query_and_subject_are_the_same ? svid : 0) + opts->node_id;
        HBN_LOG("Searching against S%s", u64_to_fixed_width_string(svid, HBN_DIGIT_WIDTH));
        if (all_vs_sj_is_mapped(opts->db_dir, 
                kBackupResultsDir,
                qvid,
                num_query_vols,
                svid,
                opts->node_id,
                opts->num_nodes)) {
            merge_all_vs_sj_results(opts->db_dir, kBackupResultsDir, qvid, num_query_vols, svid, opts->node_id, opts->num_nodes, opts, task_struct->out);
            continue;
        }
        
        hbn_task_struct_build_subject_vol_context(task_struct, svid);
        for (; qvid < num_query_vols; qvid += query_vol_stride) {
            if (qi_vs_sj_is_mapped(opts->db_dir, kBackupResultsDir, qvid, svid)) {
                merge_qi_vs_sj_results(opts->db_dir, kBackupResultsDir, qvid, svid, task_struct->subject_vol, opts, task_struct->out);
                continue;
            }
            char qibuf[64], sjbuf[64];
            u64_to_fixed_width_string_r(qvid, qibuf, HBN_DIGIT_WIDTH);
            u64_to_fixed_width_string_r(svid, sjbuf, HBN_DIGIT_WIDTH);
            sprintf(job_name, "Q%s_vs_S%s", qibuf, sjbuf);
            hbn_timing_begin(job_name);
            hbn_task_struct_build_query_vol_context(task_struct, qvid);
            hbn_align_one_volume(task_struct);
            qi_vs_sj_make_mapped(opts->db_dir, kBackupResultsDir, qvid, svid);
            hbn_timing_end(job_name);
        }
    }

    task_struct = hbn_task_struct_free(task_struct);
    free(opts);
    return 0;
}