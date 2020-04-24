#ifndef THREAD_STRUCTURE_H
#define	THREAD_STRUCTURE_H

#include <pthread.h>

#include "query_info.h"
#include "index.h"
#include "result_format.h"
#include "seq_masker.hpp"

static inline int
BindCPU(int tid)
{
#ifdef __USE_GNU
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(tid, &mask);
    return pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask);
#else
	return 0;
#endif
}

void* ThreadFunc(void* threadData)
{
    BindCPU(((SearchWorker*)threadData)->thread_id);
    SearchWorker* sw = (SearchWorker*)threadData;
    sw->Go();
    
    return NULL;
}

// This data structure will be shared by all the threads
struct ThreadCommonData
{
    Options* options;
    FMIndex* fmindex;
    DbInfo* dbinfo;
    OutputFormat* result_writter;
    
    CSymDustMasker** dust_maskers;
    CSeqMasker** window_maskers;
};

ThreadCommonData*
ThreadCommonDataNew(int argc, const char** argv)
{
    Options* opts = new Options();
    opts->ParseCmdLineArgs(argc, argv);

	cy_utility::Log::LogMsg(cy_utility::kAligner, "Loading database.");
	cy_utility::Timer timer;
	timer.start();
	
    FMIndex* index = new FMIndex(opts->input_options->db);
    DbInfo* dbinfo = index->GetDbInfo();
    dbinfo->LoadDbInfo();
    
    OutputFormat* out = new OutputFormat(opts->output_options, opts->hit_options, dbinfo);
    
    ThreadCommonData* retval = new ThreadCommonData;
    retval->options = opts;
    retval->fmindex = index;
    retval->dbinfo = dbinfo;
    retval->result_writter = out;
    
    index->RestoreBwt2();
    index->RestoreSa();
    
    retval->dust_maskers = NULL;
    retval->window_maskers = NULL;
    
    if (opts->filtering_options->windowMaskerOptions->database != NULL)
    {
		retval->window_maskers = new CSeqMasker*[opts->running_options->num_threads];
		for (int i = 0; i < opts->running_options->num_threads; ++i)
			retval->window_maskers[i] = s_BuildSeqMasker(opts->filtering_options->windowMaskerOptions->database);
	} else if (opts->filtering_options->mask_at_seeding)
    {
	retval->dust_maskers = new CSymDustMasker*[opts->running_options->num_threads];
	for (int i = 0; i < opts->running_options->num_threads; ++i)
        retval->dust_maskers[i] = 
               new CSymDustMasker (opts->filtering_options->dustOptions->kDustLevel,
                                   opts->filtering_options->dustOptions->kDustWindow,
                                   opts->filtering_options->dustOptions->kDustLinker);             
    }
	
	timer.end();
	double loading_time = timer.get_elapsed_time();
	cy_utility::Log::LogMsg(cy_utility::kAligner, "done. Time elapsed: %.2f secs.", loading_time);
	cy_utility::Log::LogMsg(NULL, "\n");
    
    return retval;
}

ThreadCommonData*
ThreadCommonDataDelete(ThreadCommonData* data)
{
    data->fmindex->Destroy();
	
	if (data->window_maskers)
	{
		for (int i = 0; i < data->options->running_options->num_threads; ++i)
			delete data->window_maskers[i];
		delete[] data->window_maskers;
		data->window_maskers = NULL;
	}

	if (data->dust_maskers)
	{
		for (int i = 0; i < data->options->running_options->num_threads; ++i)
			delete data->dust_maskers[i];
		delete[] data->dust_maskers;
		data->dust_maskers = NULL;
	}
    
    if (data == NULL)
        return NULL;
    if (data->options)
        delete data->options;
    if (data->fmindex)
        delete data->fmindex;
    if (data->result_writter)
        delete data->result_writter;

    delete data;
    data = NULL;
    
    return data;
}        

// distribute queries across CPU threads
void SetupThreadQueries(QueryInfo& qin,
                        QueryInfo& qout,
                        int qid,
                        int nthreads)
{
    ASSERT(nthreads > 0);
    ASSERT(qid >= 0);
    ASSERT(qid < nthreads);

    Int4 qstart, qend;
    Int4 total_qs = qin.num_queries;
    Int4 thread_qs = total_qs / nthreads;
    Int4 left = total_qs % nthreads;

    if (qid < left)
    {
        ++thread_qs;
        qstart = thread_qs * qid;
        qend = qstart + thread_qs - 1;
    }
    else
    {
        qstart = left * (thread_qs + 1) + (qid - left) * thread_qs;
        qend = qstart + thread_qs - 1;
    }

    qout.Clear();
    qout.query_names.destroy();
    qout.query_offsets.destroy();
    qout.org_queries.destroy();
    
    if (qend  < qstart) return;

    Int4 qs_offset = qin.query_offsets[qstart].query_offset;
    Int4 qn_offset = qin.query_offsets[qstart].header_offset;
    
    char* names = (char*)qin.query_names.get_data() + qn_offset;
    QueryOffset* offsets = (QueryOffset*)qin.query_offsets.get_data() + qstart;
    char* org_qs = (char*)qin.org_queries.get_data() + qs_offset;

    qout.query_names.set_data(names, qin.query_names.size(), qin.query_names.size());
    qout.query_offsets.set_data(offsets, thread_qs, thread_qs);
    qout.org_queries.set_data(org_qs, qin.org_queries.size(), qin.org_queries.size());

    Int4 max_len = 0;
    Int4 i;
    for (i = qstart; i <= qend; ++i)
    {
       qout.query_offsets[i - qstart].header_offset -= qn_offset;
       qout.query_offsets[i - qstart].query_offset -= qs_offset;
	   qout.query_offsets[i - qstart].result_offset = -1;
	   qout.query_offsets[i - qstart].num_alignments = 0;
       if (qout.query_offsets[i - qstart].query_length > max_len)
           max_len = qout.query_offsets[i - qstart].query_length;
    }

    qout.first_context = 0;
    qout.last_context = 2 * thread_qs - 1;
    qout.num_queries = thread_qs;
    qout.max_length = max_len;
}

void CleanUpThreadQueryInfo(QueryInfo& q)
{
    q.query_names.set_data(NULL, 0, 0, TRUE);
    q.query_offsets.set_data(NULL, 0, 0, TRUE);
    q.org_queries.set_data(NULL, 0, 0, TRUE);
}

#endif	/* THREAD_STRUCTURE_H */

