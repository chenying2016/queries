#include "search_worker.h"
#include "utility.h"
#include "thread_structure.h"

using namespace std;

/*
 * Process one query file.
 * This functions takes query file name query_file_name as input,
 * applys the 
 * seeding-->ungapped_extension-->score_only_gapped_extension-->traceback
 * procedure on each query.
 * The results are written to result_file_name.
 */
int process_one_query_file(ThreadCommonData* run_data,
						   const char* query_file_name,
						   QueryInfo* query_batch,
						   const char* result_file_name,
						   SearchWorker** sws)
{
	query_batch->Init(query_file_name);
	int nts = run_data->options->running_options->num_threads;
	Int8 num_queries = 0;
	FILE* results_file;
	if (result_file_name == NULL) results_file = NULL;
	else 
	{
		results_file = fopen(result_file_name, "w");
		if (results_file == NULL)
		{
			fprintf(stderr, "Fatal Error: Cannot open file %s for writing results.\n", result_file_name);
			exit(1);
		}
	}
	sws[0]->results->PrintProlog();
	sws[0]->results->FlushResults(results_file);

	cy_utility::Timer timer;
	cy_utility::Log::LogMsg(cy_utility::kAligner, "Processing %s.", query_file_name);
	timer.start();	
	int i;
	pthread_t tids[nts];
	
	// While more queries
    while (query_batch->GetQueryBatch(nts) > 0)
    {
		num_queries += query_batch->GetNumQueries();
		
	    cy_utility::Log::LogMsg(NULL, "\tProcessing %d queries.", query_batch->GetNumQueries());

		// To the alignment
        for (i = 0; i < nts; ++i)
        {
            pthread_create(&tids[i], NULL, ThreadFunc, sws[i]);
        }

        for (i = 0; i < nts; ++i)
        {
            pthread_join(tids[i], NULL);
        }

		// Print the results of each thread
        for (i = 0; i < nts; ++i)
        {
            sws[i]->FlushResults(results_file);
            sws[i]->CleanUp();
        }
    }
	
	sws[0]->results->PrintEpilog(num_queries, run_data->options->scoring_options);
	sws[0]->FlushResults(results_file);
	if (results_file != NULL) fclose(results_file);
	results_file = NULL;
	
	timer.end();
	double dur = timer.get_elapsed_time();
    cy_utility::Log::LogMsg(NULL, "\n\n");
    cy_utility::Log::LogMsg(cy_utility::kAligner, "done. Elpased time: %.4f secs.", dur);
	cy_utility::Log::LogMsg(cy_utility::kAligner, "%llu queries processed.", num_queries);
	
	return num_queries;
}

// The main function for the align command
int align(int argc, const char** argv)
{   	
	// The common data that shared by all the threads (if multiple threads is used)
    ThreadCommonData* run_data = ThreadCommonDataNew(argc - 1, argv + 1);
    const Int4 nts = run_data->options->running_options->num_threads;
    ASSERT(nts > 0);

	// Process a batch of queries at a time
    QueryInfo query_batch;
    // Search handler, one for each thread.
    SearchWorker** sws = new SearchWorker*[nts];
    pthread_t tids[nts];

    for (int i = 0; i < nts; ++i)
    {
        sws[i] = new SearchWorker(run_data->options,
								  &query_batch,
                                  run_data->fmindex,
                                  run_data->dbinfo,
                                  run_data->dust_maskers ? run_data->dust_maskers[i] : NULL,
                                  run_data->window_maskers ? run_data->window_maskers[i] : NULL,
								  i);
		sws[i]->thread_id = i;
    }
    
	Int8 total_queries = process_one_query_file(run_data,
											    run_data->options->input_options->query,
											    &query_batch,
											    run_data->options->output_options->output_file_name,
											    sws);

	run_data = ThreadCommonDataDelete(run_data);
    for (int i = 0; i < nts; ++i)
	{
		delete sws[i];
		sws[i] = NULL;
	}
    delete[] sws;

    return 0;
}
