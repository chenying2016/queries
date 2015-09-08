#include "search_worker.h"
#include "thread_structure.h"
#include "align.h"
#include "make_index.h"
#include <unistd.h>

using namespace std;

void print_main_help(const char* aligner)
{
	cerr << endl;
	cerr << cy_utility::kAligner << ": A nucleotide database search tool." << endl;
	cerr << "Usage:" << endl;
	cerr << "\t1) Build an FMD-index:" << endl;
	cerr << "\t  " << aligner << " index database" << endl;
	cerr << endl;
	cerr << "\t2) Search the database:" << endl;
	cerr << "\t  " << aligner << " align -h" << endl;
	cerr << "\t  or" << endl;
	cerr << "\t  " << aligner << " align -help" << endl;
	cerr << endl;
	cerr << "See README for more details." << endl;
	cerr << endl;
}

void print_command(int argc, const char** argv)
{
	cerr << "\n[CMD] ";
	for (int i = 0; i < argc; ++i)
		cerr << argv[i] << " ";
	cerr << endl;
}

int main(int argc, const char** argv)
{    	
	if (argc < 2)
	{
		print_main_help(argv[0]);
		return 1;
	}
	
	if (strcmp(argv[1], "index") == 0)
	{
		// The task is to build an index
		make_index(argc - 1, argv + 1);
	}
	else if (strcmp(argv[1], "align") == 0)
	{
		// The task is to perform alignment
		align(argc - 1, argv + 1);
	}
	else
	{
		print_main_help(argv[0]);
		return 1;
	}
	
	print_command(argc, argv);
	return 0;
}
