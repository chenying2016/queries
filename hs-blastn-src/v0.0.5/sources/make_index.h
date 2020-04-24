#include "index.h"
#include "utility.h"

void fmd_index_print_help()
{
	fprintf(stderr, "\n\n");
	fprintf(stderr, "usage:\n");
	fprintf(stderr, "./hs-blastn index database_name\n\n");
	fprintf(stderr, "\n\n");
}

static Uint8 Power(const char* arg)
{
	int n = atoi(arg);
	ASSERT(n >= 0);
	Uint8 r = 1;
	for (int i = 0; i < n; ++i)
		r <<=1;
	return r;
}

/// main function for index building

int make_index(int argc, const char** argv)
{	
	if (argc < 2) 
	{
		fmd_index_print_help();
		return 1;
	}
	
	cy_utility::Timer timer;
	timer.start();
	
	std::clog << std::endl;
	cy_utility::Log::LogMsg("IndexBuilder", "A genomic database index builder.");
	std::clog << std::endl;

	// Build the index
	FMIndex* fmindex = new FMIndex(argv[1]);
	fmindex->BuildIndex();
	delete fmindex;
	
	timer.end();
	double dur = timer.get_elapsed_time();
	fprintf(stderr, "[IndexBuilder] Time elapsed: %f secs.\n", dur);
	
	return 0;
}

