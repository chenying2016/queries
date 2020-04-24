#include "../../corelib/line_reader.h"
#include "../../corelib/khash.h"
#include "../../corelib/hbn_hit.h"
#include "../../corelib/fasta.h"
#include "../../corelib/seq_tag_report.h"
#include "../../corelib/string2hsp.h"
#include "../../corelib/ksort.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/seqdb.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>

int main(int argc, char* argv[])
{
    HbnFastaReader* reader = HbnFastaReaderNew(argv[1]);
    int i = 0;
    while (i < 100) {
        HbnFastaReaderReadOneSeq(reader);
        fprintf(stdout, ">");
        hbn_fwrite(ks_s(reader->name), 1, ks_size(reader->name), stdout);
        fprintf(stdout, "\n");
        hbn_fwrite(ks_s(reader->sequence), 1, ks_size(reader->sequence), stdout);
        fprintf(stdout, "\n");
        ++i;
    }
}