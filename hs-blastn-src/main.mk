ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libhbn.a

SOURCES      := \
	./corelib/build_db.c \
	./corelib/cstr_util.c \
	./corelib/db_format.c \
	./corelib/fasta.c \
	./corelib/gapped_candidate.c \
	./corelib/hbn_aux.c \
	./corelib/hbn_package_version.c \
	./corelib/kstring.c \
	./corelib/line_reader.c \
	./corelib/m4_record.c \
	./corelib/name2id_map.c \
	./corelib/partition_aux.c \
	./corelib/raw_reads.c \
	./corelib/seqdb_summary.c \
	./corelib/seqdb.c \
	./corelib/small_object_alloc.c \
	./algo/align.c \
	./algo/chain_dp.c \
	./algo/dalign.c \
	./algo/edlib.cpp \
	./algo/edlib_wrapper.c \
	./algo/hash_list_bucket_sort.c \
	./algo/hbn_traceback.c \
	./algo/hbn_traceback_aux.c \
	./algo/init_hit_finder.c \
	./algo/hbn_lookup_table.c \
	./algo/sort_sr_hit_seeds.cpp \
	./algo/word_finder.c \
	./ncbi_blast/c_ncbi_blast_aux.c \
	./ncbi_blast/ncbi_blast_aux.cpp \
	./ncbi_blast/cmdline_args/blast_args.cpp \
	./ncbi_blast/cmdline_args/cmdline_flags.cpp \
	./ncbi_blast/cmdline_args/format_flags.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_allow.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_desc.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_types.cpp \
	./ncbi_blast/str_util/ncbistr_util.cpp \
	./ncbi_blast/str_util/ncbistr.cpp \
	./ncbi_blast/str_util/str_cmp.cpp \
	./ncbi_blast/str_util/numeric_str_interconv.cpp \
	./ncbi_blast/str_util/str_util.cpp \
	./ncbi_blast/setup/blast_encoding.c \
	./ncbi_blast/setup/blast_hits.c \
	./ncbi_blast/setup/blast_message.c \
	./ncbi_blast/setup/blast_options.c \
	./ncbi_blast/setup/blast_stat.c \
	./ncbi_blast/setup/blast_parameters.c \
	./ncbi_blast/setup/blast_program.c \
	./ncbi_blast/setup/blast_types.cpp \
	./ncbi_blast/setup/boost_erf.c \
	./ncbi_blast/setup/hsp2string.cpp \
	./ncbi_blast/setup/ncbi_math.c \
	./ncbi_blast/setup/blast_query_info.c \
	./ncbi_blast/setup/blast_sequence_blk.c \
	./ncbi_blast/setup/gapinfo.c

SRC_INCDIRS  := ./third_party/spreadsortv2

SUBMAKEFILES := ./app/primer_map/main.mk ./app/hbnmap/main.mk