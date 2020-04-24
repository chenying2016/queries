ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := hs-blastn
SOURCES  := \
	backup_results.c \
	cmdline_args.cpp \
	find_seeding_subseqs.cpp \
	hbn_build_seqdb.c \
	hbn_extend_subseq_hit.c \
	hbn_find_subseq_hit.c \
	hbn_job_control.c \
	hbn_options_handle.c \
	hbn_task_struct.c \
	main.c \
	map_one_volume.c \
	hbn_results.c \
	search_setup.c \
	subseq_hit.cpp \
	symdust.cpp \
	tabular_format.cpp \
	traceback_stage.c \

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=