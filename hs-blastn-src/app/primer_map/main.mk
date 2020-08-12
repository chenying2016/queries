ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := hs-blastn-primer-map
SOURCES  := \
		blast_gapalign.c \
    cmdline_args.cpp \
		greedy_align.c \
		main.c \
		primer_map_chain_dp.c \
		primer_map_hit_finder.c \
		primer_map_one_volume.c \
		../hbnmap/backup_results.c \
		../hbnmap/hbn_build_seqdb.c \
		../hbnmap/hbn_options_handle.c \
		../hbnmap/hbn_results.c \
		../hbnmap/search_setup.c \
		../hbnmap/tabular_format.cpp \
		../hbnmap/traceback_stage.c

SRC_INCDIRS  := . ../../third_party/spreadsortv2

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=