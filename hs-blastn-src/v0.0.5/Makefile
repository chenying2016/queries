OBJECT_DIR      := objs

#EXTRA_CFLAGS    := -pg -g -ggdb
EXTRA_CFLAGS    := -O3
EXTRA_CFLAGS    += -Wall -Wno-unused -D__STDC_FORMAT_MACROS
LDFLAGS         :=  -pthread
src-y           := ./sources/ ./sources/commandline_options/ ./sources/mask/src/dustmask/ ./sources/mask/src/winmask/ ./sources/mask/
inc-y           := ./sources/ ./sources/commandline_options/ ./sources/mask/include/dustmask/ ./sources/mask/include/winmask/ ./sources/mask/
DEFINES         := -DFLAG1 -UFLAG2 -DPROJECT='"hs-blastn"'
TARGET_TYPE     := EXE
TARGET          := hs-blastn
SOURCE_SUFFIXES := cpp c

#==============================================================================
#          Do not modify the following text.
#==============================================================================

###############################################################################
#
# A FLEXIBLE MAKEFILE TEMPLATE
#
# The purpose of implementing this script is help quickly deploy source code
# tree during initial phase of development. It is designed to manage one whole
# project from within one single makefile and to be easily adapted to
# different directory hierarchy by simply setting user configurable variables.
# This script is expected to be used with gcc toolchains on bash-compatible 
# shell.
# 
# Author: Pan Ruochen <coderelease@163.com>
# Date:   2012/10/10
#
###############################################################################

#-----------------------------------------------------------------------------------------------------#
# User configurable variables
#
# ====================================================================================================
# GNU_TOOLCHAIN_PREFIX:   The perfix of gnu toolchain.
# ====================================================================================================
# DEFINES:        The compiler flags for macro definitions.
# EXTRA_CFLAGS:   Any other compiler flags.
# inc-y:          Header include paths.
# src-y:          Sources. The items ending with a trailing / are regarded as directories, the others
#                 are regareded as files. The files with specified suffixes in those directories will
#                 be automatically involved in compilation.
# obj-y:          Extra object file list.
# ucmd_X:         User defined command to generate targets for the prerequisites
#                 whith the specified suffix X (i.e, X could be c, cpp, etc).
#
# EXCLUDE_FILES:  The files that are not included during compilation.
# OBJECT_DIR:     The directory where object files are output.
# LD_SCRIPT:      The explicit linker script for linking.
# LIBS:           The libraries for linking.
# LDFLAGS:        All other linker flags.
# ====================================================================================================
# STRIP_UNUSED:   Remove all unreferenced functions and data during linking.
# SOURCE_SUFFIXES:The suffixes of source files.
# OBJECT_SUFFIX:  The suffix of object files.
# DEPEND_SUFFIX:  The suffix of dependency files.
# TARGET_TYPE:    The target type which can be application, shared object, archive library etc.
# TARGET:         The path name of the final target.
# IGNORE_ME:      The changes of this script will not cause remaking of any target.
# CENTRALIZED_SINGLE_DEPEND_FILE:  Use one single dependency file instead of 
#                                  generating one dependency file for each source file.
# TARGET_DEPENDS: The dependent targets by the final target.
# VERBOSE_COMMAND:Display verbose commands instead of short commands during the make process.
#-----------------------------------------------------------------------------------------------------#

#****************************************************************************#
#  PART II: FUNCTIONALITY IMPLEMENTATIONS                                    #
#****************************************************************************#

# Quiet commands
ifeq ($(VERBOSE_COMMAND),)
Q           = @
Q_compile   = @echo '  CC     $$< => $$@';
Q_link      = @echo '  LD     $@';
Q_ar        = @echo '  AR     $@';
Q_mkdir     =  echo '  MKDIR  $1';
Q_clean     = @echo '  CLEAN';
Q_distclean = @echo '  DISTCLEAN';
endif

O := $(if $(OBJECT_SUFFIX),$(OBJECT_SUFFIX),o)
D := $(if $(DEPEND_SUFFIX),$(DEPEND_SUFFIX),d)

ifndef SOURCE_SUFFIXES
SOURCE_SUFFIXES := c cpp cc cxx S s
endif

GCC    := $(GNU_TOOLCHAIN_PREFIX)gcc

src-d = $(filter %/,$(src-y))
src-f = $(foreach i,$(SOURCE_SUFFIXES),$(filter %.$i,$(src-y)))

is_equal = $(if $(filter $1,$2),$(filter $2,$1))

objdir := $(shell echo $(OBJECT_DIR)|sed -e 's:\(\./*\)*::g')
ifeq ($(objdir),)
objdir       := ./
else
objdir       := $(objdir)/
have_objdir  := y
endif

## Combine compiler flags togather.
CFLAGS   = $(foreach i,$(inc-y),-I$i) $(EXTRA_CFLAGS) $(DEFINES)

## Output file types:
##  EXE:  Application
##  AR:   static library
##  SO:   shared object
##  DLL:  dynamic link library
##  BIN:  raw binary
TARGET_TYPE := $(strip $(TARGET_TYPE))
ifeq ($(filter $(TARGET_TYPE),SO DLL AR EXE BIN),)
$(error Unknown TARGET_TYPE `$(TARGET_TYPE)')
endif

ifneq ($(filter DLL SO,$(TARGET_TYPE)),)
CFLAGS  += -shared
LDFLAGS += -shared
endif
ifneq ($(STRIP_UNUSED),)
CFLAGS  += -ffunction-sections -fdata-sections
LDFLAGS += --gc-sections
endif

ifeq ($(CENTRALIZED_SINGLE_DEPEND_FILE),)
CFLAGS += -MMD -MF $$@.$(D) -MT $$@
else
single_depend_file := $(objdir)depend
endif

g_makefile_list = $(if $(IGNORE_ME),,$(MAKEFILE_LIST))

#--------------------------------------------------#
# Exclude user-specified files from source list.   #
#  $1 -- The sources list                          #
#--------------------------------------------------#
exclude = $(filter-out $(EXCLUDE_FILES),$1)

#----------------------------------------------------------#
# List files with specified suffix inside the directory.   #
#  $1 -- The directory                                     #
#  $2 -- The suffix                                        #
#----------------------------------------------------------#
ls = $(wildcard $1*.$2)


#---------------------------------------------#
# Replace the specified suffixes with $(O).   #
#  $1 -- The file names                       #
#  $2 -- The suffixes                         #
#---------------------------------------------#
get_object_names = $(strip $(foreach i,$2,$(patsubst %.$i,%.$O,$(filter %.$i,$1))))

#---------------------------------------------#
# Get the suffix name from a file name.       #
#  $1 -- The file name                        #
#  $2 -- The favorite suffixes                #
#---------------------------------------------#
get_suffix_names = $(strip $(foreach i,$2,$(if $(filter %.$i,$1),$i)))

#-------------------------------------------------------------------#
# Replace the pattern .. with !! in the path names in order that    #
# no directories are out of the object directory                    #
#  $1 -- The path names                                             #
#-------------------------------------------------------------------#
objdir_transform = $(if $(have_objdir),$(subst ..,!!,$1),$1)


#------------------------------------------------------------------#
# Set up static pattern rules for sources with specified suffixes  #
# in specified directories.                                        #
#  $1 -- Source directories                                        #
#  $2 -- Source suffixes                                           #
#  $3 -- Equal to $(call ls $1,$2)                                 #
#------------------------------------------------------------------#
static_pattern_rules = $(if $3,$(call __static_pattern_rule,$(patsubst %.$2,$(objdir)%.$O,$3),$1,$2))


#------------------------------------#
# Command to make directory          #
#  $1 -- The directory to be made    #
#------------------------------------#
define cmd_make_directory
$(Q)if test ! -d "$1"; then $(Q_mkdir)mkdir -p "$1"; fi

endef

cmd_compile = $(Q_compile)$(if $(ucmd_$1),$(ucmd_$1),$(GCC) -I$$(dir $$<) $(CFLAGS) -c -o $$@ $$<)

#------------------------------------------------------------------#
#  Static pattern rule                                             #
#  $1 -- Targets                                                   #
#  $1 -- Source directories                                        #
#  $3 -- The source suffix                                         #
#------------------------------------------------------------------#
define __static_pattern_rule
$(call objdir_transform,$1): $(call objdir_transform,$(objdir)$2%.$(O)): $2%.$3 $(g_makefile_list)
	$(call cmd_compile,$3)

endef


#--------------------------------------------------------------#
#  Ordinary rule                                               #
#  $1 -- The prerequisite                                      #
#  $2 -- The Target                                            #
#--------------------------------------------------------------#
define ordinary_rule
$(call objdir_transform,$2): $1 $(g_makefile_list)
	$(call cmd_compile,$(call get_suffix_names,$1,$(SOURCE_SUFFIXES)))

endef

#--------------------------------------------------------#
# Make sure the default target "all" is the first target
#--------------------------------------------------------#
PHONY = all clean distclean make_sub_dirs
all: make_sub_dirs $(TARGET)

#----------------------------------------------------#
# Dynamic Targets
#----------------------------------------------------#
$(eval $(foreach i,\
    $(sort $(src-d)),\
    $(foreach j,$(SOURCE_SUFFIXES),$(call static_pattern_rules,$i,$j,$(call exclude,$(call ls,$i,$j)))))\
    $(foreach i,$(call exclude,$(sort $(src-f))),$(call ordinary_rule,$i,$(objdir)$(call get_object_names,$i,$(SOURCE_SUFFIXES)))))


#-------------------------------------#
# Get the list of all source files    #
#-------------------------------------#
srcs = $(call exclude,\
	$(foreach i,$(SOURCE_SUFFIXES),\
	$(foreach j,$(src-d),\
	$(wildcard $j*.$i)) $(filter %.$i,$(src-f))))

ifeq ($(strip $(srcs)),)
$(error Empty source list! Please check both src-y and SOURCE_SUFFIXES are correctly set.)
endif

#-------------------------------------#
# Get the list of all object files    #
#-------------------------------------#
objs = $(call objdir_transform,$(addprefix $(objdir),$(call get_object_names,$(srcs),$(SOURCE_SUFFIXES))))
objs += $(obj-y)

#----------------------------------------------------#
# Static Targets
#----------------------------------------------------#
make_sub_dirs:
	$(call cmd_make_directory,$(dir $(TARGET)))
	$(foreach i,$(call objdir_transform,$(sort $(src-d) $(dir $(src-f)))),$(call cmd_make_directory,$(objdir)$i))

ifneq ($(single_depend_file),)
$(single_depend_file): $(srcs) $(filter-out $@,$(g_makefile_list)) $(objdir)
	$(GCC) $(CFLAGS) -MM -MG $(srcs) | \
sed 's#\([^[:space:]]\+\)\.$O:\s\([^[:space:]]\+\)\.\([^[:space:].]\+\s\?\)#$(objdir)\2.$O: \2.\3#g' > $@
$(objdir): ; $(call cmd_make_directory,$(objdir))
endif

ifeq ($(TARGET_TYPE),AR)
$(TARGET): AR := $(GNU_TOOLCHAIN_PREFIX)ar
$(TARGET): $(TARGET_DEPENDS) $(objs)
	$(Q_ar)rm -f $@ && $(AR) rcvs $@ $(objs)
else

ifeq ($(TARGET_TYPE),BIN)
tmp_target   = $(basename $(TARGET)).elf
LDFLAGS     += -nodefaultlibs -nostdlibs -nostartupfiles
$(TARGET): $(tmp_target)
	$(GNU_TOOLCHAIN_PREFIX)objcopy -O binary $(tmp_target) $@
	$(GNU_TOOLCHAIN_PREFIX)objdump -d $(tmp_target) > $(basename $(@F)).lst
	$(GNU_TOOLCHAIN_PREFIX)nm $(tmp_target) | sort -k1 > $(basename $(@F)).map
else
tmp_target   = $(TARGET)
endif

$(tmp_target): LD = $(if $(foreach i,cpp cc cxx,$(filter %.$i,$(srcs))),$(GNU_TOOLCHAIN_PREFIX)g++,$(GCC))
$(tmp_target): $(TARGET_DEPENDS) $(objs) $(LD_SCRIPT)
	$(Q_link)$(LD) $(LDFLAGS) $(if $(LD_SCRIPT),-T $(LD_SCRIPT)) $(objs) $(LIBS) -o $(tmp_target)

endif

clean:
	$(Q_clean)rm -rf $(filter-out ./,$(objdir)) $(TARGET) $(filter-out $(obj-y),$(objs))
distclean: clean
	$(Q_distclean)find -name '*.$O' -o -name '*.$D' | xargs rm -f; $(if $(single_depend_file),rm -f $(single_depend_file))
print-%:
	@echo $* = $($*)

.DEFAULT_GOAL = all

sinclude $(if $(filter all,$(if $(MAKECMDGOALS),$(MAKECMDGOALS),$(.DEFAULT_GOAL))), \
$(if $(single_depend_file),$(single_depend_file),$(foreach i,$(objs),$i.$(D))))


