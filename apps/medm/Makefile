# 
# MEDM symlink builder
#
# Files are only linked from <name>App/medm/.  Majority of EPICS make options 
# are allowed to ensure compatibility with the build system.  
# 

OLD_FILES_WITH_PATH = $(shell find /usr/local/clas/devel/R3.14.8.2/EPICS/app/medm/ -name "*.adl")
ADL_FILES_WITH_PATH = $(shell find ../*App/medm/ -name "*.adl" -type f)
GIF_FILES_WITH_PATH = $(shell find ../*App/medm/ -name "*.gif" -type f)
SH_FILES_WITH_PATH  = $(shell find ../*App/medm/ -name "*.sh"  -type f)
# Special files needed by flasher scripts
DAT_FILES_WITH_PATH = $(shell find ../flasherApp/medm/ -name "*.dat" -type f)

# Used to all git tag to adl
HPS_MAIN_WITH_PATH = $(shell find ../*App/medm/ -name "HPS_EPICS.adl" -type f)
GIT_TAG = $(shell git describe --abbrev=0)


OLD_LINKS_COMMAND = $(foreach file_with_path, \
        $(OLD_FILES_WITH_PATH), \
        ln -s $(file_with_path) $(notdir $(file_with_path)) ; )


# This uses -f to overwrite old screens from clas...
SOFT_LINK_COMMAND = $(foreach file_with_path, \
	$(ADL_FILES_WITH_PATH) \
	$(GIF_FILES_WITH_PATH) \
	$(SH_FILES_WITH_PATH) \
	$(DAT_FILES_WITH_PATH), \
	ln -sf $(file_with_path) $(notdir $(file_with_path)) ; )

MEDM_HTTP_COMMAND = $(foreach file_with_path, \
	$(ADL_FILES_WITH_PATH) $(GIF_FILES_WITH_PATH) $(SH_FILES_WITH_PATH), \
	cp $(file_with_path) /u/group/clas/www/clasweb/html/medm)

all: clean
	$(OLD_LINKS_COMMAND)
	$(SOFT_LINK_COMMAND)
	$(RM) HPS_EPICS.adl
	cp $(HPS_MAIN_WITH_PATH) .
	sed -i 's/v[0-999].[0-999].[0-999]/$(GIT_TAG)/' HPS_EPICS.adl

http:
	$(RM) /u/group/clas/www/clasweb/html/medm/*.adl
	$(RM) /u/group/clas/www/clasweb/html/medm/*.gif
	$(MEDM_HTTP_COMMAND)

clean:
	$(RM) *.adl
	$(RM) *.dat
	$(RM) *.gif
	$(RM) *.sh

distclean: clean
install:   all
rebuild:   clean install
realclean: clean

# fake out other EPICS make commands
.PHONY: archclean buildInstall cvsclean help realuninstall uninstall
