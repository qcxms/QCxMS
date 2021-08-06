#!/bin/env make
# 
# created by Sebastian Ehlert in March 2018 (v 2.1)
# 2.0: recreated from scratch
# 2.1: refined the generation of dependencies, since it crashed if not
#      done in the main directory
#
##############################################################################
## this file is part of the xtb build system                                ##
##############################################################################

###------------------------------------------------------------------------###
# specifiy where you are in the source tree
MAINDIR := ..
# here resides the build system
MAKEDIR := $(MAINDIR)/MAKE

# set the source by hand
SRCDIR  := .

##############################################################################
## please edit this Files only:
include $(MAKEDIR)/qcxms.objs
include $(MAKEDIR)/Makeconfig
## here you can give optional configurations independent from Makeconfig
-include $(MAKEDIR)/qcxms.config
##############################################################################

include $(MAKEDIR)/Makerules
