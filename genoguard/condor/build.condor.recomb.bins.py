#!/usr/bin/python
#
# build.condor.dataset.py
#
# 



import os
import sys
import string

condor_header="""
####################
##
## .condor file 
##
####################


#
# your username is used in pathname
#
User	= zhihuang


# Check the Condor primer on which universe to choose
# (standard / vanilla)
Universe        = vanilla


#
# Edit these value
#
InputDir	= /home/$(User)/honeygenes


#
# Edit these value
#
OutputDir      = /home/$(User)/honeygenes/output


# The absolute path (not relative to InitialDir!) to
# your Executable
Executable      = $(InputDir)/model_test.py


# 
# Do not edit
# 
InitialDir      = $(InputDir)


Error           = $(OutputDir)/err.$(Process)
Log             = $(OutputDir)/log.$(Process)
Output          = $(OutputDir)/out.$(Process) 

# This is to be turned on.
GetEnv		= true

# Transfer-input-files is a list of all files being
# used by your job other than the Executable.
#Transfer-Input-Files = 
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = /usr/lib/libpython2.6.so.1.0, /usr/lib/libssl.so.0.9.8, /usr/lib/libcrypto.so.0.9.8, /home/zhihuang/honeygenes/GenoGuard.py

#
# End of the header
#

"""



condor_process_template="""
#
# %s
#
#Transfer-Output-Files 	= $(OutputDir)/%s
Arguments		= %s
Queue 1"""


def build_condor_file():
	print condor_header
	DatasetSize = 165
	condor_process = 1
	
	# recomb
	for i in range(1, 23):
		for j in range(10):
			comment='Condor process : %s, chromosome: %s, bin: %s' %( condor_process, i, j)
			outputfile = """$(OutputDir)/big_recomb_chr%s_bin_seq_%s_CEU.txt.txt""" % (i, j)
			arguments = """clusterCalcModelBins binSeqFileName=%s binIdx=%s modelName=RecombGenoGuard genotypeFileName=$(InputDir)/hapmap/chr%s/big_genotypes_chr%s_CEU.txt haplotypeFileName=$(InputDir)/hapmap/chr%s/big_CEU.chr%s.hap geneticDistFileName=$(InputDir)/hapmap/chr%s/big_genetic_map_chr%s_combined_b36.txt intSeqsFileName=$(InputDir)/hapmap/chr%s/big_int_seqs_chr%s_CEU.pickle""" % (outputfile, j, i, i, i, i, i, i, i, i)
			print condor_process_template % ( comment, outputfile, arguments )
			condor_process += 1


if __name__ == '__main__' :
	build_condor_file()
