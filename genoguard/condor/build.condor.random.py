#!/usr/bin/python
#
# build.condor.random.py round
#
# build.condor.random.py 0
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
Executable      = $(InputDir)/GenoGuardEval.py


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


def build_condor_file(round):
	print condor_header
	SIZE = 50
	condor_process = 1
	#Generate random samples with PubLDGenoGuard
	for i in range(SIZE):
		comment='Condor process : %s, (i=%s)' %( condor_process, i)
		outputfile = """$(OutputDir)/random-probs-publd-%s.txt""" % (round*SIZE + i)
		arguments = """generateModelSample outFileName=%s modelName=PubLDGenoGuard genotypeFileName=$(InputDir)/data/genotypes_chr22_CEU_small.txt AFFileName=$(InputDir)/data/allele_freqs_chr22_CEU_small.txt LDFileName=$(InputDir)/data/ld_chr22_CEU_small.txt sampleSize=10""" % (outputfile)
		print condor_process_template % ( comment, outputfile, arguments )
		condor_process += 1
# 	#Generate random samples with 0-order DirectCondProbGenoGuard
# 	for i in range(SIZE):
# 		comment='Condor process : %s, (i=%s)' %( condor_process, i+SIZE)
# 		outputfile = """$(OutputDir)/random-probs-imm0-direct-%s.txt""" % (round*SIZE + i)
# 		arguments = """generateModelSample outFileName=%s modelName=DirectCondProbGenoGuard genotypeFileName=$(InputDir)/data/new_genotypes_chr22_CEU.txt probFileName=$(InputDir)/data/imm0_new_chr22_CEU_ref.txt order=0 sampleSize=10""" % (outputfile)
# 		print condor_process_template % ( comment, outputfile, arguments )
# 		condor_process += 1
		
	#Generate random samples with 1-order DirectCondProbGenoGuard
	for i in range(SIZE):
		comment='Condor process : %s, (i=%s)' %( condor_process, i+SIZE*2)
		outputfile = """$(OutputDir)/random-probs-imm1-direct-%s.txt""" % (round*SIZE + i)
		arguments = """generateModelSample outFileName=%s modelName=DirectCondProbGenoGuard genotypeFileName=$(InputDir)/data/genotypes_chr22_CEU_small.txt probFileName=$(InputDir)/data/condProb1_chr22_CEU_ref.txt order=1 sampleSize=10""" % (outputfile)
		print condor_process_template % ( comment, outputfile, arguments )
		condor_process += 1
		
	#Generate random samples with 2-order DirectCondProbGenoGuard
	for i in range(SIZE):
		comment='Condor process : %s, (i=%s)' %( condor_process, i+SIZE*3)
		outputfile = """$(OutputDir)/random-probs-imm2-direct-%s.txt""" % (round*SIZE + i)
		arguments = """generateModelSample outFileName=%s modelName=DirectCondProbGenoGuard genotypeFileName=$(InputDir)/data/genotypes_chr22_CEU_small.txt probFileName=$(InputDir)/data/condProb2_chr22_CEU_ref.txt order=2 sampleSize=10""" % (outputfile)
		print condor_process_template % ( comment, outputfile, arguments )
		condor_process += 1


if __name__ == '__main__' :
	arg=sys.argv
	build_condor_file(int(arg[1]))
