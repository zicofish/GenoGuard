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
Executable      = $(InputDir)/cluster_inference.py


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
transfer_input_files = /usr/lib/libpython2.6.so.1.0, /usr/lib/libssl.so.0.9.8, /usr/lib/libcrypto.so.0.9.8

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
#   # Test Data
# 	# 10%
# 	for i in range(10):
# 		for j in range(65):
# 			comment='Condor process : %s, 10 percent, i: %s, j: %s' %( condor_process, i, j)
# 			outputfile = """$(OutputDir)/predict_recomb_chr22_CEU_10percent_%s_%s.txt""" % (i, j)
# 			arguments = """loadAndPredict predictedSNVsFileName=%s genotypeFileName=$(InputDir)/hapmap/chr22/big_genotypes_chr22_CEU.txt haplotypeFileName=$(InputDir)/hapmap/chr22/big_CEU.chr22.hap geneticDistFileName=$(InputDir)/hapmap/chr22/big_genetic_map_chr22_combined_b36.txt hiddenSNVsFileName=$(InputDir)/hapmap/chr22/infer_test/HiddenSNVs10percent_chr22_%s.txt patientIdx=%s""" % (outputfile, i, j)
# 			print condor_process_template % ( comment, outputfile, arguments )
# 			condor_process += 1
# 	# 40%
# 	for i in range(10):
# 		for j in range(65):
# 			comment='Condor process : %s, 40 percent, i: %s, j: %s' %( condor_process, i, j)
# 			outputfile = """$(OutputDir)/predict_recomb_chr22_CEU_40percent_%s_%s.txt""" % (i, j)
# 			arguments = """loadAndPredict predictedSNVsFileName=%s genotypeFileName=$(InputDir)/hapmap/chr22/big_genotypes_chr22_CEU.txt haplotypeFileName=$(InputDir)/hapmap/chr22/big_CEU.chr22.hap geneticDistFileName=$(InputDir)/hapmap/chr22/big_genetic_map_chr22_combined_b36.txt hiddenSNVsFileName=$(InputDir)/hapmap/chr22/infer_test/HiddenSNVs40percent_chr22_%s.txt patientIdx=%s""" % (outputfile, i, j)
# 			print condor_process_template % ( comment, outputfile, arguments )
# 			condor_process += 1

# 	# Training Data
# 	# 10%
# 	for i in range(10):
# 		for j in range(100):
# 			comment='Condor process : %s, 10 percent, i: %s, j: %s' %( condor_process, i, j)
# 			outputfile = """$(OutputDir)/predict_recomb_chr22_CEU_10percent_%s_%s.txt""" % (i, j)
# 			arguments = """loadAndPredict predictedSNVsFileName=%s genotypeFileName=$(InputDir)/hapmap/chr22/big_genotypes_chr22_CEU.txt haplotypeFileName=$(InputDir)/hapmap/chr22/big_CEU.chr22.hap geneticDistFileName=$(InputDir)/hapmap/chr22/big_genetic_map_chr22_combined_b36.txt hiddenSNVsFileName=$(InputDir)/hapmap/chr22/infer_training/HiddenSNVs10percent_chr22_%s.txt patientIdx=%s""" % (outputfile, i, j)
# 			print condor_process_template % ( comment, outputfile, arguments )
# 			condor_process += 1
# 	# 40%
# 	for i in range(10):
# 		for j in range(100):
# 			comment='Condor process : %s, 40 percent, i: %s, j: %s' %( condor_process, i, j)
# 			outputfile = """$(OutputDir)/predict_recomb_chr22_CEU_40percent_%s_%s.txt""" % (i, j)
# 			arguments = """loadAndPredict predictedSNVsFileName=%s genotypeFileName=$(InputDir)/hapmap/chr22/big_genotypes_chr22_CEU.txt haplotypeFileName=$(InputDir)/hapmap/chr22/big_CEU.chr22.hap geneticDistFileName=$(InputDir)/hapmap/chr22/big_genetic_map_chr22_combined_b36.txt hiddenSNVsFileName=$(InputDir)/hapmap/chr22/infer_training/HiddenSNVs40percent_chr22_%s.txt patientIdx=%s""" % (outputfile, i, j)
# 			print condor_process_template % ( comment, outputfile, arguments )
# 			condor_process += 1
			
	# Test for recombination model with varying hidden percentage from 0.05 to 1.0
	for i in range(1, 21):
		for j in range(65):
			comment='Condor process : %s, i: %s, j: %s' %( condor_process, i, j)
			outputfile = """$(OutputDir)/predict_recomb_chr22_CEU_%s_%s.txt""" % (i, j)
			arguments = """loadAndPredict predictedSNVsFileName=%s genotypeFileName=$(InputDir)/hapmap/chr22/big_genotypes_chr22_CEU.txt haplotypeFileName=$(InputDir)/hapmap/chr22/big_CEU.chr22.hap geneticDistFileName=$(InputDir)/hapmap/chr22/big_genetic_map_chr22_combined_b36.txt hiddenSNVsFileName=$(InputDir)/hapmap/chr22/infer_test/hiddenSNVs_s100_e165_chr22_%s.txt patientIdx=%s""" % (outputfile, i*0.05, j)
			print condor_process_template % ( comment, outputfile, arguments )
			condor_process += 1

if __name__ == '__main__' :
	build_condor_file()
