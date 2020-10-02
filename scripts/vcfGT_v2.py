#! /usr/local/bin python

"""
Read in VCF file and mark genotype as 'fail' if there is heterozyogsity
"""
import os
import sys

# Read in the vcf file
file = open(sys.argv[1],'r')

line = file.readline().rstrip()

while line:
	if line.startswith("#"):
		print(line)
		line = file.readline().rstrip()
	else:
		values = line.split("\t")

		sampleGenotypes = []
		for sample in (values[9:]):
			# If quoted incude next line
			#sample = sample[1:-1]
			(GT,AD,DP,GQ,PL) = sample.split(":")
			# print(sample)
			alleleDepthCounts = AD.split(",")

			if DP == "0":
				GT = "."
				
			else:
				allelicRatio = float(alleleDepthCounts[0])/float(DP)

				if allelicRatio >= 0.95:
					GT = "0"
				elif allelicRatio <= 0.05:
					GT = "1"
				else:
					GT = "."
				
			newSample = GT + ":" + AD + ":" + DP  + ":" + GQ + ":" + PL
			sampleGenotypes.append(newSample)


		print('\t'.join(values[:9]) + "\t" + '\t'.join(sampleGenotypes))	

		# Read the next line
		line = file.readline().rstrip()
