#!/usr/bin/env python

import os
import subprocess

orig_dir = '/home/exacloud/lustre1/CompBio/data/mc3_variants/'
current_dir = '/home/exacloud/lustre1/CompBio/data/mc3_variants/Data_Dec2017/'

tn_pairs = '/home/exacloud/lustre1/CompBio/projs/RThompson/immunotherapy/jobs/merged_vcfs/tn_pairs.txt'

redone = []
first_round = []
second_round = []
leftovers = []

with open(tn_pairs, 'r') as f:
	f.readline()
	for line in f:
		tokens = line.strip('\n').split('\t')
		old_dir = ''.join([orig_dir, tokens[0], '/'])
		new_dir = ''.join([current_dir, tokens[2], '/'])
		if os.path.isdir(old_dir) and not os.path.isdir(new_dir):
			first_round.append(tokens[2])
			subprocess.call(['mkdir', new_dir])
			# Copy VCFs
			muse = ''.join([old_dir, '*/root/muse/muse.vcf'])
			mutect = ''.join([old_dir, '*/root/mutect/mutect.vcf'])
			pindel = ''.join([old_dir, '*/root/pindel/pindel_somatic.vcf'])
			radia = ''.join([old_dir, '*/root/radia/radia_filtered.vcf'])
			somatic_sniper = ''.join([old_dir, 
									  '*/root/somaticsniper_fpfilter/filtered.vcf'])
			varscan = ''.join([old_dir, 
									  '*/root/varscan_fpfilter/filtered.vcf'])
			subprocess.call(['cp', muse, new_dir])
			subprocess.call(['cp', mutect, new_dir])
			subprocess.call(['cp', pindel, new_dir])
			subprocess.call(['cp', radia, new_dir])
			new_somaticsniper = ''.join([new_dir, 'somatic_sniper_fpfilter.vcf'])
			subprocess.call(['cp', somatic_sniper, new_somaticsniper])
			new_varscan = ''.join([new_dir, 'varscan_fpfilter.vcf'])
			subprocess.call(['cp', varscan, new_varscan])
			# Copy BAMs
			old_N_bam = ''.join([orig_dir, 'gatk-cocleaning-bams/', tokens[0],
								'_N.reheadered.realigned.cleaned.bam'])
			new_N_bam = ''.join([new_dir, tokens[1],
								'.reheadered.realigned.cleaned.bam'])
			old_T_bam = ''.join([orig_dir, 'gatk-cocleaning-bams/', tokens[0],
								'_T.reheadered.realigned.cleaned.bam'])
			new_T_bam = ''.join([new_dir, tokens[2],
								'.reheadered.realigned.cleaned.bam'])
			subprocess.call(['cp', old_N_bam, new_N_bam])
			subprocess.call(['cp', ''.join([old_N_bam[:-1], 'i']), 
							 ''.join([new_N_bam[-1], 'i'])])
			subprocess.call(['cp', old_T_bam, new_T_bam])
			subprocess.call(['cp', ''.join([old_T_bam[:-1], 'i']), 
							 ''.join([new_T_bam[-1], 'i'])])
		elif os.path.isdir(old_dir) and os.path.isdir(new_dir):
			redone.append(tokens[2])
		elif not os.path.isdir(old_dir) and os.path.isdir(new_dir):
			second_round.append(tokens[2])
		else:
			leftovers.append(tokens[2])

print str(len(redone)) + ' samples redone'
print str(len(first_round)) + ' samples completed first round'
print str(len(second_round)) + ' samples completed second_round'
print str(len(leftovers)) + ' samples needing work:'
for sample in leftovers:
	print sample
	