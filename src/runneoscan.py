'''
Austin is here.

'Main' Algorithm:

add argument --bam (tumor bam)
add argument --gbam (germline bam)
add argument --gvcf (germline vcf)

if args.bam == "-":
    call neoscan with args = parser.parse_args

if args.gbam == "-":
    call runhapcut2(args.vcf, args.bam)
    call neoscan with args and hapcut flag

if args.gbam != "-" and args.gvcf != "-"
    call runhapcut2(args.gvcf, args.gbam)
    call neoscan with args and hapcut flag

    call combinevcf(args.gvcf, args.vcf)
    call runhapcut2(cut.vcf, args.bam)
    call neoscan with args and hapcut flag

'''

import argparse
import os
import sys
from runhapcut2 import runhapcut2 as runhap
from combinevcf import combinevcf

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf', type=str, required=True,
        default='-',
        help='input vcf or "-" for stdin'
    )
parser.add_argument('-x', '--bowtie-index', type=str, required=True,
        help='path to Bowtie index basename'
    )
parser.add_argument('-d', '--dicts', type=str, required=True,
        help='input path to pickled dictionaries'
    )
parser.add_argument('-b', '--bam', type=str, required=False,
        default = '-', help='T/F bam is used'
    )
parser.add_argument('-gb', '--gbam', type=str, required=False,
        default = '-', help='Germline bam is used'
    )
parser.add_argument('-gv', '--gvcf', type=str, required=False,
        default = '-', help='Germline vcf is used'
    )

args = parser.parse_args()
if args.bam != '-' and args.gbam != '-' and args.gvcf != '-':
    continue
    #runhap.runhapcut2(args.gvcf, args.gbam)
    #command = "python neoscan.py --vcf " + args.gvcf + "--bowtie-index bowtie.bowtie --dicts dicts.dicts --bam " + args.gbam
    #os.system(command)
    #run mary's script and save to new output as germline vs ref...

    #combinevcf.combinevcf(args.vcf, args.gvcf)
    #runhap.runhapcut2(cutvcf.vcf, args.bam)
    #command = "python neoscan.py --vcf cut.vcf --bowtie-index bowtie.bowtie --dicts dicts.dicts --bam " + args.bam
    #os.system(command)
    #run mary's script and save to new output as tumor vs ref...

if args.bam != '-' and args.gbam == '-' and args.gvcf != '-':
    continue
    #combinevcf.combinevcf(args.vcf, args.gvcf)
    #runhap.runhapcut2(cutvcf, args.bam)
    #command = "python neoscan.py --vcf cut.vcf --bowtie-index bowtie.bowtie --dicts dicts.dicts --bam " + args.bam
    #os.system(command)
    #run mary's script and save to new output as tumor vs ref...

if (args.bam != '-' and args.gbam == '-' and args.gvcf == '-') or (args.bam != '-' and args.gbam != '-' and args.gvcf == '-'):
    continue
    #runhap.runhapcut2(args.gvcf, args.gbam)
    #command = "python neoscan.py --vcf " + args.gvcf + "--bowtie-index bowtie.bowtie --dicts dicts.dicts --bam " + args.gbam
    #os.system(command)
    #run mary's script and save to new output as tumor vs ref...

if args.bam == '-':
    continue
    #command = "python neoscan.py --vcf " + args.gvcf + "--bowtie-index bowtie.bowtie --dicts dicts.dicts"
    #os.system(command)
    #run mary's script and save to new output as tumor vs ref without phasing...



