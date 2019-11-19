#!/usr/bin/env python
# coding=utf-8
"""
file_processing.py

Part of neoepiscope
Includes functions for processing input and output files.

Licensed under the MIT license.

The MIT License (MIT)
Copyright (c) 2018 Mary A. Wood, Austin Nguyen,
                   Abhinav Nellore, and Reid Thompson

Functions parsed_md() and indels_junctions_exons_mismatches()
    are from Rail-RNA, which is copyright (c) 2015 
                    Abhinav Nellore, Leonardo Collado-Torres,
                    Andrew Jaffe, James Morton, Jacob Pritt,
                    José Alquicira-Hernández,
                    Christopher Wilks,
                    Jeffrey T. Leek, and Ben Langmead.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from __future__ import absolute_import, division, print_function
import subprocess
import warnings
import collections
import sys
import os
import re
import pickle
import datetime
from .version import version_number
from intervaltree import Interval, IntervalTree

neoepiscope_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# From https://stackoverflow.com/questions/30212413/backport-python-3-4s-regular-expression-fullmatch-to-python-2
def fullmatch(regex, string, flags=0):
    """Emulate python-3.4 re.fullmatch()."""
    return re.match("(?:" + regex + r")\Z", string, flags=flags)

def adjust_tumor_column(in_vcf, out_vcf):
    """ Swaps the sample columns in a somatic vcf

        HAPCUT2 only takes data from the first VCF sample column, so if the
            tumor sample data is in the second VCF sample column, it must be
            swapped prior to optional germline merging or running HAPCUT2

        in_vcf: input vcf that needs the tumor sample data flipped
        out_vcf: output vcf to have the correct columns

        No return value.
    """
    header_lines = []
    other_lines = []
    # Process input vcf
    with open(in_vcf, "r") as f:
        for line in f:
            # Preserve header lines with out change
            if line[0:2] == "##":
                header_lines.append(line.strip("\n"))
            # Adjust column header and variant lines
            else:
                tokens = line.strip("\n").split("\t")
                if line[0] == "#":
                    warnings.warn(
                        "".join(
                            [
                                "Reading ",
                                tokens[9],
                                " as normal tissue and ",
                                tokens[10],
                                " as tumor tissue",
                            ]
                        ),
                        Warning,
                    )
                new_line = "\t".join(
                    [
                        tokens[0],
                        tokens[1],
                        tokens[2],
                        tokens[3],
                        tokens[4],
                        tokens[5],
                        tokens[6],
                        tokens[7],
                        tokens[8],
                        tokens[10],
                        tokens[9],
                    ]
                )
                other_lines.append(new_line)
    # Write new vcf
    try:
        if out_vcf == "-":
            output_stream = sys.stdout
        else:
            output_stream = open(out_vcf, "w")
        for line in header_lines:
            print(line, file=output_stream)
        for line in other_lines:
            print(line, file=output_stream)
    finally:
        if output_stream is not sys.stdout:
            output_stream.close()


def combine_vcf(vcf1, vcf2, outfile="combined.vcf", tumor_id="TUMOR"):
    """ Combines VCFs

        vcf1: path to germline VCF file
        vcf2: path to tumor VCF file
        outfile: path to write merged VCF file
        tumor_id: identifier of tumor sample listed in VCF header

        No return value.
    """
    vcffile = open(vcf2, "r")
    temp = open(vcf2 + ".tumortemp", "w+")
    header = open(vcf2 + ".header", "w+")
    info_lines = set()
    for lines in vcffile:
        if lines[0] != "#":
            print(lines.strip(), file=temp)
        elif lines[0:2] == "##":
            if "INFO=<" in lines:
                info_lines.add(lines.strip())
            else:
                print(lines.strip(), file=header)
        else:
            tokens = lines.strip().split('\t')
            if len(tokens) == 10:
                tumor_first = True
                warnings.warn(''.join(['Only 1 sample in somatic VCF; '
                                       'treating ', tokens[9], ' column as ',
                                       'the tumor sample']))
            elif len(tokens) == 11:
                if tokens[9] == "TUMOR" and tokens[10] == "NORMAL":
                    tumor_first = True
                elif tokens[10] == "TUMOR" and tokens[9] == "NORMAL":
                    tumor_first = False
                elif tokens[9] == "PRIMARY" and tokens[10] == "NORMAL":
                    tumor_first = True
                elif tokens[10] == "PRIMARY" and tokens[9] == "NORMAL":
                    tumor_first = False
                elif tokens[9] == tumor_id:
                    tumor_first = True
                elif tokens[10] == tumor_id:
                    tumor_first = False
                elif tokens[9] == 'TUMOR':
                    warnings.warn(''.join(['Irregular sample identifiers; '
                                       'treating ', tokens[9], 'column as ',
                                       'the tumor sample']))
                    tumor_first = True
                elif tokens[10] == 'TUMOR':
                    warnings.warn(''.join(['Irregular sample identifiers; '
                                       'treating ', tokens[10], 'column as ',
                                       'the tumor sample']))
                    tumor_first = False
                elif tokens[9] == 'NORMAL':
                    warnings.warn(''.join(['Irregular sample identifiers; '
                                       'treating ', tokens[10], 'column as ',
                                       'the tumor sample']))
                    tumor_first = False
                elif tokens[10] == 'NORMAL':
                    warnings.warn(''.join(['Irregular sample identifiers; '
                                       'treating ', tokens[9], 'column as ',
                                       'the tumor sample']))
                    tumor_first = True
                else:
                    raise RuntimeError(
                                        ''.join(["Can't identify tumor sample",
                                                 " in somatic VCF; please ",
                                                 "provide tumor identifier ",
                                                 "using -t ", tokens[9],
                                                 " or -t ", tokens[10]])
                    )
            elif len(tokens) > 11:
                raise RuntimeError(
                                    "Somatic VCF contains more than two "
                                    "samples, please use a VCF that contains "
                                    "only 1 tumor and 1 normal sample or only "
                                    "1 tumor sample."
                )
            else:
                raise RuntimeError("Somatic VCF is missing sample data.")
    vcffile.close()
    temp.close()
    vcffile = open(vcf1, "r")
    temp = open(vcf2 + ".germlinetemp", "w+")
    for lines in vcffile:
        if lines[0] != "#":
            print(lines.strip(), file=temp)
        elif lines[0:2] == "##":
            if "INFO=<" in lines:
                info_lines.add(lines.strip())
            else:
                print(lines.strip(), file=header)
        else:
            tokens = lines.strip().split('\t')
            if len(tokens) > 10:
                raise RuntimeError(
                                    "Germline VCF contains more than one "
                                    "sample, please use a VCF that contains "
                                    "only  1 normal sample."
                )
            elif len(tokens) < 10:
                raise RuntimeError("Germline VCF is missing sample data.")
    vcffile.close()
    temp.close()
    for line in sorted(list(info_lines)):
        print(line, file=header)
    print('##FORMAT=<ID=VT,Number=1,Type=String,Description="Variant type, SOMATIC or GERMLINE">',
            file=header)
    print('\t'.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                     "FILTER", "INFO", "FORMAT", tumor_id]), 
                    file=header
    )
    header.close()
    markgermline = "".join(
        ["""awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9":VT\t"$10":GERMLINE"}' """, 
         vcf2, ".germlinetemp > ", vcf2, ".germline"]
    )
    if tumor_first:
        marktumor = "".join(
            ["""awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9":VT\t"$10":SOMATIC"}' """, 
             vcf2, ".tumortemp > ", vcf2, ".tumor"]
        )
    else:
        marktumor = "".join(
            ["""awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9":VT\t"$11":SOMATIC"}' """, 
             vcf2, ".tumortemp > ", vcf2, ".tumor"]
        )
    subprocess.call(markgermline, shell=True)
    subprocess.call(marktumor, shell=True)
    command = "".join(
        ["cat ", vcf2, ".germline ", vcf2, ".tumor > ", vcf2, ".combine1"]
    )
    subprocess.call(command, shell=True)
    command2 = "".join(["sort -k1,1 -k2,2n ", vcf2, ".combine1 > ", vcf2, ".sorted"])
    subprocess.call(command2, shell=True)
    command3 = "".join(
        ["cat ", vcf2, ".header ", vcf2, ".sorted > ", vcf2, ".combine2"]
    )
    subprocess.call(command3, shell=True)
    cut = "".join(["cut -f1,2,3,4,5,6,7,8,9,10 ", vcf2, ".combine2 > ", outfile])
    subprocess.call(cut, shell=True)
    for file in [
        ".tumortemp",
        ".germlinetemp",
        ".combine1",
        ".combine2",
        ".sorted",
        ".tumor",
        ".germline",
        ".header",
    ]:
        cleanup = "".join(["rm ", vcf2, file])
        subprocess.call(cleanup, shell=True)

def prep_hapcut_output(output, hapcut2_output, vcf, phased_vcf=False):
    """ Adds unphased mutations to HapCUT2 output as their own haplotypes

        output: path to output file to write adjusted haplotypes
        hapcut2_output: path to original output from HapCUT2 with only
            phased mutations, or None if using unphased mutations
        vcf: path to vcf used to generate original HapCUT2 output
        phased: vcf file is a phased vcf from GATK ReadBackedPhasing

        Return value: None
    """
    phased = collections.defaultdict(set)
    try:
        if output == "-":
            output_stream = sys.stdout
        else:
            output_stream = open(output, "w")
        if hapcut2_output is not None:
            with open(hapcut2_output) as hapcut2_stream:
                for line in hapcut2_stream:
                    if line[0] != "*" and not line.startswith("BLOCK"):
                        tokens = line.strip().split("\t")
                        if ':GERMLINE' in tokens[7]:
                            gen_end = '*'
                        else:
                            gen_end = ''
                        if tokens[1] == "-" or tokens[2] == "-":
                            continue
                        elif "," in tokens[6]:
                            alt_alleles = tokens[6].split(",")
                            try:
                                assert len(alt_alleles) == 2
                            except AssertionError:
                                warnings.warn(
                                    "".join(
                                        [
                                            "Neoepiscope does not support ",
                                            "triallellic phasing; of ",
                                            "alternate alleles ",
                                            tokens[6],
                                            " at contig ",
                                            tokens[3],
                                            " position ",
                                            tokens[4],
                                            ", only ",
                                            "the top two will be included.",
                                        ]
                                    )
                                )
                            for i in range(0, 2):
                                allele = alt_alleles[i]
                                phased[(tokens[3], int(tokens[4]))].add((tokens[5], allele))
                                if i == 0:
                                    if tokens[1] == "1":
                                        gen1 = "1"
                                        gen2 = "0"
                                    else:
                                        gen1 = "0"
                                        gen2 = "1"
                                elif i == 1:
                                    if tokens[1] == "2":
                                        gen1 = "1"
                                        gen2 = "0"
                                    else:
                                        gen1 = "0"
                                        gen2 = "1"
                                print(
                                    "\t".join(
                                        [
                                            tokens[0],
                                            gen1,
                                            gen2,
                                            tokens[3],
                                            tokens[4],
                                            tokens[5],
                                            alt_alleles[i],
                                            ''.join([tokens[7], gen_end]),
                                            tokens[8],
                                            tokens[9],
                                            tokens[10],
                                        ]
                                    ),
                                    file=output_stream,
                                )
                        else:
                            phased[(tokens[3], int(tokens[4]))].add((tokens[5], tokens[6]))
                            print(
                                    "\t".join(
                                        [
                                            tokens[0],
                                            tokens[1],
                                            tokens[2],
                                            tokens[3],
                                            tokens[4],
                                            tokens[5],
                                            tokens[6],
                                            ''.join([tokens[7], gen_end]),
                                            tokens[8],
                                            tokens[9],
                                            tokens[10],
                                        ]
                                    ),
                                    file=output_stream,
                                )
                    else:
                        print(line.strip(), file=output_stream)
            print("********", file=output_stream)
        if phased_vcf:
            haplotype_dict = collections.defaultdict(list)
            counter = 0
            with open(vcf) as vcf_stream:
                first_char = "#"
                while first_char == "#":
                    line = vcf_stream.readline().strip()
                    try:
                        first_char = line[0]
                    except IndexError:
                        first_char = "#"
                counter = 1
                while line:
                    tokens = line.strip().split("\t")
                    tokens[9] = tokens[9].strip()
                    if ':GERMLINE' in tokens[9]:
                        gen_end = '*'
                    else:
                        gen_end = ''
                    pos = int(tokens[1])
                    if 'HP' in tokens[8]:
                        hp_index = tokens[8].split(':').index('HP')
                        hap = tuple(tokens[9].split(':')[hp_index].split(','))
                        haplotype_id = (tokens[0], hap[0].split('-')[0])
                        current_haplotype = (int(hap[0].split('-')[1])-1, int(hap[1].split('-')[1])-1)
                        hap_entry = ("{vcf_line}\t{hap1}\t{hap2}\t{chrom}\t"
                                    "{pos}\t{ref}\t{alt}\t"
                                    "{genotype}\tNA\tNA\tNA"
                                    ).format(
                                        vcf_line=counter,
                                        hap1=str(current_haplotype[0]),
                                        hap2=str(current_haplotype[1]),
                                        chrom=tokens[0],
                                        pos=pos,
                                        ref=tokens[3],
                                        alt=tokens[4],
                                        genotype=''.join([tokens[9], gen_end]),
                                    )
                        haplotype_dict[haplotype_id].append(hap_entry)
                    elif '1/1' in tokens[9]:
                        print("BLOCK: unphased", file=output_stream)
                        print(
                            (
                                "{vcf_line}\t1\t1\t{chrom}\t"
                                "{pos}\t{ref}\t{alt}\t"
                                "{genotype}\tNA\tNA\tNA"
                            ).format(
                                vcf_line=counter,
                                chrom=tokens[0],
                                pos=pos,
                                ref=tokens[3],
                                alt=tokens[4],
                                genotype=''.join([tokens[9], gen_end]),
                            ),
                            file=output_stream,
                        )
                        print("********", file=output_stream)
                    else:
                        alt_alleles = tokens[4].split(",")
                        for allele in alt_alleles:
                            print("BLOCK: unphased", file=output_stream)
                            print(
                                (
                                    "{vcf_line}\t1\t0\t{chrom}\t"
                                    "{pos}\t{ref}\t{alt}\t"
                                    "{genotype}\tNA\tNA\tNA"
                                ).format(
                                    vcf_line=counter,
                                    chrom=tokens[0],
                                    pos=pos,
                                    ref=tokens[3],
                                    alt=allele,
                                    genotype=''.join([tokens[9], gen_end]),
                                ),
                                file=output_stream,
                            )
                            print("********", file=output_stream)
                    line = vcf_stream.readline().strip()
                    counter += 1
            for haplotype in haplotype_dict:
                print("BLOCK: phased", file=output_stream)
                for hap_entry in haplotype_dict[haplotype]:
                    print(hap_entry, file=output_stream)
                print("********", file=output_stream)
        else:
            with open(vcf) as vcf_stream:
                first_char = "#"
                while first_char == "#":
                    line = vcf_stream.readline().strip()
                    try:
                        first_char = line[0]
                    except IndexError:
                        first_char = "#"
                counter = 1
                while line:
                    tokens = line.strip().split("\t")
                    pos = int(tokens[1])
                    tokens[9] = tokens[9].strip()
                    if ':GERMLINE' in tokens[9]:
                        gen_end = '*'
                    else:
                        gen_end = ''
                    alt_alleles = tokens[4].split(",")
                    for allele in alt_alleles:
                        if (tokens[3], allele) not in phased[(tokens[0], pos)]:
                            print("BLOCK: unphased", file=output_stream)
                            if "1/1" in tokens[9]:
                                print(
                                    (
                                        "{vcf_line}\t1\t1\t{chrom}\t"
                                        "{pos}\t{ref}\t{alt}\t"
                                        "{genotype}\tNA\tNA\tNA"
                                    ).format(
                                        vcf_line=counter,
                                        chrom=tokens[0],
                                        pos=pos,
                                        ref=tokens[3],
                                        alt=allele,
                                        genotype=''.join([tokens[9], gen_end]),
                                    ),
                                    file=output_stream,
                                )
                            else:
                                print(
                                    (
                                        "{vcf_line}\t1\t0\t{chrom}\t"
                                        "{pos}\t{ref}\t{alt}\t"
                                        "{genotype}\tNA\tNA\tNA"
                                    ).format(
                                        vcf_line=counter,
                                        chrom=tokens[0],
                                        pos=pos,
                                        ref=tokens[3],
                                        alt=allele,
                                        genotype=''.join([tokens[9], gen_end]),
                                    ),
                                    file=output_stream,
                                )
                            print("********", file=output_stream)
                    line = vcf_stream.readline().strip()
                    counter += 1
    finally:
        if output_stream is not sys.stdout:
            output_stream.close()


def which(path):
    """ Searches for whether executable is present and returns version

        path: path to executable

        Return value: None if executable not found, else string with software
            name and version number
    """
    try:
        subprocess.check_call([path])
    except OSError as e:
        return None
    else:
        return path


def get_vaf_pos(VCF):
    """ Obtains position in VCF format/genotype fields of VAF

        VCF: path to input VCF

        Return value: None if VCF does not contain VAF,
                        otherwise position of VAF
    """
    vaf_check = False
    vaf_pos = None
    field = None
    with open(VCF) as f:
        for line in f:
            # Check header lines to see if FREQ exits in FORMAT fields
            if line[0] == "#":
                if "ID=AF" in line:
                    vaf_check = True
                    field = "AF"
                elif "ID=FREQ" in line:
                    vaf_check = True
                    field = "FREQ"
                elif "ID=FA" in line:
                    vaf_check = True
                    field = "FA"
            else:
                # Check first entry to get position of FREQ if it exists
                if vaf_check:
                    tokens = line.strip("\n").split("\t")
                    format_field = tokens[8].split(":")
                    for i in range(0, len(format_field)):
                        if format_field[i] == field:
                            vaf_pos = i
                            break
    return (vaf_pos, field)


def write_results(output_file, hla_alleles, neoepitopes, tool_dict, tx_dict):
    """ Writes predicted neoepitopes out to file

        output_file: path to output file
        hla_alleles: list of HLA alleles used for binding predictions
        neoepitopes: dictionary linking neoepitopes to their metadata
        tool_dict: dictionary storing prediction tool data
        tx_dict: dictionary linking transcript ID to list of 
                    [transcript type, gene ID, gene name]

        Return value: None.
    """
    # Load epitope to IEDB linker dicts
    with open(
            os.path.join(
                os.path.join(neoepiscope_dir, "neoepiscope", "epitopeID.pickle")
            ),
            "rb",
        ) as epitope_stream:
            epitope_to_iedb = pickle.load(epitope_stream)
    with open(
            os.path.join(
                os.path.join(neoepiscope_dir, "neoepiscope", "ambiguousEpitopeID.pickle")
            ),
            "rb",
        ) as epitope_stream:
            ambiguous_epitope_to_iedb = pickle.load(epitope_stream)
    try:
        if output_file == "-":
            output_stream = sys.stdout
        else:
            output_stream = open(output_file, "w")
        # Write file header info
        print(''.join(['# Neoepiscope version ', version_number, '; run ', 
                       str(datetime.date.today())]), 
              file=output_stream)
        headers = [
            "Neoepitope",
            "Chromosome",
            "Pos",
            "Ref",
            "Alt",
            "Mutation_type",
            "VAF",
            "Paired_normal_epitope",
            "Warnings",
            "Transcript_ID",
            "Transcript_type",
            "Gene_ID",
            "Gene_name",
            "IEDB_ID"
        ]
        for allele in hla_alleles:
            for tool in sorted(tool_dict.keys()):
                for score_method in sorted(tool_dict[tool][1]):
                    headers.append("_".join([tool, allele, score_method]))
        print("\t".join(headers), file=output_stream)
        # Write output for all epitopes
        for epitope in sorted(neoepitopes.keys()):
            # Find relevant IEDB IDs for epitope
            if epitope in epitope_to_iedb:
                iedb_id = ",".join(list(epitope_to_iedb[epitope]))
            else:
                possible_ids = set()
                for regex in ambiguous_epitope_to_iedb:
                    if fullmatch(regex, epitope) is not None:
                        possible_ids.update(ambiguous_epitope_to_iedb[regex])
                if len(possible_ids) > 0:
                    iedb_id = ",".join(list(possible_ids))
                else:
                    iedb_id = "NA"
            if len(neoepitopes[epitope]) == 1:
                # Epitope only results from 1 transcript - get variant info
                mutation = neoepitopes[epitope][0]
                if mutation[2] == "":
                    ref = "*"
                else:
                    ref = mutation[2]
                if mutation[3] == "":
                    alt = "*"
                else:
                    alt = mutation[3]
                if mutation[5] is None:
                    vaf = "NA"
                else:
                    vaf = str(mutation[5])
                # Get transcript/gene info
                tx_info = tx_dict[mutation[8]]
                out_line = [
                    epitope,
                    mutation[0],
                    str(mutation[1]),
                    ref,
                    alt,
                    mutation[4],
                    vaf,
                    mutation[6],
                    mutation[7],
                    mutation[8],
                    tx_info[0],
                    tx_info[1],
                    tx_info[2],
                    iedb_id
                ]
                for i in range(9, len(mutation)):
                    out_line.append(str(mutation[i]))
                print("\t".join(out_line), file=output_stream)
            else:
                # Epitope results from multiple transcripts
                mutation_dict = collections.defaultdict(list)
                # Get binding score info
                ep_scores = []
                for i in range(9, len(neoepitopes[epitope][0])):
                    ep_scores.append(neoepitopes[epitope][0][i])
                # Get variant info
                for mut in neoepitopes[epitope]:
                    if mut[2] == "":
                        ref = "*"
                    else:
                        ref = mut[2]
                    if mut[3] == "":
                        alt = "*"
                    else:
                        alt = mut[3]
                    if mut[5] is None:
                        vaf = "NA"
                    else:
                        vaf = str(mut[5])
                    mutation_dict[
                        (mut[0], mut[1], ref, alt, mut[4], vaf, mut[6])
                    ].append([mut[7], mut[8]])
                mutation_list = sorted(list(mutation_dict.keys()))
                # Get transcript/gene info 
                for mut in mutation_list:
                    transcripts = [str(x[1]) for x in mutation_dict[mut]]
                    tx_types = []
                    gene_ids = []
                    gene_names = []
                    for tx in transcripts:
                        tx_types.append(tx_dict[tx][0])
                        gene_ids.append(tx_dict[tx][1])
                        gene_names.append(tx_dict[tx][2])
                    out_line = [
                        epitope,
                        mut[0],
                        str(mut[1]),
                        mut[2],
                        mut[3],
                        mut[4],
                        mut[5],
                        mut[6],
                        ";".join([str(x[0]) for x in mutation_dict[mut]]),
                        ";".join(transcripts),
                        ";".join(tx_types),
                        ";".join(gene_ids),
                        ";".join(gene_names),
                        iedb_id
                    ]
                    for score in ep_scores:
                        out_line.append(str(score))
                    print("\t".join(out_line), file=output_stream)
    finally:
        if output_stream is not sys.stdout:
            output_stream.close()
