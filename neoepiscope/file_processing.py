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
import re
import datetime
from .version import version_number
from intervaltree import Interval, IntervalTree


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
    marktumor = "".join(
        ["""awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9":VT\t"$10":SOMATIC"}' """, 
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
    with open(VCF) as f:
        for line in f:
            # Check header lines to see if FREQ exits in FORMAT fields
            if line[0] == "#":
                if "FREQ" in line:
                    vaf_check = True
            else:
                # Check first entry to get position of FREQ if it exists
                if vaf_check:
                    tokens = line.strip("\n").split("\t")
                    format_field = tokens[8].split(":")
                    for i in range(0, len(format_field)):
                        if format_field[i] == "FREQ":
                            vaf_pos = i
                            break
                # Return None if VCF does not contain VAF data
                else:
                    vaf_pos = None
                    break
    return vaf_pos


def write_results(output_file, hla_alleles, neoepitopes, tool_dict):
    """ Writes predicted neoepitopes out to file

        output_file: path to output file
        hla_alleles: list of HLA alleles used for binding predictions
        neoepitopes: dictionary linking neoepitopes to their metadata
        tool_dict: dictionary storing prediction tool data

        Return value: None.
    """
    try:
        if output_file == "-":
            output_stream = sys.stdout
        else:
            output_stream = open(output_file, "w")
        print(''.join(['# Neoepiscope version ', version_number, '; run ', 
                       str(datetime.date.today())]), 
              file=output_stream)
        headers = [
            "Neoepitope",
            "Chromsome",
            "Pos",
            "Ref",
            "Alt",
            "Mutation_type",
            "VAF",
            "Paired_normal_epitope",
            "Warnings",
            "Transcript_ID",
        ]
        for allele in hla_alleles:
            for tool in sorted(tool_dict.keys()):
                for score_method in sorted(tool_dict[tool][1]):
                    headers.append("_".join([tool, allele, score_method]))
        print("\t".join(headers), file=output_stream)
        for epitope in sorted(neoepitopes.keys()):
            if len(neoepitopes[epitope]) == 1:
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
                ]
                for i in range(9, len(mutation)):
                    out_line.append(str(mutation[i]))
                print("\t".join(out_line), file=output_stream)
            else:
                mutation_dict = collections.defaultdict(list)
                ep_scores = []
                for i in range(9, len(neoepitopes[epitope][0])):
                    ep_scores.append(neoepitopes[epitope][0][i])
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
                for mut in sorted(mutation_dict.keys()):
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
                        ";".join([str(x[1]) for x in mutation_dict[mut]]),
                    ]
                    for score in ep_scores:
                        out_line.append(str(score))
                    print("\t".join(out_line), file=output_stream)
    finally:
        if output_stream is not sys.stdout:
            output_stream.close()

def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers

        md: an MD string (example: 33A^CC).

        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
        if i == 0: continue
        if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
                md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
            md_group.append(char)
    if md_group:
        md_to_parse.append(''.join(md_group))
    return [char for char in md_to_parse if char != '0']

def indels_junctions_exons_mismatches(
            cigar, md, pos, seq, drop_deletions=False, junctions_only=False
        ):
    """ Finds indels, junctions, exons, mismatches from CIGAR, MD string, POS
    
        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence
        drop_deletions: drops deletions from coverage vectors iff True
        junctions_only: does not populate mismatch list

        Return value: tuple
            (insertions, deletions, junctions, exons, mismatches).
        Insertions is a list of tuples (last genomic position before insertion, 
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). Junctions is a list
            of tuples (intron start position (inclusive),
                       intron end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)). Mismatches is a list
            of tuples (genomic position of mismatch, read base)
    """
    insertions, deletions, junctions, exons, mismatches = [], [], [], [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    seq_size = len(seq)
    cigar_chars, cigar_sizes = [], []
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            aligned_base_cap = int(cigar[cigar_index])
            aligned_bases = 0
            while True:
                try:
                    aligned_bases += int(md[md_index])
                    if aligned_bases <= aligned_base_cap:
                        md_index += 1
                except ValueError:
                    # Not an int, but should not have reached a deletion
                    assert md[md_index] != '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
                    if not junctions_only:
                        mismatches.append(
                                (pos + aligned_bases,
                                    seq[seq_index + aligned_bases])
                            )
                    correction_length = len(md[md_index])
                    m_length = aligned_base_cap - aligned_bases
                    if correction_length > m_length:
                        md[md_index] = md[md_index][:m_length]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += correction_length
                        md_index += 1
                if aligned_bases > aligned_base_cap:
                    md[md_index] = aligned_bases - aligned_base_cap
                    break
                elif aligned_bases == aligned_base_cap:
                    break
            # Add exon
            exons.append((pos, pos + aligned_base_cap))
            pos += aligned_base_cap
            seq_index += aligned_base_cap
        elif cigar[cigar_index+1] == 'N':
            skip_increment = int(cigar[cigar_index])
            # Add junction
            junctions.append((pos, pos + skip_increment,
                            seq_index, seq_size - seq_index))
            # Skip region of reference
            pos += skip_increment
        elif cigar[cigar_index+1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append(
                    (pos - 1, seq[seq_index:seq_index+insert_size])
                )
            seq_index += insert_size
        elif cigar[cigar_index+1] == 'D':
            assert md[md_index] == '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
            # Deletion
            delete_size = int(cigar[cigar_index])
            md_delete_size = len(md[md_index+1])
            assert md_delete_size >= delete_size
            deletions.append((pos, md[md_index+1][:delete_size]))
            if not drop_deletions: exons.append((pos, pos + delete_size))
            if md_delete_size > delete_size:
                # Deletion contains a junction
                md[md_index+1] = md[md_index+1][delete_size:]
            else:
                md_index += 2
            # Skip deleted part of reference
            pos += delete_size
        else:
            # Soft clip
            assert cigar[cigar_index+1] == 'S'
            # Advance seq_index
            seq_index += int(cigar[cigar_index])
        cigar_index += 2
    '''Merge exonic chunks/deletions; insertions/junctions could have chopped
    them up.'''
    new_exons = []
    last_exon = exons[0]
    for exon in exons[1:]:
        if exon[0] == last_exon[1]:
            # Merge ECs
            last_exon = (last_exon[0], exon[1])
        else:
            # Push last exon to new exon list
            new_exons.append(last_exon)
            last_exon = exon
    new_exons.append(last_exon)
    return insertions, deletions, junctions, new_exons, mismatches

def rna_support_dict(bam, transcript_to_haplotypes):
    rna_support = defaultdict(list)
    BAM_reader = pysam.AlignmentFile(bam, "rb")
    mutation_dict = {}
    for tx in transcript_to_haplotypes:
        for block in transcript_to_haplotypes[tx]:
            for i in range(0, len(block)):
                if not block[i][6].endswith('*'):
                    phased_muts = [x for x in block if x != block[i] and (x[4] == block[i][4] or x[5] == block[i][5])]
                    mutation_dict[(block[i][0], block[i][1], block[i][2], block[i][3])] = [phased_muts, []]
    contig_list = list(set([x[0] for x in mutation_dict]))
    for contig in contig_list:
        if contig in reference_index.recs.keys():
            for region in contig_dict[contig]:
                x = BAM_reader.fetch(contig=contig, start=region[0], stop=region[1])
                for alignment in x:
                    tokens = str(alignment).split('\t')
                    tags = ast.literal_eval(tokens[11])
                    md = [x[1] for x in tags if x[0] == 'MD'][0]
                    cigar = tokens[5]
                    pos = int(tokens[3])
                    seq = tokens[9]
                    jx = set()
                    insertions, deletions, junctions, exons, mismatches = indels_junctions_exons_mismatches(cigar, md, pos, seq)
                    covered_mutations = set()
                    for j in junctions:
                        jx.add((j[0], j[1]))
                    for ins in insertions:
                        if (contig, ins[0], '', ins[1]) in mutation_dict:
                            covered_mutations.add(ins)
                    for dele in deletions:
                        if (contig, dele[0], dele[1], len(dele[1])) in mutation_dict:
                            covered_mutations.add(deletions)
                    for mm in mismatches:
                        normal_base = reference_index.get_stretch(contig, mm[0]-1, len(mm[1]))
                        if (contig, mm[0], normal_base, mm[1]) in mutation_dict:
                            covered_mutations.add(mm)
                    if len(covered_mutations) > 0:
                        rna_support[tuple(covered_mutations)].append([tokens[0], jx])
    return rna_support
