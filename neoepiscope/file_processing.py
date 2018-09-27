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


def combine_vcf(vcf1, vcf2, outfile="combined.vcf"):
    """ Combines VCFs

        No return value.
    """
    vcffile = open(vcf2, "r")
    temp = open(vcf2 + ".tumortemp", "w+")
    header = open(vcf2 + ".header", "w+")
    for lines in vcffile:
        if lines[0] != "#":
            print(lines.strip(), file=temp)
        else:
            print(lines.strip(), file=header)
    vcffile.close()
    temp.close()
    header.close()
    vcffile = open(vcf1, "r")
    temp = open(vcf2 + ".germlinetemp", "w+")
    for lines in vcffile:
        if lines[0] != "#":
            print(lines.strip(), file=temp)
    vcffile.close()
    temp.close()
    markgermline = "".join(
        ["""awk '{print $0"*"}' """, vcf2, ".germlinetemp > ", vcf2, ".germline"]
    )
    marktumor = "".join(
        ["""awk '{print $0}' """, vcf2, ".tumortemp > ", vcf2, ".tumor"]
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


def prep_hapcut_output(output, hapcut2_output, vcf):
    """ Adds unphased mutations to HapCUT2 output as their own haplotypes

        output: path to output file to write adjusted haplotypes
        hapcut2_output: path to original output from HapCUT2 with only
            phased mutations, or None if using unphased mutations
        vcf: path to vcf used to generate original HapCUT2 output

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
                                            tokens[7],
                                            tokens[8],
                                            tokens[9],
                                            tokens[10],
                                        ]
                                    ),
                                    file=output_stream,
                                )
                        else:
                            phased[(tokens[3], int(tokens[4]))].add((tokens[5], tokens[6]))
                            print(line.strip(), file=output_stream)
                    else:
                        print(line.strip(), file=output_stream)
            print("********", file=output_stream)
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
                tokens = line.split("\t")
                pos = int(tokens[1])
                alt_alleles = tokens[4].split(",")
                for allele in alt_alleles:
                    if (tokens[3], allele) not in phased[(tokens[0], pos)]:
                        print("BLOCK: unphased", file=output_stream)
                        print(
                            (
                                "{vcf_line}\t1\t0\t{chrom}\t"
                                "{pos}\t{ref}\t{alt}\t"
                                "{genotype}\tNA\tNA"
                            ).format(
                                vcf_line=counter,
                                chrom=tokens[0],
                                pos=pos,
                                ref=tokens[3],
                                alt=allele,
                                genotype=tokens[9],
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
