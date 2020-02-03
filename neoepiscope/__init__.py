#!/usr/bin/env python
# coding=utf-8
"""
neoepiscope

Identifies neoepitopes from DNA-seq, VCF, GTF, and Bowtie index.

The MIT License (MIT)
Copyright (c) 2018 Mary A. Wood, Austin Nguyen,
                   Abhinav Nellore, and Reid F. Thompson

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
import argparse
from . import bowtie_index
import sys
import string
import copy
import pickle
import copy
import random
import re
import os
import collections
import tempfile
import subprocess
import warnings
from . import paths
from .transcript import (
    Transcript,
    gtf_to_cds,
    cds_to_tree,
    get_transcripts_from_tree,
    process_haplotypes,
    get_peptides_from_transcripts,
)
from .binding_scores import get_binding_tools, gather_binding_scores
from .file_processing import (
    adjust_tumor_column,
    combine_vcf,
    prep_hapcut_output,
    which,
    get_vaf_pos,
    write_results,
)
from operator import itemgetter
from intervaltree import Interval, IntervalTree

_help_intro = (
    """neoepiscope searches for neoepitopes using tumor/normal DNA-seq data."""
)


def help_formatter(prog):
    """ So formatter_class's max_help_position can be changed. """
    return argparse.HelpFormatter(prog, max_help_position=40)


def main():
    """ Entry point for neoepiscope software """
    parser = argparse.ArgumentParser(
        description=_help_intro, formatter_class=help_formatter
    )
    subparsers = parser.add_subparsers(
        help=(
            'subcommands; add "-h" or "--help" ' "after a subcommand for its parameters"
        ),
        dest="subparser_name",
    )
    index_parser = subparsers.add_parser(
        "index",
        help=(
            "produces pickled dictionaries "
            "linking transcripts to intervals and "
            " CDS lines in a GTF"
        ),
    )
    swap_parser = subparsers.add_parser(
        "swap",
        help=(
            "swaps tumor and normal columns "
            "in a somatic vcf if necessary for "
            "proper HapCUT2 results"
        ),
    )
    merge_parser = subparsers.add_parser(
        "merge",
        help=(
            "merges germline and somatic "
            "VCFS for combined mutation "
            "phasing with HAPCUT2"
        ),
    )
    download_parser = subparsers.add_parser("download", help="downloads dependencies")
    prep_parser = subparsers.add_parser(
        "prep", help=("combines HAPCUT2 output with unphased variants for call mode")
    )
    call_parser = subparsers.add_parser("call", help="calls neoepitopes")
    # Index parser options (produces pickled dictionaries for transcript data)
    index_parser.add_argument(
        "-g", "--gtf", type=str, required=True, help="input path to GTF file"
    )
    index_parser.add_argument(
        "-d",
        "--dicts",
        type=str,
        required=True,
        help="output path to pickled CDS dictionary directory",
    )
    # Swap parser options (swaps columns in somatic VCF)
    swap_parser.add_argument(
        "-i", "--input", type=str, required=True, help="input path to somatic VCF"
    )
    swap_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="-",
        help="output path to column-swapped VCF; use - for stdout",
    )
    # Merger parser options (merges somatic and germline VCFs)
    merge_parser.add_argument(
        "-g", "--germline", type=str, required=True, help="input path to germline VCF"
    )
    merge_parser.add_argument(
        "-s", "--somatic", type=str, required=True, help="input path to somatic VCF"
    )
    merge_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="-",
        help="output path to combined VCF; use - for stdout",
    )
    merge_parser.add_argument(
        "-t",
        "--tumor-id",
        type=str,
        required=False,
        default="TUMOR",
        help="tumor ID (matching the sample in your tumor BAM file "
             "if using GATK ReadBackedPhasing)",
    )
    # Prep parser options (adds unphased mutations as their own haplotype)
    prep_parser.add_argument("-v", "--vcf", type=str, required=True, help="input VCF")
    prep_parser.add_argument(
        "-c",
        "--hapcut2-output",
        type=str,
        required=False,
        help="path to output file of HAPCUT2 run on input VCF",
    )
    prep_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="-",
        help="path to output file to be input to call mode; use - for stdout",
    )
    prep_parser.add_argument(
        "-p",
        "--phased",
        required=False,
        action="store_true",
        help="indicates that input VCF is phased using GATK ReadBackedPhasing",
    )
    # Call parser options (calls neoepitopes)
    call_parser.add_argument(
        "-x",
        "--bowtie-index",
        type=str,
        required=False,
        help="path to Bowtie index basename",
    )
    call_parser.add_argument(
        "-v", "--vcf", type=str, required=False, help="input path to somatic VCF"
    )
    call_parser.add_argument(
        "-d",
        "--dicts",
        type=str,
        required=False,
        help="input path to pickled CDS dictionary directory",
    )
    call_parser.add_argument(
        "-c",
        "--merged-hapcut2-output",
        type=str,
        required=False,
        default="-",
        help="path to output of prep subcommand; use - for stdin",
    )
    call_parser.add_argument(
        "-k",
        "--kmer-size",
        type=str,
        required=False,
        default="8,11",
        help="kmer size for epitope calculation",
    )
    call_parser.add_argument(
        "-p",
        "--affinity-predictor",
        type=str,
        nargs=3,
        required=False,
        action="append",
        default=[["mhcflurry", "1", "affinity,rank"]],
        help="binding affinity prediction software,"
        "associated version number, and scoring method(s) "
        "(e.g. -p netMHCpan 4 rank,affinity); "
        "for multiple programs, repeat the argument; "
        "see documentation for details",
    )
    call_parser.add_argument(
        "-n",
        "--no-affinity",
        required=False,
        action="store_true",
        help="do not run binding affinity predictions; overrides any "
        "binding affinity prediction tools specified via "
        "--affinity-predictor option",
    )
    call_parser.add_argument(
        "-a",
        "--alleles",
        type=str,
        required=False,
        help="comma separated list of alleles; "
        "see documentation online for more information",
    )
    call_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="-",
        help="path to output file; use - for stdout",
    )
    call_parser.add_argument(
        "-f",
        "--fasta",
        required=False,
        action="store_true",
        help="produce additional fasta output; see documentation",
    )
    call_parser.add_argument(
        "-u",
        "--upstream-atgs",
        type=str,
        required=False,
        default="none",
        help="how to handle upstream start codons, see "
        "documentation online for more information",
    )
    call_parser.add_argument(
        "-g",
        "--germline",
        type=str,
        required=False,
        default="background",
        help="how to handle germline mutations in "
        "neoepitope enumeration; documentation online for more information",
    )
    call_parser.add_argument(
        "-s",
        "--somatic",
        type=str,
        required=False,
        default="include",
        help="how to handle somatic mutations in "
        "neoepitope enumeration; documentation online for more information",
    )
    call_parser.add_argument(
        "-b",
        "--build",
        type=str,
        required=False,
        help="which default genome build to use (human hg19 or GRCh38, or mouse mm9); "
        "must have used download.py script to install these",
    )
    call_parser.add_argument(
        "-i",
        "--isolate",
        required=False,
        action="store_true",
        help="isolate mutations - do not use phasing information to "
        "combine nearby mutations in the same neoepitope",
    )
    call_parser.add_argument(
        "-r",
        "--rna-bam",
        type=str,
        required=False,
        help="path to tumor RNA-seq BAM alignment file")
    call_parser.add_argument(
        "--nmd",
        required=False,
        action="store_true",
        default=False,
        help="enumerate neoepitopes from nonsense mediated decay transcripts",
    )
    call_parser.add_argument(
        "--pp",
        required=False,
        action="store_true",
        default=False,
        help="enumerate neoepitopes from polymorphic pseudogene transcripts",
    )
    call_parser.add_argument(
        "--igv",
        required=False,
        action="store_true",
        default=False,
        help="enumerate neoepitopes IGV transcripts",
    )
    call_parser.add_argument(
        "--trv",
        required=False,
        action="store_true",
        default=False,
        help="enumerate neoepitopes from TRV transcripts",
    )
    call_parser.add_argument(
        "--allow-nonstart",
        required=False,
        action="store_true",
        default=False,
        help="enumerate neoepitopes from transcripts without annotated start codons",
    )
    call_parser.add_argument(
        "--allow-nonstop",
        required=False,
        action="store_true",
        default=False,
        help="enumerate neoepitopes from transcripts without annotated stop codons",
    )
    args = parser.parse_args()
    if args.subparser_name == "download":
        from .download import NeoepiscopeDownloader

        downloader = NeoepiscopeDownloader()
        downloader.run()
    elif args.subparser_name == "index":
        cds_dict, tx_data_dict = gtf_to_cds(args.gtf, args.dicts)
        tree = cds_to_tree(cds_dict, args.dicts)
    elif args.subparser_name == "swap":
        adjust_tumor_column(args.input, args.output)
    elif args.subparser_name == "merge":
        combine_vcf(args.germline, args.somatic, outfile=args.output, 
                    tumor_id=args.tumor_id)
    elif args.subparser_name == "prep":
        prep_hapcut_output(args.output, args.hapcut2_output, args.vcf, args.phased)
    elif args.subparser_name == "call":
        # Check that output options are compatible
        if args.fasta and args.output == "-":
            sys.exit(
                "Cannot write fasta results when writing output to standard out; "
                "please specify an output file using the -o/--output option when "
                "using the -f/--fasta flag"
            )
        # Load pickled dictionaries and prepare bowtie index
        if args.build is not None:
            if (
                args.build == "GRCh38"
                and paths.gencode_v29 is not None
                and paths.bowtie_grch38 is not None
            ):
                with open(
                    os.path.join(paths.gencode_v29, "intervals_to_transcript.pickle"),
                    "rb",
                ) as interval_stream:
                    interval_dict = pickle.load(interval_stream)
                with open(
                    os.path.join(paths.gencode_v29, "transcript_to_CDS.pickle"), "rb"
                ) as cds_stream:
                    cds_dict = pickle.load(cds_stream)
                with open(
                    os.path.join(paths.gencode_v29, "transcript_to_gene_info.pickle"), "rb"
                ) as info_stream:
                    info_dict = pickle.load(info_stream)
                reference_index = bowtie_index.BowtieIndexReference(paths.bowtie_grch38)
            elif (
                args.build == "hg19"
                and paths.gencode_v19 is not None
                and paths.bowtie_hg19 is not None
            ):
                with open(
                    os.path.join(paths.gencode_v19, "intervals_to_transcript.pickle"),
                    "rb",
                ) as interval_stream:
                    interval_dict = pickle.load(interval_stream)
                with open(
                    os.path.join(paths.gencode_v19, "transcript_to_CDS.pickle"), "rb"
                ) as cds_stream:
                    cds_dict = pickle.load(cds_stream)
                with open(
                    os.path.join(paths.gencode_v19, "transcript_to_gene_info.pickle"), "rb"
                ) as info_stream:
                    info_dict = pickle.load(info_stream)
                reference_index = bowtie_index.BowtieIndexReference(paths.bowtie_hg19)
            elif (
                args.build == "mm9"
                and paths.gencode_vM1 is not None
                and paths.bowtie_mm9 is not None
            ):
                with open(
                    os.path.join(paths.gencode_vM1, "intervals_to_transcript.pickle"),
                    "rb",
                ) as interval_stream:
                    interval_dict = pickle.load(interval_stream)
                with open(
                    os.path.join(paths.gencode_vM1, "transcript_to_CDS.pickle"), "rb"
                ) as cds_stream:
                    cds_dict = pickle.load(cds_stream)
                with open(
                    os.path.join(paths.gencode_vM1, "transcript_to_gene_info.pickle"), "rb"
                ) as info_stream:
                    info_dict = pickle.load(info_stream)
                reference_index = bowtie_index.BowtieIndexReference(paths.bowtie_mm9)
            else:
                raise RuntimeError(
                    "".join(
                        [
                            args.build,
                            " is not an available genome build. Please "
                            "check that you have run neoepiscope download and are "
                            "using 'hg19', 'GRCh38', or 'mm9' for this argument.",
                        ]
                    )
                )
        else:
            if args.bowtie_index is not None and args.dicts is not None:
                intervals_path = os.path.join(
                    args.dicts, "intervals_to_transcript.pickle"
                )
                if os.path.isfile(intervals_path):
                    with open(intervals_path, "rb") as interval_stream:
                        interval_dict = pickle.load(interval_stream)
                else:
                    raise RuntimeError(
                        "".join(
                            [
                                "Cannot find ",
                                intervals_path,
                                "; have you indexed your GTF with neoepiscope index?",
                            ]
                        )
                    )
                cds_path = os.path.join(args.dicts, "transcript_to_CDS.pickle")
                if os.path.isfile(cds_path):
                    with open(cds_path, "rb") as cds_stream:
                        cds_dict = pickle.load(cds_stream)
                else:
                    raise RuntimeError(
                        "".join(
                            [
                                "Cannot find ",
                                cds_path,
                                "; have you indexed your GTF with neoepiscope index?",
                            ]
                        )
                    )
                info_path = os.path.join(args.dicts, "transcript_to_gene_info.pickle")
                if os.path.isfile(info_path):
                    with open(info_path, "rb") as info_stream:
                        info_dict = pickle.load(cds_stream)
                else:
                    raise RuntimeError(
                        "".join(
                            [
                                "Cannot find ",
                                info_path,
                                "; have you indexed your GTF with neoepiscope index?",
                            ]
                        )
                    )
                bowtie_files = [
                    "".join([args.bowtie_index, ".", str(x), ".ebwt"]) for x in range(1, 5)
                ]
                if list(set([os.path.isfile(x) for x in bowtie_files])) == [True]:
                    reference_index = bowtie_index.BowtieIndexReference(
                        args.bowtie_index
                    )
                else:
                    raise RuntimeError("Cannot find specified bowtie index")
            else:
                raise RuntimeError(
                    "User must specify either --build OR "
                    "--bowtie_index and --dicts options"
                )
        # Check affinity predictor(s)
        if args.no_affinity:
            args.affinity_predictor = None
            tool_dict = {}
        if args.affinity_predictor is not None:
            tool_dict = get_binding_tools(args.affinity_predictor)
        if not tool_dict:
            warnings.warn(
                "No binding prediction tools specified; "
                "will proceed without binding predictions",
                Warning,
            )
            hla_alleles = []
        else:
            if args.alleles:
                hla_alleles = sorted(args.alleles.split(","))
            else:
                raise RuntimeError(
                    "To perform binding affinity predictions, "
                    "user must specify at least one allele "
                    "via the --alleles option"
                )
        # Obtain peptide sizes for kmerizing peptides
        if "," in args.kmer_size:
            size_list = args.kmer_size.split(",")
            size_list.sort(reverse=True)
            for i in range(0, len(size_list)):
                size_list[i] = int(size_list[i])
        elif "-" in args.kmer_size:
            size_list = args.kmer_size.split("-")
            size_list.sort(reverse=True)
            for i in range(0, len(size_list)):
                size_list[i] = int(size_list[i])
        else:
            size_list = [int(args.kmer_size)]
        # Establish handling of ATGs
        if args.upstream_atgs == "none":
            only_novel_upstream = False
            only_downstream = True
            only_reference = False
        elif args.upstream_atgs == "novel":
            only_novel_upstream = True
            only_downstream = False
            only_reference = False
        elif args.upstream_atgs == "all":
            only_novel_upstream = False
            only_downstream = False
            only_reference = False
        elif args.upstream_atgs == "reference":
            only_novel_upstream = False
            only_downstream = False
            only_reference = True
        else:
            raise RuntimeError(
                "--upstream-atgs must be one of "
                '{"novel", "all", "none", "reference"}'
            )
        # Establish VAF position as None (can be modified)
        vaf_pos = None
        # Establish handling of germline mutations:
        if args.germline == "background":
            include_germline = 2
        elif args.germline == "include":
            include_germline = 1
        elif args.germline == "exclude":
            include_germline = 0
        else:
            raise RuntimeError(
                "--germline must be one of " '{"background", "include", "exclude"}'
            )
        # Establish handling of somatic mutations:
        if args.somatic == "include":
            include_somatic = 1
            # If VCF is given, search for VAF position
            if args.vcf:
                vaf_pos = get_vaf_pos(args.vcf)
        elif args.somatic == "background":
            include_somatic = 2
        elif args.somatic == "exclude":
            include_somatic = 0
        else:
            raise RuntimeError(
                "--somatic must be one of " '{"background", "include", "exclude"}'
            )
        # Warn if somatic and germline are both excluded
        if include_somatic == 0 and include_germline == 0:
            warnings.warn(
                "Germline and somatic mutations are both being "
                "excluded, no epitopes will be returned",
                Warning,
            )
        if not args.isolate:
            phase_mutations = True
        else:
            phase_mutations = False
        # Find transcripts that haplotypes overlap
        relevant_transcripts, homozygous_variants = process_haplotypes(
            args.merged_hapcut2_output, interval_dict, phase_mutations
        )
        # Apply mutations to transcripts and get neoepitopes
        neoepitopes, fasta = get_peptides_from_transcripts(
            relevant_transcripts,
            homozygous_variants,
            vaf_pos,
            cds_dict,
            only_novel_upstream,
            only_downstream,
            only_reference,
            reference_index,
            size_list,
            args.nmd,
            args.pp,
            args.igv,
            args.trv,
            args.allow_nonstart,
            args.allow_nonstop,
            include_germline,
            include_somatic,
            protein_fasta=args.fasta,
        )
        # If neoepitopes are found, get binding scores and write results
        if len(neoepitopes) > 0:
            full_neoepitopes = gather_binding_scores(
                neoepitopes, tool_dict, hla_alleles, size_list
            )
            write_results(args.output, hla_alleles, full_neoepitopes, tool_dict, info_dict)
            if args.fasta:
                fasta_file = "".join([args.output, ".fasta"])
                with open(fasta_file, "w") as f:
                    for tx in fasta:
                        proteins = sorted(list(fasta[tx]))
                        for i in range(0, len(proteins)):
                            identifier = "".join([">", tx, "_v", str(i)])
                            print(identifier, file=f)
                            print(proteins[i], file=f)
        else:
            print("No neoepitopes found", file=sys.stderr)
    else:
        parser.print_usage()


if __name__ == "__main__":
    main()
