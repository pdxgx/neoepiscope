#!/usr/bin/env python
"""
neoepiscope

Identifies neoepitopes from DNA-seq, VCF, GTF, and Bowtie index.
"""
import argparse
import bowtie_index
import sys
import string
import copy
import pickle
import copy
import os
import random
import re
import collections
import tempfile
import subprocess
import string
import warnings
#import exe_paths
import exe_paths_mary as exe_paths
from transcript import Transcript, gtf_to_cds, cds_to_tree, get_transcripts_from_tree
from operator import itemgetter
from intervaltree import Interval, IntervalTree
#import Hapcut2interpreter as hap

_revcomp_translation_table = string.maketrans('ATCG', 'TAGC')

# X below denotes a stop codon
_codon_table = {
        'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
        'TAT':'Y', 'TAC':'Y', 'TAA':'X', 'TAG':'X',
        'TGT':'C', 'TGC':'C', 'TGA':'X', 'TGG':'W',
        'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
        'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
        'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
    }
_complement_table = string.maketrans('ATCG', 'TAGC')

_help_intro = '''neoepiscope searches for neoepitopes in seq data.'''
def help_formatter(prog):
    """ So formatter_class's max_help_position can be changed. """
    return argparse.HelpFormatter(prog, max_help_position=40)

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
    with open(in_vcf, 'r') as f:
        for line in f:
            # Preserve header lines with out change
            if line[0:2] == '##':
                header_lines.append(line.strip('\n'))
            # Adjust column header and variant lines
            else:
                tokens = line.strip('\n').split('\t')
                if line[0] == '#':
                    warnings.warn(''.join(['Reading ', tokens[9], 
                                           'as normal tissue and ', tokens[10],
                                           'as tumor tissue']), 
                                  Warning)
                new_line = '\t'.join([tokens[0], tokens[1], tokens[2], 
                                        tokens[3], tokens[4], tokens[5], 
                                        tokens[6], tokens[7], tokens[8], 
                                        tokens[10], tokens[9]])
                other_lines.append(new_line)
    # Write new vcf
    with open(out_vcf, 'w') as f:
        for line in header_lines:
            f.write(line + '\n')
        for line in other_lines:
            f.write(line + '\n')

def combinevcf(vcf1, vcf2, outfile='Combined.vcf'):
    """ Combines VCFs

        No return value.
    """
    vcffile = open(vcf2, 'r')
    temp = open(vcf2 + '.tumortemp', 'w+');
    header = open(vcf2 + '.header', 'w+');
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
        else:
            header.write(lines)
    vcffile.close()
    temp.close()
    header.close()
    vcffile = open(vcf1, 'r')
    temp = open(vcf2 + '.germlinetemp', 'w+');
    for lines in vcffile:
        if (lines[0] != '#'):
            temp.write(lines)
    vcffile.close()
    temp.close()    
    markgermline = ''.join(['''awk '{print $0"*"}' ''', vcf2, 
                            ".germlinetemp > ", vcf2, '.germline'])
    marktumor    = ''.join(['''awk '{print $0}' ''', vcf2, 
                            '.tumortemp > ', vcf2, '.tumor'])
    subprocess.call(markgermline, shell=True)
    subprocess.call(marktumor, shell=True)
    command = ''.join(['cat ', vcf2, '.germline ', vcf2, '.tumor > ', 
                        vcf2, '.combine1'])
    subprocess.call(command, shell=True)
    command2 = ''.join(['sort -k1,1 -k2,2n ', vcf2, '.combine1 > ', 
                        vcf2, '.sorted'])
    subprocess.call(command2, shell=True)
    command3 = ''.join(['cat ', vcf2, '.header ', vcf2, '.sorted > ', 
                        vcf2, '.combine2'])
    subprocess.call(command3, shell=True)
    cut = ''.join(['cut -f1,2,3,4,5,6,7,8,9,10 ', vcf2, 
                    '.combine2 > ', outfile])
    subprocess.call(cut, shell=True)
    for file in ['.tumortemp', '.germlinetemp', '.combine1', '.combine2', 
                    '.sorted', '.tumor', '.germline', '.header']:
        cleanup = ''.join(['rm ', vcf2, file])
        subprocess.call(cleanup, shell=True)

def prep_hapcut_output(output, hapcut2_output, vcf):
    """ Adds unphased mutations to HapCUT2 output as their own haplotypes
        
        output: path to output file to write adjusted haplotypes
        hapcut2_output: path to original output from HapCUT2 with only 
            phased mutations
        vcf: path to vcf used to generate original HapCUT2 output

        Return value: None
    """
    phased = collections.defaultdict(set)
    with open(output, 'w') as output_stream:
        with open(hapcut2_output) as hapcut2_stream:
            for line in hapcut2_stream:
                if line[0] != '*' and not line.startswith('BLOCK'):
                    tokens = line.strip().split('\t')
                    phased[(tokens[3], int(tokens[4]))].add(
                                                    (tokens[5], tokens[6])
                                                )
                print >>output_stream, line,
        print >>output_stream, '********'
        with open(vcf) as vcf_stream:
            first_char = '#'
            while first_char == '#':
                line = vcf_stream.readline().strip()
                try:
                    first_char = line[0]
                except IndexError:
                    first_char = '#'
            counter = 1
            while line:
                tokens = line.split('\t')
                pos = int(tokens[1])
                alt_alleles = tokens[4].split(',')
                for allele in alt_alleles:
                    if (tokens[3], allele) not in phased[
                                                (tokens[0], pos)
                                            ]:
                        print >>output_stream, 'BLOCK: unphased'
                        print >>output_stream, ('{vcf_line}\t1\t0\t{chrom}\t'
                                               '{pos}\t{ref}\t{alt}\t'
                                               '{genotype}\tNA\tNA').format(
                                                    vcf_line=counter,
                                                    chrom=tokens[0],
                                                    pos=pos,
                                                    ref=tokens[3],
                                                    alt=tokens[4],
                                                    genotype=tokens[9]
                                                )
                        print >>output_stream, '********' 
                line = vcf_stream.readline().strip()
                counter += 1

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

def get_VAF_pos(VCF):
    """ Obtains position in VCF format/genotype fields of VAF

        VCF: path to input VCF

        Return value: None if VCF does not contain VAF, 
                        otherwise position of VAF
    """
    VAF_check = False
    with open(VCF) as f:
        for line in f:
            # Check header lines to see if FREQ exits in FORMAT fields
            if line[0] == '#':
                if 'FREQ' in line:
                    VAF_check = True
            else:
                # Check first entry to get position of FREQ if it exists
                if VAF_check:
                    tokens = line.strip('\n').split('\t')
                    format_field = tokens[8].split(':')
                    for i in range(0,len(format_field)):
                        if format_field[i] == 'FREQ':
                            VAF_pos = i
                            break
                # Return None if VCF does not contain VAF data
                else:
                    VAF_pos = None
                    break
    return VAF_pos

def process_haplotypes(hapcut_output, interval_dict):
    """ Stores all haplotypes relevant to different transcripts as a dictionary

        hapcut_output: output from HAPCUT2, adjusted to include unphased 
                        mutations as their own haplotypes (performed in 
                        software's prep mode)
        interval_dict: dictionary linking genomic intervals to transcripts

        Return value: dictinoary linking haplotypes to transcripts
    """
    affected_transcripts = collections.defaultdict(list)
    with open(hapcut_output, 'r') as f:
        block_transcripts = collections.defaultdict(list)
        for line in f:
            if line.startswith('BLOCK'):
                # Skip block header lines
                continue
            elif line[0] == '*':
                # Process all transcripts for the block
                for transcript_ID in block_transcripts:
                    block_transcripts[transcript_ID].sort(key=itemgetter(1))
                    haplotype = []
                    for mut in block_transcripts[transcript_ID]:
                        haplotype.append(mut)
                    affected_transcripts[transcript_ID].append(haplotype)
                # Reset transcript dictionary
                block_transcripts = collections.defaultdict(list)
            else:
                # Add mutation to transcript dictionary for the block
                tokens = line.strip("\n").split()
                if len(tokens[5]) == len(tokens[6]):
                    mutation_type = 'V'
                    pos = int(tokens[4])
                    ref = tokens[5]
                    alt = tokens[6]
                    mut_size = len(tokens[5])
                    end = pos + mut_size
                elif len(tokens[5]) > len(tokens[6]):
                    mutation_type = 'D'
                    deletion_size = len(tokens[5]) - len(tokens[6])
                    pos = int(tokens[4]) + (len(tokens[5]) - deletion_size)
                    ref = tokens[5]
                    alt = deletion_size
                    end = pos + deletion_size
                elif len(tokens[5]) < len(tokens[6]):
                    mutation_type = 'I'
                    insertion_size = len(tokens[6]) - len(tokens[5])
                    pos = int(tokens[4])
                    ref = tokens[5]
                    alt = tokens[6][len(ref):]
                    end = pos + 1
                overlapping_transcripts = get_transcripts_from_tree(tokens[3], 
                                                                pos, 
                                                                end,
                                                                interval_dict)
                # For each overlapping transcript, add mutation entry
                # Contains chromosome, position, reference, alternate, allele
                #   A, allele B, genotype line from VCF
                for transcript in overlapping_transcripts:
                    block_transcripts[transcript].append([tokens[3], pos, 
                                                          ref, alt, 
                                                          tokens[1], tokens[2], 
                                                          tokens[7], 
                                                          mutation_type])
    return affected_transcripts

def get_peptides_from_transcripts(relevant_transcripts, VAF_pos, cds_dict,
                                  only_novel_upstream, only_downstream, 
                                  only_reference, reference_index, size_list):
    """ For transcripts that are affected by a mutation, mutations are applied
        and neoepitopes resulting from mutations are called
        
        relevant_transcripts: dictionary linking haplotypes to transcripts;
            output from process_haplotypes()
        VAF_pos: position of VAF in VCF mutation data from HapCUT2
        cds_dict: dictionary linking transcript IDs, to lists of
            relevant CDS/stop codon data; output from gtf_to_cds()
        only_novel_upstream: whether to start translation from novel upstream
            start codons (boolean)
        only_downstream: whether to start translation from only downstream of
            a disrupted canonical start codon (boolean)
        only_reference: whether to start translation only from the canonical
            start codon for a transcript
        reference_index: BowtieIndexReference object for retrieving
            reference genome sequence
        size_list: list of peptide sizes for neoepitope enumeration

        return value: dictionary linking neoepitopes to their associated 
            metadata

        """
    neoepitopes = collections.defaultdict(list)
    for affected_transcript in relevant_transcripts:
        # Create transcript object
        transcriptA = Transcript(reference_index, 
                        [[str(chrom), 'blah', seq_type, str(start), 
                          str(end), '.', strand] for (chrom, seq_type, 
                                                      start, end, strand) 
                          in cds_dict[affected_transcript]], affected_transcript
                        )
        transcriptB = Transcript(reference_index, 
                        [[str(chrom), 'blah', seq_type, str(start), 
                          str(end), '.', strand] for (chrom, seq_type, 
                                                      start, end, strand) 
                          in cds_dict[affected_transcript]], affected_transcript
                        )
        # Iterate over haplotypes associated with this transcript
        haplotypes = relevant_transcripts[affected_transcript]
        for ht in haplotypes:
            # Make edits for each mutation
            for mutation in ht:
                # Determine if mutation is somatic or germline
                if mutation[6][-1] == '*':
                    mutation_class = 'G'
                else:
                    mutation_class = 'S'
                # Determine VAF if available
                if VAF_pos is not None:
                    VAF = float(mutation[6].strip('*').split(':')[VAF_pos].strip('%'))
                else:
                    VAF = None
                # Determine which copies variant exists on & make edits
                if mutation[4] == '1':
                    transcriptA.edit(mutation[3], mutation[1], 
                                mutation_type=mutation[7], 
                                mutation_class=mutation_class,
                                vaf=VAF)
                if mutation[5] == '1':
                    transcriptB.edit(mutation[3], mutation[1], 
                                mutation_type=mutation[7], 
                                mutation_class=mutation_class,
                                vaf=VAF)
            # Extract neoepitopes
            A_peptides = transcriptA.neopeptides(
                                        min_size=size_list[0], 
                                        max_size=size_list[-1],
                                        include_somatic=1,
                                        include_germline=2, 
                                        only_novel_upstream=only_novel_upstream,
                                        only_downstream=only_downstream, 
                                        only_reference=only_reference
                                        )
            B_peptides = transcriptB.neopeptides(
                                        min_size=size_list[0], 
                                        max_size=size_list[-1],
                                        include_somatic=1,
                                        include_germline=2, 
                                        only_novel_upstream=only_novel_upstream,
                                        only_downstream=only_downstream, 
                                        only_reference=only_reference)
            # Store neoepitopes and their metadata
            for pep in A_peptides:
                for meta_data in A_peptides[pep]:
                    adj_meta_data = meta_data + (transcriptA.transcript_id,)
                    if adj_meta_data not in neoepitopes[pep]:
                        neoepitopes[pep].append(adj_meta_data)
            for pep in B_peptides:
                for meta_data in B_peptides[pep]:
                    adj_meta_data = meta_data + (transcriptB.transcript_id,)
                    if adj_meta_data not in neoepitopes[pep]:
                        neoepitopes[pep].append(adj_meta_data)
            transcriptA.reset(reference=True)
            transcriptB.reset(reference=True)
    return neoepitopes

def get_affinity_netMHCIIpan(peptides, allele, netmhciipan, scores,
                                            remove_files=True):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: Allele to use for binding affinity (string)
        netmhciipan: path to netMHCIIpan executable
        scores: list of scoring methods
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities 
                        as strings)
    """
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(os.path.join(
                            os.path.dirname(
                                            os.path.realpath(__file__)
                                ), 'availableAlleles.pickle'
                        ), 'rb'
                ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        # Homogenize format
        allele = allele.replace('HLA-', '')
        if allele not in avail_alleles['netMHCIIpan']:
            warnings.warn(' '.join([allele, 
                                    'is not a valid allele for netMHCIIpan']),
                          Warning)
            if len(scores) == 1:
                return [('NA',) for i in range(0,len(peptides))]
            else:
                return [('NA', 'NA') for i in range(0,len(peptides))]
        allele = allele.replace('*', '_').replace(':', '')
        # Establish return list and sample id
        sample_id = '.'.join([peptides[0],
                                str(len(peptides)), allele,
                                'netmhciipan'])
        affinities = []
        # Write one peptide per line to a temporary file for 
        #   input if peptide length is at least 9
        # Count instances of smaller peptides
        na_count = 0
        peptide_file = tempfile.mkstemp(
                        suffix='.peptides', prefix='id.', text=True)[1]
        files_to_remove.append(peptide_file)
        with open(peptide_file, 'w') as f:
            for sequence in peptides:
                if len(sequence) >= 9:
                    print >>f, sequence
                else:
                    na_count += 1
        if na_count > 0:
            warnings.warn(' '.join([str(na_count),
                                    'peptides not compatible with',
                                    'netMHCIIpan will not receive score']),
                            Warning)
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(suffix='.netMHCIIpan.out', 
                                    prefix='id.', text=True)[1]
        files_to_remove.append(mhc_out)
        # Run netMHCIIpan
        subprocess.check_call(
                        [netmhciipan, '-a', allele, '-inptype', '1', 
                         '-xls', '-xlsfile', mhc_out, '-f', peptide_file]
                    )
        # Retrieve scores for valid peptides
        score_dict = {}
        with open(mhc_out, 'r') as f:
            # Skip headers
            f.readline()
            f.readline()
            for line in f:
                # token 1 is peptide; token 4 is affinity; token[5] is rank
                tokens = line.split('\t')
                if sorted(scores) == ['affinity', 'rank']:
                    score_dict[tokens[1]] = (tokens[4], tokens[5])
                elif sorted(scores) == ['affinity']:
                    score_dict[tokens[1]] = (tokens[4],)
                elif (scores) == ['rank']:
                    score_dict[tokens[1]] = (tokens[5],)
        # Produce list of scores for valid peptides
        # Invalid peptides receive "NA" score
        for sequence in peptides:
            if sequence in score_dict:
                nM = score_dict[sequence]
            else:
                if len(scores) == 1:
                    nM = ('NA',)
                else:
                    nM = ('NA', 'NA')
            affinities.append(nM)
        return affinities
    finally:
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)

def get_affinity_netMHCpan(peptides, allele, netmhcpan, version, scores,
                                remove_files=True):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: allele to use for binding affinity 
                    (string, format HLA-A02:01)
        netmhcpan: path to netMHCpan executable
        version: version of netMHCpan software
        scores: list of scoring methods
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities 
                        as strings)
    """
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(os.path.join(
                            os.path.dirname(
                                            os.path.realpath(__file__)
                                ), 'availableAlleles.pickle'
                        ), 'rb'
                ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        allele = allele.replace('*', '')
        if allele not in avail_alleles['netMHCpan']:
            warnings.warn(' '.join([allele, 
                                    'is not a valid allele for netMHCpan']),
                            Warning)
            if len(scores) == 1:
                return [('NA',) for i in range(0,len(peptides))]
            else:
                return [('NA', 'NA') for i in range(0,len(peptides))]
        # Establish return list and sample id
        sample_id = '.'.join([peptides[0], str(len(peptides)), 
                                allele, 'netmhcpan'])
        affinities = []
        # Write one peptide per line to a temporary file for input
        peptide_file = tempfile.mkstemp(suffix='.peptides', 
                                        prefix=''.join([sample_id, 
                                                        '.']), 
                                        text=True)[1]
        files_to_remove.append(peptide_file)
        with open(peptide_file, 'w') as f:
            for sequence in peptides:
                print >>f, sequence
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(suffix='.netMHCpan.out', 
                                    prefix=''.join([sample_id, 
                                                    '.']), 
                                    text=True)[1]
        files_to_remove.append(mhc_out)
        # Run netMHCpan
        if version == '3':
            subprocess.check_call(
                [netmhcpan, '-a', allele, '-inptype', '1', '-p', '-xls', 
                    '-xlsfile', mhc_out, peptide_file])
        elif version == '4':
            subprocess.check_call(
                [netmhcpan, '-BA', '-a', allele, '-inptype', '1', '-p', '-xls', 
                    '-xlsfile', mhc_out, peptide_file])
        with open(mhc_out, 'r') as f:
            f.readline()
            f.readline()
            for line in f:
                tokens = line.strip('\n').split('\t')
                # for v3, tokens[5] is affinity, tokens[6] is rank
                # for v4, tokens[6] is affinty, tokens[7] is rank
                if sorted(scores) == ['affinity', 'rank']:
                    if version == '3':
                        nM = (tokens[5], tokens[6])
                    elif version == '4':
                        nM = (tokens[6], tokens[7])
                elif sorted(scores) == ['affinity']:
                    if version == '3':
                        nM = (tokens[5], )
                    elif version == '4':
                        nM = (tokens[6], )
                elif sorted(scores) == ['rank']:
                    if version == '3':
                        nM = (tokens[6], )
                    elif version == '4':
                        nM = (tokens[7], )
                affinities.append(nM)
        return affinities
    finally:
        # Remove temporary files
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)

def gather_binding_scores(neoepitopes, tool_dict, hla_alleles):
    """ Adds binding scores from desired programs to neoepitope metadata
        
        neoepitopes: dictionary linking neoepitopes to their metadata
        tool_dict: dictionary storing prediction tool data
        hla_alleles: list of HLA alleles used for binding predictions
        
        Return value: dictionary linking neoepitopes to their metadata,
            which now includes binding scores
    """
    for allele in hla_alleles:
        for tool in sorted(tool_dict.keys()):
            if tool == 'netMHCIIpan3':
                binding_scores = get_affinity_netMHCIIpan(
                                                sorted(neoepitopes.keys()), 
                                                allele, tool_dict[tool][0], 
                                                tool_dict[tool][1],
                                                remove_files=True
                                                )
            elif tool == 'netMHCpan3':
                binding_scores = get_affinity_netMHCpan(
                                                sorted(neoepitopes.keys()), 
                                                allele, tool_dict[tool][0], 
                                                '3', tool_dict[tool][1],
                                                remove_files=True
                                                )
            elif tool == 'netMHCpan4':
                binding_scores = get_affinity_netMHCpan(
                                                sorted(neoepitopes.keys()), 
                                                allele, tool_dict[tool][0], 
                                                '4', tool_dict[tool][1],
                                                remove_files=True
                                                )
            for i in range(0, len(neoepitopes.keys())):
                meta_data = neoepitopes[sorted(neoepitopes.keys())[i]]
                for j in range(0, len(meta_data)):
                    neoepitopes[sorted(neoepitopes.keys())[i]][j] = meta_data[j] + binding_scores[i]
    return neoepitopes

def write_results(output_file, hla_alleles, neoepitopes, tool_dict):
    """ Writes predicted neoepitopes out to file
        
        output_file: path to output file
        hla_alleles: list of HLA alleles used for binding predictions
        neoepitopes: dictionary linking neoepitopes to their metadata
        tool_dict: dictionary storing prediction tool data

        Return value: None.   
    """
    with open(output_file, 'w') as o:
        headers = ['Neoepitope', 'Chromsome', 'Pos', 'Ref', 'Alt', 
                   'Mutation_type', 'VAF', 'Transcript_ID']
        for allele in hla_alleles:
            for tool in sorted(tool_dict.keys()):
                for score_method in sorted(tool_dict[tool][1]):
                    headers.append('_'.join([tool, allele, score_method]))
        o.write('\t'.join(headers) + '\n')
        for epitope in sorted(neoepitopes.keys()):
            if len(neoepitopes[epitope]) == 1:
                mutation = neoepitopes[epitope][0]
                if mutation[2] == '':
                    ref = "''"
                else:
                    ref = mutation[2]
                if mutation[3] == '':
                    alt = "''"
                else:
                    alt = mutation[2]
                out_line = [epitope, mutation[0], str(mutation[1]), ref, alt,
                            mutation[4], str(mutation[5]), mutation[6]]
                for i in range(7,len(mutation)):
                    out_line.append(str(mutation[i]))
                o.write('\t'.join(out_line) + '\n')
            else:
                mutation_dict = collections.defaultdict(list)
                ep_scores = []
                for i in range(7, len(neoepitopes[epitope][0])):
                    ep_scores.append(neoepitopes[epitope][0][i])
                for mut in neoepitopes[epitope]:
                    if mut[2] == '':
                        ref = "''"
                    else:
                        ref = mut[2]
                    if mut[3] == '':
                        alt = "''"
                    else:
                        alt = mut[3]
                    if mut[4] == 'V':
                        mut_type = 'SNV'
                    elif mut[4] == 'I':
                        mut_type = 'Insertion'
                    elif mut[4] == 'D':
                        mut_type = 'Deletion'
                    mutation_dict[(mut[0], mut[1], ref, alt, mut_type)].append(
                                                                [str(mut[5]), 
                                                                 mut[6]]
                                                                 )
                for mut in sorted(mutation_dict.keys()):
                    out_line = [epitope, mut[0], str(mut[1]), mut[2], mut[3],
                                mut[4],
                                ';'.join([x[0] for x in mutation_dict[mut]]),
                                ';'.join([x[1] for x in mutation_dict[mut]])]
                    for score in ep_scores:
                        out_line.append(str(score))
                    o.write('\t'.join(out_line) + '\n')

def test():
    """ Runs all unit tests """
    del sys.argv[1:] # Don't choke on extra command-line parameters
    import unittest
    import filecmp
    class TestGTFprocessing(unittest.TestCase):
        """Tests proper creation of dictionaries store GTF data"""
        def setUp(self):
            """Sets up gtf file and creates dictionaries for tests"""
            self.gtf = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Ychrom.gtf'
                            )
            self.Ycds = gtf_to_cds(self.gtf, 'NA', pickle_it=False)
            self.Ytree = cds_to_tree(self.Ycds, 'NA', pickle_it=False)
        def test_transcript_to_CDS(self):
            """Fails if dictionary was built incorrectly"""
            self.assertEqual(len(self.Ycds.keys()), 164)
        def test_CDS_tree(self):
            """Fails if dictionary was built incorrectly"""
            self.assertEqual(len(self.Ytree.keys()), 1)
            self.assertEqual(len(self.Ytree['Y']), 2174)
        def test_transcript_extraction(self):
            """Fails if incorrect transcripts are pulled"""
            self.assertEqual(len(get_transcripts_from_tree(
                                                            'Y', 150860, 
                                                            150861,
                                                            self.Ytree)), 
                                                            3
                                                            )
            self.coordinate_search = list(self.Ytree['Y'].search(150860,
                                                                 150861))
            self.transcripts = []
            for interval in self.coordinate_search:
                self.transcripts.append(interval[2])
            self.transcripts.sort()
            self.assertEqual(
                             self.transcripts, ['ENST00000381657.7_3_PAR_Y',
                                                'ENST00000381663.8_3_PAR_Y',
                                                'ENST00000399012.6_3_PAR_Y']
                            )
    class TestVCFmerging(unittest.TestCase):
        """Tests proper merging of somatic and germline VCFS"""
        def setUp(self):
            """Sets up files to use for tests"""
            self.varscan = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Ychrom.varscan.vcf'
                            )                
            self.germline = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Ychrom.germline.vcf'
                            )   
            self.precombined = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Ychrom.combined.vcf'
                            )  
            self.outvcf = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Ychrom.testcombine.vcf'
                            )  
            combinevcf(self.germline, self.varscan, self.outvcf)
        def test_merge(self):
            """Fails if VCFs were merged improperly"""
            self.assertTrue(filecmp.cmp(self.outvcf, self.precombined))
        def tearDown(self):
            """Removes test file"""
            os.remove(self.outvcf)
    class TestPrepHapCUT(unittest.TestCase):
        """Tests addition of unphased mutations to HapCUT2 output"""
        def setUp(self):
            """Sets up hapcut and vcf files to use for tests"""
            self.hapcut = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'test.hapcut.out'
                            )
            self.vcf = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'test.vcf'
                            )
            self.complete_hapcut = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'complete_hapcut.out'
                            )
            self.test_hapcut = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'test_complete_hapcut.out'
                            )
        def test_haplotype_prep(self):
            """Tests that output of haplotype prep is correct"""
            prep_hapcut_output(self.test_hapcut, self.hapcut, self.vcf)
            self.assertTrue(filecmp.cmp(self.test_hapcut,
                                        self.complete_hapcut))
        def tearDown(self):
            """Removes test file"""
            os.remove(self.test_hapcut)
    class TestVAFpos(unittest.TestCase):
        """Tests fetching of VAF position from VCF file"""
        def setUp(self):
            """Sets up vcf files to use for tests"""
            self.varscan = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'test', 'Ychrom.varscan.vcf'
                            )
            self.mutect = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Ychrom.mutect.vcf'
                            )
        def test_position(self):
            """Fails if incorrect positions are returned"""
            self.assertEqual(get_VAF_pos(self.varscan), 5)
            self.assertEqual(get_VAF_pos(self.mutect), None)
    class TestHaplotypeProcessing(unittest.TestCase):
        """Tests proper processing of HAPCUT2 files"""
        def setUp(self):
            """Sets up input files and dictionaries to use for tests"""
            self.ref_prefix = os.path.join(
                        os.path.dirname(
                                os.path.dirname(
                                        os.path.realpath(__file__)
                                    )
                            ), 'tests', 'Chr11.ref'
                    )
            self.reference_index = bowtie_index.BowtieIndexReference(
                                                            self.ref_prefix)
            self.Chr11gtf = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Chr11.gtf'
                            )
            self.Chr11cds = gtf_to_cds(self.Chr11gtf, 'NA', pickle_it=False)
            self.Chr11tree = cds_to_tree(self.Chr11cds, 'NA',
                                            pickle_it=False)
            self.Chr11hapcut = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'Chr11.hapcut.out'
                            )
        def test_hap_processing(self):
            """Fails if file is processed incorrectly"""
            Chr11_txs = process_haplotypes(self.Chr11hapcut, self.Chr11tree)
            self.assertEqual(sorted(Chr11_txs.keys()), 
                                    ['ENST00000398531.2_2'])
            self.assertEqual(Chr11_txs['ENST00000398531.2_2'],
                            [[['11', 71276862, 'TGT', 2, '0', '0', 
                               '0/0:.:53:52:0:0%:22,30,0,0:.:2', 'D'], 
                              ['11', 71276900, 'C', 'G', '0', '0',
                               '0/0:.:35:34:0:0%:19,15,0,0:.:2', 'V'],
                              ['11', 71277000, 'G', 'AA', '0', '0',
                               '0/0:.:35:34:0:0%:19,15,0,0:.:2', 'I']]]
                )
        def test_peptide_gathering(self):
            Chr11_txs = {'ENST00000398531.2_2': [
                                [['11', 71276651, 'CTC', 3, '0', '1', 
                                  '0/0:.:53:52:0:3.0%:22,30,0,0:.:2', 'D'], 
                                 ['11', 71277229, 'A', 'C', '0', '1',
                                  '0/0:.:35:34:0:15.7%:19,15,0,0:.:2', 'V'],
                                 ['11', 71277056, 'G', 'AAA', '1', '1',
                                  '0/0:.:35:34:0:0.1%:19,15,0,0:.:2', 'I']]
                                 ]
                        }
            transcript_blocks = self.Chr11cds['ENST00000398531.2_2']
            neoepitopes = get_peptides_from_transcripts(Chr11_txs, 5, 
                                                        self.Chr11cds,
                                                        True, False, False, 
                                                        self.reference_index,
                                                        [8,9,10,11])
            self.assertEqual(len(neoepitopes.keys()), 70)
            self.assertEqual(neoepitopes['CGCSQKCN'], [('11', 71277056, '',
                                                        'AAA', 'I', 0.1, 
                                                        'ENST00000398531.2_2')])
            self.assertEqual(neoepitopes['PVCCPCKI'], [('11', 71277229, 
                                                        'A', 'C', 'V', 15.7, 
                                                        'ENST00000398531.2_2')])
            self.assertEqual(sorted(neoepitopes.keys())[0], 'CCGCGGCG')
            self.assertEqual(sorted(neoepitopes.keys())[-1], 'VPVCCPCKI')
    class TestBindingPrediction(unittest.TestCase):
        """Tests binding prediction functions"""
        def setUp(self):
            """"""
            self.neoepitopes = {'CGCSQKCN': [('11', 71277056, '', 'AAA', 
                                              'I', 0.1, 
                                              'ENST00000398531.2_2')],
                                'PVCCPCKI': [('11', 71277229, 'A', 'C', 'V', 
                                              15.7, 'ENST00000398531.2_2')]}
            self.tools = {'netMHCpan4': [exe_paths.netMHCpan4, 
                                         ['rank', 'affinity']],
                          'netMHCpan3': [exe_paths.netMHCpan3, 
                                         ['rank', 'affinity']]}
            self.alleles = ['HLA-A*02:01', 'HLA-B*07:02']
            self.neoepitopesII = {'AHWTEARIMLDNINM': [('11', 71277056, '', 'AAA', 
                                              'I', 0.1, 
                                              'ENST00000398531.2_2')]}
            self.toolsII = {'netMHCIIpan3': [exe_paths.netMHCIIpan3, 
                                           ['rank', 'affinity']]}
            self.allelesII = ['HLA-DRB1*01:01']
        def test_binding_scores(self):
            new_neoepitopes = gather_binding_scores(self.neoepitopes, 
                                                    self.tools, 
                                                    self.alleles)
            new_neoepitopesII = gather_binding_scores(self.neoepitopesII, 
                                                    self.toolsII, 
                                                    self.allelesII)
            self.assertEqual(new_neoepitopes['CGCSQKCN'], [('11', 71277056, 
                                                            '', 'AAA', 
                                                            'I', 0.1, 
                                                            'ENST00000398531.2_2',
                                                            '48126.8320', 
                                                            '100.0000', 
                                                            '47206.2734', 
                                                            '96.4206', 
                                                            '48380.5781', 
                                                            '100.0000',
                                                            '47300.3477',
                                                            '97.0230'
                                                            )])
            self.assertEqual(new_neoepitopes['PVCCPCKI'], [('11', 71277229, 
                                                            'A', 'C', 'V', 
                                                            15.7, 
                                                            'ENST00000398531.2_2',
                                                            '41484.8242', 
                                                            '70.0000', 
                                                            '38923.1328', 
                                                            '66.8880',
                                                            '22375.3965', 
                                                            '18.0000',
                                                            '29030.9199',
                                                            '25.8751'
                                                            )])
            self.assertEqual(new_neoepitopesII['AHWTEARIMLDNINM'], 
                             [('11', 71277056, '', 'AAA', 'I', 0.1, 
                               'ENST00000398531.2_2', '295.78', '49.00')])
            self.assertEqual(15, len(new_neoepitopes['PVCCPCKI'][0]))
            self.assertEqual(15, len(new_neoepitopes['CGCSQKCN'][0]))
            self.assertEqual(9, len(new_neoepitopesII['AHWTEARIMLDNINM'][0]))
    class TestOutput(unittest.TestCase):
        """Tests function to write output"""
        def setUp(self):
            """Sets up paths and dictionaries"""
            self.out_file = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'neoepiscope.out'
                            )
            self.correct_out = os.path.join(
                                os.path.dirname(
                                        os.path.dirname(
                                                os.path.realpath(__file__)
                                            )
                                    ), 'tests', 'expected.neoepiscope.out'
                            )
            self.tools = {'netMHCpan4': ['netMHCpan', ['rank', 'affinity']],
                          'netMHCIIpan3': ['netMHCIIpan', ['rank']]}
            self.HLA_alleles = ['HLA*A01:01', 'HLA*A02:01']
            self.neoepitopes = {'CGCSQKCN': [('11', 71277056, '',
                                              'AAA', 'I', 0.1, 
                                              'ENST00000398531.2_2',
                                              5, 10000.0, 1, 5, 150, 4),
                                             ('4', 300000, 'A',
                                              'T', 'V', 10.2, 
                                              'ENST00000398554.1_1',
                                              0.5, 100, 0.5, 3, 150, 4)],
                                'PVCCPCKI': [('11', 71277229, 'A', 'C', 'V', 
                                              15.7, 'ENST00000398531.2_2',
                                              10, 50.57, 1.2, 3, 10.1, 7),
                                             ('11', 71277229, 'A', 'C', 'V', 
                                              20.3, 'ENST00000398200.3_1',
                                              10, 50.57, 1.2, 3, 10.1, 7)]}
        def testwrite(self):
            """Tests that output file is written correctly"""
            write_results(self.out_file, self.HLA_alleles, self.neoepitopes,
                          self.tools)
            self.assertTrue(filecmp.cmp(self.out_file, self.correct_out))
        def tearDown(self):
            """Removes test file"""
            os.remove(self.out_file)
    unittest.main()

def main():
    """ Entry point for neoepiscope software """
    parser = argparse.ArgumentParser(description=_help_intro, 
                                        formatter_class=help_formatter)
    subparsers = parser.add_subparsers(help=(
                                    'subcommands; add "-h" or "--help" '
                                    'after a subcommand for its parameters'
                                ), dest='subparser_name')
    test_parser = subparsers.add_parser('test', help='performs unit tests')
    index_parser = subparsers.add_parser('index',
                                        help=('produces pickled dictionaries '
                                        'linking transcripts to intervals and '
                                        ' CDS lines in a GTF'))
    swap_parser = subparsers.add_parser('swap',
                                        help=('swaps tumor and normal columns '
                                        'in a somatic vcf if necessary for '
                                        'proper HapCUT2 results'))
    merge_parser = subparsers.add_parser('merge',
                                         help=('merges germline and somatic '
                                               'VCFS for combined mutation '
                                               'phasing with HAPCUT2'))
    prep_parser = subparsers.add_parser('prep',
                                        help=('combines HAPCUT2 output with '
                                              'unphased variants for call '
                                              'mode'))
    call_parser = subparsers.add_parser('call', help='calls neoepitopes')
    # Index parser options (produces pickled dictionaries for transcript data)
    index_parser.add_argument('-g', '--gtf', type=str, required=True,
            help='input path to GTF file'
        )  
    index_parser.add_argument('-d', '--dicts', type=str, required=True,
            help='output path to pickled CDS dictionary'
        )
    # Swap parser options (swaps columns in somatic VCF)
    swap_parser.add_argument('-i', '--input', type=str, required=True,
            help='input path to somatic VCF'
        )
    swap_parser.add_argument('-o', '--output', type=str, required=False,
            help='output path to column-swapped VCF'
        )
    # Merger parser options (merges somatic and germline VCFs)
    merge_parser.add_argument('-g', '--germline', type=str, required=True,
            help='input path to germline VCF'
        )
    merge_parser.add_argument('-s', '--somatic', type=str, required=True,
            help='input path to somatic VCF'
        )
    merge_parser.add_argument('-o', '--output', type=str, required=False,
            help='output path to combined VCF'
        )
    # Prep parser options (adds unphased mutations as their own haplotype)
    prep_parser.add_argument('-v', '--vcf', type=str, required=True,
            help='input VCF'
        )
    prep_parser.add_argument('-c', '--hapcut2-output', type=str, required=True,
            help='path to output file of HAPCUT2 run on input VCF'
        )
    prep_parser.add_argument('-o', '--output', type=str, required=True,
            help='path to output file to be input to call mode'
        )
    # Call parser options (calls neoepitopes)
    call_parser.add_argument('-x', '--bowtie-index', type=str, required=True,
            help='path to Bowtie index basename'
        )
    call_parser.add_argument('-v', '--vcf', type=str, required=True,
            help='input path to VCF'
        )
    call_parser.add_argument('-d', '--dicts', type=str, required=True,
            help='input path to pickled CDS dictionary'
        )
    call_parser.add_argument('-c', '--merged-hapcut2-output', type=str,
            required=True,
            help='path to output of prep subcommand'
        )
    call_parser.add_argument('-k', '--kmer-size', type=str, required=False,
            default='8,11', help='kmer size for epitope calculation'
        )
    call_parser.add_argument('-p', '--affinity-predictor', type=str, 
            nargs=3, required=False, action='append',
            help='binding affinity prediction software,'
                'associated version number, and scoring method(s) '
                '(e.g. -p netMHCpan 4 rank,affinity); '
                'for multiple softwares, repeat the argument;'
                'see documentation for details'
        )
    call_parser.add_argument('-a', '--alleles', type=str, required=True,
            help='comma separated list of alleles; '
                 'see documentation online for more information'
        )
    call_parser.add_argument('-o', '--output_file', type=str, required=True,
            help='path to output file'
        )
    call_parser.add_argument('-u', '--upstream_atgs', type=str, required=False,
            default='novel', help='how to handle upstream start codons, see '
            'documentation online for more information'
        )
    args = parser.parse_args()
    if args.subparser_name == 'test':
        test()
    elif args.subparser_name == 'index':
        cds_dict = gtf_to_cds(args.gtf, args.dicts)
        tree = cds_to_tree(cds_dict, args.dicts)
    elif args.subparser_name == 'swap':
        adjust_tumor_column(args.input, args.output)
    elif args.subparser_name == 'merge':
        combinevcf(args.germline, args.somatic, outfile=args.output)
    elif args.subparser_name == 'prep':
        prep_hapcut_output(args.output, args.hapcut2_output, args.vcf)
    elif args.subparser_name == 'call':
        # Load pickled dictionaries
        interval_dict = pickle.load(open(args.dicts + ''.join([dictdir, 
                                '/', 'intervals_to_transcript.pickle']), 'rb'))
        cds_dict = pickle.load(open(args.dicts + ''.join([dictdir, 
                                '/', 'transcript_to_CDS.pickle']), 'rb'))
        # Check affinity predictor
        tool_dict = {}
        if args.affinity_predictor is not None:
            for tool in args.affinity_predictor:
                program = tool[0]
                version = tool[1]
                scoring = tool[2].split(',')
                if 'netMHCIIpan' in program:
                    if version == '3' and 'netMHCIIpan3' not in tool_dict:
                        program = exe_paths.netMHCIIpan3
                        if program is None:
                            program = which('netMHCIIpan3')
                        else:
                            program = which(program)
                        if program is None:
                            warnings.warn(' '.join(['No valid install of', 
                                            'netMHCIIpan available']),
                                            Warning)
                            continue
                        acceptable_scoring = ['rank', 'affinity']
                        for method in scoring:
                            if method not in acceptable_scoring:
                                warnings.warn(' '.join([method, 
                                            'not compatible with netMHCIIpan']),
                                            Warning)
                                scoring.remove(method)
                        if len(scoring) > 0:
                            tool_dict['netMHCIIpan3'] = [program, sorted(scoring)]
                    elif 'netMHCIIpan3' in tool_dict:
                        raise RuntimeError('Conflicting or repetitive installs'
                                            'of netMHCIIpan given')
                    else:
                        raise NotImplementedError(
                            ' '.join(['Neoepiscope does not support version', 
                                      version, 'of netMHCIIpan'])
                            )
                elif 'netMHCpan' in program:
                    if (('netMHCpan3' not in tool_dict and version == '3') or 
                                    ('netMHCpan4' not in tool_dict and 
                                        version == '4')):
                        if version == '3':
                            program = exe_paths.netMHCpan3
                            if program is None:
                                program = which('netMHCpan3')
                            else:
                                program = which(program)
                        elif version == '4':
                            if program is None:
                                program = which('netMHCpan4')
                            else:
                                program = which(program)
                        if program is None:
                            warnings.warn(' '.join(['No valid install of ', 
                                            'netMHCIIpan available']),
                                            Warning)
                            continue
                        if program is None:
                            warnings.warn(' '.join(['No valid install of', 
                                            'netMHCpan version', version, 
                                            'available']),  Warning)
                            continue
                        acceptable_scoring = ['rank', 'affinity']
                        for method in scoring:
                            if method not in acceptable_scoring:
                                warnings.warn(' '.join([method, 
                                            'not compatible with netMHCpan']),
                                        Warning)
                                scoring.remove(method)
                        if len(scoring) > 0:
                            if version == '3':
                                name = 'netMHCpan3'
                            elif version == '4':
                                name = 'netMHCpan4'
                            tool_dict[name] = [program, sorted(scoring)]
                    elif (('netMHCpan3' in tool_dict and version == '3') or 
                                ('netMHCpan4' in tool_dict and version == '4')):
                        raise RuntimeError('Conflicting or repetitive installs'
                                            'of netMHCpan given')
                    else:
                        raise NotImplementedError(
                            ' '.join(['Neoepiscope does not support version', 
                                      version, 'of netMHCpan'])
                            )
                else:
                    raise NotImplementedError(
                                    ' '.join(['Neoepiscope does not support', 
                                              program, 
                                              'for binding predictions'])
                                    )
        if len(tool_dict.keys()) == 0:
            warnings.warn(' '.join(['No binding prediction tools given,',
                                    'will proceed without binding predictions']),
                            Warning)
        # Obtain VAF frequency VCF position
        VAF_pos = get_VAF_pos(args.vcf)
        # Obtain peptide sizes for kmerizing peptides
        if ',' in args.kmer_size:
            size_list = args.kmer_size.split(',')
            size_list.sort()
            for i in range(0, len(size_list)):
                size_list[i] = int(size_list[i])
        hla_alleles = sorted(args.alleles.split(','))
        # For retrieving genome sequence
        reference_index = bowtie_index.BowtieIndexReference(args.bowtie_index)
        # Find transcripts that haplotypes overlap 
        relevant_transcripts = process_haplotypes(args.merged_hapcut2_output, 
                                                    interval_dict)
        # Establish handling of ATGs
        if arg.upstream_atgs == 'novel':
            only_novel_upstream = True
            only_downstream = False
            only_reference = False
        elif args.upstream_atgs == 'all':
            only_novel_upstream = False
            only_downstream = False
            only_reference = False
        elif args.upstream_atgs == 'none':
            only_novel_upstream = False
            only_downstream = True
            only_reference = False
        elif args.upstream_atgs == 'reference':
            only_novel_upstream = False
            only_downstream = False
            only_reference = True
        else:
            raise RuntimeError('--upstream_atgs must be one of {"novel", "all", '
                                '"none", "reference"}')
        # Apply mutations to transcripts and get neoepitopes
        neoepitopes = get_peptides_from_transcripts(relevant_transcripts, 
                                                    VAF_pos, cds_dict,
                                                    only_novel_upstream,
                                                    only_downstream, 
                                                    only_reference,
                                                    reference_index,
                                                    size_list)
        full_neoepitopes = gather_binding_scores(neoepitopes, tool_dict, 
                                                 hla_alleles)
        write_results(args.output_file, hla_alleles, full_neoepitopes, tool_dict)
    else:
        raise RuntimeError(''.join([args.subparser_name, 
                            ' is not a valid software mode']))

if __name__ == '__main__':
    main()
