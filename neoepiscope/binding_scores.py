#!/usr/bin/env python
# coding=utf-8
"""
binding_scores.py

Part of neoepiscope
Includes functions for grabbing MHC binding scores from third-party software.

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
from . import paths
from .file_processing import which
import os
import warnings
import tempfile
import pickle
import subprocess

neoepiscope_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get_binding_tools(binding_tool_list):
    """ Processes user-specified binding tools to ensure usability

        binding_tool_list: nested list of binding affinity tool info
            where inner lists are [name of tool, version of tool, 
            comma separated scoring methods]

        Return value: dictionary with binding tool IDs as keys and
            lists as values, where lists are [program executable,
            sorted list of scoring methods]
    """
    tool_dict = {}
    if len(binding_tool_list) > 1:
        binding_tool_list.remove(["mhcflurry", "1", "affinity,rank"])
    for tool in binding_tool_list:
        program = tool[0]
        version = tool[1]
        scoring = tool[2].split(",")
        if "mhcflurry" in program.lower():
            if version == "1" and "mhcflurry1" not in tool_dict:
                program = "mhcflurry-predict"
                acceptable_scoring = ["rank", "affinity", "high", "low"]
                for method in scoring:
                    if method not in acceptable_scoring:
                        warnings.warn(
                            " ".join([method, "not compatible with mhcflurry"]),
                            Warning,
                        )
                        scoring.remove(method)
                if len(scoring) > 0:
                    tool_dict["mhcflurry1"] = [program, sorted(scoring)]
                else:
                    warnings.warn(
                            " ".join(
                                [
                                    "No compatible scoring methods given",
                                    "for mhcflurry version", version,
                                    "- will not use this tool for",
                                    "binding predictions"
                                ]
                            ),
                            Warning,
                        )
            elif "mhcflurry1" in tool_dict:
                raise RuntimeError(
                    "Conflicting or repetitive installs of mhcflurry given"
                )
            else:
                raise NotImplementedError(
                    " ".join(
                        [
                            "neoepiscope does not support version",
                            version,
                            "of mhcflurry",
                        ]
                    )
                )
        elif "mhcnuggets" in program.lower():
            if version == "2" and "mhcnuggets2" not in tool_dict:
                program = "NA"
                acceptable_scoring = ["affinity"]
                for method in scoring:
                    if method not in acceptable_scoring:
                        warnings.warn(
                            " ".join(
                                [method, "not compatible with mhcnuggets"]
                            ),
                            Warning,
                        )
                        scoring.remove(method)
                if len(scoring) > 0:
                    tool_dict["mhcnuggets2"] = [program, sorted(scoring)]
                else:
                    warnings.warn(
                            " ".join(
                                [
                                    "No compatible scoring methods given",
                                    "for mhcnuggets version", version,
                                    "- will not use this tool for",
                                    "binding predictions"
                                ]
                            ),
                            Warning,
                        )
            elif "mhcnuggets2" in tool_dict:
                raise RuntimeError(
                    "Conflicting or repetitive installs of mhcnuggets given"
                )
            else:
                raise NotImplementedError(
                    " ".join(
                        [
                            "neoepiscope does not support version",
                            version,
                            "of mhcnuggets",
                        ]
                    )
                )
        elif "netmhciipan" in program.lower():
            if version == "3" and "netMHCIIpan3" not in tool_dict:
                program = paths.netMHCIIpan3
                if program is None:
                    program = which("netMHCIIpan3")
                else:
                    program = which(program)
                if program is None:
                    warnings.warn(
                        " ".join(
                            ["No valid install of", "netMHCIIpan available"]
                        ),
                        Warning,
                    )
                    continue
                acceptable_scoring = ["rank", "affinity"]
                for method in scoring:
                    if method not in acceptable_scoring:
                        warnings.warn(
                            " ".join(
                                [method, "not compatible with netMHCIIpan"]
                            ),
                            Warning,
                        )
                        scoring.remove(method)
                if len(scoring) > 0:
                    tool_dict["netMHCIIpan3"] = [program, sorted(scoring)]
                else:
                    warnings.warn(
                            " ".join(
                                [
                                    "No compatible scoring methods given",
                                      "for netMHCIIpan version", version,
                                      "- will not use this tool for",
                                      "binding predictions"
                                ]
                            ),
                            Warning,
                        )
            elif "netMHCIIpan3" in tool_dict:
                raise RuntimeError(
                    "Conflicting or repetitive installs of netMHCIIpan given"
                )
            else:
                raise NotImplementedError(
                    " ".join(
                        [
                            "neoepiscope does not support version",
                            version,
                            "of netMHCIIpan",
                        ]
                    )
                )
        elif "netmhcpan" in program.lower():
            if ("netMHCpan3" not in tool_dict and version == "3") or (
                "netMHCpan4" not in tool_dict and version == "4"
            ):
                if version == "3":
                    program = paths.netMHCpan3
                    if program is None:
                        program = which("netMHCpan3")
                    else:
                        program = which(program)
                elif version == "4":
                    program = paths.netMHCpan4
                    if program is None:
                        program = which("netMHCpan4")
                    else:
                        program = which(program)
                if program is None:
                    warnings.warn(
                        " ".join(
                            [
                                "No valid install of",
                                "netMHCpan version",
                                version,
                                "available",
                            ]
                        ),
                        Warning,
                    )
                    continue
                acceptable_scoring = ["rank", "affinity"]
                for method in scoring:
                    if method not in acceptable_scoring:
                        warnings.warn(
                            " ".join([method, "not compatible with netMHCpan"]),
                            Warning,
                        )
                        scoring.remove(method)
                if len(scoring) > 0:
                    if version == "3":
                        name = "netMHCpan3"
                    elif version == "4":
                        name = "netMHCpan4"
                    tool_dict[name] = [program, sorted(scoring)]
                else:
                    warnings.warn(
                            " ".join(
                                [
                                    "No compatible scoring methods given",
                                      "for NetMHCpan version", version,
                                      "- will not use this tool for",
                                      "binding predictions"
                                ]
                            ),
                            Warning,
                        )
            elif ("netMHCpan3" in tool_dict and version == "3") or (
                "netMHCpan4" in tool_dict and version == "4"
            ):
                raise RuntimeError(
                    "Conflicting or repetitive installs of netMHCpan given"
                )
            else:
                raise NotImplementedError(
                    " ".join(
                        [
                            "neoepiscope does not support version",
                            version,
                            "of netMHCpan",
                        ]
                    )
                )
        elif "netmhccons" in program.lower():
            raise NotImplementedError("Binding predictions with netMHCcons "
                                      "not currently supported due to instability"
                                      " of predictions - support for netMHCcons "
                                      "may be included in future releases.")
            if "netMHCcons1" not in tool_dict and version == "1":
                program = paths.netMHCcons1
                if program is None:
                    program = which("netMHCcons1")
                else:
                    program = which(program)
                if program is None:
                    warnings.warn(
                        " ".join(
                            [
                                "No valid install of",
                                "netMHCcons version",
                                version,
                                "available",
                            ]
                        ),
                        Warning,
                    )
                    continue
                acceptable_scoring = ["rank", "affinity"]
                for method in scoring:
                    if method not in acceptable_scoring:
                        warnings.warn(
                            " ".join([method, "not compatible with netMHCcons"]),
                            Warning,
                        )
                        scoring.remove(method)
                if len(scoring) > 0:
                    if version == "1":
                        name = "netMHCcons1"
                    tool_dict[name] = [program, sorted(scoring)]
                else:
                    warnings.warn(
                            " ".join(
                                [
                                    "No compatible scoring methods given",
                                      "for NetMHCcons version", version,
                                      "- will not use this tool for",
                                      "binding predictions"
                                ]
                            ),
                            Warning,
                        )
            elif "netMHCcons1" in tool_dict and version == "1":
                raise RuntimeError(
                    "Conflicting or repetitive installs of netMHC given"
                )
            else:
                raise NotImplementedError(
                    " ".join(
                        [
                            "neoepiscope does not support version",
                            version,
                            "of netMHCcons",
                        ]
                    )
                )
        elif "netmhc" in program.lower():
            if "netMHC4" not in tool_dict and version == "4":
                program = paths.netMHC4
                if program is None:
                    program = which("netMHC4")
                else:
                    program = which(program)
                if program is None:
                    warnings.warn(
                        " ".join(
                            [
                                "No valid install of",
                                "netMHC version",
                                version,
                                "available",
                            ]
                        ),
                        Warning,
                    )
                    continue
                acceptable_scoring = ["rank", "affinity"]
                for method in scoring:
                    if method not in acceptable_scoring:
                        warnings.warn(
                            " ".join([method, "not compatible with netMHC"]),
                            Warning,
                        )
                        scoring.remove(method)
                if len(scoring) > 0:
                    if version == "4":
                        name = "netMHC4"
                    tool_dict[name] = [program, sorted(scoring)]
                else:
                    warnings.warn(
                            " ".join(
                                [
                                    "No compatible scoring methods given",
                                      "for NetMHC version", version,
                                      "- will not use this tool for",
                                      "binding predictions"
                                ]
                            ),
                            Warning,
                        )
            elif "netMHC4" in tool_dict and version == "4":
                raise RuntimeError(
                    "Conflicting or repetitive installs of netMHC given"
                )
            else:
                raise NotImplementedError(
                    " ".join(
                        [
                            "neoepiscope does not support version",
                            version,
                            "of netMHC",
                        ]
                    )
                )
        elif "pickpocket" in program.lower():
            if "pickpocket1" not in tool_dict and version == "1":
                program = paths.PickPocket1
                if program is None:
                    program = which("pickpocket1")
                else:
                    program = which(program)
                if program is None:
                    warnings.warn(
                        " ".join(
                            [
                                "No valid install of",
                                "PickPocket version",
                                version,
                                "available",
                            ]
                        ),
                        Warning,
                    )
                    continue
                acceptable_scoring = ["affinity"]
                for method in scoring:
                    if method not in acceptable_scoring:
                        warnings.warn(
                            " ".join([method, "not compatible with netMHC"]),
                            Warning,
                        )
                        scoring.remove(method)
                if len(scoring) > 0:
                    if version == "1":
                        name = "pickpocket1"
                    tool_dict[name] = [program, sorted(scoring)]
                else:
                    warnings.warn(
                            " ".join(
                                [
                                    "No compatible scoring methods given",
                                      "for PickPocket version", version,
                                      "- will not use this tool for",
                                      "binding predictions"
                                ]
                            ),
                            Warning,
                        )
            elif "pickpocket1" in tool_dict and version == "1":
                raise RuntimeError(
                    "Conflicting or repetitive installs of PickPocket given"
                )
            else:
                raise NotImplementedError(
                    " ".join(
                        [
                            "neoepiscope does not support version",
                            version,
                            "of PickPocket",
                        ]
                    )
                )
        else:
            raise NotImplementedError(
                " ".join(
                    [
                        "neoepiscope does not support",
                        program,
                        "for binding predictions",
                    ]
                )
            )
    return tool_dict


def get_affinity_netMHCIIpan(
    peptides, allele, netmhciipan, version, scores, remove_files=True
):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: Allele to use for binding affinity (string)
        netmhciipan: path to netMHCIIpan executable
        version: version of netMHCIIpan
        scores: list of scoring methods
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(
            os.path.join(
                os.path.join(neoepiscope_dir, "neoepiscope", "availableAlleles.pickle")
            ),
            "rb",
        ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        # Homogenize format
        if "DRB" in allele:
            allele = allele.replace("HLA-", "").replace(":", "").replace("*", "_")
        elif "DP" in allele or "DQ" in allele:
            allele = allele.replace(":", "").replace("*", "")
        if allele not in avail_alleles["".join(["netMHCIIpan", str(version)])]:
            warnings.warn(
                " ".join([allele, "is not a valid allele for netMHCIIpan"]), Warning
            )
            if len(scores) == 1:
                return [(peptides[i], "NA") for i in range(0, len(peptides))]
            else:
                return [(peptides[i], "NA", "NA") for i in range(0, len(peptides))]
        allele = allele.replace("*", "_").replace(":", "")
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "netmhciipan", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for
        #   input if peptide length is at least 9
        # Count instances of smaller peptides
        na_count = 0
        peptide_file = tempfile.mkstemp(
            suffix=".peptides", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(peptide_file)
        with open(peptide_file, "w") as f:
            for sequence in peptides:
                if len(sequence) >= 9:
                    print(sequence, file=f)
                else:
                    na_count += 1
        if na_count > 0:
            warnings.warn(
                " ".join(
                    [
                        str(na_count),
                        "peptides not compatible with",
                        "netMHCIIpan will not receive score",
                    ]
                ),
                Warning,
            )
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".netMHCIIpan.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        # Run netMHCIIpan
        subprocess.check_call(
            [
                netmhciipan,
                "-a",
                allele,
                "-inptype",
                "1",
                "-xls",
                "-xlsfile",
                mhc_out,
                "-f",
                peptide_file,
            ]
        )
        # Retrieve scores for valid peptides
        score_dict = {}
        with open(mhc_out, "r") as f:
            # Skip headers
            f.readline()
            f.readline()
            for line in f:
                # token 1 is peptide; token 4 is affinity; token[5] is rank
                tokens = line.strip("\n").split("\t")
                if sorted(scores) == ["affinity", "rank"]:
                    score_dict[tokens[1]] = (tokens[4], tokens[5])
                elif sorted(scores) == ["affinity"]:
                    score_dict[tokens[1]] = (tokens[4],)
                elif sorted(scores) == ["rank"]:
                    score_dict[tokens[1]] = (tokens[5],)
        # Produce list of scores for valid peptides
        # Invalid peptides receive "NA" score
        for sequence in peptides:
            if sequence in score_dict:
                nM = (sequence,) + score_dict[sequence]
            else:
                if len(scores) == 1:
                    nM = (sequence, "NA")
                else:
                    nM = (sequence, "NA", "NA")
            affinities.append(nM)
        return affinities
    finally:
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)


def get_affinity_mhcflurry(peptides, allele, scores, version, remove_files=True):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: Allele to use for binding affinity (string)
        scores: list of scoring methods
        version: version of mhc-flurry
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(
            os.path.join(neoepiscope_dir, "neoepiscope", "availableAlleles.pickle"),
            "rb",
        ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        if allele not in avail_alleles["mhcflurry"]:
            warnings.warn(
                " ".join([allele, "is not a valid allele for mhcflurry"]), Warning
            )
            score_form = tuple(["NA" for i in range(0, len(scores))])
            return [(peptides[i],) + score_form for i in range(0, len(peptides))]
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "mhcflurry", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for
        #   input if peptide length is at least 9
        # Count instances of smaller peptides
        # Establish temporary file to hold output
        peptide_file = tempfile.mkstemp(
            suffix=".csv", prefix="".join([sample_id, "."]), text=True
        )[1]
        na_count = 0
        with open(peptide_file, "w") as f:
            f.write("allele,peptide\n")
            for sequence in peptides:
                if len(sequence) < 8 or len(sequence) > 15:
                    na_count += 1
                else:
                    print("".join([allele, ",", sequence, "\n"]), file=f)
        if na_count > 0:
            warnings.warn(
                " ".join(
                    [
                        str(na_count),
                        "peptides not compatible with",
                        "mhcflurry will not receive score",
                    ]
                ),
                Warning,
            )
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".mhcflurry.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        # Run netMHCIIpan
        command = ["mhcflurry-predict", "--out", mhc_out, peptide_file]
        subprocess.check_call(command)
        # Retrieve scores for valid peptides
        score_dict = {}
        with open(mhc_out, "r") as f:
            # Skip headers
            f.readline()
            for line in f:
                # token 1 is peptide; token 4 is affinity; token[5] is rank
                tokens = line.strip("\n").split(",")
                if sorted(scores) == ["affinity", "rank"]:
                    score_dict[tokens[1]] = (tokens[2], tokens[5])
                elif sorted(scores) == ["affinity", "high", "low", "rank"]:
                    score_dict[tokens[1]] = (tokens[2], tokens[4], tokens[3], tokens[5])
                elif sorted(scores) == ["affinity", "high", "low"]:
                    score_dict[tokens[1]] = (tokens[2], tokens[4], tokens[3])
                elif sorted(scores) == ["affinity", "high", "rank"]:
                    score_dict[tokens[1]] = (tokens[2], tokens[4], tokens[5])
                elif sorted(scores) == ["affinity", "low", "rank"]:
                    score_dict[tokens[1]] = (tokens[2], tokens[3], tokens[5])
                elif sorted(scores) == ["high", "low", "rank"]:
                    score_dict[tokens[1]] = (tokens[4], tokens[3], tokens[3])
                elif sorted(scores) == ["affinity", "high"]:
                    score_dict[tokens[1]] = (tokens[2], tokens[4])
                elif sorted(scores) == ["affinity", "low"]:
                    score_dict[tokens[1]] = (tokens[2], tokens[3])
                elif sorted(scores) == ["high", "low"]:
                    score_dict[tokens[1]] = (tokens[4], tokens[3])
                elif sorted(scores) == ["high", "rank"]:
                    score_dict[tokens[1]] = (tokens[4], tokens[5])
                elif sorted(scores) == ["low", "rank"]:
                    score_dict[tokens[1]] = (tokens[3], tokens[5])
                elif sorted(scores) == ["affinity"]:
                    score_dict[tokens[1]] = (tokens[2],)
                elif sorted(scores) == ["rank"]:
                    score_dict[tokens[1]] = (tokens[5],)
                elif sorted(scores) == ["high"]:
                    score_dict[tokens[1]] = (tokens[4],)
                elif sorted(scores) == ["low"]:
                    score_dict[tokens[1]] = (tokens[3],)
        # Produce list of scores for valid peptides
        # Invalid peptides receive "NA" score
        for sequence in peptides:
            if sequence in score_dict:
                nM = (sequence,) + score_dict[sequence]
            else:
                nM = (sequence,) + tuple(["NA" for i in range(0, len(scores))])
            affinities.append(nM)
        return affinities
    finally:
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)


def get_affinity_netMHC(
    peptides, allele, netmhc, version, scores, remove_files=True
):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: allele to use for binding affinity
                    (string, format HLA-A02:01)
        netmhc: path to netMHC executable
        version: version of netMHC software
        scores: list of scoring methods
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(
            os.path.join(neoepiscope_dir, "neoepiscope", "availableAlleles.pickle"),
            "rb",
        ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        allele = allele.replace("*", "").replace(':', '')
        if allele not in avail_alleles["".join(["netMHC", str(version)])]:
            warnings.warn(
                " ".join([allele, "is not a valid allele for netMHC"]), Warning
            )
            if len(scores) == 1:
                return [(peptides[i], "NA") for i in range(0, len(peptides))]
            else:
                return [(peptides[i], "NA", "NA") for i in range(0, len(peptides))]
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "netmhc", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for input
        peptide_file = tempfile.mkstemp(
            suffix=".peptides", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(peptide_file)
        with open(peptide_file, "w") as f:
            for sequence in peptides:
                print(sequence, file=f)
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".netMHC.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        err_file = tempfile.mkstemp(
            suffix=".netMHC.err", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(err_file)
        with open(err_file, "w") as e:
            # Run netMHC
            subprocess.check_call(
                [
                    netmhc,
                    "-a",
                    allele,
                    "-inptype",
                    "1",
                    "-p",
                    "-xls",
                    "-xlsfile",
                    mhc_out,
                    peptide_file,
                ],
                stderr=e,
            )
        with open(mhc_out, "r") as f:
            f.readline()
            f.readline()
            for i in range(0, len(peptides)):
                tokens = f.readline().strip("\n").split("\t")
                # for v4, tokens[3] is affinty, tokens[4] is rank
                if sorted(scores) == ["affinity", "rank"]:
                    nM = (peptides[i], tokens[3], tokens[4])
                elif sorted(scores) == ["affinity"]:
                    nM = (peptides[i], tokens[3])
                elif sorted(scores) == ["rank"]:
                    nM = (peptides[i], tokens[4])
                affinities.append(nM)
        return affinities
    finally:
        # Remove temporary files
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)


def get_affinity_pickpocket(
    peptides, allele, pickpocket, version, scores, remove_files=True
):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: allele to use for binding affinity
                    (string, format HLA-A02:01)
        pickpocket: path to pickpocket executable
        version: version of pickpocket software
        scores: list of scoring methods
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    import math
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(
            os.path.join(neoepiscope_dir, "neoepiscope", "availableAlleles.pickle"),
            "rb",
        ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        allele = allele.replace("*", "")
        if allele not in avail_alleles["".join(["pickpocket", str(version)])]:
            warnings.warn(
                " ".join([allele, "is not a valid allele for PickPocket"]), Warning
            )
            return [(peptides[i], "NA") for i in range(0, len(peptides))]
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "pickpocket", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for input
        peptide_file = tempfile.mkstemp(
            suffix=".peptides", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(peptide_file)
        with open(peptide_file, "w") as f:
            for sequence in peptides:
                print(sequence, file=f)
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".pickpocket.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        err_file = tempfile.mkstemp(
            suffix=".pickpocket.err", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(err_file)
        with open(err_file, "w") as e:
            with open(mhc_out, "w") as o:
                # Run PickPocket
                subprocess.check_call(
                    [
                        pickpocket,
                        "-a",
                        allele,
                        "-inptype",
                        "1",
                        "-p",
                        peptide_file
                    ],
                    stderr=e, stdout=o
                )
        with open(mhc_out, "r") as f:
            first_char = "#"
            while first_char == "#":
                line = f.readline().strip()
                try:
                    first_char = line[0]
                except IndexError:
                    first_char = "#"
            for i in range(2):
                f.readline()
            for i in range(0, len(peptides)):
                tokens = line.strip().split()
                aff = 50000.0**((-1*float(tokens[4]))+1)
                nM = (peptides[i], str(aff))
                affinities.append(nM)
        return affinities
    finally:
        # Remove temporary files
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)


def get_affinity_netMHCcons(
    peptides, allele, netmhccons, version, scores, size_list, remove_files=True
):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: allele to use for binding affinity
                    (string, format HLA-A02:01)
        netmhccons: path to netMHCcons executable
        version: version of netMHCcons software
        scores: list of scoring methods
        size_list: list of [min size, ..., max size] of peptide sizes
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(
            os.path.join(neoepiscope_dir, "neoepiscope", "availableAlleles.pickle"),
            "rb",
        ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        allele = allele.replace("*", "")
        if allele not in avail_alleles["".join(["netMHCcons", str(version)])]:
            warnings.warn(
                " ".join([allele, "is not a valid allele for netMHCcons"]), Warning
            )
            if len(scores) == 1:
                return [(peptides[i], "NA") for i in range(0, len(peptides))]
            else:
                return [(peptides[i], "NA", "NA") for i in range(0, len(peptides))]
        # Establish score dict and return list
        score_dict = {}
        affinities = []
        # Get scores for peptides of each size
        for i in range(size_list[0], size_list[-1]+1):
            # Sample id
            sample_id = ".".join(
                [peptides[0], str(len(peptides)), allele, "netmhccons", version, str(i)]
            )
            # Write one peptide per line to a temporary file for input
            peptide_file = tempfile.mkstemp(
                suffix=".peptides", prefix="".join([sample_id, "."]), text=True
            )[1]
            files_to_remove.append(peptide_file)
            # Get peptides of correct size
            sized_peps = [x for x in peptides if len(x) == i]
            if len(sized_peps) == 0:
                continue
            with open(peptide_file, "w") as f:
                for sequence in sized_peps:
                    print(sequence, file=f)
            # Establish temporary file to hold output
            mhc_out = tempfile.mkstemp(
                suffix=".netMHCcons.out", prefix="".join([sample_id, "."]), text=True
            )[1]
            files_to_remove.append(mhc_out)
            err_file = tempfile.mkstemp(
                suffix=".netMHCcons.err", prefix="".join([sample_id, "."]), text=True
            )[1]
            files_to_remove.append(err_file)
            with open(err_file, "w") as e:
                # Run netMHC
                subprocess.check_call(
                    [
                        netmhccons,
                        "-a",
                        allele,
                        "-inptype",
                        "1",
                        "-xls",
                        "-xlsfile",
                        mhc_out,
                        "-f",
                        peptide_file,
                    ],
                    stderr=e,
                )
            with open(mhc_out, "r") as f:
                f.readline()
                f.readline()
                for j in range(0, len(sized_peps)):
                    tokens = f.readline().strip("\n").split("\t")
                    # for v4, tokens[4] is affinty, tokens[5] is rank
                    if sorted(scores) == ["affinity", "rank"]:
                        score_dict[peptides[j]] = (tokens[4], tokens[5])
                    elif sorted(scores) == ["affinity"]:
                        score_dict[peptides[j]] = (tokens[4],)
                    elif sorted(scores) == ["rank"]:
                        score_dict[peptides[j]] = (tokens[5],)
        for i in range(0, len(peptides)):
            if peptides[i] in score_dict:
                nM = (peptides[i], ) + score_dict[peptides[i]]
            else:
                if len(scores) == 2:
                    nM = (peptides[i], 'NA', 'NA')
                else:
                    nM = (peptides[i], 'NA')
            affinities.append(nM)
        return affinities
    finally:
        # Remove temporary files
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)


def get_affinity_netMHCpan(
    peptides, allele, netmhcpan, version, scores, remove_files=True
):
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
        with open(
            os.path.join(neoepiscope_dir, "neoepiscope", "availableAlleles.pickle"),
            "rb",
        ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        allele = allele.replace("*", "")
        if allele not in avail_alleles["".join(["netMHCpan", str(version)])]:
            warnings.warn(
                " ".join([allele, "is not a valid allele for netMHCpan"]), Warning
            )
            if len(scores) == 1:
                return [(peptides[i], "NA") for i in range(0, len(peptides))]
            else:
                return [(peptides[i], "NA", "NA") for i in range(0, len(peptides))]
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "netmhcpan", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for input
        peptide_file = tempfile.mkstemp(
            suffix=".peptides", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(peptide_file)
        with open(peptide_file, "w") as f:
            for sequence in peptides:
                print(sequence, file=f)
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".netMHCpan.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        err_file = tempfile.mkstemp(
            suffix=".netMHCpan.err", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(err_file)
        with open(err_file, "w") as e:
            # Run netMHCpan
            if version == "3":
                subprocess.check_call(
                    [
                        netmhcpan,
                        "-a",
                        allele,
                        "-inptype",
                        "1",
                        "-p",
                        "-xls",
                        "-xlsfile",
                        mhc_out,
                        peptide_file,
                    ],
                    stderr=e,
                )
            elif version == "4":
                subprocess.check_call(
                    [
                        netmhcpan,
                        "-BA",
                        "-a",
                        allele,
                        "-inptype",
                        "1",
                        "-p",
                        "-xls",
                        "-xlsfile",
                        mhc_out,
                        peptide_file,
                    ],
                    stderr=e,
                )
        with open(mhc_out, "r") as f:
            f.readline()
            f.readline()
            for i in range(0, len(peptides)):
                tokens = f.readline().strip("\n").split("\t")
                # for v3, tokens[5] is affinity, tokens[6] is rank
                # for v4, tokens[6] is affinty, tokens[7] is rank
                if sorted(scores) == ["affinity", "rank"]:
                    if version == "3":
                        nM = (peptides[i], tokens[5], tokens[6])
                    elif version == "4":
                        nM = (peptides[i], tokens[6], tokens[7])
                elif sorted(scores) == ["affinity"]:
                    if version == "3":
                        nM = (peptides[i], tokens[5])
                    elif version == "4":
                        nM = (peptides[i], tokens[6])
                elif sorted(scores) == ["rank"]:
                    if version == "3":
                        nM = (peptides[i], tokens[6])
                    elif version == "4":
                        nM = (peptides[i], tokens[7])
                affinities.append(nM)
        return affinities
    finally:
        # Remove temporary files
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)


def get_affinity_mhcnuggets(peptides, allele, version, remove_files=True):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: Allele to use for binding affinity (string)
        scores: list of scoring methods
        version: version of mhcnuggets
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    from mhcnuggets.src.predict import predict

    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(
            os.path.join(neoepiscope_dir, "neoepiscope", "availableAlleles.pickle"),
            "rb",
        ) as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        # Check that allele is valid for method
        allele = allele.replace("*", "")
        if allele in avail_alleles["mhcnuggets_mhcI"]:
            allele_class = "I"
            max_length = 15
        elif allele in avail_alleles["mhcnuggets_mhcII"]:
            allele_class = "II"
            max_length = 30
        else:
            warnings.warn(
                " ".join([allele, "is not a valid allele for mhcnuggets"]), Warning
            )
            return [(peptides[i], "NA") for i in range(0, len(peptides))]
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "mhcnuggets", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for
        #   input if peptide length is at least 9
        # Count instances of smaller peptides
        # Establish temporary file to hold output
        peptide_file = tempfile.mkstemp(
            suffix=".txt", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(peptide_file)
        na_count = 0
        with open(peptide_file, "w") as f:
            for sequence in peptides:
                if len(sequence) > max_length:
                    na_count += 1
                else:
                    print(sequence, file=f)
        if na_count > 0:
            warnings.warn(
                " ".join(
                    [
                        str(na_count),
                        "peptides not compatible with",
                        "mhcnuggets will not receive score",
                    ]
                ),
                Warning,
            )
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".mhcnuggets.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        # Run mhcnuggets
        predict(
            class_=allele_class, peptides_path=peptide_file, mhc=allele, output=mhc_out
        )
        # Retrieve scores for valid peptides
        score_dict = {}
        with open(mhc_out, "r") as f:
            # Skip headers
            f.readline()
            for line in f:
                tokens = line.strip("\n").split(",")
                score_dict[tokens[0]] = tokens[1]
        # Produce list of scores for valid peptides
        # Invalid peptides receive "NA" score
        for sequence in peptides:
            if sequence in score_dict:
                nM = (sequence, score_dict[sequence])
            else:
                nM = (sequence, "NA")
            affinities.append(nM)
        return affinities
    finally:
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)


def gather_binding_scores(neoepitopes, tool_dict, hla_alleles, size_list):
    """ Adds binding scores from desired programs to neoepitope metadata

        neoepitopes: dictionary linking neoepitopes to their metadata
        tool_dict: dictionary storing prediction tool data
        hla_alleles: list of HLA alleles used for binding predictions
        size_list: list of [min size, ..., max size] of peptide sizes

        Return value: dictionary linking neoepitopes to their metadata,
            which now includes binding scores
    """
    for allele in hla_alleles:
        for tool in sorted(tool_dict.keys()):
            if tool == "mhcflurry1":
                binding_scores = get_affinity_mhcflurry(
                    list(neoepitopes.keys()),
                    allele,
                    tool_dict[tool][1],
                    "1",
                    remove_files=True,
                )
            elif tool == "mhcnuggets2":
                binding_scores = get_affinity_mhcnuggets(
                    list(neoepitopes.keys()), allele, "2", remove_files=True
                )
            elif tool == "netMHCIIpan3":
                binding_scores = get_affinity_netMHCIIpan(
                    list(neoepitopes.keys()),
                    allele,
                    tool_dict[tool][0],
                    "3",
                    tool_dict[tool][1],
                    remove_files=True,
                )
            elif tool == "netMHCpan3":
                binding_scores = get_affinity_netMHCpan(
                    list(neoepitopes.keys()),
                    allele,
                    tool_dict[tool][0],
                    "3",
                    tool_dict[tool][1],
                    remove_files=True,
                )
            elif tool == "netMHCpan4":
                binding_scores = get_affinity_netMHCpan(
                    list(neoepitopes.keys()),
                    allele,
                    tool_dict[tool][0],
                    "4",
                    tool_dict[tool][1],
                    remove_files=True,
                )
            elif tool == "netMHC4":
                binding_scores = get_affinity_netMHC(
                    list(neoepitopes.keys()),
                    allele,
                    tool_dict[tool][0],
                    "4",
                    tool_dict[tool][1],
                    remove_files=True,
                )
            elif tool == "pickpocket1":
                binding_scores = get_affinity_pickpocket(
                    list(neoepitopes.keys()),
                    allele,
                    tool_dict[tool][0],
                    "1",
                    tool_dict[tool][1],
                    remove_files=True,
                )
            elif tool == "netMHCcons1":
                binding_scores = get_affinity_netMHCcons(
                    list(neoepitopes.keys()),
                    allele,
                    tool_dict[tool][0],
                    "1",
                    tool_dict[tool][1],
                    size_list,
                    remove_files=True,
                )
            for score in binding_scores:
                meta_data = neoepitopes[score[0]]
                for i in range(0, len(meta_data)):
                    neoepitopes[score[0]][i] = meta_data[i] + score[1:]
    return neoepitopes
