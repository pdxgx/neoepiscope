# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.

"""Code to work with the sprotXX.dat file from SwissProt.

https://web.expasy.org/docs/userman.html

Classes:
 - Record             Holds SwissProt data.

Functions:
 - read               Read one SwissProt record
 - parse              Read multiple SwissProt records

"""

import io
import re
import os
import pickle
from collections import defaultdict

neoepiscope_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

try:
    with open(
        os.path.join(os.path.join(neoepiscope_dir, "neoepiscope", "uniprotPTMreference.pickle")),
    ) as ptm_ref_stream:
        ptm_reference_dict = pickle.load(ptm_ref_stream)
except IOError as e:
    raise

class SwissProtParserError(ValueError):
    """An error occurred while parsing a SwissProt file."""

    def __init__(self, *args, line=None):
        """Create a SwissProtParserError object with the offending line."""
        super().__init__(*args)
        self.line = line


class Record:
    """Holds information from a SwissProt record.

    Attributes:
     - entry_name        Name of this entry, e.g. RL1_ECOLI
     - accession         Primary accession number (first in AC line) as string e.g. P00321
     - organism          The source of the sequence
     - comments          List of strings
     - isoforms          List of tuples (isoform_id, [vsp1, vsp2, ...])
     - iso_seqs          Defaultdict with:
                            key = isoform_id
                            value = List of variant sequence IDs
     - cross_references  List of tuples (database, id1[, id2][, id3])
     - tx_ids            List of tuples (transcript_id, isoform_id)
     - features          List of tuples (key name, from, to, description)
       from and to can be either integers for the residue
       numbers, '<', '>', or '?'
     - ptms              Defaultdict with:
                            key = isoform_id
                            value = List of tuples
                                (1-based location,
                                 MOD_RES or CARBOHYD,
                                 String of description)
     - var_seqs          List with:
                            [tuple of 1-based (from_location(inclusive), to_location(exclusive)),
                            {'note': String of alt sequence description,
                             'evidence': String of ECO codes|PubMed IDs,
                             'id': String of VSP ID}]
     - sequence          String of the canonical amino acid sequence

    Example:
    --------
    >>> example_filename = "SwissProt/P68308.txt"
    >>> species = 'BALPH'
    >>> with open(example_filename) as handle:
    ...     records = parse(handle, species)
    ...     for record in records:
    ...         print(record.entry_name)
    ...         print(record.accession)
    ...         print(record.organism)
    ...         print(record.sequence[:20] + "...")
    ...
    NU3M_BALPH
    P68308
    Balaenoptera physalus (Fin whale) (Balaena physalus).
    MNLLLTLLTNTTLALLLVFI...
    """

    def __init__(self):
        """Initialize the class."""
        self.entry_name = None
        self.accession = None
        self.organism = []
        self.comments = []
        self.isoforms = []
        self.cross_references = []
        self.tx_ids = []
        self.features = []
        self.ptms = None
        self.var_seqs = []
        self.sequence = ""

def parse(source, species):
    """Read multiple SwissProt records from file.

    Argument source is a file-like object or a path to a file.

    Returns a generator object which yields Bio.SwissProt.Record() objects.
    """
    handle = _open(source)
    try:
        while True:
            record = _read(handle, species)
            if not record:
                return
            yield record
    finally:
        if handle is not source:
            handle.close()

def read(source):
    """Read one SwissProt record from file.

    Argument source is a file-like object or a path to a file.

    Returns a Record() object.
    """
    handle = _open(source)
    try:
        record = _read(handle)
        if not record:
            raise ValueError("No SwissProt record found")
        # We should have reached the end of the record by now.
        # Try to read one more line to be sure:
        try:
            next(handle)
        except StopIteration:
            return record
        raise ValueError("More than one SwissProt record found")
    finally:
        if handle is not source:
            handle.close()

# Everything below is considered private

def _open(source):
    try:
        handle = open(source)
        return handle
    except TypeError:
        handle = source
        if handle.read(0) == "":
            # handle is text; assume the encoding is compatible with ASCII
            return handle
        # handle is binary; SwissProt encoding is always ASCII
        return io.TextIOWrapper(handle, encoding="ASCII")


def _read(handle, species):
    record = None
    unread = ""
    try:
        line = next(handle)
    except StopIteration:
        return record
    key, value = line[:2], line[5:].rstrip()
    if key != "ID":
        raise SwissProtParserError("Failed to find ID in first line", line=line)
    # Skip entries for unspecified organisms
    while species not in value:
        try:
            line = next(handle)
        except StopIteration:
            return record
        key, value = line[:2], line[5:].rstrip()
        if key == "ID" and species in value:
            break
    record = Record()
    #_read_id(record, line)
    _sequence_lines = []
    for line in handle:
        key, value = line[:2], line[5:].rstrip()
        if unread:
            value = unread + " " + value
            unread = ""
        if key == "AC":
            record.accession = value.rstrip(";").split("; ")[0]
        # Possible inclusion of record date/time info in future
        #elif key == "DT":
            #_read_dt(record, line)
        elif key == "OS":
            record.organism.append(value)
        elif key == "CC":
            _read_cc(record, line)
        elif key == "DR":
            _read_dr(record, value)
        elif key == "FT":
            read_ft(record, line)
        elif key == "SQ":
            cols = value.split()
            assert len(cols) == 7, f"I don't understand SQ line {line}"
            # Do more checking here?
            record.seqinfo = int(cols[1]), int(cols[3]), cols[5]
        elif key == "  ":
            _sequence_lines.append(value.replace(" ", "").rstrip())
        elif key == "//":
            # Join multiline data into one string
            record.organism = " ".join(record.organism)
            record.sequence = "".join(_sequence_lines)
            # Link isoform IDs to VSP IDs
            alts = [c for c in record.comments if 'ALTERNATIVE PRODUCTS:' in c]
            for c in alts:
                tmp_ids = re.findall('(?<=IsoId=)\w{6}-\d+', c)
                tmp_vsps = re.findall('Sequence=(.*?);', c)
                for iso_id, VSPs in zip(tmp_ids, tmp_vsps):
                    vsp_list = [VSP.strip() for VSP in re.split(",", VSPs)]
                    record.isoforms.append((iso_id, vsp_list))
            # Link ENSEMBL transcript IDs to isoform IDs
            txs = [i for i in record.cross_references if 'Ensembl' in i]
            for i in txs:
                try:
                    iso = re.search('(?<=\[)\w{6}-\d+', i[3]).group()
                # Outdated DR format with ENSEMBL lines missing canonical isoform ID
                # As of March 2023 only an issue with non-human entries
                except AttributeError:
                    try:
                        if record.isoforms[0][1] == ['Displayed']:
                            iso = record.isoforms[0][0]
                    except IndexError:
                        iso = None
                record.tx_ids.append((i[1], iso))
            # Keep PTMs for modified residues and glycosylation
            ptms = [p for p in record.features if p[0] == 'MOD_RES' or p[0] == 'CARBOHYD']
            if ptms:
                record.ptms = defaultdict(list)
            for p in ptms:
                # Add isoform ID for canonical sequence where possible otherwise iso_id = None
                if p[1] == None:
                    try:
                        p[1] = record.isoforms[0][0]
                    except IndexError:
                        pass
                descrip = p[3]['note']
                simple_descrip = descrip.split(';')[0]
                # Skipping PTMs from pathogens
                if simple_descrip.startswith('(Microbial infection)'):
                    pass
                else:
                    # PTMs format: {iso_id: [(1-based location, [type, base, PTM ID]), ...]}
                    record.ptms[p[1]].append((p[2][0], ptm_reference_dict[simple_descrip]))
                
            # Retain variant sequence information, see Record description for formatting
            record.var_seqs = [(v[2], v[3]) for v in record.features if v[0] == 'VAR_SEQ']
            return record
        elif key == "**":
            # Do this one last, as it will almost never occur.
            # See Bug 2353, some files from the EBI have extra lines
            # starting "**" (two asterisks/stars).  They appear
            # to be unofficial automated annotations. e.g.
            # **
            # **   #################    INTERNAL SECTION    ##################
            # **HA SAM; Annotated by PicoHamap 1.88; MF_01138.1; 09-NOV-2003.
            pass
        #else:
            #raise SwissProtParserError(f"Unknown keyword {key!r} found", line=line)
    if record:
        raise ValueError("Unexpected end of stream.")

def _read_id(record, line):
    cols = line[5:].split()
    # Prior to release 51, included with MoleculeType:
    # ID   EntryName DataClass; MoleculeType; SequenceLength AA.
    #
    # Newer files lack the MoleculeType:
    # ID   EntryName DataClass; SequenceLength AA.
    if len(cols) == 5:
        record.entry_name = cols[0]
        record.data_class = cols[1].rstrip(";")
        record.molecule_type = cols[2].rstrip(";")
        record.sequence_length = int(cols[3])
    elif len(cols) == 4:
        record.entry_name = cols[0]
        record.data_class = cols[1].rstrip(";")
        record.molecule_type = None
        record.sequence_length = int(cols[2])
    else:
        raise SwissProtParserError("ID line has unrecognised format", line=line)
    # check if the data class is one of the allowed values
    allowed = ("STANDARD", "PRELIMINARY", "IPI", "Reviewed", "Unreviewed")
    if record.data_class not in allowed:
        message = f"Unrecognized data class {record.data_class!r}"
        raise SwissProtParserError(message, line=line)

    # molecule_type should be 'PRT' for PRoTein
    # Note that has been removed in recent releases (set to None)
    if record.molecule_type not in (None, "PRT"):
        message = f"Unrecognized molecule type {record.molecule_type!r}"
        raise SwissProtParserError(message, line=line)


def _read_dt(record, line):
    value = line[5:]
    uprline = value.upper()
    cols = value.rstrip().split()
    if (
        "CREATED" in uprline
        or "LAST SEQUENCE UPDATE" in uprline
        or "LAST ANNOTATION UPDATE" in uprline
    ):
        # Old style DT line
        # =================
        # e.g.
        # DT   01-FEB-1995 (Rel. 31, Created)
        # DT   01-FEB-1995 (Rel. 31, Last sequence update)
        # DT   01-OCT-2000 (Rel. 40, Last annotation update)
        #
        # or:
        # DT   08-JAN-2002 (IPI Human rel. 2.3, Created)
        # ...

        # find where the version information will be located
        # This is needed for when you have cases like IPI where
        # the release version is in a different spot:
        # DT   08-JAN-2002 (IPI Human rel. 2.3, Created)
        uprcols = uprline.split()
        rel_index = -1
        for index in range(len(uprcols)):
            if "REL." in uprcols[index]:
                rel_index = index
        assert rel_index >= 0, f"Could not find Rel. in DT line: {line}"
        version_index = rel_index + 1
        # get the version information
        str_version = cols[version_index].rstrip(",")
        # no version number
        if str_version == "":
            version = 0
        # dot versioned
        elif "." in str_version:
            version = str_version
        # integer versioned
        else:
            version = int(str_version)
        date = cols[0]

        if "CREATED" in uprline:
            record.created = date, version
        elif "LAST SEQUENCE UPDATE" in uprline:
            record.sequence_update = date, version
        elif "LAST ANNOTATION UPDATE" in uprline:
            record.annotation_update = date, version
        else:
            raise SwissProtParserError("Unrecognised DT (DaTe) line", line=line)
    elif (
        "INTEGRATED INTO" in uprline
        or "SEQUENCE VERSION" in uprline
        or "ENTRY VERSION" in uprline
    ):
        # New style DT line
        # =================
        # As of UniProt Knowledgebase release 7.0 (including
        # Swiss-Prot release 49.0 and TrEMBL release 32.0) the
        # format of the DT lines and the version information
        # in them was changed - the release number was dropped.
        #
        # For more information see bug 1948 and
        # http://ca.expasy.org/sprot/relnotes/sp_news.html#rel7.0
        #
        # e.g.
        # DT   01-JAN-1998, integrated into UniProtKB/Swiss-Prot.
        # DT   15-OCT-2001, sequence version 3.
        # DT   01-APR-2004, entry version 14.
        #
        # This is a new style DT line...

        # The date should be in string cols[1]
        # Get the version number if there is one.
        # For the three DT lines above: 0, 3, 14
        try:
            version = 0
            for s in cols[-1].split("."):
                if s.isdigit():
                    version = int(s)
        except ValueError:
            version = 0
        date = cols[0].rstrip(",")

        # Re-use the historical property names, even though
        # the meaning has changed slightly:
        if "INTEGRATED" in uprline:
            record.created = date, version
        elif "SEQUENCE VERSION" in uprline:
            record.sequence_update = date, version
        elif "ENTRY VERSION" in uprline:
            record.annotation_update = date, version
        else:
            raise SwissProtParserError("Unrecognised DT (DaTe) line", line=line)
    else:
        raise SwissProtParserError("Failed to parse DT (DaTe) line", line=line)

def _read_cc(record, line):
    key, value = line[5:8], line[9:].rstrip()
    if key == "-!-":  # Make a new comment
        record.comments.append(value)
    elif key == "   ":  # add to the previous comment
        if not record.comments:
            # TCMO_STRGA in Release 37 has comment with no topic
            record.comments.append(value)
        else:
            record.comments[-1] += " " + value

def _read_dr(record, value):
    cols = value.rstrip(".").split("; ")
    record.cross_references.append(tuple(cols))

def read_ft(record, line):
    name = line[5:13].rstrip()
    if name:
        # originally 1-based
        location = line[21:80].rstrip()
        try:
            isoform_id, location = location.split(':')
        except ValueError:
            isoform_id = None
        # strip location, keep 1-based
        try:
            from_res, to_res = map(int, location.split('..'))
        except ValueError:
            try:
                from_res = int(location) 
                to_res = int(location) + 1 
            except ValueError:
                # if location is inexact we can't map from protein to genome
                from_res = '?'
                to_res = '?'
        # {} is placeholder for qualifiers
        feature = [name, isoform_id, (from_res, to_res), {}]
        record.features.append(feature)
        return

    # continuation of the previous feature
    elif line[5:21] == "                ":
        feature = record.features[-1]
        value = line[21:].rstrip()
        match = re.match(r"^/([a-z_]+)=", value)
        if match:
            qualifier_type = match.group(1)
            value = value[len(match.group(0)) :]
            if not value.startswith('"'):
                raise ValueError("Missing starting quote in feature")
            if qualifier_type == "id":
                if not value.endswith('"'):
                    raise ValueError("Missing closing quote for id")
                # index 3 of feature is the qualifier dict
                feature[3][qualifier_type] = value[1:-1]
            else:
                if value.endswith('"'):
                    value = value[1:-1]
                else:  # continues on the next line
                    value = value[1:]
                if qualifier_type in feature[3].keys():
                    raise ValueError(
                        f"Feature qualifier {qualifier_type!r} already exists for feature"
                    )
                feature[3][qualifier_type] = value
            return
        
        # this line is a continuation of the description of the previous feature
        keys = list(feature[3].keys())
        key = keys[-1]
        description = value.rstrip('"')
        old_description = feature[3][key]
        if key == "evidence" or old_description.endswith("-"):
            description = f"{old_description}{description}"
        else:
            description = f"{old_description} {description}"
        # index 0 of feature is name
        if feature[0] == "VAR_SEQ":
            try:
                first_seq, second_seq = description.split(" -> ")
            except ValueError:
                pass
            else:
                extra_info = ""
                # we might have more information at the end of the
                # second sequence, which should be in parenthesis
                extra_info_pos = second_seq.find(" (")
                if extra_info_pos != -1:
                    extra_info = second_seq[extra_info_pos:]
                    second_seq = second_seq[:extra_info_pos]
                # now clean spaces out of the first and second string
                first_seq = first_seq.replace(" ", "")
                second_seq = second_seq.replace(" ", "")
                # reassemble the description
                description = first_seq + " -> " + second_seq + extra_info
        feature[3][key] = description

if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
