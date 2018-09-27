#!/usr/bin/env python
# coding=utf-8
"""
download.py

Part of neoepiscope
Downloads dependencies according to user input to facilitate installation.

Licensed under the MIT license.

The MIT License (MIT)
Copyright (c) 2018 Mary A. Wood, Austin Nguyen,
                   Abhinav Nellore, and Reid Thompson

Adapted from Rail-RNA, which is copyright (c) 2015 
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
from six.moves import input
from .version import version_number
from .file_processing import which
import signal
import shutil
import tempfile
import sys
import os
import subprocess
from .transcript import gtf_to_cds
from .transcript import cds_to_tree
from distutils.core import Command

download = {
    "Gencode v27 annotation": [
        "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human"
        "/release_27/gencode.v27.annotation.gtf.gz"
    ],
    "Gencode v19 annotation": [
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/"
        "Gencode_human/release_19/"
        "gencode.v19.annotation.gtf.gz"
    ],
    "Bowtie NCBI GRCh38 index": [
        "ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip",
        "http://verve.webfactional.com/GRCh38_no_alt.zip",
    ],
    "Bowtie UCSC hg19 index": [
        "ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip",
        "http://verve.webfactional.com/hg19.ebwt.zip",
    ],
}


def is_exe(fpath):
    """ Tests whether a file is executable.
        fpath: path to file
        Return value: True iff file exists and is executable.
    """
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)


def remove_temporary_directories(temp_dir_paths):
    """ Deletes temporary directory.

        temp_dir_paths: iterable of paths of temporary directories

        No return value.
    """
    for temp_dir_path in temp_dir_paths:
        try:
            shutil.rmtree(temp_dir_path, ignore_errors=True)
        except Exception as e:
            # Don't know what's up, but forge on
            pass


def sig_handler(signum, frame):
    """ Helper function for register_cleanup that's called on signal. """
    import sys

    sys.exit(0)


def register_cleanup(handler, *args, **kwargs):
    """ Registers cleanup on normal and signal-induced program termination.

        Executes previously registered handler as well as new handler.

        handler: function to execute on program termination
        args: named arguments of handler
        kwargs includes keyword args of handler as well as:
            signals_to_handle: list of signals to handle, e.g. [signal.SIGTERM,
                signal.SIGHUP]

        No return value.
    """
    if "signals_to_handle" in kwargs:
        signals_to_handle = kwargs["signals_to_handle"]
        del kwargs["signals_to_handle"]
    else:
        signals_to_handle = [signal.SIGTERM, signal.SIGHUP]
    from atexit import register

    register(handler, *args, **kwargs)
    old_handlers = [
        signal.signal(a_signal, sig_handler) for a_signal in signals_to_handle
    ]
    for i, old_handler in enumerate(old_handlers):
        if (old_handler != signal.SIG_DFL) and (old_handler != sig_handler):

            def new_handler(signum, frame):
                try:
                    sig_handler(signum, frame)
                finally:
                    old_handler(signum, frame)

        else:
            new_handler = sig_handler
        signal.signal(signals_to_handle[i], new_handler)


def print_to_screen(message, newline=True, carriage_return=False):
    """ Prints message to stdout as well as stderr if stderr is redirected.

        message: message to print
        newline: True iff newline should be printed
        carriage_return: True iff carriage return should be printed; also
            clears line with ANSI escape code

        No return value.
    """
    full_message = (
        "\x1b[K"
        + message
        + ("\r" if carriage_return else "")
        + ("\n" if newline else "")
    )
    try:
        sys.stderr.write(full_message)
        if sys.stderr.isatty():
            sys.stderr.flush()
        else:
            try:
                # So the user sees it too
                sys.stdout.write(full_message)
                sys.stdout.flush()
            except UnicodeEncodeError:
                sys.stdout.write(
                    unicodedata.normalize("NFKD", full_message).encode(
                        "ascii", "ignore"
                    )
                )
                sys.stdout.flush()
    except UnicodeEncodeError:
        sys.stderr.write(
            unicodedata.normalize("NFKD", full_message).encode("ascii", "ignore")
        )
        sys.stderr.flush()


class NeoepiscopeDownloader(object):
    """ Convenience class for downloading files so neoepiscope is ready to use.

        Init vars
        -------------
        curl_exe: path to cURL executable; if None, use 'curl'
        download_dir: where to put files; if None, use home directory
            + 'neoepiscope.data'
        yes: if True, answer yes to all questions automatically
        print_log_on_error: print log on error
    """

    def __init__(
        self,
        curl_exe=None,
        download_dir=None,
        yes=False,
        print_log_on_error=False,
        version_number=version_number,
    ):
        print_to_screen("""Configuring neoepiscope v{0} ...""".format(version_number))
        self.download_dir = download_dir
        self.curl_exe = curl_exe
        log_dir = tempfile.mkdtemp()
        self.log_file = os.path.join(log_dir, "neoepiscope_config.err")
        self.log_stream = open(self.log_file, "w")
        self.finished = False
        register_cleanup(remove_temporary_directories, [log_dir])
        self.print_log_on_error = print_log_on_error
        self.yes = yes

    def __enter__(self):
        return self

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def _print_to_screen_and_log(self, message, **kwargs):
        print(message, file=self.log_stream)
        print_to_screen(message, **kwargs)

    def _bail(self):
        """ Copy log to some temporary dir and GTFO. """
        new_log_file = os.path.join(tempfile.mkdtemp(), "neoepiscope_config.err")
        shutil.copyfile(self.log_file, new_log_file)
        if self.print_log_on_error:
            print_to_screen("Log (also at %s):" % new_log_file)
            print_to_screen("===")
            with open(new_log_file, "r") as log_fh:
                for ln in log_fh:
                    print_to_screen(ln.rstrip())
            print_to_screen("===")
        else:
            print_to_screen("Installation log may be found at %s." % new_log_file)
        sys.exit(1)

    def _yes_no_query(self, question, answer=None):
        """ Gets a yes/no answer from the user if self.yes is not True.

            question: string with question to be printed to console
            answer: boolean that overrides self.yes if an answer should
                be forced

            Return value: boolean
        """
        from distutils.util import strtobool

        if answer is None:
            if self.yes:
                print("\x1b[K%s [y/n]: y" % question)
                return True
        elif answer:
            print("\x1b[K%s [y/n]: y" % question)
            return True
        else:
            print("\x1b[K%s [y/n]: n" % question)
            return False
        while True:
            sys.stdout.write("\x1b[K%s [y/n]: " % question)
            try:
                try:
                    return strtobool(input().lower())
                except KeyboardInterrupt:
                    sys.stdout.write("\n")
                    sys.exit(0)
            except ValueError:
                sys.stdout.write("Please enter 'y' or 'n'.\n")

    def _request_path(self, request, program="", use_which=True):
        """ Gets a path from a user to an installed software

            request: string with request to be printed to console
            use_which: whether to use the which function to verify a path

            Return value: string with response to request
        """
        if use_which:
            which_function = lambda x: which(x)
        else:
            which_function = lambda x: x
        while True:
            sys.stdout.write("%s: " % request)
            try:
                path = which_function(input())
                if path is not None:
                    return path
                else:
                    raise ValueError("".join(["Not a valid install of ", program]))
            except KeyboardInterrupt:
                sys.stdout.write("\n")
                sys.exit(0)

    def check_exe(self, program):
        """ Tests whether an executable is in PATH.
            program: executable to search for
            Return value: path to executable or None if the executable is not
                found.
        """

        def ext_candidates(fpath):
            yield fpath
            for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
                yield fpath + ext

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                for candidate in ext_candidates(exe_file):
                    if is_exe(candidate):
                        return candidate
        return None

    def _grab_and_explode(self, urls, name, explode=True):
        """ Special method for grabbing and exploding a package, if necessary.

            Does not verify URLs since these are preverified. Package is
            downloaded to current directory.

            url: list of urls to grab
            name: name of download
            explode: True iff downloaded file should be decompressed

            No return value
        """
        from collections import deque

        self._print_to_screen_and_log(
            "[Configuring] Downloading %s..." % name,
            newline=False,
            carriage_return=True,
        )
        url_deque = deque(urls)
        while url_deque:
            url = url_deque.popleft()
            """Follow redirects (-L), write to file (-O), respect
            content disposition"""
            command = [self.curl_exe, "-L", "-O", url]
            filename = url.rpartition("/")[2]
            print("[Configuring] Downloading {} from {}...".format(name, url))
            try:
                subprocess.check_output(command, stderr=self.log_stream)
            except subprocess.CalledProcessError as e:
                if not url_deque:
                    self._print_to_screen_and_log(
                        (
                            "Error encountered downloading file %s; exit "
                            'code was %d; command invoked was "%s".'
                        )
                        % (url, e.returncode, " ".join(command))
                    )
                    self._print_to_screen_and_log("Make sure web is accessible.")
                    self._bail()
                else:
                    self._print_to_screen_and_log(
                        "[Configuring] Download failed; "
                        "trying alternate URL for %s..." % name,
                        newline=False,
                        carriage_return=True,
                    )
            else:
                if explode:
                    explode_command = None
                    if url[-8:] == ".tar.bz2":
                        explode_command = ["tar", "xvjf", filename]
                    elif url[-7:] == ".tar.gz" or url[-4:] == ".tgz":
                        explode_command = ["tar", "xvzf", filename]
                    elif url[-3:] == ".gz":
                        explode_command = ["gunzip", filename]
                    elif url[-4:] == ".zip":
                        explode_command = ["unzip", filename]
                    if explode_command is not None:
                        self._print_to_screen_and_log(
                            "[Configuring] Extracting %s..." % name,
                            newline=False,
                            carriage_return=True,
                        )
                        try:
                            subprocess.check_output(explode_command)
                        except subprocess.CalledProcessError as e:
                            if not url_deque:
                                self._print_to_screen_and_log(
                                    (
                                        "Error encountered exploding file %s; exit "
                                        'code was %d; command invoked was "%s".'
                                    )
                                    % (
                                        filename,
                                        e.returncode,
                                        " ".join(explode_command),
                                    )
                                )
                                self._bail()
                            else:
                                self._print_to_screen_and_log(
                                    "[Configuring] Extraction failed; "
                                    "trying alternate URL for %s..." % name,
                                    newline=False,
                                    carriage_return=True,
                                )
                                continue
                        finally:
                            try:
                                os.remove(filename)
                            except OSError:
                                break
                break

    def _quote(self, path=None):
        """ Decks path with single quotes if it's not None

            path: path or None

            Return value: None if path is None, else path decked by single
                quotes
        """
        if path is None:
            return "None"
        return "'%s'" % path

    def run(self):
        """ Downloads neoepiscope dependencies. """
        import multiprocessing

        if self.curl_exe is None:
            self.curl_exe = self.check_exe("curl")
            if self.curl_exe is None:
                print_to_screen(
                    "Configuring neoepiscope requires Curl if "
                    "dependencies are to be installed. "
                    "Download it at "
                    "http://curl.haxx.se/download.html and use "
                    "--curl to specify its path, or "
                    "disable installing dependencies with "
                    "--no-dependencies."
                )
                sys.exit(1)
            else:
                self.curl_exe = self.curl_exe.replace(" --help", "")
        self._print_to_screen_and_log("[Configuring] Downloading mhcflurry data...")
        subprocess.call(["mhcflurry-downloads", "fetch"])
        if self.download_dir is None:
            self.download_dir = os.path.join(
                os.path.expanduser("~"), "neoepiscope.data"
            )
            if not self._yes_no_query(
                "Install neoepiscope data in {}?".format(self.download_dir)
            ):
                self.download_dir = os.path.abspath(
                    os.path.expanduser(
                        self._request_path(
                            "Please enter a directory for storing neoepiscope data",
                            use_which=False,
                        )
                    )
                )
        # Download to a temporary directory first, then move to final dest
        temp_install_dir = tempfile.mkdtemp()
        register_cleanup(remove_temporary_directories, [temp_install_dir])
        if os.path.exists(self.download_dir):
            if self._yes_no_query(
                (
                    "The download directory path {dir} already exists.\n    "
                    "* Overwrite {dir}?"
                ).format(dir=self.download_dir)
            ):
                try:
                    shutil.rmtree(self.download_dir)
                except OSError:
                    # Handle this later if directory creation fails
                    pass
                try:
                    os.remove(self.download_dir)
                except OSError:
                    pass
            else:
                print_to_screen("Specify a different download directory")
                sys.exit(0)
        try:
            os.makedirs(self.download_dir)
        except OSError as e:
            self._print_to_screen_and_log(
                (
                    "Problem encountered trying to create "
                    "directory %s for installation. May need "
                    "sudo permissions."
                )
                % self.download_dir
            )
            self._bail()
        else:
            # So it's possible to move temp installation dir there
            os.rmdir(self.download_dir)
            pass
        os.chdir(temp_install_dir)
        if self._yes_no_query("Download Gencode v27 gtf annotation file?"):
            self._grab_and_explode(
                download["Gencode v27 annotation"],
                "Gencode v27 annotation",
                explode=False,
            )
            gencode_v27_temp = os.path.join(temp_install_dir, "gencode_v27")
            gencode_v27 = os.path.join(self.download_dir, "gencode_v27")
            gencode_v27_gtf = os.path.join(
                temp_install_dir,
                os.path.basename(download["Gencode v27 annotation"][0]),
            )
            try:
                os.makedirs(gencode_v27_temp)
            except OSError as e:
                self._print_to_screen_and_log(
                    (
                        "Problem encountered trying to create "
                        "directory %s for installation. May need "
                        "sudo permissions."
                    )
                    % gencode_v27_temp
                )
                self._bail()
            self._print_to_screen_and_log("[Configuring] Indexing Gencode v27...")
            cds_dict = gtf_to_cds(gencode_v27_gtf, gencode_v27_temp)
            cds_to_tree(cds_dict, gencode_v27_temp)
        else:
            gencode_v27 = None
        if self._yes_no_query("Download Gencode v19 gtf annotation file?"):
            self._grab_and_explode(
                download["Gencode v19 annotation"],
                "Gencode v19 annotation",
                explode=False,
            )
            gencode_v19_temp = os.path.join(temp_install_dir, "gencode_v19")
            gencode_v19 = os.path.join(self.download_dir, "gencode_v19")
            gencode_v19_gtf = os.path.join(
                temp_install_dir,
                os.path.basename(download["Gencode v19 annotation"][0]),
            )
            try:
                os.makedirs(gencode_v19_temp)
            except OSError as e:
                self._print_to_screen_and_log(
                    (
                        "Problem encountered trying to create "
                        "directory %s for installation. May need "
                        "sudo permissions."
                    )
                    % gencode_v19_temp
                )
                self._bail()
            self._print_to_screen_and_log("[Configuring] Indexing Gencode v19...")
            cds_dict = gtf_to_cds(gencode_v19_gtf, gencode_v19_temp)
            cds_to_tree(cds_dict, gencode_v19_temp)
        else:
            gencode_v19 = None
        if self._yes_no_query("Download Bowtie NCBI GRCh38 index?"):
            self._grab_and_explode(
                download["Bowtie NCBI GRCh38 index"], "Bowtie NCBI GRCh38 index"
            )
            bowtie_grch38 = os.path.join(
                self.download_dir, "GCA_000001405.15_GRCh38_no_alt_analysis_set"
            )
        else:
            bowtie_grch38 = None
        if self._yes_no_query("Download Bowtie UCSC hg19 index?"):
            self._grab_and_explode(
                download["Bowtie UCSC hg19 index"], "Bowtie UCSC hg19 index"
            )
            bowtie_hg19 = os.path.join(self.download_dir, "hg19")
        else:
            bowtie_hg19 = None
        programs = []
        for program in ["netMHCIIpan v3", "netMHCpan v3", "netMHCpan v4"]:
            if self._yes_no_query(
                (
                    "Do you have an install of {} "
                    "you would like to use for "
                    "binding score predictions with "
                    "neoepiscope?"
                ).format(program)
            ):
                exe = self._request_path(
                    "Please enter the path to your {} executable".format(program),
                    program,
                )
                if is_exe(exe):
                    programs.append(exe)
                else:
                    self._print_to_screen_and_log(
                        (
                            "The path to %s is not an executable. Please "
                            "verify that this is the correct path."
                        )
                        % program
                    )
                    self._bail()
            else:
                programs.append(None)
        if self._yes_no_query(
            "Do you have an install of HapCUT2 "
            "you would like to integrate with neoepiscope?"
        ):
            hapcut_build_dir = self._request_path(
                "Please enter the path to your HapCUT2 build directory", use_which=False
            )
            if is_exe(os.path.join(hapcut_build_dir, "extractHAIRS")) and is_exe(
                os.path.join(hapcut_build_dir, "HAPCUT2")
            ):
                programs.append(os.path.join(hapcut_build_dir, "extractHAIRS"))
                programs.append(os.path.join(hapcut_build_dir, "HAPCUT2"))
            else:
                self._print_to_screen_and_log(
                    "This build directory does not contain the "
                    "relevant HapCUT2 script - please check "
                    "that you have entered the correct path."
                )
                self._bail()
        else:
            for i in range(0, 2):
                programs.append(None)
        # Write paths to exe_paths
        with open(os.path.join(temp_install_dir, "paths.py"), "w") as paths_stream:
            print(
                (
                    """\"""
paths.py
Part of neoepiscope
Defines default paths of neoepiscope's dependencies
None indicates the user didn't install the tool or data
\"""

gencode_v27 = {gencode_v27}
gencode_v19 = {gencode_v19}
bowtie_grch38 = {bowtie_grch38}
bowtie_hg19 = {bowtie_hg19}
netMHCIIpan3 = {netMHCIIpan3}
netMHCpan3 = {netMHCpan3}
netMHCpan4 = {netMHCpan4}
hapcut2_hairs = {hapcut2_hairs}
hapcut2 = {hapcut2}
"""
                ).format(
                    gencode_v27=(
                        "None" if gencode_v27 is None else self._quote(gencode_v27)
                    ),
                    gencode_v19=(
                        "None" if gencode_v19 is None else self._quote(gencode_v19)
                    ),
                    bowtie_grch38=(
                        "None" if bowtie_grch38 is None else self._quote(bowtie_grch38)
                    ),
                    bowtie_hg19=(
                        "None" if bowtie_hg19 is None else self._quote(bowtie_hg19)
                    ),
                    netMHCIIpan3=(
                        "None" if programs[0] is None else self._quote(programs[0])
                    ),
                    netMHCpan3=(
                        "None" if programs[1] is None else self._quote(programs[1])
                    ),
                    netMHCpan4=(
                        "None" if programs[2] is None else self._quote(programs[2])
                    ),
                    hapcut2_hairs=(
                        "None" if programs[3] is None else self._quote(programs[3])
                    ),
                    hapcut2=(
                        "None" if programs[4] is None else self._quote(programs[4])
                    ),
                ),
                file=paths_stream,
            )
        # Move to final directory
        try:
            if os.path.isdir(self.download_dir):
                shutil.rmtree(self.download_dir)
            shutil.copytree(temp_install_dir, self.download_dir)
            shutil.rmtree(temp_install_dir)
            if os.path.isfile(
                os.path.join(
                    sys.prefix,
                    "lib",
                    "".join(
                        [
                            "python",
                            str(sys.version_info.major),
                            ".",
                            str(sys.version_info.minor),
                        ]
                    ),
                    "site-packages",
                    "neoepiscope",
                    "paths.py",
                )
            ):
                os.remove(
                    os.path.join(
                        sys.prefix,
                        "lib",
                        "".join(
                            [
                                "python",
                                str(sys.version_info.major),
                                ".",
                                str(sys.version_info.minor),
                            ]
                        ),
                        "site-packages",
                        "neoepiscope",
                        "paths.py",
                    )
                )
            shutil.copy2(
                os.path.join(self.download_dir, "paths.py"),
                os.path.join(
                    sys.prefix,
                    "lib",
                    "".join(
                        [
                            "python",
                            str(sys.version_info.major),
                            ".",
                            str(sys.version_info.minor),
                        ]
                    ),
                    "site-packages",
                    "neoepiscope",
                ),
            )
            os.remove(os.path.join(self.download_dir, "paths.py"))
        except Exception as e:
            self._print_to_screen_and_log(
                (
                    'Problem "%s" encountered moving '
                    "temporary installation directory "
                    "%s to final destination %s."
                )
                % (e, temp_install_dir, self.download_dir)
            )
            self._bail()
        self._print_to_screen_and_log("Configured neoepiscope.")
        self.finished = True
