from version import version_number
from file_processing import which
import signal
import shutil
import dependency_urls
import tempfile

download = {
            'Gencode v27 annotation': [
                            'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human'
                            '/release_27/gencode.v27.annotation.gtf.gz'
                            ],
            'Gencode v19 annotation': [
                            'ftp://ftp.ebi.ac.uk/pub/databases/gencode/'
                            'Gencode_human/release_19/'
                            'gencode.v19.annotation.gtf.gz'
                            ],
            'Bowtie NCBI GRCh38 index': [
                            'ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/'
                            'GRCh38_no_alt.zip'
                            ],
            'Bowtie UCSC hg19 index': [
                            'ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/'
                            'hg19.ebwt.zip'
                            ],
            'HapCUT2': [
                            'https://github.com/vibansal/HapCUT2/archive/'
                            'cf8ce6d70be9412e96a8a9364cb1aa4556db4c92.zip'
                            ]
            }

def remove_temporary_directories(temp_dir_paths):
    """ Deletes temporary directory.

        temp_dir_paths: iterable of paths of temporary directories

        No return value.
    """
    for temp_dir_path in temp_dir_paths:
        try:
            shutil.rmtree(temp_dir_path,
                            ignore_errors=True)
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
    if 'signals_to_handle' in kwargs:
        signals_to_handle = kwargs['signals_to_handle']
        del kwargs['signals_to_handle']
    else:
        signals_to_handle = [signal.SIGTERM, signal.SIGHUP]
    from atexit import register
    register(handler, *args, **kwargs)
    old_handlers = [signal.signal(a_signal, sig_handler)
                    for a_signal in signals_to_handle]
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
    full_message = ('\x1b[K' + message + ('\r' if carriage_return else '')
                        + ('\n' if newline else ''))
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
                                unicodedata.normalize(
                                        'NFKD', full_message
                                    ).encode('ascii', 'ignore')
                            )
                sys.stdout.flush()
    except UnicodeEncodeError:
        sys.stderr.write(
                        unicodedata.normalize(
                                'NFKD', full_message
                            ).encode('ascii', 'ignore')
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

    def __init__(self, curl_exe=None, download_dir=None,
                    yes=False, print_log_on_error=False):
        print_to_screen(u"""neoepiscope v{1} configuration""".format(
                                        version_number)
                                    )
        self.download_dir = download_dir
        self.curl_exe = curl_exe
        log_dir = tempfile.mkdtemp()
        self.log_file = os.path.join(log_dir, 'neoepiscope_config.err')
        self.log_stream = open(self.log_file, 'w')
        self.finished = False
        register_cleanup(remove_temporary_directories, [log_dir])
        self.print_log_on_error = print_log_on_error

    def __enter__(self):
        return self

    def _print_to_screen_and_log(self, message, **kwargs):
        print >>self.log_stream, message
        print_to_screen(message, **kwargs)

    def _bail(self):
        """ Copy log to some temporary dir and GTFO. """
        new_log_file = os.path.join(tempfile.mkdtemp(),
                                            'neoepiscope_config.err')
        shutil.copyfile(self.log_file, new_log_file)
        if self.print_log_on_error:
            print_to_screen('Log (also at %s):' % new_log_file)
            print_to_screen('===')
            with open(new_log_file, 'r') as log_fh:
                for ln in log_fh:
                    print_to_screen(ln.rstrip())
            print_to_screen('===')
        else:
            print_to_screen('Installation log may be found at %s.'
                            % new_log_file)
        sys.exit(1)

    def _yes_no_query(self, question, answer=None):
        """ Gets a yes/no answer from the user if self.yes is not True.

            question: string with question to be printed to console
            answer: boolean that overrides self.yes if an answer should
                be forced

            Return value: boolean
        """
        if answer is None:
            if self.yes:
                print '%s [y/n]: y' % question
                return True
        elif answer:
            print '%s [y/n]: y' % question
            return True
        else:
            print '%s [y/n]: n' % question
            return False
        while True:
            sys.stdout.write('%s [y/n]: ' % question)
            try:
                try:
                    return strtobool(raw_input().lower())
                except KeyboardInterrupt:
                    sys.stdout.write('\n')
                    sys.exit(0)
            except ValueError:
                sys.stdout.write('Please enter \'y\' or \'n\'.\n')

    def _request_path(self, request, software_name):
        """ Gets a path from a user to an installed software

            request: string with request to be printed to console

            Return value: string with response to request
        """
        while True:
            sys.stdout.write('%s: ' % question)
            try:
                path = which(raw_input())
                if path is not None:
                    return path
                else:
                    raise ValueError(''.join(['Not a valid install of ',
                                                software_name]))
            except KeyboardInterrupt:
                sys.stdout.write('\n')
                sys.exit(0)

    def _grab_and_explode(self, urls, name):
        """ Special method for grabbing and exploding a package, if necessary.

            Does not verify URLs since these are preverified. Package is
            downloaded to current directory.

            url: list of urls to grab
            name: name of download

            No return value
        """
        self._print_to_screen_and_log('[Configuring] Downloading %s...' % name,
                                        newline=False,
                                        carriage_return=True)
        url_deque = deque(urls)
        while url_deque:
            url = url_deque.popleft()
            '''Follow redirects (-L), write to file (-O), respect
            content disposition'''
            command = [self.curl_exe, '-L', '-O', url]
            filename = url.rpartition('/')[2]
            print >>self.log_stream, ('[Configuring] '
                                      'Downloading {} from {}...'.format(
                                                                    name, url
                                                                ))
            try:
                subprocess.check_output(command, stderr=self.log_stream)
            except subprocess.CalledProcessError as e:
                if not url_deque:
                    self._print_to_screen_and_log(
                            ('Error encountered downloading file %s; exit '
                             'code was %d; command invoked was "%s".') %
                                (url, e.returncode, ' '.join(command))
                        )
                    self._print_to_screen_and_log(
                            'Make sure web is accessible.'
                        )
                    self._bail()
                else:
                    self._print_to_screen_and_log(
                        '[Configuring] Download failed; '
                        'trying alternate URL for %s...' % name,
                        newline=False,
                        carriage_return=True
                    )
            else:
                # Explode
                explode_command = None
                if url[-8:] == '.tar.bz2':
                    explode_command = ['tar', 'xvjf', filename]
                elif url[-7:] == '.tar.gz' or url[-4:] == '.tgz':
                    explode_command = ['tar', 'xvzf', filename]
                elif url[-4:] == '.zip':
                    self._print_to_screen_and_log(
                            'Extracting %s...' % name,
                            newline=False,
                            carriage_return=True)
                    try:
                        with zipfile.ZipFile(filename) as zip_object:
                            zip_object.extractall()
                    except Exception as e:
                        if not url_deque:
                            self._print_to_screen_and_log(
                                    'Error encountered exploding %s.'
                                        % filename
                                )
                            self._bail()
                        else:
                            self._print_to_screen_and_log(
                                '[Installing] Extraction failed; '
                                'trying alternate URL for %s...' % name,
                                newline=False,
                                carriage_return=True
                            )
                            continue
                    finally:
                        os.remove(filename)
                if explode_command is not None:
                    self._print_to_screen_and_log(
                            '[Configuring] Extracting %s...' % name,
                            newline=False,
                            carriage_return=True)
                    try:
                        subprocess.check_output(explode_command,
                                                stderr=self.log_stream)
                    except subprocess.CalledProcessError as e:
                        if not url_deque:
                            self._print_to_screen_and_log(
                                ('Error encountered exploding file %s; exit '
                                 'code was %d; command invoked was "%s".') %
                                (filename, e.returncode,
                                    ' '.join(explode_command))
                            )
                            self._bail()
                        else:
                            self._print_to_screen_and_log(
                                '[Configuring] Extraction failed; '
                                'trying alternate URL for %s...' % name,
                                newline=False,
                                carriage_return=True
                            )
                            continue
                    finally:
                        os.remove(filename)
                break

    def _quote(self, path=None):
        """ Decks path with single quotes if it's not None

            path: path or None

            Return value: None if path is None, else path decked by single
                quotes
        """
        if path is None:
            return 'None'
        return "'%s'" % path

    def configure(self):
        """ Downloads neoepiscope dependencies. """
        if not self.no_dependencies and self.curl_exe is None:
            self.curl_exe = which('curl')
            if self.curl_exe is None:
                print_to_screen('Configuring neoepiscope requires Curl if '
                                'dependencies are to be installed. '
                                'Download it at '
                                'http://curl.haxx.se/download.html and use '
                                '--curl to specify its path, or '
                                'disable installing dependencies with '
                                '--no-dependencies.')
                sys.exit(1)
        if self.download_dir is None:
            self.download_dir = os.path.join(os.path.expanduser('~'),
                                                'neoepiscope.data')
        if self._yes_no_query('Download Gencode v27 gtf annotation file?'):
            self._grab_and_explode(download['Gencode v27 annotation'], 
                                   'Gencode v27 annotation')
        if self._yes_no_query('Download Gencode v19 gtf annotation file?'):
            self._grab_and_explode(download['Gencode v19 annotation'], 
                                   'Gencode v19 annotation')
        if self._yes_no_query('Download Bowtie NCBI GRCh38 index?'):
            self._grab_and_explode(download['Bowtie NCBI GRCh38 index'], 
                                   'Bowtie NCBI GRCh38 index')
        if self._yes_no_query('Download Bowtie UCSC hg19 index?'):
            self._grab_and_explode(download['Bowtie UCSC hg19 index'], 
                                   'Bowtie UCSC hg19 index')
        self._grab_and_explode(download['HapCUT2'], 
                                   'HapCUT2')
        # Download to a temporary directory first, then move to final dest
        temp_install_dir = tempfile.mkdtemp()
        register_cleanup(remove_temporary_directories, [temp_install_dir])

        if os.path.exists(self.final_install_dir):
            if self._yes_no_query(
                    ('The installation path {dir} already exists.\n    '
                     '* Overwrite {dir}?').format(dir=self.final_install_dir)
                ):
                try:
                    shutil.rmtree(self.final_install_dir)
                except OSError:
                    # Handle this later if directory creation fails
                    pass
                try:
                    os.remove(self.final_install_dir)
                except OSError:
                    pass
            else:
                print_to_screen(
                        'Specify a different installation directory with '
                        '--install-dir.'
                    )
                sys.exit(0)
        self._print_to_screen_and_log('[Installing] Extracting neoepiscope...',
                                        newline=False,
                                        carriage_return=True)
        try:
            os.makedirs(self.final_install_dir)
        except OSError as e:
            self._print_to_screen_and_log(
                            ('Problem encountered trying to create '
                             'directory %s for installation. May need '
                             'sudo permissions.') % self.final_install_dir
                        )
            self._bail()
        else:
            # So it's possible to move temp installation dir there
            os.rmdir(self.final_install_dir)
            pass
        with cd(temp_install_dir):
            # Have to make HapCUT2
            hapcut2_dir = os.path.join(temp_install_dir,
                    download['HapCUT2'][0].split('/')[4]
                )
            with cd(hapcut2_dir):
                # Make on all but one cylinder
                thread_count = max(1, multiprocessing.cpu_count() - 1)
                hapcut2_command = ['make', '-j%d' % thread_count]
                self._print_to_screen_and_log(
                            '[Installing] Making HapCUT2...',
                            newline=False,
                            carriage_return=True
                        )
                try:
                    subprocess.check_output(hapcut2_command,
                                                stderr=self.log_stream)
                except subprocess.CalledProcessError as e:
                    self._print_to_screen_and_log(
                            ('Error encountered making HapCUT2; exit '
                             'code was %d; command invoked was "%s".') %
                                (e.returncode, ' '.join(hapcut2_command))
                        )
                    self._bail()
            hapcut2 = os.path.join(self.final_install_dir,
                                   download['HapCUT2'][0].split('/')[4],
                                   'HAPCUT2')
            hapcut2_hairs = os.path.join(self.final_install_dir,
                                        download['HapCUT2'][0].split('/')[4],
                                        'extractHAIRS')

            if self._yes_no_query('Do you have an install of netMHCIIpan '
                                  'version 3 that you would like to use for '
                                  'binding score predictions with neoepiscope?'
                                  ):
                netMHCIIpan3 = self._request_path('Please enter the path to '
                                                  'your netMHCIIpan v3 '
                                                  'executable', 
                                                  'netMHCIIpan v3')
            else:
                netMHCIIpan3 = 'None'

            if self._yes_no_query('Do you have an install of netMHCpan '
                                  'version 3 that you would like to use for '
                                  'binding score predictions with neoepiscope?'
                                  ):
                netMHCpan3 = self._request_path('Please enter the path to '
                                                'your netMHCpan v3 '
                                                'executable', 'netMHCpan v3')
            else:
                netMHCpan3 = 'None'

            if self._yes_no_query('Do you have an install of netMHCpan '
                                  'version 4 that you would like to use for '
                                  'binding score predictions with neoepiscope?'
                                  ):
                netMHCpan4 = self._request_path('Please enter the path to'
                                                'your netMHCpan v4 '
                                                'executable', 'netMHCpan v4')
            else:
                netMHCIIpan3 = 'None'
    
            # Write paths to exe_paths
            with open(
                            os.path.join(temp_install_dir, 'neoepiscope',
                                            'exe_paths.py'), 'w'
                        ) as exe_paths_stream:
                print >>exe_paths_stream, (
"""\"""
exe_paths.py
Part of neoepiscope

Defines default paths of neoepiscope's executable dependencies. Set a given
variable equal to None if the default path should be in PATH.
\"""

hapcut2_hairs = {hapcut2_hairs}
hapcut2 = {hapcut2}
netMHCIIpan3 = {netMHCIIpan3}
netMHCpan3 = {netMHCpan3}
netMHCpan4 = {netMHCpan4}
"""
                ).format(hapcut2_hairs=self._quote(hapcut2_hairs), 
                            hapcut2=self._quote(hapcut2),
                            netMHCIIpan3=self._quote(netMHCIIpan3),
                            netMHCpan3=self._quote(netMHCpan3),
                            netMHCpan4=self._quote(netMHCpan4))
        # Move to final directory
        try:
            shutil.move(temp_install_dir, self.final_install_dir)
        except Exception as e:
            self._print_to_screen_and_log(('Problem "%s" encountered moving '
                                           'temporary installation directory '
                                           '%s to final destination %s.') % (
                                                e,
                                                temp_install_dir,
                                                self.final_install_dir
                                            ))
            self._bail()
        # Create shell-script executable
        try:
            os.makedirs(bin_dir)
        except Exception as e:
            if not os.path.isdir(bin_dir):
                self._print_to_screen_and_log(('Problem "%s" encountered '
                                               'creating directory %s.') % (
                                                    e,
                                                    bin_dir
                                                )
                                            )
                self._bail()
        self._print_to_screen_and_log('Configured neoepiscope')
        self.finished = True
