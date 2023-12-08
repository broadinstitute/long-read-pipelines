import glob
import logging
import subprocess
from pathlib import Path, PurePosixPath, PurePath

Logger = logging.getLogger(__name__)


def get_absolute_path(path: Path or str) -> Path:
    """
    Get the absolute path for a path
    @param path: path to get absolute path for
    @return: absolute path
    """

    Logger.debug(f'Getting absolute path for {path}...')
    if PurePosixPath(path).is_absolute():
        return path
    else:
        return Path(path).resolve()


def set_logging_level(args):
    """
    Set the logging level
    @param args: command-line arguments with a debug flag
    @return:
    """
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)


def get_basename(path: str) -> str:
    """
    Get the basename of a path
    @param path: path to get basename of
    @return: basename
    """
    Logger.debug(f'Getting basename of {path}...')
    return PurePath(path).name.split('.')[0]


def check_dir_exists(directory: Path or str) -> bool:
    """
    Check if a directory exists
    @param directory:  directory to check
    @return: True if directory exists, False otherwise
    """
    Logger.debug(f'Checking if {directory} exists...')
    return Path(directory).exists()


def get_all_files_with_extension(directory: Path, ext: str) -> list:
    """
    Get all files with a given extension in a directory
    @param directory: directory to search
    @param ext: extension (e.g. 'md', 'dot.png')
    @return: list of files with extension
    """
    check_dir_exists(directory)

    Logger.debug(f'Getting all files with extension {ext} in {directory}...')
    return glob.glob(f'{directory}/**/*.{ext}', recursive=True)


def run_command(command: list, log_output=True, ) -> None:
    """
    Run a shell command and wait for it to complete.

    :param command: A list representing the shell command to execute.
    :param log_output: Whether to log the command's output (default is True).
    :return: None

    This function runs the specified shell command and waits for it to complete.
    It logs the command before execution and raises an exception if the command fails.
    """
    cmd_str = ' '.join(command)  # Convert the command list to a string for logging
    logging.debug(f'Running command: {cmd_str}...')

    try:
        result = subprocess.run(command, check=True, shell=False,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, encoding='utf-8')

        if log_output:
            if result.stdout:
                logging.debug(f'Command output (stdout):\n{result.stdout}')
            if result.stderr:
                logging.debug(f'Command output (stderr):\n{result.stderr}')
    except subprocess.CalledProcessError as e:
        logging.error(f'Command failed with error: {e}')
        raise
    except Exception as e:
        logging.error(f'An unexpected error occurred: {e}')
        raise

