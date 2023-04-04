import glob
import logging
from pathlib import Path, PurePosixPath, PurePath

Logger = logging.getLogger(__name__)


def get_absolute_path(path) -> Path:
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
