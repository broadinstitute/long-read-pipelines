import errno
import glob

import os
import argparse

from pathlib import Path, PurePosixPath, PurePath

import logging

Logger = logging.getLogger(__name__)
logging.basicConfig()


def main():
    parser = argparse.ArgumentParser(
        description="Add dot header and dot.png link to workflow markdown files"
    )
    parser.add_argument(
        "--md_dir", metavar="MARKDOWN_DIRECTORY", nargs='+', required=True,
        help="One of or more directories containing markdown "
             "documentation files for workflow"
    )
    parser.add_argument(
        "--dot_dir", metavar="DOT_DIRECTORY", nargs='+', required=True,
        help="One of or more directories containing dot.png "
             "files for workflows"
    )
    parser.add_argument("--debug", action="store_true", help="verbose logging")

    args = parser.parse_args()

    set_logging_level(args)

    md_dir_names = args.md_dir
    dot_png_dir_names = args.dot_dir

    dot_png_paths = get_all_dot_png_in_directory(dot_png_dir_names)

    if len(dot_png_paths) == 0:
        Logger.error(f'Could not find any dot.png files in {dot_png_dir_names}!')
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), dot_png_paths)

    for md_dir_name in md_dir_names:

        Logger.debug(f'Processing {md_dir_name}...')

        md_dir = get_absolute_path(md_dir_name)

        md_paths = get_all_files_with_extension(directory=md_dir, ext='md')

        for md in md_paths:
            Logger.debug(f'Processing {md}...')
            append_png_to_md(md_path=md, dot_paths=dot_png_paths)

    print("Done!")


def append_png_to_md(md_path: str, dot_paths: list) -> None:
    """
    Append dot.png link to markdown file

    @param md_path:
    @param dot_paths:
    @return:
    """

    md_basename = get_basename(md_path)

    dot_png_matches = []

    for dot_path in dot_paths:
        if md_basename == get_basename(dot_path):
            dot_png_matches.append(dot_path)

    if len(dot_png_matches) == 1:
        with open(md_path, 'a') as md_file:
            md_file.write(f'### Dot Diagram\n')
            md_file.write(f'![{md_basename}]({dot_png_matches[0]})\n')
    elif len(dot_png_matches) > 1:
        Logger.error(f'Found more than one dot.png for {md_basename}!')


def get_absolute_path(path) -> Path:
    """
    Get the absolute path for a path
    @param path:
    @return:
    """

    Logger.debug(f'Getting absolute path for {path}...')
    if PurePosixPath(path).is_absolute():
        return path
    else:
        return Path(path).resolve()


def get_all_files_with_extension(directory: Path, ext: str) -> list:
    """
    Get all files with a given extension in a directory
    @param directory: directory to search
    @param ext: extension (e.g. 'md', 'dot.png')
    @return:
    """
    check_dir_exists(directory)

    Logger.debug(f'Getting all files with extension {ext} in {directory}...')
    return glob.glob(f'{directory}/**/*.{ext}', recursive=True)


def check_dir_exists(directory: Path or str) -> bool:
    """
    Check if a directory exists
    @param directory:
    @return:
    """
    Logger.debug(f'Checking if {directory} exists...')
    return Path(directory).exists()


def get_basename(path: str) -> str:
    """
    Get the basename of a path
    @param path:
    @return:
    """
    Logger.debug(f'Getting basename of {path}...')
    return PurePath(path).name.split('.')[0]


def set_logging_level(args):
    """
    Set the logging level
    @param args:
    @return:
    """
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)


def get_all_dot_png_in_directory(directory_names: list) -> list:
    """
    Get all dot.png files in a list of directories
    @param directory_names:
    @return:
    """

    dot_png_paths = []
    for dot_png_dir_name in directory_names:
        dot_png_dir = get_absolute_path(dot_png_dir_name)

        dot_png_paths = (
            dot_png_paths + get_all_files_with_extension(
                directory=dot_png_dir, ext='dot.png'
            )
        )
    return dot_png_paths


if __name__ == "__main__":
    main()

