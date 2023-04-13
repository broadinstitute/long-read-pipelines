import errno

import os
import argparse
from pathlib import Path
import logging

import utility as util

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

    util.set_logging_level(args)

    md_dir_names = args.md_dir
    dot_png_dir_names = args.dot_dir

    dot_png_paths = get_all_dot_png_in_directory(dot_png_dir_names)

    if len(dot_png_paths) == 0:
        Logger.error(f'Could not find any dot.png files in {dot_png_dir_names}!')
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), dot_png_paths)

    for md_dir_name in md_dir_names:

        Logger.debug(f'Processing {md_dir_name}...')

        md_dir = util.get_absolute_path(md_dir_name)

        md_paths = util.get_all_files_with_extension(directory=md_dir, ext='md')

        for md in md_paths:
            Logger.debug(f'Processing {md}...')
            append_png_to_md(md_path=md, dot_paths=dot_png_paths)

    print("Done!")


def append_png_to_md(md_path: str, dot_paths: list) -> None:
    """
    Append dot.png link to markdown file

    @param md_path: path to markdown file
    @param dot_paths: list of paths to dot.png files
    @return:
    """

    md_basename = util.get_basename(md_path)

    dot_png_matches = []

    for dot_path in dot_paths:
        if md_basename == util.get_basename(dot_path):
            dot_png_matches.append(dot_path)

    if len(dot_png_matches) == 1:
        relative_d_p_path = get_relative_path(dot_png_matches[0], md_path)
        with open(md_path, 'a') as md_file:
            md_file.write(f'### Dot Diagram\n')
            md_file.write(f'![{md_basename}](../{relative_d_p_path})\n')
    elif len(dot_png_matches) > 1:
        Logger.error(f'Found more than one dot.png for {md_basename}!')


def get_all_dot_png_in_directory(directory_names: list) -> list:
    """
    Get all dot.png files in a list of directories
    @param directory_names: list of directories to search
    @return: list of paths to dot.png files
    """

    dot_png_paths = []
    for dot_png_dir_name in directory_names:
        dot_png_dir = util.get_absolute_path(dot_png_dir_name)

        dot_png_paths = (
                dot_png_paths + util.get_all_files_with_extension(
                    directory=dot_png_dir, ext='dot.png'
                )
        )
    return dot_png_paths


def get_relative_path(path, relative_to_path) -> Path:
    """
    Get the relative path for a path
    @param relative_to_path:  Path to get relative path to
    @param path: Path to get relative path for
    @return:
    """
    p = Path(path)
    parent_dir = Path(relative_to_path).parent
    return p.relative_to(parent_dir.parent)


if __name__ == "__main__":
    main()

