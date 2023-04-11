import argparse
from pathlib import Path
import logging
from typing import Union

import utility as util

Logger = logging.getLogger(__name__)
logging.basicConfig()


def main():
    parser = argparse.ArgumentParser(
        description="Performs post processing on markdown files. - Adds title header to markdown files if it does not "
    )
    parser.add_argument(
        "--md_dir", metavar="MARKDOWN_DIRECTORY", required=True,
        help="Directory containing markdown documentation files for workflows"
    )
    parser.add_argument("--debug", action="store_true", help="verbose logging")

    args = parser.parse_args()

    util.set_logging_level(args)

    md_dir_name = args.md_dir

    md_dir = util.get_absolute_path(md_dir_name)

    md_paths = util.get_all_files_with_extension(directory=md_dir, ext='md')
    Logger.debug(f'Found {len(md_paths)} markdown files')

    for md in md_paths:
        Logger.debug(f'Processing {md}...')
        if not check_if_md_file_has_title_header(md_path=md):
            Logger.debug(f'Adding title header to {md}')
            add_title_header(md_path=md)

        if md_path_in_workflows(md_path=md) and workflow_md_has_multiple_subheaders(md_path=md):
            Logger.debug(f'Removing workflow subheaders from {md}')
            remove_subheaders_after_first(md_path=md)

    Logger.debug('Finished processing all markdown files')


def add_title_header(md_path: str) -> None:
    """
    Adds a title header to the markdown file
    :param md_path: Path to markdown file
    :return: None
    """
    with open(md_path, 'r+') as md_file:
        md_lines = md_file.read()
        header = f"# {Path(md_path).stem}\n\n"

        md_file.seek(0, 0)
        md_file.write(header + md_lines)

    Logger.debug(f'Added title header to {md_path}')


def check_if_md_file_has_title_header(md_path: str) -> bool:
    """
    Checks if the markdown file has a title header
    :param md_path: Path to markdown file
    :return: True if the markdown file has a title header
    """

    with open(md_path, 'r') as md_file:
        md_lines = [s.strip() for s in md_file.readlines()]

    if len(md_lines) == 0:
        return False

    for line in md_lines:
        if line.startswith('# '):
            return True

    return False


def md_path_in_workflows(md_path: str) -> bool:
    """
    Checks if the markdown file is a workflow by checking whether path contains
    the word 'workflows'
    :param md_path: Path to markdown file
    :return: True if the markdown file path contains 'workflows'
    """
    return "/workflows/" in md_path.lower()


def workflow_md_has_multiple_subheaders(md_path: str) -> bool:
    """
    Checks if the markdown file has multiple subheaders
    :param md_path: Path to markdown file
    :return: True if the markdown file has multiple subheaders
    """
    with open(md_path, 'r') as md_file:
        md_lines = [s.strip() for s in md_file.readlines()]

    if len(md_lines) == 0:
        return False

    subheader_count = 0
    for line in md_lines:
        if line.startswith('## '):
            subheader_count += 1

    return subheader_count > 1


def remove_subheaders_after_first(md_path: Union[str, Path]) -> None:
    """
    Removes the subheader from the markdown file after the first subheader.
    :param md_path: Path to markdown file
    :return: None
    """
    try:
        with Path(md_path).open('r+') as md_file:
            lines = md_file.readlines()

            second_subhead_index = None
            count = 0
            for i, line in enumerate(lines):
                if line.startswith("## "):
                    count += 1
                    if count == 2:
                        second_subhead_index = i
                        break
            if second_subhead_index is None:
                return
            md_file.seek(0)
            md_file.truncate()
            md_file.writelines(lines[:second_subhead_index])

    except (FileNotFoundError, IsADirectoryError):
        Logger.debug(f"Unable to open file at {md_path}")


if __name__ == "__main__":
    main()
