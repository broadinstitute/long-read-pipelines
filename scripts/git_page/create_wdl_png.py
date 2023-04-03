import errno
import glob

import os
import argparse
import wdlviz

from pathlib import Path, PurePosixPath

import logging

Logger = logging.getLogger(__name__)
logging.basicConfig()


def main():
    # read command-line arguments
    parser = argparse.ArgumentParser(
        description="Visualize a WDL workflow using miniwdl and graphviz"
    )
    parser.add_argument(
        "wdl_dir", metavar="WDL_DIRECTORY", nargs='+',
        help="One of or more directories containing WDL workflow files"
    )
    parser.add_argument("--output_path", help="Path to output directory")
    parser.add_argument("--debug", action="store_true", help="verbose logging")

    args = parser.parse_args()

    logging.getLogger().setLevel(logging.DEBUG) if args.debug else logging.getLogger().setLevel(logging.INFO)

    dir_names = args.wdl_dir
    output_path = str(get_absolute_path(args.output_path)) if args.output_path else None

    for dir_name in dir_names:

        Logger.debug(f'Processing {dir_name}...')

        wdl_dir = get_absolute_path(dir_name)
        if not wdl_dir.exists():
            Logger.error(f'Could not find {wdl_dir}!')
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), wdl_dir)
        wdl_paths = wdl_dir / "**/*.wdl"

        for wdl in glob.glob(str(wdl_paths), recursive=True):
            Logger.debug(f'Processing {wdl}...')
            print(f'{wdl}:')
            run_wdlviz(wdl_path=wdl, output_path=output_path)


def run_wdlviz(wdl_path: str, output_path: str = None) -> None:
    """
    Run wdlviz on a WDL file to create dot and png files
    @param wdl_path:
    @param output_path:
    @return:
    """
    wdlviz_args = [wdl_path]
    if output_path is not None:
        wdlviz_args.append("--output_path")
        wdlviz_args.append(output_path)

    Logger.debug(f'Running wdlviz with args: {wdlviz_args}...')
    wdlviz.main(args=wdlviz_args)


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


if __name__ == "__main__":
    main()

