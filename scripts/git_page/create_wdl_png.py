import glob

import argparse
import wdlviz

import logging

from scripts.git_page.utility import set_logging_level, \
    get_all_files_with_extension, get_absolute_path, check_dir_exists

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

    set_logging_level(args)

    dir_names = args.wdl_dir
    output_path = str(get_absolute_path(args.output_path)) if args.output_path else None

    for dir_name in dir_names:

        Logger.debug(f'Processing {dir_name}...')

        wdl_dir = get_absolute_path(dir_name)

        check_dir_exists(wdl_dir)

        wdl_paths = get_all_files_with_extension(directory=wdl_dir, ext='wdl')

        for wdl in glob.glob(str(wdl_paths), recursive=True):
            Logger.debug(f'Processing {wdl}...')

            run_wdlviz(wdl_path=wdl, output_path=output_path)


def run_wdlviz(wdl_path: str, output_path: str = None) -> None:
    """
    Run wdlviz on a WDL file to create dot and png files
    @param wdl_path: path to WDL file
    @param output_path: path to output directory
    @return:
    """
    wdlviz_args = [wdl_path]
    if output_path is not None:
        wdlviz_args.append("--output_path")
        wdlviz_args.append(output_path)

    Logger.debug(f'Running wdlviz with args: {wdlviz_args}...')
    wdlviz.main(args=wdlviz_args)


if __name__ == "__main__":
    main()

