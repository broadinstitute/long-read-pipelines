import argparse
import logging
from pathlib import Path

import utility as util
import wdlviz
import WDL

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

    util.set_logging_level(args)

    dir_names = args.wdl_dir
    output_path = str(util.get_absolute_path(args.output_path)) if args.output_path else None

    for dir_name in dir_names:

        Logger.debug(f'Processing {dir_name}...')

        wdl_dir = util.get_absolute_path(dir_name)

        wdl_paths = util.get_all_files_with_extension(directory=wdl_dir, ext='wdl')

        for wdl in wdl_paths:
            Logger.debug(f'Processing {wdl}...')

            if WDL.load(wdl).workflow:  # check if wdl is a workflow
                run_wdlviz(wdl_path=wdl, output_path=output_path)
            else:
                Logger.debug(f"{wdl} is not a WDL workflow (does not contain a "
                             "workflow block). Skipping...")


def run_wdlviz(
    wdl_path: str, use_basename: bool = True, output_path: str = None
) -> None:
    """
    Run wdlviz on a WDL file to create dot and png files
    @param use_basename: name output files using the basename of the WDL file
    @param wdl_path: path to WDL file
    @param output_path: path to output directory
    @return:
    """
    wdlviz_args = [wdl_path]

    if output_path is not None:
        wdlviz_args.append("--output_path")
        wdlviz_args.append(output_path)
    if use_basename:
        wdlviz_args.append("--dot-basename")
        wdlviz_args.append(Path(wdl_path).stem)

    Logger.debug(f'Running wdlviz with args: {wdlviz_args}...')
    wdlviz.main(args=wdlviz_args)


if __name__ == "__main__":
    main()

