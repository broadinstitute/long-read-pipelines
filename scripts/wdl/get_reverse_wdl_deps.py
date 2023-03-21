import glob
import json
import os
import sys


def get_wdls_that_import(wdl_basename, dir_path_to_check) -> dict:
    """
    Recursively searches for wdls that import the provided wdl_basename
    @param wdl_basename: name of wdl to search for
    @param dir_path_to_check: path to directory to search for imports
    @return:
    """
    wild_wdl_path = os.path.join(dir_path_to_check, "**/*.wdl")
    wdls_that_import = {}

    for wdl in glob.glob(wild_wdl_path, recursive=True):
        with open(wdl) as f:
            lines = f.readlines()
            for line in lines:
                if "import" and "/" + wdl_basename in line or "import" and '\"' + wdl_basename in line:
                    more_wdls = get_wdls_that_import(os.path.basename(wdl),
                                                     dir_path_to_check)
                    wdls_that_import[wdl] = more_wdls

    return wdls_that_import


def main():
    if sys.argv[1] in ["-h", "--help"]:
        print("This script searches for the reverse dependency of a "
              "provided wdl in a given directory.\n")
        print("Usage:\n\tget_reverse_wdl_deps [wdl name] [directory to search]")
        print(
            "Example: \n\tget_reverse_wdl_deps VariantCaller.wdl "
            "/Users/bob/longreadpipes/"
        )
        print("Warning: \n\tScript cannot distinguish WDLs with same names.")
        exit(0)

    if len(sys.argv) < 2:
        print("Error: Missing Arguments!!")
        print(
            "Please provide the name of wdl file and the directory to search for "
            "imports"
        )
        print(
            "Example:\n\tget_reverse_wdl_deps VariantCaller.wdl "
            "/Users/bob/longreadpipes/"
        )
        exit(0)

    input_wdl_name = sys.argv[1]
    input_dir_path = sys.argv[2]

    print(f"WDL name: {input_wdl_name}")
    print(f"Directory to search: {input_dir_path}")

    if not os.path.isdir(input_dir_path):
        raise Exception(f"The provided directory {input_dir_path} does not exist.")

    reverse_dep_wdl: dict = get_wdls_that_import(
        wdl_basename=input_wdl_name,
        dir_path_to_check=input_dir_path
    )

    print("-----------------------------------")
    print(json.dumps(reverse_dep_wdl, indent=4))


if __name__ == "__main__":
    main()
