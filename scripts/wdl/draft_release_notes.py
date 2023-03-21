import os
import sys
from pathlib import Path
import subprocess
import get_reverse_wdl_deps

ROOT_REPO_PATH = Path(__file__).resolve().parents[2]

if len(sys.argv) > 1:
    if sys.argv[1] in ["-h", "--help"]:
        print("""
        Compares two branches/tags to retrieve a list of edited wdl files,
        and a list of affected wdl pipelines. By default compares 'main' branch and
        latest tag release.
        Usage:
        \tdraft_release_notes 
        or..
        \tdraft_release_notes [branch/tag 1] [branch/tag 2]
        
        Example: 
        \tdraft_release_notes main 3.0.58
        """
              )
        exit(0)

if len(sys.argv) == 3:
    input_branch1 = sys.argv[1]
    input_branch2 = sys.argv[2]
else:
    input_branch1 = "main"
    input_branch2 = subprocess.run(
        ["git", "describe", "--tags", "--abbrev=0"], capture_output=True
    ).stdout.decode("utf-8").split("\n")[0]

print(f"Comparing branches {input_branch1} and {input_branch2}")


def get_edited_files_in_repo(branch1, branch2) -> list:
    """
    Returns a list of edited files between two branches/tags
    @param branch1: name of branch/tag
    @param branch2: name of branch/tag
    @return:
    """
    ret = subprocess.run(["git", "diff", "--name-only", branch1, branch2],
                         capture_output=True)
    path_2_edited_files = ret.stdout.decode("utf-8").split("\n")

    # remove empty indexes before sending
    return list(filter(None, path_2_edited_files))


def parse_for_wdl_files_in_list(files: list) -> list:
    """
    Returns a list of wdl files from a list of files
    @param files: list of files
    @return:
    """
    wdl_files = []
    for f in files:
        if ".wdl" in f:
            wdl_files.append(f)

    return wdl_files


def get_basename_for_files_in_list(files: list) -> list:
    """
    Returns a list of basenames for a list of files
    @param files: list of files
    @return:
    """
    basenames = []
    for f in files:
        basenames.append(os.path.basename(f))

    return basenames


def convert_dep_wdl_map_dict_to_list(dep_wdl_map: dict) -> list:
    """
    Converts a dict of wdl files to a list of wdl files
    @param dep_wdl_map: dict of wdl files
    @return:
    """
    # dict is turned into list of tuples
    some_list = list(dep_wdl_map.items())

    flattened_list = list(sum(some_list, ()))

    # Script below flattens the list of tuples to just a list
    # It will also automatically not include duplicates items
    affected_wdls = []
    for file in flattened_list:
        if type(file) is not dict:
            affected_wdls.append(file)

    return affected_wdls


def parse_pipeline_wdls_from_list(files: list) -> list:
    """
    Returns a list of pipeline wdls from a list of files
    @param files: list of files
    @return:
    """
    pipeline_wdls = []
    for file in files:
        if "/pipelines/" in file:
            pipeline_wdls.append(file)

    return pipeline_wdls


def print_list_as_markdown_bullet_point(items: list) -> None:
    """
    Prints a list of items as a markdown bullet point list
    @param items:
    @return:
    """
    for item in items:
        print("- " + item)


def remove_edited_wdls_from_list(wdls_affected: list, wdl_edited_names: list) -> list:
    """
    Removes edited wdls from the list of affected wdls
    @param wdls_affected:
    @param wdl_edited_names:
    @return:
    """
    for name in wdl_edited_names:
        for wa in wdls_affected:
            if name in wa:
                wdls_affected.remove(wa)
    return wdls_affected


edited_files: list = get_edited_files_in_repo(
    branch1=input_branch1,
    branch2=input_branch2,
)

edited_wdl_files = parse_for_wdl_files_in_list(edited_files)

if len(edited_wdl_files) == 0:
    print("No edited WDLs found.")
else:
    print("The following wdls haved been edited")
    print_list_as_markdown_bullet_point(edited_wdl_files)

    edited_wdl_names: list = get_basename_for_files_in_list(edited_wdl_files)

    all_rev_dep_wdl = {}
    for wdl in edited_wdl_names:
        reverse_dep_wdl = get_reverse_wdl_deps.get_wdls_that_import(
            wdl_basename=wdl,
            dir_path_to_check=str(ROOT_REPO_PATH),
        )

        if reverse_dep_wdl:
            all_rev_dep_wdl.update(reverse_dep_wdl)

    affected_wdls: list = convert_dep_wdl_map_dict_to_list(all_rev_dep_wdl)
    nonedited_affected_wdls: list = remove_edited_wdls_from_list(
        wdls_affected=affected_wdls, wdl_edited_names=edited_wdl_names
    )
    affected_pipeline_wdls: list = parse_pipeline_wdls_from_list(
        files=nonedited_affected_wdls
    )

    if affected_pipeline_wdls:
        print("The following pipeline wdls may have been affected")
        print_list_as_markdown_bullet_point(affected_pipeline_wdls)
    else:
        print("No affected pipeline WDL found (or is already in the above edited list)")
