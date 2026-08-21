import argparse
import os
import re
import subprocess
import urllib.request
from urllib.error import HTTPError, URLError
import json
import logging

logging.basicConfig(level=logging.INFO)


# TODO: Future suggestion: have the results be generated for main branch for each merge

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
Collects docker usage summary from wdl files.

Output: 
    dockers.in_use.tsv
Notes: 
    - This script is not perfect. It will not be able to detect dockers that are
    imported from other wdl files. It will only detect dockers that are
    explicitly defined in the wdl file.
    - The script assumes it is executed from the scripts/docker directory, and the
    wdl files are in ../../wdl directory.
    - The script will retrieve the "latest" tag by date, so if an unofficial tag was
    created, after the official tag was created, the script will retrieve the
    unofficial tag as the latest tag. (It tries to avoid this by filtering out
    tags with no digits).
    - The script occasionally uses gcloud to retrieve the latest tag. Its suggested
    to have gcloud installed.
        ''',
    )
    parser.parse_args()

    current_dir = os.path.abspath(os.path.dirname(__file__))
    repo_dir = os.path.abspath(os.path.join(current_dir, os.pardir, os.pardir))

    print("COLLECTING DOCKERS IN USE...")
    WDLS_DIR = os.path.abspath(os.path.join(current_dir, "../../wdl"))
    OUT_SUMMARY_TSV = os.path.join(current_dir, "dockers_in_use.tsv")

    if os.path.exists(OUT_SUMMARY_TSV):
        os.rename(OUT_SUMMARY_TSV, OUT_SUMMARY_TSV + ".bak")

    wdl_files = get_wdl_files(dir_to_wdls=WDLS_DIR)

    pattern = re.compile(r'^\s*docker:\s*"')
    global_docker_info = process_wdl_files(
        wdl_files=wdl_files, pattern=pattern, repo_dir=repo_dir
    )

    # Remove empty elements in list
    non_empty_docker_info = [x for x in global_docker_info if x]

    write_docker_info_to_tsv(
        output_summary_tsv=OUT_SUMMARY_TSV,
        docker_info=non_empty_docker_info
    )

    print(f"DONE. PLEASE CHECKOUT TSV FILE: {OUT_SUMMARY_TSV}")


def get_wdl_files(dir_to_wdls: str) -> list:
    """
    Returns a list of wdl files
    @return:
    """
    wdl_files = []
    for root, _, files in os.walk(dir_to_wdls):
        for filename in files:
            if filename.endswith(".wdl"):
                wdl_path = os.path.join(root, filename)
                wdl_files.append(wdl_path)

    return wdl_files


def get_docker_info_from_string(wdl_lines: [tuple], wdl_path: str, repo_dir: str) -> list:
    """
    Returns a list of docker info
    @param wdl_path: path to wdl file
    @param wdl_lines: (line_number, line_content)
    @return: list of docker info e.g. [" docker_name\tlatest_tag\tused_tag\tline_num\twdl_path, ..."]
    """
    docker_detail = []

    # Get the path after /wdl/ for better readability
    wdl_path_sum = wdl_path[wdl_path.find("/wdl/"):]

    for line_num, line_content in wdl_lines:
        docker_image = re.findall(r'docker:\s*"(\S*?)"', line_content)
        if docker_image:
            docker_name_and_version = docker_image[0]
            used_tag = os.path.basename(docker_name_and_version).split(":")[1]
            docker_name = docker_name_and_version.split(":")[0]

            # Get latest tag from list of docker details if it was already retrieved
            latest_tag = get_tag_from_docker_details(
                docker_detail=docker_detail, docker_name=docker_name
            )
            # Get latest tag from local docker if it was not retrieved from list of docker details
            latest_tag = get_latest_local_docker_tag(
                docker_name=docker_name,  repo_dir=repo_dir
            ) if latest_tag == "NA" else latest_tag

            # If the latest tag is not found locally, try to get it from remote
            latest_tag = get_latest_remote_docker_tag(
                docker_name) if latest_tag == "NA" else latest_tag
            docker_detail.append(
                f"{docker_name}\t{latest_tag}\t{used_tag}\t{line_num}\t{wdl_path_sum}")

    return docker_detail


def get_tag_from_docker_details(docker_detail: list, docker_name: str) -> str:
    """
    Returns the latest tag of a docker from a list of docker details
    @param docker_detail: list of docker details e.g. ["docker_name\tlatest_tag\tused_tag\tline_num\twdl_path, ..."]
    @param docker_name: docker name
    @return: latest tag
    """
    latest_tag = "NA"
    for docker_info in docker_detail:
        if docker_name in docker_info:
            latest_tag = docker_info.split("\t")[1]
            break
    return latest_tag


def get_latest_remote_docker_tag(docker_path: str) -> str:
    """
    Returns the latest tag of a docker from gcr, quay or dockerhub
    @param docker_path:
    @return:
    """
    if "gcr" in docker_path or "ghcr" in docker_path:
        latest_tag = get_latest_tag_from_gcr(docker_path)
        if latest_tag in ["NA", "None", "latest"] or latest_tag is None:
            latest_tag = get_gcr_tag_with_gcloud(docker_path)
    elif "quay.io" in docker_path:
        latest_tag = get_latest_tag_from_quay(docker_path)
    else:
        latest_tag = get_latest_tag_from_dockerhub(docker_path)
    return latest_tag


def get_latest_tag_from_dockerhub(docker_path: str) -> str:
    """
    Returns the latest tag of a docker from dockerhub using the dockerhub API
    @param docker_path:
    @return:
    """

    image_name = docker_path
    registry_url = f"https://registry.hub.docker.com/v2/repositories/{image_name}/tags/?page_size=1&ordering=last_updated"
    try:
        with urllib.request.urlopen(registry_url) as response:
            data = response.read().decode("utf-8")
            json_data = json.loads(data)
            tags = json_data.get("results")
            if tags:
                latest_tag = tags[0].get("name")
                return latest_tag
            else:
                return "NA"
    except urllib.error.HTTPError as e:
        # print(f"Error: {e.code} - {e.reason}")
        pass
    except urllib.error.URLError as e:
        # print(f"Error: Failed to reach the server - {e.reason}")
        pass


def get_latest_tag_from_gcr(docker_path: str) -> str:
    """
    Returns the latest tag of a docker from GCR using the Container Registry API
    @param docker_path:
    @return:
    """

    # Split the image string into project ID and image name
    # us.gcr.io/broad-dsp-lrma/lr-transcript_utils:latest
    parts = docker_path.split("/")
    gcr_repo = parts[0]  # Example: us.gcr.io
    project_id = parts[1]  # Example: broad-dsp-lrma
    image_name = "/".join(parts[2:])  # Example: lr-transcript_utils

    # Construct the URL for retrieving tags
    # https://cloud.google.com/artifact-registry/docs/reference/docker-api
    registry_url = f"https://{gcr_repo}/v2/{project_id}/{image_name}/tags/list"

    try:
        # Send the GET request to the Container Registry API
        with urllib.request.urlopen(registry_url) as response:
            return extract_latest_tag_from_registry_response(response)
    except urllib.error.HTTPError as e:
        # print(f"Error: {e.code} - {e.reason}")
        pass
    except urllib.error.URLError as e:
        # print(f"Error: Failed to reach the server - {e.reason}")
        pass


def extract_latest_tag_from_registry_response(response) -> str:
    """
    Extracts the latest tag from the response of a registry API call
    @param response:
    @return:
    """

    response_content = response.read().decode("utf-8")
    latest_tag = "NA"

    if response_content:
        json_data = json.loads(response_content)
        tags = json_data.get("tags")
        manifest = json_data.get("manifest")

        # The manifest is a dict of dicts for each version of an image, each dict
        # has a key called "tag", which is a list of tags for that version.
        # The image version having "latest" as part of its tag is returned.
        if manifest:
            for sha_key in manifest:
                sha_key_tags = manifest[sha_key].get("tag")
                if "latest" in sha_key_tags:
                    latest_tag = sha_key_tags[0]
                    return latest_tag

        # If the image doesn't have a "latest" tag, return the tag with the
        # highest version number.
        tags_str_removed = [item for item in tags if
                            any(char.isdigit() for char in item)]
        if tags_str_removed:
            latest_tag = max(tags_str_removed)
            return latest_tag
        else:
            return latest_tag  # If no numerical version tags are found, return NA
    else:
        return latest_tag  # If manifest is empty, return NA


def get_gcr_tag_with_gcloud(docker_path: str) -> str or None:
    """
    Returns the latest tag of a docker using gcloud
    @param docker_path:
    @return:
    """

    # Split the image string into project ID and image name

    if is_gcloud_installed():

        command = [
            "gcloud",
            "container",
            "images",
            "list-tags",
            docker_path,
            "--format=json",
            "--limit=1",
            "--sort-by=~timestamp.datetime",
            "--filter=tags:latest",
        ]

        gc_container_results = subprocess.run(command, capture_output=True, text=True)
        if gc_container_results and gc_container_results.returncode == 0:
            gc_container_results_json = json.loads(gc_container_results.stdout)
            try:
                latest_tag = gc_container_results_json[0].get("tags")[0]
                return latest_tag if latest_tag is not None else "NA"
            except IndexError:

                logging.warning(
                    f"Gcloud Container obtain empty tag for : {gc_container_results_json} - {docker_path}")
                return "NA"
        else:
            # Error handling
            error_message = gc_container_results.stderr.strip() if gc_container_results.stderr else gc_container_results.stdout.strip()
            # print(f"Error: {error_message}")
    else:
        return "NA"


def is_gcloud_installed() -> bool:
    """
    Checks if gcloud is installed
    @return:
    """

    command = ["gcloud", "--version"]

    try:
        subprocess.run(command, check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        return False


def get_latest_tag_from_quay(docker_path: str) -> str:
    """
    Returns the latest tag of a docker from quay.io
    @param docker_path:
    @return:
    """
    # Split the image string into project ID and image name
    parts = docker_path.split("/")
    quayio_repo = parts[0]
    project_id = parts[1]
    image_name = "/".join(parts[2:])
    # Construct the URL for retrieving tags
    registry_url = f"https://{quayio_repo}/v2/{project_id}/{image_name}/tags/list"

    try:
        # Send the GET request to the Container Registry API
        with urllib.request.urlopen(registry_url) as response:
            return extract_latest_tag_from_registry_response(response)
    except urllib.error.HTTPError as e:
        # print(f"Error: {e.code} - {e.reason}")
        pass
    except urllib.error.URLError as e:
        # print(f"Error: Failed to reach the server - {e.reason}")
        pass


def get_latest_local_docker_tag(docker_name: str, repo_dir: str) -> str:
    """
    Returns the latest tag of a docker from the local docker directory
    @param docker_name: name of the docker e.g. "gatk"
    @return:
    """
    docker_name = os.path.basename(docker_name)
    docker_dir = os.path.join(repo_dir, "docker")
    latest_tag = "NA"

    for docker_im_dir in os.listdir(docker_dir):
        if docker_im_dir == docker_name:
            docker_dir_path = os.path.join(docker_dir, docker_im_dir)
            for makefile in os.listdir(docker_dir_path):
                if not makefile.endswith("Makefile"):
                    continue

                with open(os.path.join(docker_dir_path, makefile)) as f:
                    for makefile_line in f:
                        if "VERSION =" in makefile_line:
                            latest_tag = makefile_line.split("=")[1].strip()

    return latest_tag


def process_wdl_files(wdl_files: list[str], pattern, repo_dir: str):
    """
    Returns a list of docker info
    @param wdl_files: list of wdl files
    @param pattern: pattern to search for
    @param repo_dir: directory of the repo
    @return:
    """
    global_docker_info = []
    total_files = len(wdl_files)

    for index, wdl_path in enumerate(wdl_files, start=1):
        wdl_name = wdl_path
        matched_lines = []

        with open(wdl_path, "r") as file_content:
            lines = file_content.readlines()
            for line_number, line in enumerate(lines, start=1):
                if pattern.search(line):
                    matched_lines.append((line_number, line.strip()))

        docker_info: list[str] = get_docker_info_from_string(
            wdl_lines=matched_lines, wdl_path=wdl_name, repo_dir=repo_dir
        )

        sorted_info: list = sorted(docker_info, reverse=False)
        global_docker_info.append(sorted_info)

        # Visual progression to show percentage of files processed
        progress: float = (index + 1) / total_files * 100

        # Clear the previous line and print the progress
        print(f"Progress: {progress:.2f}%\r", end="")

    return global_docker_info


def write_docker_info_to_tsv(
    output_summary_tsv: str, docker_info: list[list]
):
    """
    Writes docker info to tsv file

    @param output_summary_tsv:
    @param docker_info:
    @return:
    """

    with open(output_summary_tsv, "w") as tsv_file:
        # Add header
        tsv_file.write("DOCKER_NAME\tLATEST_TAG\tUSED_TAG\tFILE_LINE\tWDL_PATH\n")
        # Add content
        for docker_info_line in sorted(docker_info):
            delimiter = "\n"
            tsv_file.write(delimiter.join(docker_info_line) + "\n")



if __name__ == "__main__":
    main()
