import os
import re
import urllib.request
from urllib.error import HTTPError, URLError
import json


# A script to collect which dockers are in use and which latest dockers are available

def main():
    dir = os.path.abspath(os.path.dirname(__file__))
    # os.chdir(dir)

    print("COLLECTING DOCKERS IN USE...")
    # os.chdir("../../wdl")
    wdls_dir = os.path.join(dir, "../../wdl")
    if os.path.exists("dockers.in_use.tsv"):
        os.remove("dockers.in_use.tsv")

    wdl_files = get_wdl_files(dir_to_wdls=wdls_dir)
    global_docker_info = []

    for wdl_path in wdl_files:

        wdl_name = wdl_path

        with open(wdl_path, "r") as file:
            content = file.read()
            pattern = re.compile(r'.*docker.*"')
            if pattern.search(content):
                matched_lines = []
                file.seek(0)
                lines = file.readlines()

                for line_number, line in enumerate(lines, start=1):
                    if pattern.search(line):
                        matched_lines.append((line_number, line.strip()))

                docker_info: list[str] = get_docker_info_from_string(
                    wdl_lines=matched_lines, wdl_name=wdl_name
                )

                sorted_info: list = sorted(docker_info, reverse=False)

                global_docker_info.append(sorted_info)

    with open("dockers.in_use.tsv", "a") as tsv_file:
        tsv_file.write(f"name\tused_tag\tlatest_tag\tline\twdl\n")
        for line in sorted(global_docker_info):
            tsv_file.write("\n".join(line) + "\n")

    print("DONE. PLEASE CHECKOUT TSV FILE: dockers.in_use.tsv")
    # os.chdir(dir)
    # os.rename("../../wdl/dockers.in_use.tsv", "dockers.in_use.tsv")


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


def get_docker_info_from_string(wdl_lines: [tuple], wdl_name: str) -> list:
    """
    Returns a list of docker info
    @param wdl_name:
    @param wdl_lines: (line_number, line_content)
    @return:
    """
    docker_detail = []

    for line_num, line_content in wdl_lines:
        docker_names = re.findall(r'docker.*"(\S*?)"', line_content)
        if docker_names:
            docker_name = docker_names[0]
            used_tag = os.path.basename(docker_name).split(":")[1]
            docker_path = docker_name.split(":")[0]
            latest_tag = get_latest_local_docker_tag(docker_path)
            latest_tag = get_latest_remote_docker_tag(
                docker_path) if latest_tag == "NA" else latest_tag
            docker_detail.append(
                f"{docker_path}\t{used_tag}\t{latest_tag}\t{line_num}\t{wdl_name}")
        else:
            pass

    return docker_detail


def get_latest_remote_docker_tag(docker_path: str) -> str:
    """
    Returns the latest tag of a docker
    @param docker_path:
    @return:
    """
    if "gcr" in docker_path:
        latest_tag = get_latest_tag_from_gcr(docker_path)
    elif "quay.io" in docker_path:
        latest_tag = get_latest_tag_from_quayio(docker_path)
    else:
        latest_tag = get_latest_tag_from_duckerhub(docker_path)
    return latest_tag


def get_latest_tag_from_duckerhub(docker_path: str) -> str:
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
    # Split the image string into project ID and image name
    parts = docker_path.split("/")
    gcr_repo = parts[0]
    project_id = parts[1]
    image_name = "/".join(parts[2:])
    # Construct the URL for retrieving tags
    registry_url = f"https://{gcr_repo}/v2/{project_id}/{image_name}/tags/list"

    try:
        # Send the GET request to the Container Registry API
        with urllib.request.urlopen(registry_url) as response:
            data = response.read().decode("utf-8")
            json_data = json.loads(data)
            tags = json_data.get("tags")
            if tags:
                latest_tag = max(tags)
                return latest_tag
            else:
                return "NA"
    except urllib.error.HTTPError as e:
        # print(f"Error: {e.code} - {e.reason}")
        pass
    except urllib.error.URLError as e:
        # print(f"Error: Failed to reach the server - {e.reason}")
        pass


def get_latest_tag_from_quayio(docker_path: str) -> str:
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
            data = response.read().decode("utf-8")
            json_data = json.loads(data)
            tags = json_data.get("tags")
            if tags:
                latest_tag = max(tags)
                return latest_tag
            else:
                return "NA"
    except urllib.error.HTTPError as e:
        # print(f"Error: {e.code} - {e.reason}")
        pass
    except urllib.error.URLError as e:
        # print(f"Error: Failed to reach the server - {e.reason}")
        pass


def get_latest_local_docker_tag(docker_path: str) -> str:
    """
    Returns the latest tag of a docker
    @param docker_path:
    @return:
    """
    docker_name = os.path.basename(docker_path)
    docker_dir = "../docker"
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


if __name__ == "__main__":
    main()
