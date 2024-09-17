#!/usr/bin/env python

'''
This script will load IGV in a virtual X window, load all supplied input files
as tracks, and take snapshots at the coordinates listed in the BED formatted
region file.
'''

# ~~~~ LOAD PACKAGES ~~~~~~ #
import sys
import os
import subprocess as sp
import argparse
import datetime
from pathlib import Path

# ~~~~ DIRECTORY AND DEFAULTS ~~~~~~ #
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
SNAPSHOT_DIR = "/cromwell_root/output/IGV_Snapshots"  # Snapshot output directory
default_igv_sh = os.path.join(THIS_DIR, 'igv.sh')
default_regions_bed = os.path.join(THIS_DIR, 'regions.bed')

# ~~~~ CUSTOM FUNCTIONS ~~~~~~ #
def file_exists(myfile, kill=False):
    '''
    Checks to make sure a file exists, optionally kills the script if file is missing.
    '''
    if not os.path.isfile(myfile):
        print(f"ERROR: File '{myfile}' does not exist!")
        if kill:
            print("Exiting...")
            sys.exit()

def check_for_fai(fasta_file):
    '''
    Check to make sure a .fai index file exists for the FASTA file.
    If not, alert the user.
    '''
    fai_file = fasta_file + ".fai"
    if not os.path.isfile(fai_file):
        print(f"ERROR: FASTA index file '{fai_file}' is missing!")
        print("Please generate it with 'samtools faidx' before running the script.")
        sys.exit(1)

def check_for_bai(bam_file):
    '''
    Check to make sure a .bam.bai file is present in the same directory as the .bam file.
    '''
    bai_file = bam_file + ".bai"
    if not os.path.isfile(bai_file):
        print(f"ERROR: BAM index file '{bai_file}' is missing!")
        sys.exit(1)

def verify_input_files_list(files_list):
    '''
    Check to make sure input files meet criteria.
    Add more criteria as issues are found.
    '''
    for file in files_list:
        # Check if the file exists
        if not os.path.isfile(file):
            print(f"ERROR: Input file '{file}' does not exist!")
            sys.exit(1)

        # For BAM files, ensure the corresponding BAI file exists
        if file.endswith(".bam"):
            check_for_bai(file)

def subprocess_cmd(command):
    '''
    Runs a terminal command with stdout piping enabled.
    '''
    process = sp.Popen(command, stdout=sp.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)

def make_chrom_region_list(region_file):
    '''
    Creates a list of tuples representing the regions from the BED file [(chrom, start, stop), ...].
    '''
    region_list = []
    with open(region_file) as f:
        for line in f:
            chrom, start, stop = line.split()[0:3]
            region_list.append((chrom, start, stop))
    return region_list

def make_IGV_chrom_loc(region):
    '''
    Return a chrom location string in IGV format.
    '''
    chrom, start, stop = region[0:3]
    return f'{chrom}:{start}-{stop}'

def make_snapshot_filename(region, height, sample_name=None):
    '''
    Formats a filename for the IGV snapshot.
    Adds more useful context to filenames by including region information.
    '''
    chrom, start, stop = region[0:3]
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    if sample_name:
        return os.path.join(SNAPSHOT_DIR, f"{sample_name}_{chrom}_{start}_{stop}_region_h{height}.png")
    return os.path.join(SNAPSHOT_DIR, f"{chrom}_{start}_{stop}_region_h{height}.png")

def mkdir_p(path):
    '''
    Recursively create a directory using pathlib.
    '''
    Path(path).mkdir(parents=True, exist_ok=True)

def get_open_X_server():
    '''
    Search for an open Xvfb port to render into.
    '''
    x_serv_command = '''
    for serv_num in $(seq 1 1000); do
        if ! (xdpyinfo -display :${serv_num}) &>/dev/null; then
            echo "$serv_num" && break
        fi
    done
    '''
    process = sp.Popen(x_serv_command, stdout=sp.PIPE, shell=True)
    output = process.communicate()[0].strip().decode('utf-8')

    # Handle if xdpyinfo returns unexpected output
    try:
        x_serv_port = int(output.split('\n')[0].strip())  # Take only the port number
    except ValueError:
        print(f"Unexpected xdpyinfo output: {output}")
        sys.exit(1)

    return x_serv_port

def write_IGV_script(input_files, region_file, IGV_batchscript_file, IGV_snapshot_dir, fasta_file, image_height, sample_name):
    '''
    Write out a batchscript for IGV.
    '''
    with open(IGV_batchscript_file, "w") as f:
        # Initialize IGV
        f.write("new\n")
        f.write(f"genome {fasta_file}\n")
        f.write(f"snapshotDirectory {IGV_snapshot_dir}\n")
        f.write(f"maxPanelHeight {image_height}\n")

        # Load BAM files
        for file in input_files:
            f.write(f"load {file}\n")

        # Write regions and snapshots
        region_list = make_chrom_region_list(region_file)
        for region in region_list:
            chrom_loc = make_IGV_chrom_loc(region)
            svsize = int(region[2]) - int(region[1])

            # For large regions, split snapshots into two: start and end
            if svsize > 10000:
                f.write(f"goto {region[0]}:{int(region[1]) - 1000}-{int(region[1]) + 500}\n")
                f.write(f"snapshot {sample_name}_{region[0]}_{region[1]}_start.png\n")
                f.write(f"goto {region[0]}:{int(region[2]) - 500}-{int(region[2]) + 1000}\n")
                f.write(f"snapshot {sample_name}_{region[0]}_{region[2]}_end.png\n")
            else:
                f.write(f"goto {chrom_loc}\n")
                f.write(f"snapshot {sample_name}_{region[0]}_{region[1]}_{region[2]}_region.png\n")

        f.write("exit\n")

def run_IGV_script(igv_script, igv_sh, memMB):
    '''
    Run an IGV batch script and ensure snapshots are generated.
    '''
    # Ensure the output directory exists
    mkdir_p(SNAPSHOT_DIR)

    # Get an open Xvfb port
    x_serv_port = get_open_X_server()
    print(f"Open Xvfb port found on: {x_serv_port}")

    # Build and run IGV command using igv.sh
    igv_command = f"xvfb-run --auto-servernum --server-num=1 bash {igv_sh} -b {igv_script}"
    print(f"IGV command: {igv_command}")

    # Record start time
    startTime = datetime.datetime.now()
    print(f"Started at: {startTime}")

    # Run the IGV command
    subprocess_cmd(igv_command)

    # Check if snapshots were generated
    snapshot_files = os.listdir(SNAPSHOT_DIR)
    if len(snapshot_files) == 0:
        print("ERROR: No snapshot files were generated.")
    else:
        print(f"Generated {len(snapshot_files)} snapshot files.")

    elapsed_time = datetime.datetime.now() - startTime
    print(f"Elapsed time: {elapsed_time}")

def main(input_files, region_file, fasta_file, image_height, igv_sh_bin, igv_mem, sample_name):
    '''
    Main control function for the script.
    '''
    batchscript_file = os.path.join(SNAPSHOT_DIR, "IGV_snapshots.bat")

    # Check if input files, regions, and IGV script exist
    file_exists(region_file, kill=True)
    file_exists(igv_sh_bin, kill=True)
    verify_input_files_list(input_files)

    # Check if the reference FASTA file and its index exist
    file_exists(fasta_file, kill=True)
    check_for_fai(fasta_file)

    print(f"\n~~~ IGV SNAPSHOT AUTOMATOR ~~~\n")
    print(f"Reference FASTA: {fasta_file}")
    print(f"Track height: {image_height}")
    print(f"IGV script file: {igv_sh_bin}")
    print(f"Batchscript file: {batchscript_file}")
    print(f"Region file: {region_file}")

    # Create output directory
    mkdir_p(SNAPSHOT_DIR)

    # Write the IGV batch script
    write_IGV_script(input_files=input_files, region_file=region_file, IGV_batchscript_file=batchscript_file, IGV_snapshot_dir=SNAPSHOT_DIR, fasta_file=fasta_file, image_height=image_height, sample_name=sample_name)

    # Run the IGV batch script
    run_IGV_script(igv_script=batchscript_file, igv_sh=igv_sh_bin, memMB=igv_mem)

def run():
    '''
    Parse script args to run the script.
    '''
    parser = argparse.ArgumentParser(description='IGV snapshot automator')
    parser.add_argument("input_files", nargs='+', help="Paths to the files to create snapshots from (e.g., .bam files).")
    parser.add_argument("-r", default=default_regions_bed, type=str, dest='region_file', help="BED file with regions to create snapshots over.")
    parser.add_argument("-f", "--fasta_file", required=True, help="Reference FASTA file to use.")
    parser.add_argument("-ht", default='500', type=str, dest='image_height', help="Height for the IGV tracks.")
    parser.add_argument("-bin", default=default_igv_sh, type=str, dest='igv_sh_bin', help="Path to the IGV sh binary to run.")
    parser.add_argument("-mem", default="4000", type=str, dest='igv_mem', help="Amount of memory to allocate to IGV, in Megabytes (MB).")
    parser.add_argument("--sample_name", required=True, help="Sample name to include in snapshot filenames.")

    args = parser.parse_args()

    main(input_files=args.input_files, region_file=args.region_file, fasta_file=args.fasta_file, image_height=args.image_height, igv_sh_bin=args.igv_sh_bin, igv_mem=args.igv_mem, sample_name=args.sample_name)

if __name__ == "__main__":
    run()

