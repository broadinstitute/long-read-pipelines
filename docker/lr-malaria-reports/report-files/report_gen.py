import os

# Data manipulation
import numpy as np
import pandas as pd
import random as rnd
import tabulate
import re
from enum import Enum
import glob

# Plotting
import matplotlib.pyplot as plt
import matplotlib
import plotly.express as px
import folium
import plotly.graph_objects as go
import math

# HTML
from markupsafe import Markup
import jinja2
import io
from io import StringIO
import base64
import urllib.parse

# Argument Parsing
import argparse

''' 
matplotlib render settings
'''
# Make big figures:
gFIG_SIZE_in = [14, 10]

# Set plotting defaults:
gPLOT_PARAMS = {
    "legend.fontsize": "x-large",
    "figure.figsize": gFIG_SIZE_in,
    "axes.labelsize": "x-large",
    "axes.titlesize": "x-large",
    "xtick.labelsize": "x-large",
    "ytick.labelsize": "x-large"
}
matplotlib.rcParams.update(gPLOT_PARAMS)

# Some single-place definitions of sizes for plots / figures:
gFONT_SIZE_UNITS = "pt"
gTITLE_FONT_SIZE = 36
gAXIS_LABEL_FONT_SIZE = 24
gTICK_LABEL_FONT_SIZE = 16
gTEXT_FONT_SIZE = 16

# To track global figure number creation:
gFIG_NUM = 0

def fix_plot_visuals(fig,
                     titlesize=gTITLE_FONT_SIZE,
                     labelsize=gAXIS_LABEL_FONT_SIZE,
                     ticklabelsize=gTICK_LABEL_FONT_SIZE,
                     textsize=gTEXT_FONT_SIZE,
                     tight_rect=None):
    """Fix the plot elements to be appropriate sizes for a slide / presentation."""

    if not textsize:
        textsize = ticklabelsize

    for ax in fig.get_axes():

        for ticklabel in (ax.get_xticklabels()):
            ticklabel.set_fontsize(ticklabelsize)
        for ticklabel in (ax.get_yticklabels()):
            ticklabel.set_fontsize(ticklabelsize)
        for c in ax.get_children():
            if c.__class__ == matplotlib.text.Text:
                c.set_fontsize(textsize)

        ax.xaxis.get_label().set_fontsize(labelsize)
        ax.yaxis.get_label().set_fontsize(labelsize)
        ax.title.set_fontsize(titlesize)

    for c in fig.get_children():
        if c.__class__ == matplotlib.legend.Legend:
            c.prop.set_size(ticklabelsize)
            c.get_title().set_size(ticklabelsize)

    if tight_rect:
        fig.tight_layout(rect=tight_rect)
    else:
        fig.tight_layout()
    
    if fig._suptitle:
        sup_title = fig._suptitle.get_text()
        fig.suptitle(sup_title, fontsize=titlesize)
    
    # Make it so we can actually see what's happening on the plots with "dark mode":
    fig.patch.set_facecolor("white")


'''
Coverage Plot
'''
def plot_coverage(directory, sample_name, bin_width=500):
    
    # only looking for .bed.gz files
    ext = (".bed.gz")

    # Set up data
    beds = pd.DataFrame(columns=["stop", "depth"])
    sorted_file_list = list()
    
    # Reading and sorting files if coverage data is available
    if os.path.exists(directory):
        file_list = os.listdir(directory)
        sorted_file_list = sorted(file_list,  key=lambda item: item.split(".")[-4].split("_")[1] if item.endswith(".bed.gz") else "00")
        print(f"Files to be plotted: {sorted_file_list}")
    else:
        return None
    
    # Plot setup
    fig, ax = plt.subplots(1,1, figsize=(15, 9))
    plt.title(f"Coverage Plot of Sample {sample_name}", pad=12, fontsize=12)
    plt.xlabel("Contig (bp)", labelpad=10, fontsize=12)
    plt.ylabel("Depth", labelpad=10, fontsize=12)
    color = "#2494C7"
    tick_positions = []
    contigs = []
    bin_max = 0
 
    # Iterate over each bed.gz in test_data folder
    for idx, file in enumerate(sorted_file_list):

        if file.endswith(ext):
            # Reading current .bed.gz file
            f = os.path.join(directory, file)
            #print(f"\n{idx} Current file: {f}")

            # Create DataFrame per .bed.gz
            bed_df = pd.read_csv(f, delimiter="\\t", names=["contig", "start", "stop", "depth"], engine="python")
            
            # Get bins
            start_min = bed_df["start"].min()
            stop_max = bed_df["stop"].max()
            bins = np.arange(start_min, stop_max+1, bin_width)
            values = np.zeros(len(bins)) 

            # Iterrate through each DataFrame and populate bins
            #print("Populating bins...")
            for _, row in bed_df.iterrows():
                avg = (row["stop"]+row["start"])/2
                index = int(np.floor(avg/bin_width))
                values[index]+=row["depth"]

            # Append new data to DF
            #print("Plotting data...")
            if(idx == 0):
                ax.plot(bins, values, ".", c=color)
                tick_positions.append(math.floor(stop_max/2)+bin_max)
                bin_max = max(bed_df["stop"])
            else:
                ax.plot(bins+bin_max, values, ".", c=color)
                tick_positions.append(math.floor(stop_max/2)+bin_max)
                bin_max = bin_max + max(bed_df["stop"])                      
                               
            # Saving xtick data
            contigs.append(bed_df["contig"][0])
                
    # Setting xticks
    ax.set_xticks(ticks=tick_positions)
    ax.set_xticklabels(labels=contigs, rotation=90)
    ax.set_yscale("log")
    fix_plot_visuals(fig)
            
    return fig
    
'''
Drug Resistance Table
'''
def create_drug_table(file):
    if not file:
        resistances = ["UNDETERMINED","UNDETERMINED","UNDETERMINED","UNDETERMINED","UNDETERMINED","UNDETERMINED"]
    else:
        data = open(file, 'r').read()
        resistances = list(get_drug_resistance(data, None, None, do_print=True))
    
    resistances_tbl = pd.DataFrame(columns = ["Chloroquine", "Pyrimethamine", "Sulfadoxine", "Mefloquine", "Artemisinin", "Piperaquine"])
    resistances = map(str, resistances)
    resistances = [s.replace('DrugSensitivity.','') for s in list(resistances)]
    
    return resistances

# Set up some functions to determine drug sensitivity:
# vars for markers in drug_resistance_report.txt
present = '+'
absent = '-'

pchange_regex = re.compile(r"""p\.([A-z]+)([0-9]+)([A-z]+)""")
AA_3_2 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

DrugSensitivity = Enum("DrugSensitivity", ["UNDETERMINED", "SENSITIVE", "RESISTANT"])

def parse_pchange(pchange):
    old, pos, new = pchange_regex.match(pchange).groups()
    return AA_3_2[old.upper()], int(pos), AA_3_2[new.upper()]

def get_chloroquine_sensitivity(dr_report):
# Locus utilized: PF3D7_0709000 (crt)
# Codon: 76
# Workflow:
# Step Genetic change Interpretation Classification
# 1 76 K/T heterozygote Heterozygous mutant Undetermined
# 2 76 missing Missing Undetermined
# 3 K76 Wild type Sensitive
# 4 76T Mutant Resistant
# 5 otherwise Unknown mutant Undetermined

    for line in StringIO(dr_report):
        line = ' '.join(line.split())
        # We only care about this locus for this drug:
        if line.startswith("pfcrt PF3D7_0709000"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            if pos == 76:
                if (old == "L") and (new == "T") and (marker == absent):
                    return DrugSensitivity.SENSITIVE
                elif (new == "T") and (marker == present):
                    return DrugSensitivity.RESISTANT
                
    return DrugSensitivity.UNDETERMINED
                
def get_pyrimethamine_sensitivity(dr_report):
# Locus utilized: PF3D7_0417200 (dhfr)
# Codon: 108
# Workflow:
# Step Genetic change Interpretation Classification
# 1 108 S/N heterozygote Heterozygous mutant Undetermined
# 2 108 missing Missing Undetermined
# 3 S108 Wild type Sensitive
# 4 108N Mutant Resistant
# 5 otherwise Unknown mutant Undetermined

    for line in StringIO(dr_report):
        line = ' '.join(line.split())
        # We only care about this locus for this drug:
        if line.startswith("pfdhfr PF3D7_0417200"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            if pos == 108:
                if (old == "S") and (new == "N") and (marker == absent):
                    return DrugSensitivity.SENSITIVE
                elif (old == "S") and (new == "N") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif (new == "N") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif marker == absent:
                    return DrugSensitivity.SENSITIVE
                
    return DrugSensitivity.UNDETERMINED

def get_sulfadoxine_sensitivity(dr_report):
# Locus utilized: PF3D7_0810800 (dhps)
# Codon: 437
# Workflow:
# Step Genetic change Interpretation Classification
# 1 437 A/G heterozygote Heterozygous mutant Undetermined
# 2 437 missing Missing Undetermined
# 3 A437 Wild type Sensitive
# 4 437G Mutant Resistant
# 5 otherwise Unknown mutant Undetermined

    for line in StringIO(dr_report):
        line = ' '.join(line.split())

        # We only care about this locus for this drug:
        if line.startswith("pfdhps PF3D7_0810800"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            if pos == 437:
                if (old == "A") and (new == "G") and (marker == absent):
                    return DrugSensitivity.SENSITIVE
                elif (old == "A") and (new == "G") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif (new == "G") and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif marker == absent:
                    return DrugSensitivity.SENSITIVE
                
    return DrugSensitivity.UNDETERMINED

def get_mefloquine_sensitivity(dr_report):
# Locus utilized: PF3D7_0523000 (mdr1)
# Codons: Amplification status of whole gene
# Workflow:
# Step Genetic change Interpretation Classification
# 1 Missing Missing Undetermined
# 2 Heterozygous duplication Heterozygous mutant Undetermined
# 3 Single copy Wild type Sensitive
# 4 Multiple copies Mutant Resistant

    # Currently we can't determine this.
    # We need to get CNV calling working first.
                
    return DrugSensitivity.UNDETERMINED

def get_artemisinin_sensitivity(dr_report):
# Locus utilized: PF3D7_1343700 (kelch13)
# Codons: 349-726 (BTB/POZ and propeller domains)
# Workflow:
# Step Genetic change Interpretation Classification
# 1 Homozygous non-synonymous mutations in the kelch13 BTB/POZ and propeller
# domain classified by the World Health Organisation as associated with delayed
# parasite clearance 
# Mutant – associated with delayed clearance Resistant
# 2 Heterozygous non-synonymous mutations in the kelch13 BTB/POZ and
# propeller domain classified by the World Health Organisation as associated
# with delayed parasite clearance
# Mutant - heterozygous Undetermined
# 3 578S as homozygous Mutant - not associated Sensitive
# 4 Any missing call in amino acids 349-726 Missing Undetermined
# 5 No non-reference calls in amino acids 349-726 Wild-type Sensitive
# 6 otherwise Mutant - not in WHO list Undetermined

    has_variants = False

    for line in StringIO(dr_report):
        line = ' '.join(line.split())
        # We only care about this locus for this drug:
        if line.startswith("pfkelch13 PF3D7_1343700"):
            gene, locus, pchange, marker = line.strip().split(" ")
            old, pos, new = parse_pchange(pchange)
            
            has_non_ref = False
            has_variants = False
            if 349 <= pos <= 726:
                if (old != new) and (marker == present):
                    return DrugSensitivity.RESISTANT
                elif (new == "S") and (marker == present):
                    return DrugSensitivity.SENSITIVE
                elif marker == present:
                    has_non_ref = True
                has_variants = True
    if (has_variants) and (not has_non_ref):
        return DrugSensitivity.SENSITIVE
    
    return DrugSensitivity.UNDETERMINED

def get_piperaquine_sensitivity(dr_report):
# Loci utilized: PF3D7_1408000 (plasmepsin 2) and PF3D7_1408100 (plasmepsin 3)
# Codons: Amplification status of both genes
# Workflow:
# Step Genetic change Interpretation Classification
# 1 Missing Missing Undetermined
# 2 Heterozygous duplication Heterozygous mutant Undetermined
# 3 Single copy Wild type Sensitive
# 4 Multiple copies Mutant Resistant

    # Currently we can't determine this.
    # We need to get CNV calling working first.
                
    return DrugSensitivity.UNDETERMINED

def get_drug_resistance(data, sample_id, sample_df, do_print=False):
    
    chloroquine = get_chloroquine_sensitivity(data)
    pyrimethamine = get_pyrimethamine_sensitivity(data)
    sulfadoxine = get_sulfadoxine_sensitivity(data)
    mefloquine = get_mefloquine_sensitivity(data)
    artemisinin = get_artemisinin_sensitivity(data)
    piperaquine = get_piperaquine_sensitivity(data)
        
    return chloroquine, pyrimethamine, sulfadoxine, mefloquine, artemisinin, piperaquine
 
'''
Map
'''
def create_map(coordinates, sample_name):
    m = folium.Map(location=coordinates, zoom_start = 5)

    if coordinates != [0,0]:
        folium.Marker(location=coordinates, popup = ('Sample: '+sample_name), icon=folium.Icon(color='red',prefix='fa',icon='circle'), parse_html=True).add_to(m)

    m.get_root().width = "473px"
    m.get_root().height = "397px"
    map_html = m.get_root()._repr_html_()
    
    #with open('/templates/map.html', mode='w+') as f:
    #    f.write(map_html)
    
    return map_html

'''
Quality Report (FastQC or ONT QC)
'''
def read_qc_report(directory):
    '''
    Function to read HTML from FastQC or ONT QC report as string to pass to Jinja template IFrame
    '''
    if(os.path.exists(directory)):
        # There should only be one file in the directory (report-files/data/quality_report)
        for file in os.listdir(directory):
            print(f"{directory} + {file}")
            with open(os.path.join(directory, file), "r", encoding="utf-8") as f:
                html = f.read()
    else:
        html = None
    
    return html

''' 
Snapshots
'''
def get_res_loci_names(bed_file):
    '''
    Input: Bed file containing drug resistance loci and their names.
    Output: A list of the fourth column of the bed file containing the names of all the drug resistance regions.
    '''
    column_names = ["chromosome", "start", "end", "name"]
    bed = pd.read_csv(bed_file, sep="\s+", header=None, names=column_names) 
    print(bed.name)
    return list(bed.name) 

'''
Utility
'''
def plot_to_b64(plot):
    '''
    Input: A matplotlib plot.
    Output: The b64 string conversion of the input plot.
    '''
    # Create I/O buffer to save image in bytes
    stringIObytes = io.BytesIO()
    plot.savefig(stringIObytes, format="jpeg")
    
    # Retrieve byte-string and encode it in base64
    stringIObytes.seek(0)
    plot_b64 = base64.b64encode(stringIObytes.read())
    
    return plot_b64

def img_to_b64(img):    
    with open(img, "rb") as image_file:
        return base64.b64encode(image_file.read())

def format_dates(date_string):
    '''
    Input: A string of numbers in the form of YYYYMMDD, such as 20240528 (May 28th, 2024).
    Output: A formatted version of the input string using slashes to separate year, month, and day or "N/A".
    '''
    if date_string in ["", None, "N/A"]:
        return "N/A"
    else:
        return "/".join([date_string[:4], date_string[4:6], date_string[6:]])

def check_na(data):
    '''
    Input: A variable of any type to check if it is valid for passing to the report.
    Output: The variable or "N/A".
    '''
    return "N/A" if data in [None, "None", 0] else data

'''
Report Generation & Templating
'''
class Sample:
    '''
    This class holds all variables that will be used on the summary page of the report.
    
    Additionally, this class holds variables that are used on both the summary and analysis pages,
    such as the sample name.
    '''
    
    def __init__(self, sample_name, hrp2, hrp3, drug_res, info, _map, location_info, qc_pass):
        '''This function defines the class variables and retrieves them from their respective functions.'''
        
        self.sample_name = sample_name
        self.hrp2 = hrp2
        self.hrp3 = hrp3
        self.drug_res = drug_res
        self.info = info
        self.map = _map
        self.location_info = location_info
        self.qc_pass = qc_pass

class Analysis:
    '''
    This class holds all variables used on the analysis page of the report.
    '''
    
    def __init__(self, sequencing_summary, qscore_reads, qscore_scores, barcode, coverage_plot, snapshots, res_loci_names):
        '''This function defines the variables used and retrieves them from their respective functions.'''
        
        self.sequencing_summary = sequencing_summary
        self.scores = qscore_scores
        self.reads = qscore_reads
        self.barcode = barcode
        self.coverage_plot = coverage_plot
        self.snapshots = snapshots
        self.res_loci_names = res_loci_names

class QCReport:
    def __init__(self, qc_report_html):
        self.qc_report_html = qc_report_html

def prepare_summary_data(arg_dict):
    '''
    Gathers all data needed for summary page from inputs.
    '''
    
    sample_name = arg_dict['sample_name']

    upload_date = arg_dict['upload_date'][0]
    collection_date = arg_dict['collection_date']
    sequencing_date = arg_dict['sequencing_date']
    species = ' '.join(arg_dict['species'])
    
    info = [upload_date, format_dates(collection_date), format_dates(sequencing_date), species, round(arg_dict['aligned_coverage'], 2), check_na(arg_dict['aligned_read_length']), 
            check_na(arg_dict['pct_properly_paired_reads']), arg_dict['read_qual_median'], round(arg_dict['read_qual_mean'], 2)]
    processed_info = [item if item not in ["", None] else "N/A" for item in info]
    
    qc_pass = arg_dict["qc_pass"]
    if (qc_pass == "true"):
        qc_pass = "PASS"
    elif (qc_pass == "false"):
        qc_pass = "FAIL"

    resistances = create_drug_table(None if arg_dict["drug_resistance_text"] in [None, "None", ""] else arg_dict["drug_resistance_text"])

    location_info = [round(arg_dict['latitude'], 2), round(arg_dict['longitude'], 2), arg_dict['location']]
    coordinates = [arg_dict['latitude'], arg_dict['longitude']]
    _map = create_map(coordinates, sample_name)

    HRP2 = arg_dict['HRP2']
    HRP3 = arg_dict['HRP3']

    return Sample(sample_name, HRP2, HRP3, resistances, processed_info, _map, location_info, qc_pass)

def prepare_analysis_data(arg_dict):
    '''
    Gathers all data needed for analysis page from inputs.
    '''
    
    frac_bases = check_na(arg_dict["fraction_aligned_bases"])
        
    sequencing_summary = [arg_dict['sample_type'], arg_dict['analysis_success'], arg_dict['aligned_bases'], arg_dict['aligned_reads'], 
                          round(frac_bases, 2), round(arg_dict['average_identity'], 2)]

    barcode = arg_dict['barcode']
        
    qscorex = [5, 7, 10, 12, 15] # available q-score measures are predetermined
    qscorey = [arg_dict['num_reads_q5'], arg_dict['num_reads_q7'], arg_dict['num_reads_q10'], arg_dict['num_reads_q12'], arg_dict['num_reads_q15']]

    # Create coverage plot and convert it to base64
    coverage_bin_size = arg_dict["coverage_bin_size"]
    coverage_plot = plot_coverage("/report-files/data/coverage", arg_dict["sample_name"], coverage_bin_size) # default bin size = 750
    coverage_b64 = plot_to_b64(coverage_plot)
    
    # IGV Snapshots
    snapshots = str(arg_dict["snapshots"])
    
    snapshots = snapshots.split(",")
    snapshots_b64 = [img_to_b64(image) for image in snapshots]
    
    bed_file = open(arg_dict["regions_bed"], "r")
    
    return Analysis(sequencing_summary, qscorey, qscorex, barcode, coverage_b64, snapshots_b64, get_res_loci_names(bed_file))

def create_report(sample, analysis, qc_report):
    '''
    This function handles the final steps of report generation.
    
    It pools all variables and data needed for the report and creates two objects:
    analysis and summary. These objects hold the variables and are passed to their respective
    page templates (analysis or summary).
    '''
    
    print("path: "+os.path.dirname(os.path.realpath(__file__)))
    print("cwd: "+os.getcwd())

    # creating summary page
    templateLoader = jinja2.FileSystemLoader(searchpath='/report-files/templates/')
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = 'report.html' # may need to change if file is moved
    template = templateEnv.get_template(TEMPLATE_FILE)
    output = template.render(sample=sample, analysis=analysis, qc_report = qc_report)

    file_name = os.path.join(os.getcwd(),sample.sample_name+'_lrma_report.html')
    print(file_name)

    with open(file_name, 'w') as fh:
        fh.write(output)  
        
    print('Report generated!')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    ''' Summary Page '''
    parser.add_argument("--barcode", help="barcode of the sample", default="N/A")
    
    # Sample Info
    parser.add_argument("--sample_name", help="name of sequenced sample", required=True)
    parser.add_argument("--upload_date", help="date sample was uploaded", nargs='+', required=True)
    parser.add_argument("--collection_date", help="date sample was collected", required=True, default="N/A", nargs="?")
    parser.add_argument("--sequencing_date", help="date sample was sequenced", default="N/A", nargs="?")
    parser.add_argument("--species", help="species of sample", nargs='+', default="P. falciparum")
    parser.add_argument("--aligned_coverage", help="number of times the bases in the sequenced reads cover the target genome", required=True, type=float)
    parser.add_argument("--aligned_read_length", help="number at which 50%\ of the read lengths are longer than this value", type=float)
    parser.add_argument("--pct_properly_paired_reads", help="median read length", type=float)
    parser.add_argument("--read_qual_median", help="median measure of the uncertainty of base calls", required=True, type=float, default="N/A", nargs="?")
    parser.add_argument("--read_qual_mean", help="mean measure of the uncertainty of base calls", required=True, type=float, default="N/A", nargs="?")

    # Drug Resistance
    parser.add_argument("--drug_resistance_text", help="path of text file used for determining and displaying drug resistances", default=None)
    parser.add_argument("--HRP2", help="value denoting whether the HRP2 marker is present or not -- true or false", default="N/A")
    parser.add_argument("--HRP3", help="value denoting whether the HRP3 marker is present or not -- true or false", default="N/A")

    # Map
    parser.add_argument("--longitude", help="longitude value of where the sample was collected", type=float, default=0)
    parser.add_argument("--latitude", help="latitude value of where the sample collected", type=float, default=0)
    parser.add_argument("--location", help="location where the sample was collected", default="N/A")
    
    # QC Status
    parser.add_argument("--qc_pass", help="status to determine whether or not the sequencing run passes quality control standards", required=True)

    ''' Analysis Page '''
    # Q-Scores Plot
    parser.add_argument("--num_reads_q5", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 3", required=True)
    parser.add_argument("--num_reads_q7", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 5", required=True)
    parser.add_argument("--num_reads_q10", help="the number of reads where the probability of a given base call being wrong is 1 in 10", required=True)
    parser.add_argument("--num_reads_q12", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 16", required=True)
    parser.add_argument("--num_reads_q15", help="the number of reads where the probability of a given base call being wrong is approximately 1 in 32", required=True)

    # Sequencing Summary
    parser.add_argument("--sample_type", help="type of sample", default="N/A", nargs="?")
    parser.add_argument("--analysis_success", help="whether the analysis process completed successfully", required=True)
    parser.add_argument("--aligned_bases", help="total number of bases aligned to the reference genome", required=True)
    parser.add_argument("--aligned_reads", help="total number of reads aligned to the reference genome", required=True)
    parser.add_argument("--fraction_aligned_bases", help="number of bases aligned out of all bases sequenced", required=True, type=float)
    parser.add_argument("--average_identity", help="", required=True, type=float) # check
    
    # Coverage Plot
    parser.add_argument("--fastqc_path", help="location of fastqc_report file; used to locate BAM files for coverage report generation")
    parser.add_argument("--coverage_bin_size", help="number to use as size of bins for coverage plot generation; default is 1500", default=500, required=False, type=int)
    
    # Snapshots
    parser.add_argument("--snapshots", required=True)
    parser.add_argument("--regions_bed")
    
    # parse given arguments
    args = parser.parse_args()
    arg_dict = vars(args)

    qc_report_html = read_qc_report("/report-files/data/quality_report")
    
    summary = prepare_summary_data(arg_dict)
    analysis = prepare_analysis_data(arg_dict)
    qc_report = QCReport(qc_report_html)

    create_report(summary, analysis, qc_report)




