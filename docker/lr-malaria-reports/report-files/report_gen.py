# Data manipulation
import numpy as np
import pandas as pd
import random as rnd
import tabulate
import re
from enum import Enum
from io import StringIO

# Plotting
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import plotly.express as px
import folium
import plotly.graph_objects as go

# HTML
from markupsafe import Markup
import jinja2


'''
Coverage Plot
'''
def get_coverage_pos():
    #load in phage genome fasta
    fasta_file = open("Plasmodium_falciparum.ASM276v2.dna.chromosome.1.fa")
    ref_seq = fasta_file.read().split('\n', 1)[1] #remove the header
    ref_seq = ref_seq.rstrip('\n') #remove trailing newline
    ref_seq = ref_seq.replace('\n','') #remove middle newlines
    ref_seq_list = list(ref_seq)
    len(ref_seq_list)

    #generate 5000 150bp reads - sequence positions are 1-based as in SAM/BAM files
    rnd.seed(42)
    start_pos_all = np.array(rnd.choices(range(1,len(ref_seq_list)-150), k = 5000)) #sample with replacement
    end_pos_all = start_pos_all + 150 #add 150bp to sstart to get where end of the read maps

    print(f'total reads: {len(start_pos_all)}')
    print(f'covered sequence: {np.min(start_pos_all)} - {np.max(end_pos_all)}bp')

def get_sequence_coverage1(start_pos_all:np.array, end_pos_all:np.array, start_region:int, end_region:int):
    '''
    Get coverage over a region of the sequence by jumping between 
    start and end positions of reads, incrementing coverage as you go. 
    Infers that until you encounter the start or end of a read, you
    do not change the coverage value, and therefore you avoid looping 
    over all bp of interest, saving time.
    
    start_pos_all: array of all sequence start positions for reads 
    end_pos_all: array of all sequence end positions for reads. Must be greater than start_pos_all for a single read
    start_region: starting bp of region of interest for plotting (int). min value = 1
    end_region: end bp of region of interest for plotting (int). Greater than start_region
    
    Requires numpy (import as np)
    '''
    #start and end values must be sorted low to high
    start_pos = start_pos_all[~(end_pos_all < start_region)] #ignore all reads that end before region of interest
    end_pos = end_pos_all[~(end_pos_all < start_region)]

    start_pos.sort()
    end_pos.sort()
    
    start_pos[start_pos < start_region] = start_region #truncate reads that start before region of interest
    
    #moving pointer for sequence position
    seq_pos = start_region
    #running coverage value
    cov = 0
    #initialize empty array to store coverage values
    coverage = np.zeros(end_region-(start_region-1))
    #moving pointer for index in coverage array
    idx = 0
    #initialize indices at 0
    i=0
    j=0
    
    #while within the index range of the reads
    while (j < len(start_pos) and i < len(end_pos)): 
        #if you haven't hit the start of the last read and encounter the start position of a read first
        if (j < (len(start_pos) - 1) and start_pos[j] < end_pos[i]):
            seq_pos_jump = start_pos[j] #jump to start_pos[j]
            j += 1 #increment index
            cov_increment = 1 #increment coverage value
        #if you encounter the end position of a read first or you've reached the start of the last read
        elif (start_pos[j] > end_pos[i] or j == len(start_pos)-1):
            seq_pos_jump = end_pos[i] #jump to end_pos[i]
            i += 1 #increment index
            cov_increment = -1 #decrement coverage value    
        #if you encounter start and end position of a read at the same position   
        else:
            seq_pos_jump = start_pos[j]+1 #jump by 1bp
            cov_increment = 1 #increment coverage value
            j+=1
        #store the constant coverage value from seq_pos to seq_pos_jump, correcting for 0-indexing
        jump = idx + (seq_pos_jump - seq_pos)  
        if jump <= len(coverage):
            coverage[idx:jump] = [cov] * (seq_pos_jump - seq_pos)
        else:
            coverage[idx:] = [cov] * (len(coverage) - idx) #if you run out of room in the array
            break

        idx = jump #update index
        cov += cov_increment #update coverage value
        seq_pos = seq_pos_jump #update seq position
    
    return coverage

def get_sequence_coverage2(start_pos_all:np.array, end_pos_all:np.array, start_region:int, end_region:int):
    '''
    Get coverage over a region of the sequence by iterating through
    the reads and incrementing coverage in the region of the genome
    that the read covers. 

    start_pos_all: array of all sequence start positions for reads 
    end_pos_all: array of all sequence end positions for reads. Must be greater than start_pos_all for a single read
    start_region: starting bp of region of interest for plotting (int). min value = 1
    end_region: end bp of region of interest for plotting (int). Greater than start_region
    
    Requires numpy (import as np)
    '''
    #initialize empty array to store coverage values
    genome_coverage = np.zeros(end_region-(start_region-1))

    #pick reads that fall within the region of interest
    #ignore all reads that end before region of interest and start after region of interest
    start_pos = start_pos_all[~((end_pos_all < start_region) & (start_pos_all > end_region))]
    end_pos = end_pos_all[~((end_pos_all < start_region) & (start_pos_all > end_region))]

    #iterate through reads and increment coverage
    for i in range(0, len(start_pos)):
        genome_coverage[start_pos[i]:end_pos[i]] += 1
    
    return genome_coverage

def plot_coverage(start_pos_all:np.array, end_pos_all:np.array, start_region:int, end_region:int, version = 1):
    
    '''
    Plots output of get_sequence_coverage(). Requires matplotlib.pyplot (plt)
    version: 1 uses get_sequence_coverage1, 2 uses get_sequence_coverage2
    '''
    if version == 2:
        coverage = get_sequence_coverage2(start_pos_all, end_pos_all, start_region, end_region)
    else:
        coverage = get_sequence_coverage1(start_pos_all, end_pos_all, start_region, end_region)
        
    #plotting
    fig, ax = plt.subplots(1,1, figsize=(20,6))
    ax.bar(range(start_region, end_region+1), coverage, width=1, \
           color = 'silver', edgecolor = 'dimgray',linewidth = 1,align='center')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
    ax.set_xlabel('Sequence Position', fontsize = 18)
    ax.set_ylabel('Coverage', fontsize = 18)
    ax.tick_params(labelsize = 16)
    plt.margins(0,0.05)



'''
Drug Resistance Table
'''
def create_drug_table(file):
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
    # Get the GSURI to the drug resistance report:
    #row = reportable_samples.loc[reportable_samples["entity:sample_id"] == sample_id]
    #dr_report_gs_url = row.iloc[0]["drug_res_report"]

    #blob = storage.Blob.from_string(dr_report_gs_url)
    
    # Download the file contents:
    #dr_report_contents = blob.download_as_text(client=my_storage_client)
    
    chloroquine = get_chloroquine_sensitivity(data)
    pyrimethamine = get_pyrimethamine_sensitivity(data)
    sulfadoxine = get_sulfadoxine_sensitivity(data)
    mefloquine = get_mefloquine_sensitivity(data)
    artemisinin = get_artemisinin_sensitivity(data)
    piperaquine = get_piperaquine_sensitivity(data)

   # if do_print:
   #     print("\t".join([
   #         "Sample\t", "Chloroquine", "Pyrimethamine", "Sulfadoxine", "Mefloquine", "Artemisinin", "Piperaquine"
   #     ]))
   #     print(f"{sample_id}\t", end='')
   #     for r in [chloroquine, pyrimethamine, sulfadoxine, mefloquine, artemisinin, piperaquine]:
   #         print(f"{r.name}\t", end='')
   #     print()
        
    return chloroquine, pyrimethamine, sulfadoxine, mefloquine, artemisinin, piperaquine
 
 
'''
Map
'''
def create_map(coordinates, sample_name):
    m = folium.Map(location=coordinates, zoom_start = 5)
    folium.Marker(location=coordinates, popup = ('Sample: '+sample_name), icon=folium.Icon(color='red',prefix='fa',icon='circle'), parse_html=True).add_to(m)
    m.get_root().width = "473px"
    m.get_root().height = "397px"
    map_html = m.get_root()._repr_html_()
    
    with open('github/lr-malaria-automation-report/report_files/templates/map.html', mode='w') as f:
        f.write(map_html)
    
    return map_html


'''
Q-Score Plot
'''
def create_qscore_plot():
    '''
    Function to create a q-score index vs. reads plot
    
    Plot is returned as an HTML file. Currently, this plot is written
    in Javascript. This function is unused in report generation.
    '''
    
    x=[5, 7, 10, 12, 15] # index
    y=[3413205, 3413204, 3413073, 3062945, 2120402] # q scores
    
    df = pd.DataFrame({
        'Q Score': x,
        'Number of Reads': y
    })
    
    fig1 = px.line(df,x="Q Score",y="Number of Reads", width=450, height=400)
    
    outpath="github/lr-malaria-automation-report/report_files/templates/qscore_plot.html"
    return fig1.write_html(outpath,
                full_html=False,
                include_plotlyjs=True)


'''
Report Generation & Templating
'''
class Sample:
    '''
    This class holds all variables that will be used on the summary page of the report.
    
    Additionally, this class holds variables that are used on both the summary and analysis pages,
    such as the sample name.
    '''
    
    def __init__(self, sample_name, hrp2, hrp3, qc_status, drug_res, info, _map, location_info):
        '''This function defines the class variables and retrieves them from their respective functions.'''
        
        self.sample_name = sample_name
        self.hrp2 = hrp2
        self.hrp3 = hrp3
        self.qc_status = qc_status
        self.drug_res = drug_res
        self.info = info
        self.map = _map
        self.location_info = location_info

# %%
class Analysis:
    '''
    This class holds all variables used on the analysis page of the report.
    
    Additionally, it passes in variables needed for plot generation in Javascript.
    '''
    
    def __init__(self, sequencing_summary, qscore_reads, qscore_scores, reference_info):
        '''This function defines the variables used and retrieves them from their respective functions.'''
        
        self.sequencing_summary = sequencing_summary
        self.scores = qscore_scores
        self.reads = qscore_reads
        self.reference_info = reference_info


def create_report(sample, analysis):
    '''
    This function handles the final steps of report generation.
    
    It pools all variables and data needed for the report and creates two objects:
    analysis and summary. These objects hold the variables and are passed to their respective
    page templates (analysis or summary).

    These values will later be taken from the Terra workspace. 
    QC_status will be determined via an undetermined method.
    '''

    ''' Summary Page Variables '''
    hb3_info = ['2021-07-16', 'P. falciparum', '210', '5,394'+' bp', '1,294'+' bp','17.1']
    hb3_resistances = create_drug_table('drug_resistance_report.txt')
    coordinates = [14.5,14.5]
    _map = create_map(coordinates, 'HB3')
    coordinates.append('Doundodji, Senegal')
    hb3_summary = Sample('hb3', 'Positive', 'Negative', 'PASS', hb3_resistances, hb3_info, _map, coordinates)

    '''Analysis Page Variables'''
    hb3_sequencing_summary = ['sWGA','Pass','15,650,871,313 bp','13,381,441 bp', '0.3', '92.6']
    with open('github/lr-malaria-automation-report/report_files/templates/qscore_plot.html', 'r') as file:
        graph_html = file.read()
    hb3_qscorex = [5, 7, 10, 12, 15] # index
    hb3_qscorey = [3413205, 3413204, 3413073, 3062945, 2120402] # q scores
    hb3_reference_info = ['P. falciparum', 'Genes', 'BAM']
    hb3_analysis = Analysis(hb3_sequencing_summary, hb3_qscorey, hb3_qscorex, hb3_reference_info)
        
    
    
    # creating summary page
    templateLoader = jinja2.FileSystemLoader(searchpath='./')
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = './summary.html' # may need to change if file is moved
    template = templateEnv.get_template(TEMPLATE_FILE)
    output = template.render(sample=sample)
    
    with open('./lrma_report_summary.html', 'w') as fh:
        fh.write(output)        
        templateLoader = jinja2.FileSystemLoader(searchpath='./')
        
    templateEnv = jinja2.Environment(loader=templateLoader)
    TEMPLATE_FILE = './analysis.html' # may need to change if file is moved
    template = templateEnv.get_template(TEMPLATE_FILE)
    output = template.render(sample=sample, analysis=analysis)
    
    with open('./lrma_report_analysis.html', 'w') as fh:
        fh.write(output)
        
    print('report generated')

if __name__ == '__main__':
    '''Reports will be generated in the current working directory.'''
    create_report(hb3_summary, hb3_analysis) 




