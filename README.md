# RDpkkd
a stubborn and ambitious child with curiosity in science and talent in spacing out to the imaginary possibilities of beyondness.
for the nanopore qualified sequences, quality check is done via NanoPlot software which i installed into my HPC environment env_0 by following below mentioned steps;
```
module load bioconda
conda activate env_0

pip install NanoPlot

pip install NanoPlot --upgrade


```
once installed then the following command is used to qc the samples
```
#!/bin/bash
#SBATCH --job-name=nanoplot
#SBATCH --account=da494
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00

### Load conda module and initialize
module load bioconda  # Adjust module name if necessary
conda init bash  # Initialize conda for bash shell

### Activate conda environment
source activate env_0  # Ensure correct environment activation

### Navigate to where the data is
cd /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs

NanoPlot --fastq FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.fastq.gz -o MSCprojectwork2024gutmicrobiomealopeciaareata
```
(-o 'directory_location' is to assign an output directory to the results genrated)  the above command will work if the samples are properly formatted as gzip compressed files, in my case the files were not in proper gzip format so i had to reformat them by using the following command;
```
mv FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.fastq.gz FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.fastq

```
the above command means the samples were gunzipped and renamed one by one and then again Gzipped by following the below command;
```
gzip *.fastq*

```
using the *.fastq* means to collectively apply same command for all the sample files ending with the same extension.
the NanoPlot files were successfully generated with outputs ending with .txt, .html, .png, these files included the information such as :
Mean read length:             
Mean read quality:               
Median read length:           
Median read quality:             
Number of reads:                
Read length N50:              
STDEV read length:            
Total bases:

Number, percentage and megabases of reads above quality cutoffs falling into a scale of :
>Q10:	
>Q15:	
>Q20:	
>Q25:	
>Q30:	

Top 5 highest mean basecall quality scores and their read lengths
and,
Top 5 longest reads and their mean basecall quality score

quality check fastqc was also done for the samples using the following command which was provided by Dr. Tedder (my project supervisor):
```
#!/bin/sh
#SBATCH --job-name=fastQC
#SBATCH --account=da494
#SBATCH --partition=compute
#
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00


### load conda environment

module load bioconda
conda activate env_0

### navigate to where the data is

cd /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs

### run fastQC

fastqc FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.fastq.gz -o /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs/MSCprojectwork2024gutmicrobiomealopeciaareata/FAY_scripts
```
i did a little change to the command by assigning an output directory to save the fastqc files in a different directory as the samples so to avoid mixing them and getting confused.
when reading and finding relevant options regarding the command generation i searched through google and github repositories of several people including Dr. Tedder and also saw some youtube videos, from where i learnt several small changes that can be done to either activate the conda shell or assigning specific directories to save results into. below is another command that i generated with some minor changes in it just to see what it does but have actually used the command provided by Dr. Tedder directly;
```
#!/bin/sh 
#SBATCH --job-name=fastQC
#SBATCH --account=da494
#SBATCH --partition=compute
#
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00

### Initialize conda for the current shell
module load conda
conda init bash  # Make sure the correct shell is used, replace 'bash' with the appropriate shell i$
source ~/.bashrc # Source the bashrc to apply the conda initialization

### Load conda environment
conda activate env_0

### Navigate to where the data is
cd /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs

### Check if the file is a valid GZIP file
gzip -t "input_file"
```
both commands can be used to generate the fastqc results for the samples depending upon if the file formatting is done properly or not and if the conda environment for bash shell is properly activated or not. the prior one is for the files properly Gzipped and the later one is for the samples not properly Gzipped and needs to be checked and to be settled accordingly, the -t in gzip -t denotes test, literally meaning test if the samples are gzipped properly.
```
the fastqc result files are saved into directory path:
cd /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs/MSCprojectwork2024gutmicrobiomealopeciaareata/FAY_scripts
```
then the samples are gotten through NanoFilt for trimming and filtering to generate higher quality data using the below command
```
#!/bin/bash
#SBATCH --job-name=nanofilt
#SBATCH --account=da494
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00

### Load conda module and initialize
module load bioconda  # Adjust module name if necessary
conda init bash  # Initialize conda for bash shell

### Activate conda environment
source activate env_0  # Ensure correct environment activation

### Navigate to where the data is
cd /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs

### Execute the command
gunzip -c FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.fastq.gz |
NanoFilt -q 12 --headcrop 10 | gzip > FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.nanofilt_trimmed.fastq.gz

echo "NanoFilt filtering complete."
```

nanoplot script for nanofilt samples:
```

#!/bin/bash
#SBATCH --job-name=nanoplot
#SBATCH --account=da494
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00

### Load conda module and initialize
module load bioconda  # Adjust module name if necessary
conda init bash  # Initialize conda for bash shell

### Activate conda environment
source activate env_0  # Ensure correct environment activation

### Navigate to where the data is
cd /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs

### Execute the command
NanoPlot --fastq FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.nanofilt_trimmed.fastq.gz
 -o /storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs/MSCprojectwork2024gutmicrobiomealopeciaareata/nanofilt_nanoplot
```

centrifuge script:
```
#!/bin/sh
#SBATCH --job-name=centrifuge
#SBATCH --account=da494
#SBATCH --partition=compute
#
#SBATCH --cpus-per-task=20
#SBATCH --time=2:00:00

# Load the centrifuge module
module load centrifuge

# Define input and output directories
INPUT=/storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs
OUTPUT=/storage/pkdosanj/alopecia_nanopore_RAW/nonHuman_fastqs/MSCprojectwork2024gutmicrobiomealopeciaareata/centrifuge_files

# Run centrifuge command
centrifuge -p $SLURM_CPUS_PER_TASK \
-x /storage02/data/centrifuge-dbs/park-et-al-2020/hpvc \
-U ${INPUT}/FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.nanofilt_trimmed.fastq.gz \
-S ${OUTPUT}/FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.nanofilt_trimmed.fastq.gz_centrifugeOutputs.txt \
--report-file ${OUTPUT}/FAY24736_pass_barcode01_4fde6033_100a1825_COMBINED_nonHuman.nanofilt_trimmed.fastq.gz_centrifugeReport.txt
```
decontamination of centrifuged samples
```
!python "C:/Users/LENOVO/Downloads/centrifuge_env_decontam.py" "C:/Users/LENOVO/Downloads/centrifugefiles_todo_decon_barcode01-13" "C:/Users/LENOVO/Downloads/contaminantbarcode13" "C:/Users/LENOVO/Downloads/metadata_barcode01-13.txt" species
```
pie-chart generation of centrifuged decontaminated samples
```
!python "C:/Users/LENOVO/Downloads/cent_out_2_pie_chart.py" "C:/Users/LENOVO/Downloads/decon_files_todo_abundance_filter" pdf
```
abundance filtering of decontaminated samples
```
!python "C:/Users/LENOVO/Downloads/genus_level_read_count_abundance.py" "C:/Users/LENOVO/Downloads/decon_files_todo_abundance_filter" 0.01 species
```
PCA-plot generation
```
!python "C:/Users/LENOVO/Downloads/abundance_PCA_3D_variance.py" "C:/Users/LENOVO/Downloads/abundancefiltered_barcode01-13_todo_stats" "C:/Users/LENOVO/Downloads/pca_plot_output_barcode01-13" "C:\Users\LENOVO\Downloads\metadata_pca_barcode01-13.txt" 3D show_variance
```
estimations of alpha and beta diversities
```
!python "C:/Users/LENOVO/Downloads/alpha_beta_diversity.py" "C:\Users\LENOVO\Downloads\abundancefiltered_barcode01-13_todo_stats"
```
heatmap generation
```
import glob
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
 
directory = "C:/Users/LENOVO/Downloads/abundancefiltered_barcode01-13_todo_stats"
 
def abundance_list_prep(folder, top_hits):
    '''Prepares a list of unique OTUs collected from all abundance files.
    User can specify how many "top_hit" OTUs per sample.'''
    OTU_list = []
    for file_path in glob.glob(os.path.join(folder, "*.txt")):
        tmp_dict = {}
        print(f"Processing file: {file_path}")
        with open(file_path, "r") as file:
            for line in file:
                species, reads, abun = line.strip().split(',')
                tmp_dict[float(abun)] = species
        tmp_dict = dict(reversed(sorted(tmp_dict.items())))
        count = 0
        for item in tmp_dict.values():
            if count < int(top_hits):
                if item not in OTU_list:
                    OTU_list.append(item)
                    count += 1
                else:
                    count += 1
    print(f"OTU list: {OTU_list}")
    return OTU_list
 
# Main code block
OTU_dict = {}
sample_list = []
for taxa in abundance_list_prep(directory, 100):
    OTU_dict[taxa] = []
 
for file_path in glob.glob(os.path.join(directory, "*.txt")):
    sample_list.append(os.path.basename(file_path).split('_')[0])
    with open(file_path, "r") as file:
        tmp_dict = {}
        for line in file:
            species, reads, abun = line.strip().split(',')
            if species in OTU_dict:
                tmp_dict[species] = float(abun.strip())
    for taxa in OTU_dict:
        if taxa in tmp_dict:
            OTU_dict[taxa].append(tmp_dict[taxa])
        else:
            OTU_dict[taxa].append(0)
 
# Debug print statements
print(f"Sample list: {sample_list}")
for taxa, abundances in OTU_dict.items():
    print(f"{taxa}: {abundances}")
 
# Create a list of taxa
taxa = list(OTU_dict.keys())
 
# Convert OTU_dict to a numpy array
data = np.array([OTU_dict[taxon] for taxon in taxa])
 
# Check if the data is non-empty
if data.size == 0:
    raise ValueError("The input data is empty. Please check your input files.")
 
# Perform hierarchical clustering
col_linkage = hierarchy.linkage(data.T, method='ward')
 
# Plot heatmap with integrated dendrogram
plt.figure(figsize=(12, 8))
sns.clustermap(data, cmap='coolwarm', linewidths=0.5, linecolor='black', 
               xticklabels=sample_list, yticklabels=taxa, cbar_kws={'label': 'Abundance'}, 
               col_cluster=True, col_linkage=col_linkage)
 

 
# Show the plot
plt.show()
```
statistical analysis
-scipy installation
```
pip install scipy
```
installation of Shapiro-Wilk
```
from scipy.stats import shapiro
```
script for normality test
```
# Function to test normality using Shapiro-Wilk test
def test_normality(data):
 """
 Perform Shapiro-Wilk test for normality.
 Parameters:
 data (list or array-like): "C:Users/LENOVO/Downloads/diversityFAY24736/"
 
 Returns:
 float: Test statistic.
 float: p-value.
 """
 stat, p = shapiro(data)
 return stat, p

# Example data 
alpha_diversity = [0.373898274159285, 0.168729041954119, 0.135735044079991, 
    0.52797315192924, 0.232995277642313, 0.264033240197922, 
    0.436631761952671, 0.986728211044696, 0.313193232563141, 
    0.542805854346446, 0.208898638095607, 0]

beta_diversity = [0.961448598130841, 0.966717479674797, 0.183800623052959, 
    0.856201975850713, 0.487279843444227, 0.967266366816592, 
    0.743137254901961, 0.87588820464235, 0.901094752736882, 
    0.741362290227048, 0.682926829268293, 0.961448598130841, 
    0.084325671559125, 0.947773103148121, 0.615439599714081, 
    0.897501058873359, 0.091794731064764, 0.775553213909378, 
    0.582402395485431, 0.469471857716083, 0.780619644034278, 
    0.992280285035629, 0.084325671559125, 0.954528230390299, 
    0.657676789127203, 0.91044776119403, 0.019615384615385, 
    0.802924791086351, 0.625630987946842, 0.522498779892631, 
    0.807410849111395, 0.993298969072165, 0.183800623052959, 
    0.947773103148121, 0.954528230390299, 0.807177289769684, 
    0.352517985611511, 0.955273947074171, 0.661971830985916, 
    0.833024118738404, 0.866369710467706, 0.659735349716446, 
    0.751196172248804, 0.856201975850713, 0.615439599714081, 
    0.657676789127203, 0.807177289769684, 0.647058823529412, 
    0.660523560209424, 0.335931410756041, 0.113480995351381, 
    0.233373063170441, 0.35052754982415, 0.969590643274854, 
    0.487279843444227, 0.897501058873359, 0.91044776119403, 
    0.352517985611511, 0.647058823529412, 0.911882510013351, 
    0.421513944223108, 0.686274509803922, 0.748266296809986, 
    0.418269230769231, 0.869674185463659, 0.967266366816592, 
    0.091794731064764, 0.019615384615385, 0.955273947074171, 
    0.660523560209424, 0.911882510013351, 0.805669867398262, 
    0.627807704034963, 0.528765539173171, 0.807802310948404, 
    0.99341104916371, 0.743137254901961, 0.775553213909378, 
    0.802924791086351, 0.661971830985916, 0.335931410756041, 
    0.421513944223108, 0.805669867398262, 0.404553415061296, 
    0.499557913351017, 0.057484348321002, 0.942731277533039, 
    0.87588820464235, 0.582402395485431, 0.625630987946842, 
    0.833024118738404, 0.113480995351381, 0.686274509803922, 
    0.627807704034963, 0.404553415061296, 0.193577163247101, 
    0.415730337078652, 0.973986993496748, 0.901094752736882, 
    0.469471857716083, 0.522498779892631, 0.866369710467706, 
    0.233373063170441, 0.748266296809986, 0.528765539173171, 
    0.499557913351017, 0.193577163247101, 0.511518015357354, 
    0.979503350413875, 0.741362290227048, 0.780619644034278, 
    0.807410849111395, 0.659735349716446, 0.35052754982415, 
    0.418269230769231, 0.807802310948404, 0.057484348321002, 
    0.415730337078652, 0.511518015357354, 0.942286348501665, 
    0.682926829268293, 0.992280285035629, 0.993298969072165, 
    0.751196172248804, 0.969590643274854, 0.869674185463659, 
    0.99341104916371, 0.942731277533039, 0.973986993496748, 
    0.979503350413875, 0.942286348501665]

# Test normality for list1
statistic1, p_value1 = test_normality(list1)
print("List 1:")
print("Test statistic:", statistic1)
print("p-value:", p_value1)

if p_value1 > 0.05:
 print("Data is normally distributed (fail to reject H0)")
else:
 print("Data is not normally distributed (reject H0)")

# Test normality for list2
statistic2, p_value2 = test_normality(list2)
print("\nList 2:")
print("Test statistic:", statistic2)
print("p-value:", p_value2)

if p_value2 > 0.05:
 print("Data is normally distributed (fail to reject H0)")
else:
 print("Data is not normally distributed (reject H0)")
```
