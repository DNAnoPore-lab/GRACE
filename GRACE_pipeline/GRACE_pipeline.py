import sys
import joblib
import shutil
import subprocess
import numpy as np
import pandas as pd


# all file will store in current working directory

########## pre-processing ##########
def require_executable(cmd: str) -> str:
    path = shutil.which(cmd)
    if path is None:
        print(f"\n>>> ERROR: {cmd} is not installed or not in PATH.\n")
        sys.exit(1)
    print(f">>> find {cmd} in {path} <<<")


def run_cmd(command: str, action, capture_output=False) -> str:
    try:
        print(f">>> start {action} <<<\n")
        result = subprocess.run(command, shell=True, check=True, capture_output=capture_output, text=capture_output)
        print(f">>> end {action} <<<\n")
        if capture_output:
            return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print("\n>>> something went wrong, please check error info <<<\n")
        sys.exit(1)


def filter() -> None:
    require_executable("samtools")
    command = f"samtools view aligned.bam -h -q 20 -F 0x4 -F 0x100 -F 0x800 | perl Fragmentomics/Scripts/samCigarToTlen.pl | awk '( $9 < 700 || $1 ~ /^@/ )' | samtools view -b -o filtered.bam"
    run_cmd(command, "samtools view")


def sort():
    command = "samtools sort filtered.bam -o filtered.sorted.bam"
    run_cmd(command, "samtools sort")


def index():
    command = "samtools index filtered.sorted.bam"
    run_cmd(command, "samtools index")


########## Fragmentomics ###########
#### short_fragment_ratio ####
def short_fragment_ratio() -> tuple[float, float, float]:
    require_executable("samtools")

    short_mononucleosome_count_command = "samtools view filtered.sorted.bam | awk 'length($10) >= 100 && length($10) <= 150 {print $0}' | wc -l"
    short_mononucleosome_count = run_cmd(short_mononucleosome_count_command, "samtools view",capture_output=True)

    all_mononucleosome_count_command = "samtools view filtered.sorted.bam | awk 'length($10) >= 100 && length($10) <= 220 {print $0}' | wc -l"
    all_mononucleosome_count = run_cmd(all_mononucleosome_count_command, "samtools view",capture_output=True)

    short_dinucleosome_count_command = "samtools view filtered.sorted.bam | awk 'length($10) >= 275 && length($10) <= 325 {print $0}' | wc -l"
    short_dinucleosome_count = run_cmd(short_dinucleosome_count_command, "samtools view",capture_output=True)

    all_dinucleosome_count_command = "samtools view filtered.sorted.bam | awk 'length($10) >= 275 && length($10) <= 400 {print $0}' | wc -l"
    all_dinucleosome_count = run_cmd(all_dinucleosome_count_command, "samtools view",capture_output=True)

    length_180_220bp_count_command = "samtools view filtered.sorted.bam | awk 'length($10) >= 180 && length($10) <= 220 {print $0}' | wc -l"
    length_180_220bp_count = run_cmd(length_180_220bp_count_command, "samtools view",capture_output=True)

    all_filtered_reads_count_command = "samtools flagstat filtered.sorted.bam | awk 'NR==1{print $1}'"
    all_filtered_reads_count = run_cmd(all_filtered_reads_count_command, "samtools flagstat",capture_output=True)

    counts = {
        'short_mono': short_mononucleosome_count,
        'all_mono'  : all_mononucleosome_count,
        'short_di'  : short_dinucleosome_count,
        'all_di'    : all_dinucleosome_count,
        '180-220bp' : length_180_220bp_count,
        'all_reads' : all_filtered_reads_count,
    }
    if not all(counts.values()):
        raise ValueError(f"Count cannot be empty; check BAM file. \nEmpty fields: "
                         f"{[k for k, v in counts.items() if not v]}\n"
                         "please check filtered.sorted.bam")
    elif all_mononucleosome_count == 0 or all_dinucleosome_count ==0 or all_filtered_reads_count ==0:
        raise ValueError("all_mononucleosome_count/all_dinucleosome_count/all_filtered_reads_count "
                         "can't be 0,please check filtered.sorted.bam")
    else:
        short_mononucleosome_ratio = float(short_mononucleosome_count) / float(all_mononucleosome_count)
        short_dinucleosome_ratio = float(short_dinucleosome_count) / float(all_dinucleosome_count)
        proportion_180_220bp = float(length_180_220bp_count) / float(all_filtered_reads_count)
        return short_mononucleosome_ratio,short_dinucleosome_ratio,proportion_180_220bp

#### end_motif ####
def create_stats_file():
    require_executable("samtools")
    require_executable("perl")
    command = ("samtools view filtered.sorted.bam | perl Fragmentomics/Scripts/General/stats_maker_0basedstart.pl > "
               "Fragmentomics/data/STATS/rapid.stats")
    run_cmd(command,"samtools view")

def create_motif_file():
    require_executable("bedtools")
    command_1 = "awk '{print $1\":\"$2\"-\"$3}' Fragmentomics/data/STATS/rapid.stats > Fragmentomics/data/A-rapid.txt"
    command_2 = ("bedtools getfasta -s -name -tab -fi "
                 "Homo_sapiens_assembly38.fasta "
                 "-bed Fragmentomics/data/STATS/rapid.stats | "
                 "awk '{print substr($1,1,length($1)-3), substr($1,length($1)-2,length($1)),substr($2,1,4)}' "
                 "> Fragmentomics/data/B-rapid.txt")
    command_3 = "awk 'NR==FNR {a[NR]=$1; next} {print $1,a[FNR]$2,$3}' Fragmentomics/data/A-rapid.txt Fragmentomics/data/B-rapid.txt > Fragmentomics/data/STATS/MOTIF/rapid.motif"
    run_cmd(command_1,"awk")
    run_cmd(command_2,"bedtools")
    run_cmd(command_3,"awk")

def create_motif_R():
    require_executable("Rscript")
    command = ("Rscript "
               "Fragmentomics/Scripts/Motifs/count_motif.R "
               "Fragmentomics/Utility/chr_list_hg38.txt "
               "Fragmentomics/data/STATS "
               "Fragmentomics/data/STATS/MOTIF/ "
               "Fragmentomics/data/STATS/MOTIF/MOTIF_COUNTS/ "
               "rapid")
    run_cmd(command,"Rscript")

def extract_motif():
    require_executable("Rscript")
    command = ("Rscript "
               "Fragmentomics/Scripts/Motifs/motif_extraction_rapid.R "
               "Fragmentomics/Utility/sample_info_motif_rapid.txt "
               "Fragmentomics/Utility/endmotifs.txt "
               "Fragmentomics/motif_rapid.csv")
    run_cmd(command,"Rscript")

def end_motif():
    create_stats_file()
    create_motif_file()
    create_motif_R()
    extract_motif()


############## Epigenetic ###################
def bam2bed():
    require_executable("modkit")
    command = f"modkit pileup filtered.sorted.bam traditional.pileup.bed --ref Homo_sapiens_assembly38.fasta --preset traditional"
    run_cmd(command, "modkit pileup")

def bed2bedgraph():
    command = "cat traditional.pileup.bed | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$11}' > traditional.pileup.bedgraph"
    run_cmd(command, "bed2bedgraph")

#### genome-level_methylation ####
def gene_level_analysis():
    # intersect and tag
    require_executable("bedtools")
    require_executable("Rscript")
    command = f"bedtools intersect -a Epigenetics/gene_analysis/gencode.v45.annotation.genetype.sorted.bed -b traditional.pileup.bedgraph -wa -wb > rapid.gene_summary.txt"
    run_cmd(command, "bed2bedgraph")

    # use Rscript to extract genome
    command = f"Rscript Epigenetics/gene_analysis/genome_extraction_rapid.R Epigenetics/gene_analysis/genome_rapid.csv"
    run_cmd(command, "Rscript")

#### bin-level_methylation ####
def bins_10Mb():
    # intersect and tag
    require_executable("bedtools")
    require_executable("Rscript")
    command = "bedtools intersect -a Epigenetics/10mbin_analysis/10mb-bin_pos.bed -b traditional.pileup.bedgraph -wa -wb > rapid.10mbin_summary.txt"
    run_cmd(command, "bedtools")

    command = "Rscript Epigenetics/10mbin_analysis/bin_extraction_rapid.R Epigenetics/10mbin_analysis/bin_rapid.csv"
    run_cmd(command, "Rscript")



############# Machine Learning #############
def run_diagnosis_model(features: np.ndarray):
    model_path = "GRACEmodel.pkl"
    model = joblib.load(model_path)
    y_scores = model.predict_proba(features)[:, 1]
    return y_scores


if __name__ == "__main__":
    filter()
    sort()
    index()


    ########## Fragmentomics #########

    # short_fragment_ratio
    short_mononucleosome_ratio,short_dinucleosome_ratio,proportion_180_220bp = short_fragment_ratio()
    short_fragment_ratio = np.array([short_mononucleosome_ratio,short_dinucleosome_ratio,proportion_180_220bp]).reshape(1,3)
    # end_motif
    end_motif()
    end_motif_feature = pd.read_csv("Fragmentomics/motif_rapid.csv",header=None).values[1,:]
    end_motif_feature = np.array(end_motif_feature).reshape(1,-1)


    ########## Epigenetics ############
    bam2bed()
    bed2bedgraph()
    gene_level_analysis()
    bins_10Mb()
    gene_level_feature = pd.read_csv("Epigenetics/gene_analysis/genome_rapid.csv",header=None).values
    bin_level_feature = pd.read_csv("Epigenetics/10mbin_analysis/bin_rapid.csv",header=None).values[1,:]
    bin_level_feature = np.array(bin_level_feature).reshape(1,-1)


    input_feature = np.concat([short_fragment_ratio,end_motif_feature,gene_level_feature,bin_level_feature],axis=1)
    print(f"short_fragment_ratio: {short_fragment_ratio}")
    print(f"end_motif_feature: {end_motif_feature}")
    print(f"gene_level_feature: {gene_level_feature}")
    print(f"bin_level_feature: {bin_level_feature}")
    cancer_score = run_diagnosis_model(input_feature)
    print(f"G-score: {cancer_score}")