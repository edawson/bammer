import sys
import multiprocessing
import argparse
import subprocess

fb = "freebayes"
samtools = "samtools"
bgzip = "bgzip"
tabix = "tabix"
vcf_sort = "vcfsort-parallel"
pylauncher = "launcher.py"

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", "--first", help="First bam file for comparison", required=True)
    parser.add_argument("-2", "--second", help="Second bam file for comparison", required=True)
    parser.add_argument("-v", "--vcf", help="VCF file of sites for comparison", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads to run for variant calling", required = False)

    return parser.parse_args()

def index_bam(bam):
    outname = basename(bam, "") + ".bai"
    subprocess.check_call(samtools + " index " + bam)
    return outname

def basename(fi, suffix):
    return fi.split("/")[-1].strip(suffix)

def sort_bam(bam):
    outname = basename(bam, ".bam")
    subprocess.call(samtools + " sort " + bam + " > " + outname)
    return outname

def index_vcf(vcf):
    subprocess.call(
    return None

def sort_vcf(vcf, outdir="."):
    subprocess.call(vcf_sort + " " + vcf " > " + outdir + "/" + basename(vcf, ".vcf") + ".sorted.vcf", shell = True)
    return basename(vcf, ".vcf") + ".sorted.vcf"

def call_variants_freebayes(bam, bed):
    return None

def check_sex(bam):
    return None

def calculate_concordance(first_bed, second_bed):
    return None

def check_sex_concordance(first_ratio, second_ratio):
    return None

def compare_gt(g1, g2):
    return g1 == g2

def gt_to_text(gt):
    if gt == "./." or gt == ".|.":
        return "missing"
    elif gt == "0/0" or gt == "0|0":
        return "homref"
    elif gt == "1/1" or gt == "1|1" or gt == "2/2" or gt == "2|2":
        return "homalt"
    else:
        return "het"

def parse_vcf_line_to_bed(line):
    tokens = line.strip().split("\t")
    chrom = tokens[0]
    pos = int(tokens[1])
    end = pos + len(tokens[3]) - 1

    ## Extract genotype
    fmt = tokens[7].split(";")
    samp = tokens[8].split(";")
    g_index = 0
    for i in range(0, len(fmt)):
        if fmt[i] == "GT":
            g_index = i
            break
    gt = samp[i]
    gt_text = gt_to_text(gt)
    return "\t".join([str(i) for i in [chrom, pos, end, ";".join([gt_text, gt])]]), gt, gt_text

def compare_positions(chr1, pos1, chr2, pos2):
    return chr1 == chr2 and pos1 == pos2

def compare_vcfs(vcf1, vcf2):
    cp = compare_positions

    same_count = 0
    diff_count = 0
    het_hom_count = 0
    hom_het_count = 0
    missing_first_count = 0
    missing_second_count = 0
    with open(vcf1, "r") as first and \
    open(vcf2, "r") as second:
        for line in  





if __name__ == "__main__":


