library(R.utils)
library(glue)

genome_version <- snakemake@wildcards[['genome_version']]

input_tumor_bam <- snakemake@input[['Tum']]
input_tumor_bam <- getAbsolutePath(input_tumor_bam)

input_normal_bam <- snakemake@input[['NC']]
input_normal_bam <- getAbsolutePath(input_normal_bam)

output_rdata <- snakemake@output[['rdata']]
output_rdata <- getAbsolutePath(output_rdata)

threads = as.numeric(snakemake@threads)

## ASCAT config
alleles.prefix <- snakemake@params[['allelesprefix']]
loci.prefix <- snakemake@params[['lociprefix']]
GCcontentfile <- snakemake@params[['GCcontentfile']]
rp_file <- snakemake@params[['replictimingfile']]


# gender <- snakemake@params[['gender']]
gender = 'XX'
work_dir <- snakemake@params[['wd']]
sample_index <- snakemake@params[['sample_name']]

if(dir.exists(work_dir)){
} else {
  dir.creat(work_dir)
}
setwd(work_dir)

if(genome_version == 'T2T'){
  rp_file <- NULL
  f_gv = 'hg38'
} else if (genome_version == 'b37'){
  genome_version = 'hg19'
  f_gv = 'hg19'
} else {
  f_gv = genome_version
}

library(ASCAT)
ascat.prepareHTS(
  tumourseqfile = input_tumor_bam,
  normalseqfile = input_normal_bam,
  tumourname = "Tumor",
  normalname = "Germline",
  loci.prefix = loci.prefix,
  alleles.prefix = alleles.prefix,
    allelecounter_exe = 'alleleCounter',
  gender = gender,
  genomeVersion = f_gv,
  nthreads = 8,
  chrom_names = c(1:22),
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt")
  

ascat.bc = ascat.loadData(
  Tumor_LogR_file = "Tumor_LogR.txt", 
  Tumor_BAF_file = "Tumor_BAF.txt", 
  Germline_LogR_file = "Germline_LogR.txt", 
  Germline_BAF_file = "Germline_BAF.txt", 
  gender = gender, genomeVersion = f_gv, isTargetedSeq=F)


ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")

ascat.bc = ascat.correctLogR(ascat.bc, 
GCcontentfile = GCcontentfile, 
replictimingfile = rp_file)

ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = output_rdata)
