# ALIGN GUIDES TO REFERENCE GENOME

# Additional requirements: the bowtie and samtools command line tools 
# (see https://github.com/cancerdatasci/ceres installation instructions)

# This script is adapted from CERES scripts:
# Alignment code is adapted from map_guide_to_locus() / guideAlignments() functions
# Mapping guide to locus is adapted from load_ccds_genes() + get_gene_annotations() functions

# Genome: hg38

# Inputs:
#   - Guide sequences from DepMap log-fold change file - sgrna_sequences.csv
#   - Bowtie indexes for hg38
#   - CCDS gene annotations for GRCh38.p12: release 22, 06/14/2018

# Output: File with guide alignments and corresponding gene mappings - guide_alignments_hg19.csv

# Set location of bowtie indices and CCDS gene annotations
third_party_dir <- "/Users/barbaradekegel/Google Drive/3rd_party_data"
local_dir <- "../../local_data/processed/depmap20Q2"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.9")

BiocManager::install("Biostrings")
BiocManager::install("Rsamtools")
BiocManager::install("GenomeInfoDb")
BiocManager::install("BSgenome")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("GenomicRanges")

library(tidyverse)

# Inputs
bowtie_indexes <- file.path(third_party_dir, "bowtie_indexes", "hg38")
file_sgrna_sequences <-  file.path(local_dir, "sgrna_sequences.csv")
file_gene_annot <- file.path(third_party_dir, "ccds_gene_annotations", "CCDS.20180614.txt")

# Output
file_alignments <- file.path(local_dir, "guide_alignments_hg38.csv")

# Temp files
dir.create('temp_data', showWarnings = FALSE)
sam_file <- file.path("temp_data", "guides.sam", fsep=.Platform$file.sep)
bam_file <- file.path("temp_data", "guides.bam", fsep=.Platform$file.sep)
fasta_file <- file.path("temp_data", "guides.fa", fsep=.Platform$file.sep)


# Read in guide sequences
guide_sequences <- read_csv(file_sgrna_sequences)
guide_sequences <- pull(guide_sequences, 'sgrna')

# Write guide sequences to FASTA file, set names of the guides to be the sequences themselves
guide_sequences %>% 
  set_names(.,.) %>%
  Biostrings::DNAStringSet() %>%  # DNAStringSet knows the semantics of DNA base pairs
  Biostrings::writeXStringSet(fasta_file)

# Align guides with bowtie, to the specified reference genome - output is a SAM file
# -v: alignments may have no more than v mismatches, where v is a number between 0-3.
# -f: indicates that the query input files are FASTA files
# -S: print alignments in SAM format
Sys.setenv(BOWTIE_INDEXES = bowtie_indexes)
bowtie_cmd <- paste('bowtie', "-t -a -v 2 -f -S", "hg38", fasta_file, sam_file)
print('Doing alignment with bowtie...')
print(system(bowtie_cmd)) # 0 = success

# Convert SAM to BAM using samtools
samtools_cmd <- paste('samtools', "view -bS -o", bam_file, sam_file)
print('Converting SAM to BAM...')
system(samtools_cmd)

# Read in BAM files
# tag: options, MD:Z is mismatches
params <- Rsamtools::ScanBamParam(tag=c("MD"), 
                                  what=c("qname", "flag", "rname", "strand", "pos", "seq"))
bam <- Rsamtools::scanBam(bam_file, param=params)

# Sex chromosomes should be filtered out
chromosomes <- paste0("chr", c(as.character(1:22)))

print('Converting BAM file to data frame...')
# Bam file to data frame, count number of alignments
# Assumes guides are all the same length - indicated by guide_length variable
guide_length <- 20
bamdf <- as_tibble(bam[[1]][1:5]) %>% 
  rename(start_pos = pos) %>%
  mutate(end_pos = start_pos + guide_length,
         ref_seq = as.character(bam[[1]]$seq),
         aligned = !is.na(rname), 
         alignment = bam[[1]]$tag$MD) %>%
  filter(aligned & rname %in% chromosomes) %>%
  group_by(qname) %>% mutate(num_alns = sum(aligned)) %>% ungroup()

# Calculate PAM start and end positions
# Cas9 cuts 3 bps upstream from the PAM site
bamdf.pam <- bamdf %>% mutate(cut_pos = ifelse(strand=="+", end_pos-4, start_pos+2),
                              PAM_start = ifelse(strand=="+", end_pos, start_pos-3),
                              PAM_end = ifelse(strand=="+", end_pos+2, start_pos-1))

# Find PAM sequences
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
bamdf.pam$PAM <- BSgenome::getSeq(genome,
                                  names=bamdf.pam$rname,
                                  start=bamdf.pam$PAM_start,
                                  end=bamdf.pam$PAM_end,
                                  strand=bamdf.pam$strand,
                                  as.character=TRUE)

print('Filter on PAM sequences...')
# Check that PAM sequences match NGG | GGN
pam_pattern <- "[ACGTN]GG|GG[ACGTN]"
# Filter out guides that don't have correct PAM and calculate the number of alignments
bamdf.pam <- bamdf.pam %>% filter(str_detect(PAM, pam_pattern)) %>% 
  group_by(qname) %>% mutate(num_alns = sum(aligned)) %>% ungroup()


# Read in + clean CCDS data, gene annoations are used to overlap guide sequence 
# locations with gene locations
ccds <- read_tsv(file_gene_annot, col_types=cols("#chromosome" = col_character(), 
                                                 "cds_from" = col_integer(), 
                                                 "cds_to" = col_integer()))
ccds <- ccds %>%
  dplyr::rename(chromosome = "#chromosome") %>%
  mutate(chromosome = str_c("chr", chromosome)) %>%
  filter(ccds_status %in% c("Public", "Reviewed, update pending", "Under review, update")) %>%
  filter(chromosome %in% chromosomes, !is.na(cds_from), !is.na(cds_to))

# Transform gene annotations to format expected by makeGRangesFromDataFrame
# One row per exon location, gene start and end spans the full exon interval.
ccds_exon <- ccds %>%
  mutate(cds_interval = str_replace_all(cds_locations, "[\\[\\]]", "") %>% str_split("\\s*,\\s*")) %>%
  tidyr::unnest(cds_interval) %>%
  group_by(gene, gene_id, cds_locations) %>% 
  mutate(exon_code = ifelse(cds_strand=="+", 1:dplyr::n(), dplyr::n():1)) %>% 
  ungroup() %>%
  dplyr::mutate(cds_start = str_extract(cds_interval, "^[0-9]+") %>% as.integer,
                cds_end = str_extract(cds_interval, "[0-9]+$") %>% as.integer) %>%
  dplyr::select(gene, gene_id, chromosome, start=cds_start, end=cds_end, strand=cds_strand,
                gene_start=cds_from, gene_end=cds_to, exon_code)

genome_info <- GenomeInfoDb::Seqinfo(genome="hg38")[chromosomes]
gene_annot_granges <- ccds_exon %>% 
  GenomicRanges::makeGRangesFromDataFrame(seqinfo=genome_info, keep.extra.columns=T)

# Transform aligned_sequences to format expected by makeGRangesFromDataFrame
# Annotation is based on whether the cut position of an sgRNA overlaps with an exon of the gene
guide_aln_granges <- bamdf.pam %>%
  select(Guide=qname, Chr=rname, Start=cut_pos, strand, alignment) %>%
  mutate(End = Start) %>%
  GenomicRanges::makeGRangesFromDataFrame(seqinfo=genome_info, keep.extra.columns=T)



print('Use CCDS to annotate regions with genes...')
# Find overlapping genomics regions
# Annotations is based on whether the cut position of an sgRNA overlaps with an exon of the gene
hits <- GenomicRanges::findOverlaps(guide_aln_granges, gene_annot_granges, ignore.strand=T) %>% as_tibble

gene_df <- hits %>%
  transmute(sgrna = guide_aln_granges$Guide[queryHits],
            chr = GenomicRanges::seqnames(guide_aln_granges)[queryHits] %>% as.character(),
            cut_pos = GenomicRanges::start(guide_aln_granges)[queryHits] %>% as.integer(),
            strand = GenomicRanges::strand(guide_aln_granges)[queryHits] %>% as.character(),
            alignment = guide_aln_granges$alignment[queryHits],
            gene = gene_annot_granges$gene[subjectHits],
            geneID = gene_annot_granges$gene_id[subjectHits],
            CDS_strand = GenomicRanges::strand(gene_annot_granges)[subjectHits] %>% as.character(),
            CDS_start = GenomicRanges::start(gene_annot_granges)[subjectHits] %>% as.integer(),
            CDS_end = GenomicRanges::end(gene_annot_granges)[subjectHits] %>% as.integer()) %>% 
  distinct()

# Save out list of guides mapped to genes, including single and double mismatch guides
write_csv(gene_df, file_alignments)
