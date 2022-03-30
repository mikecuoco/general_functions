#!/usr/bin/env Rscript
# Author: Mike Cuoco
# Date: Mar 29, 2022
# Usage: Functions to import RNA-seq data into R

require(tidyverse)
require(rtracklayer)
require(tximeta)
require(GenomicFeatures)
require(DESeq2)

#' Create a DESeq2 or RangedSummarizedExperiment object from
#' RSEM *.genes.results files, sample metadata, and the GENCODE annotation
#' 
#' @param coldata, a data.frame with the following columns: files, names, ...
#' @param se.obj, a logical to indicate whether to return data in RangedSummarizedExperiment format
#' @param deseq.obj, a logical to indicate whether to return data in DESeq2 format
read_rsem_genes <- function(coldata, se.obj = TRUE, deseq.obj = FALSE){
    
    # catch argument error
    if (se.obj & deseq.obj){
        stop("Must specify either se.obj or deseq.obj, but not both")
    }

    # read gtf file
    gtf = "/raidixshare_logg01/mcuoco/references/Homo_sapiens/GENCODE/gencode.v39.annotation.gtf.gz"
    id2name = rtracklayer::readGFF(gtf) %>%
        tibble::as_tibble() %>% 
        dplyr::filter(type == "gene", gene_type %in% c("protein_coding","lncRNA")) %>%
        dplyr::select(gene_id, gene_name, gene_type) %>%
        dplyr::distinct() 
    gencode = GenomicFeatures::makeTxDbFromGFF(gtf) %>%
        GenomicFeatures::genes() 
    gencode = gencode[id2name$gene_id] 
    gencode$gene_name = id2name$gene_name
    gencode$gene_type = id2name$gene_type
    
    # read rsem counts with tximeta, saving to RangedSummarizedExperiment 
    se = tximeta::tximeta(coldata, type = "rsem", txIn = F, txOut = F)
    se = se[rownames(assays(se)[["counts"]]) %in% gencode$gene_id,]
    rowRanges(se) <- gencode

    # Remove zero-length genes and return as DESeq2 object
    if (deseq.obj) {
        zero_length_genes = se@rowRanges$gene_id[apply(assays(se)[['length']], 1, min) == 0]
        se = se[!se@rowRanges$gene_id %in% zero_length_genes,]
        dds = DESeq2::DESeqDataSet(se, design = ~ 1)
        return(dds)
    }

    # return as RangedSummarizedExperiment
    if (se.obj) {
        return(se)
    }
   
}

#' Read multiqc_stats.txt file and return a cleaned up data frame, removing rows with R1 and R2
read_qc_stats <- function(path){
    readr::read_delim(path) %>%
        dplyr::rename_with(.fn = function(x) stringr::str_remove(x, ".*-generalstats-")) %>%
        dplyr::filter(!grepl("_[0-9]$",Sample)) %>%
        dplyr::select(where(function(x) !all(is.na(x))))
}

#' Collect paths to RSEM *.genes.results files and return a data.frame with the following columns: files, names
prepare_rsem <- function(path){
    data.frame(files = list.files(path, pattern = "genes.results", full.names = T)) %>%
        dplyr::mutate(names = stringr::str_remove(basename(files), ".genes.results")) 
}