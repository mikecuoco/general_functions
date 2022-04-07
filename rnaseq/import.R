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
read_rsem_genes <- function(coldata, deseq.obj = TRUE, filter.genes.by = NULL){
    
    # read gtf file, slow! 
    # TODO make this faster
    gtf = "/raidixshare_logg01/mcuoco/references/Homo_sapiens/GENCODE/gencode.v39.basic.annotation.gtf.gz"
    id2name = rtracklayer::readGFF(gtf) %>%
        tibble::as_tibble() %>% 
        dplyr::filter(type == "gene") %>%
        dplyr::select(gene_id, gene_name, gene_type) %>%
        dplyr::distinct() 
    if (!is.null(filter.genes.by)){
        message(paste0("Filtering genes by: ", filter.genes.by))
        id2name = id2name %>% 
            dplyr::filter(gene_type %in% filter.genes.by)
    }
    gencode = GenomicFeatures::makeTxDbFromGFF(gtf) %>%
        GenomicFeatures::genes() 
    gencode = gencode[id2name$gene_id] # sort and filter
    gencode$gene_name = id2name$gene_name
    gencode$gene_type = id2name$gene_type
    
    # read rsem counts with tximeta, saving to RangedSummarizedExperiment 
    coldata[sapply(coldata, is.character)] = lapply(coldata[sapply(coldata, is.character)], as.factor)
    se = tximeta::tximeta(coldata, type = "rsem", txIn = F, txOut = F, ignoreTxVersion = TRUE, lengthCol = NA)
    se = se[rownames(assays(se)[["counts"]]) %in% gencode$gene_id,] # filter
    se = se[match(rownames(assays(se)[["counts"]]), gencode$gene_id),] # sort
    rowRanges(se) <- gencode

    # compute library size
    se$lib.size = colSums(assays(se)[["counts"]])

    if (deseq.obj) { # Remove zero-length genes and return as DESeq2 object
        zero_length_genes = se@rowRanges$gene_id[apply(assays(se)[['length']], 1, min) == 0]
        se = se[!se@rowRanges$gene_id %in% zero_length_genes,]
        dds = DESeq2::DESeqDataSet(se, design = ~ 1)
        return(dds)
    } else { # return as RangedSummarizedExperiment
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
