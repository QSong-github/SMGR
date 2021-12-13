#' Convert a peak matrix to a gene activity matrix
#'
#' This function will take in a peak matrix and an annotation file (gtf) and collapse the peak
#' matrix to a gene activity matrix. It makes the simplifying assumption that all counts in the gene
#' body plus X kb up and or downstream should be attributed to that gene.
#'
#' @param peak.matrix Matrix of peak counts
#' @param annotation.file Path to GTF annotation file
#' @param seq.levels Which seqlevels to keep (corresponds to chromosomes usually)
#' @param include.body Include the gene body?
#' @param upstream Number of bases upstream to consider
#' @param downstream Number of bases downstream to consider
#' @param keep.sparse Leave the matrix as a sparse matrix. Setting this option to
#' TRUE will take much longer but will use less memory. This can be useful if
#' you have a very large matrix that cannot fit into memory when converted to
#' a dense form.
#' @importFrom future nbrOfWorkers
#' @export feature mapped annotations and activity matrix

GeneActivityMatrix <- function(
                               peak.matrix,
                               annotation.file,
                               seq.levels = c(1:22, "X", "Y"),
                               include.body = TRUE,
                                     upstream = 2000,
                               downstream = 0,
                               keep.sparse = FALSE,
                               verbose = FALSE
                               ) {
    
    #' convert peak matrix to GRanges object
    peak.df <- rownames(x = peak.matrix)

    peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))

    peak.df <- as.data.frame(x = peak.df)
    colnames(x = peak.df) <- c("chromosome", 'start', 'end')
    peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)

    #' if any peaks start at 0, change to 1
    #' otherwise GenomicRanges::distanceToNearest will not work
    BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1

    #' get annotation file, select genes
    gtf <- rtracklayer::import(con = annotation.file)
    gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')

    #' change seqlevelsStyle if not the same
    if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
        GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
    }

    gtf.genes <- gtf[gtf$type == 'gene']

    #' Extend definition up/downstream
    gtf.body_prom <- Seurat:::Extend(x = gtf.genes, upstream = upstream, downstream = downstream)

    gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
    keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
    peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
    gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]

    #' Some GTF rows will not have gene_name attribute
    #' Replace it by gene_id attribute
    gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]

    peak.ids$gene.name <- gene.ids$gene_name
    peak.ids <- as.data.frame(x = peak.ids)
    peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
    annotations <- peak.ids[, c('peak', 'gene.name')]
    colnames(annotations) <- c('feature', 'new_feature')

    #' collapse into expression matrix
    if (!keep.sparse) {
        peak.matrix <- as(object = peak.matrix, Class = 'matrix')
    }
    all.features <- unique(annotations$new_feature)

    library(future); library(pbapply)
    if (nbrOfWorkers() > 1) {
        mysapply <- future_sapply
    } else {
        mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
    }

    newmat.list <- mysapply(1:length(all.features), FUN = function(x){
        features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
        submat <- peak.matrix[features.use, ]
        if (length(features.use) > 1) {
            submat <- Matrix::colSums(submat)
        }
        if (keep.sparse) {
            return(as(object = as.matrix(submat), Class = 'dgCMatrix'))
        } else {
            return(as.matrix(submat))
        }
    }, simplify = FALSE)

    newmat = do.call(cbind,newmat.list)
    newmat <- t(newmat)
    rownames(newmat) <- all.features

    out <- list(peak_id=peak.ids,feature_map=annotations,activity_matrix=as(object = newmat, Class = 'dgCMatrix'))

    return (out)
}


