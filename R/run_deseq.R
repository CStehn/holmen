#' Function for DESeq2 differential expression analysis workflow
#'
#' This function will output results and plots for each comparison within each variable presented in the samples table
#'
#' @param counts a tibble with each row representing counts from a gene of interest, the 1st column being the gene names, and each column representing a sample
#' @param samples a tibble with each row representing samples in the counts tibble with variables of interest included
#' @param design a formula specifying the variables of interest to be included in the analysis in the form ~ v1 + v2 + ... + vn. Defaults to ~ trt
#' @param human a logical value specifying whether to use a human annotation table. Defaults to TRUE. FALSE defaults to mouse
#' @param filt_threshold value specifying the minimum amount of counts needed by a sample within a gene to pass to DESeq
#' @param top either an integer specifying the number of top ranked genes (by p-value) to create heatmap, or a character vector of specific genes of interest
#' @param blind a logical value specifying whether the rlog or vst transformations are blind to experimental design. Defaults to TRUE
#' @param font_size integer specifying the font size for rows of the heatmap
#' @param volcano_cols a length two character vector specifying two colors (either base R colors or hexadecimal) to be used for non-significant and significant genes in the volcano, respectively
#'
#' @return Returns PCA for each variable along with results table, volcano, heatmap, and sample correlation heatmap using the top genes between each comparison
#'
#' @author Chris Stehn
#'
#' @examples
#' \dontrun{
#' run_deseq(counts, samples, design = ~ gender + height, human = TRUE, filt_threshold = 10, top = 200, blind = TRUE, font_size = 3)
#' }
#' @export

run_deseq <- function(counts, samples, design = ~ trt, human = TRUE, filt_threshold = 1, top = 40, blind = TRUE, font_size = 4, volcano_cols = c('black','red')) {
  ## Check if inputs are tibble
  if(!tibble::is_tibble(counts)) {stop('counts object must be a tibble')}
  if(!tibble::is_tibble(samples)) {stop('samples object must be a tibble')}

  ## Get variable names from formula
  vars <- formula.tools::get.vars(design)

  ## Pre-process data with filtering defined by filt_threshold
  message('Filtering low-signal genes')
  counts <- suppressMessages(hciR::filter_counts(counts, n = filt_threshold)) %>% hciR::as_matrix()

  ## Round counts if not already rounded
  if(any(round(counts[1:20,]) != counts[1:20,])) {
    message('Rounding counts')
    counts <- round(counts, 0)
  }
  ## Create DESeqDataSet and run differential expression analysis
  message('Running DESeq2')
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(counts, samples, design = design))

  dds <- suppressMessages(DESeq2::DESeq(dds))
  ## Create annotation table using read_biomart function of hciR
  message('Downloading annotation table')
  if(human) {bio <- suppressMessages(hciR::read_biomart('human'))
  } else {
    bio <- suppressMessages(hciR::read_biomart('mouse'))
  }


  message('Transforming counts')
  if(ncol(counts) <= 50) {
    rld <- DESeq2::rlog(dds, blind = blind)
  } else {
    rld <- DESeq2::vst(dds, blind = blind)
  }

  for (i in 1:length(vars)) {
    ## Allows for comparison using multiple variables
    htmlwidgets::saveWidget(hciR::plot_pca(rld, vars[[i]]), file = paste0('pca_',vars[[i]],'.html'))
    # Extract results from dds using biomart table
    message(paste('Extracting results and plotting for variable', vars[[i]]))
    res <- suppressMessages(hciR::results_all(dds, bio, trt = vars[[i]]))

    ## Filter NA p-values and create Result column denoting statistical significance

    if (!is.data.frame(res)) {
      ## For more than one comparison
      for (j in 1:choose(length(levels(dds[[vars[[i]]]])), 2)) {
        names <- gsub('__','_', gsub('[ .]', '_', names(res)))
        res1 <- res[[j]]
        res1 <- dplyr::filter(res, !is.na(padj)) %>% dplyr::mutate(Result = ifelse(padj < 0.05, 'Significant','Not Significant'))

        volcano <- plotly::ggplotly(ggplot2::ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), group = gene_name)) + geom_point(aes(color = Result), size = 0.5) + xlab('log2 Fold Change') + ylab('-log10 p-value') + scale_color_manual(values = volcano_cols, name = 'Result'))
        htmlwidgets::saveWidget(volcano, file = paste0('volcano_', vars[[i]],'_',names[[j]], '.html'))
        ## Arrange results by adjusted p-value, then write to csv
        res1 <- dplyr::arrange(res1, padj)
        write.csv(res1, row.names = FALSE, file = paste0('DESeq_results_', vars[[i]],'_',names[[j]], '.csv'))

        ## Store all filtered results in another list of data frames
        if (is.character(top)) {
          resfilt <- filter(res1, gene_name %in% top)
        } else {
          resfilt <- res1[1:top,]
        }

        ## Convert each filtered result table into normalized count tables for plotting
        x <- hciR::top_counts(resfilt, rld, filter = FALSE) %>% hciR::as_matrix()
        x <- x - rowMeans(x)
        n <- max(abs(range(x)))
        brks <- seq(-n, n, length = 101)
        ## Create annotation table from colData, convert to data frame and create heatmap
        annot <- attr(dds, 'colData')[, vars, drop = FALSE] %>% as.data.frame()
        pheatmap::pheatmap(x, annotation_col = annot, breaks = brks, border_color = NA, fontsize_row = font_size, color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdBu')))(100), filename = paste0('heatmap_', vars[[i]],'_',names[[j]], '.png'))
        ## Create correlation matrix between samples and create heatmap

        cor <- cor(x, method = 'pearson')
        pheatmap::pheatmap(cor, breaks = seq(-1, 1, length = 101), annotation_col = annot, annotation_row = annot, fontsize_row = font_size, border_color = NA, show_colnames = FALSE, annotation_names_col = FALSE, color = colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(100), filename = paste0('corheatmap_', vars[[i]],'_',names[[j]], '.png'))
      }
    } else {
      ## For one comparison only
      res <- filter(res, !is.na(padj)) %>% mutate(Result = ifelse(padj < 0.05, 'Significant','Not Significant'))

      volcano <- plotly::ggplotly(ggplot2::ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), group = gene_name)) + geom_point(aes(color = Result), size = 0.5) + xlab('log2 Fold Change') + ylab('-log10 p-value') + scale_color_manual(values = volcano_cols, name = 'Result'))
      htmlwidgets::saveWidget(volcano, file = paste0('volcano_', vars[[i]],'.html'))
      ## Arrange results by adjusted p-value, then write to csv
      res <- arrange(res, padj)
      write.csv(res, row.names = FALSE, file = paste0('DESeq_results_', vars[[i]],'.csv'))

      ## Store all filtered results in another list of data frames
      if (is.character(top)) {
        res1 <- filter(res, gene_name %in% top)
      } else {
        res1 <- res[1:top,]
      }

      ## Convert each filtered result table into normalized count tables for plotting
      x <- hciR::top_counts(res1, rld, filter = FALSE) %>% hciR::as_matrix()
      x <- x - rowMeans(x)
      n <- max(abs(range(x)))
      ## Create annotation table from colData, convert to data frame and create heatmap
      annot <- attr(dds, 'colData')[, vars, drop = FALSE] %>% as.data.frame()
      pheatmap::pheatmap(x, annotation_col = annot, breaks = seq(-n, n, length = 101), border_color = NA, fontsize_row = font_size, color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdBu')))(100), filename = paste0('heatmap_', vars[[i]],'.png'))
      ## Create correlation matrix between samples and create heatmap

      cor <- cor(x, method = 'pearson')
      pheatmap::pheatmap(cor, breaks = seq(-1, 1, length = 101), annotation_col = annot, annotation_row = annot, fontsize_row = font_size, border_color = NA, show_colnames = FALSE, annotation_names_col = FALSE, color = colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(100), filename = paste0('corheatmap_', vars[[i]],'.png'))
    }
  }
  return(message('Finished'))
}
