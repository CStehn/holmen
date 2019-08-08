#' Function for converting between Ensembl IDs and gene symbols (canonical names)
#'
#' This function will output the converted names of each gene as a data frame variable, either replacing the names in a given data frame or simply returning a column.
#'
#' @param IDs either a data frame of counts in which the first column is a vector of gene symbols or ensembl IDs or, simply, a vector
#' @param human a logical value specifying if the IDs are human. Defaults to TRUE. FALSE specifies mouse
#' @param mart a biomart data frame returned from hciR::read_biomart, if available. If set to NULL, the function will read its own biomart
#'
#' @return Returns either the original counts table with converted gene IDs, or a data frame of one column with the converted IDs
#'
#' @author Chris Stehn
#'
#' @examples
#' \dontrun{
#' ensembl_to_symbol(counts, mart = mmu, human = FALSE)
#' }
#' @export

ensembl_to_symbol <- function(IDs, mart = NULL, human = TRUE) {
  ## Read biomart
  if (missing(IDs)) {stop('No IDs have been supplied')}
  
  if (is.data.frame(IDs)) {
    if (!grepl('ENS', IDs[1, 1])) {
      stop('Ensembl IDs have not been supplied. Consider using symbol_to_ensembl instead')
    }
  } else {
    if (!grepl('ENS', IDs[1])) {
      stop('Ensembl IDs have not been supplied. Consider using symbol_to_ensembl instead')
    }
  }
  
  if (!is.null(mart)) {
    names <- dplyr::select(mart, id, gene_name) %>% dplyr::filter(!duplicated(gene_name))
  } else {
    if (human) {
      mart <- data(hsa, package = 'holmen')
    } else {
      mart <- data(mmu, package = 'holmen')
    }
    names <- dplyr::select(mart, id, gene_name) %>% dplyr::filter(!duplicated(gene_name))
  }

  ## Convert ensembl IDs to gene symbols using mart
  if (is.data.frame(IDs)) {
    temp <- as.data.frame(IDs[[1]])
    ## Set colname of temp so match can be made consistently
    colnames(temp)[1] <- 'ensembl'
    ## Perform left_join of gene symbols to ensembl ids
    symbol_match <- suppressMessages(dplyr::left_join(temp, names, by = c('ensembl' = 'id')))
    ## Set column of ensembl ids equal to corresponding gene symbols
    IDs[[1]] <- symbol_match[[2]]
    colnames(IDs)[1] <- 'gene_name'
  } else {
    temp <- as.data.frame(IDs)
    ## Set colname of temp so match can be made consistently
    colnames(temp)[1] <- 'ensembl'
    ## Perform left_join of gene symbols to ensembl ids
    symbol_match <- suppressMessages(dplyr::left_join(temp, names, by = c('ensembl' = 'id')))
    ## Set column of ensembl ids equal to corresponding gene symbols
    IDs <- as.data.frame(symbol_match[[2]])
    colnames(IDs)[1] <- 'gene_name'
  }
  return(IDs)
}

#' Function for converting between Ensembl IDs and gene symbols (canonical names)
#'
#' This function will output the converted names of each gene as a data frame variable, either replacing the names in a given data frame or simply returning a column.
#'
#' @inheritParams ensembl_to_symbol
#'
#' @return Returns either the original counts table with converted gene IDs, or a data frame of one column with the converted IDs
#'
#' @author Chris Stehn
#'
#' @examples
#' \dontrun{
#' ensembl_to_symbol(counts, mart = mmu, human = FALSE)
#' }
#' @export

symbol_to_ensembl <- function(IDs, mart = NULL, human = TRUE) {
  ## Read biomart
  if (missing(IDs)) {stop('No IDs have been supplied')}
  
  if (is.data.frame(IDs)) {
    if (grepl('ENS', IDs[1, 1])) {
      stop('Ensembl IDs have been supplied. Consider using ensembl_to_symbol instead')
    }
  } else {
    if (grepl('ENS', IDs[1])) {
      stop('Ensembl IDs have been supplied. Consider using ensembl_to_symbol instead')
    }
  }
  
  if (!is.null(mart)) {
    names <- dplyr::select(mart, id, gene_name) %>% dplyr::filter(!duplicated(gene_name))
  } else {
    if (human) {
      mart <- data(hsa, package = 'holmen')
    } else {
      mart <- data(mmu, package = 'holmen')
    }
    names <- dplyr::select(mart, id, gene_name) %>% dplyr::filter(!duplicated(gene_name))
  }

  ## Convert ensembl IDs to gene symbols using mart
  if (is.data.frame(IDs)) {
    temp <- as.data.frame(IDs[[1]])
    ## Set colname of temp so match can be made consistently
    colnames(temp)[1] <- 'symbol'
    ## Perform left_join of ensembl ids to gene symbols
    symbol_match <- suppressMessages(dplyr::left_join(temp, names, by = c('symbol' = 'gene_name')))
    ## Set column of ensembl ids equal to corresponding gene symbols
    IDs[[1]] <- symbol_match[[2]]
    colnames(IDs)[1] <- 'ensembl'
  } else {
    temp <- as.data.frame(IDs)
    ## Set colname of temp so match can be made consistently
    colnames(temp)[1] <- 'symbol'
    ## Perform left_join of ensembl ids to gene symbols
    symbol_match <- suppressMessages(dplyr::left_join(temp, names, by = c('symbol' = 'gene_name')))
    ## Set column of ensembl ids equal to corresponding gene symbols
    IDs <- as.data.frame(symbol_match[[2]])
    colnames(IDs)[1] <- 'ensembl'
  }
  return(IDs)
}
