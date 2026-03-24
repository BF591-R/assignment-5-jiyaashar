library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
    
    # Read the full counts matrix — genes are rows, samples are columns
    # row.names=1 makes the first column (gene IDs) into rownames
    counts <- read.table(counts_csv, header = TRUE, sep = '\t', row.names = 1)
    
    # Read the sample metadata
    meta <- read.csv(metafile_csv, header = TRUE)
    
    # Keep only rows for the timepoints we want, and only the two columns we need
    meta_subset <- meta %>%
      filter(timepoint %in% selected_times) %>%
      dplyr::select(samplename, timepoint)
    
    # Set vP0 as the reference factor level — without this R would pick
    # alphabetically (vAd first), which would flip all fold change directions
    meta_subset$timepoint <- factor(
      meta_subset$timepoint,
      levels = c('vP0', selected_times[selected_times != 'vP0'])
    )
    
    # Reorder counts columns to match the row order of our colData
    # This is required — SummarizedExperiment will error if they don't align
    counts_subset <- counts[, meta_subset$samplename]
    
    # DESeq2 requires a matrix, not a data.frame
    counts_matrix <- as.matrix(counts_subset)
    
    # Build the SummarizedExperiment object
    # The assay MUST be named "counts" — the test checks this explicitly
    se <- SummarizedExperiment(
      assays = list(counts = counts_matrix),
      colData = meta_subset
    )
    
    return(se)
  }

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
    # Wrap the SummarizedExperiment in a DESeqDataSet with our model formula
    dds <- DESeqDataSet(se, design = design)
    
    # DESeq() runs three steps internally:
    #   1. estimateSizeFactors() — median-of-ratios normalization
    #   2. estimateDispersions() — per-gene dispersion estimates
    #   3. nbinomWaldTest() — negative binomial GLM + Wald test
    dds <- DESeq(dds)
    
    # Extract results as a dataframe
    # Compares last factor level (vAd) vs reference level (vP0)
    # so positive log2FC = higher in adult
    res <- results(dds) %>% as.data.frame()
    
    # Return BOTH objects — dds is needed later to extract normalized counts
    return(list(dds = dds, results = res))
  }


#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
    labeled <- deseq2_res %>%
      # rownames='genes' promotes the row names (Ensembl IDs) to a proper column
      as_tibble(rownames = 'genes') %>%
      mutate(volc_plot_status = case_when(
        # Check !is.na(padj) FIRST — NA < 0.10 returns NA, not FALSE
        !is.na(padj) & padj < padj_threshold & log2FoldChange > 0 ~ 'UP',
        !is.na(padj) & padj < padj_threshold & log2FoldChange < 0 ~ 'DOWN',
        TRUE ~ 'NS'  # catch-all: covers NA padj and non-significant genes
      ))
    
    return(labeled)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
    p <- labeled_results %>%
      filter(!is.na(pvalue)) %>%
      ggplot(aes(x = pvalue)) +
      geom_histogram(bins = 50, fill = 'steelblue', color = 'white') +
      labs(
        title = 'Distribution of Raw P-values',
        x = 'P-value (unadjusted)',
        y = 'Count'
      ) +
      theme_bw()
    
    return(p)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
    p <- labeled_results %>%
      filter(!is.na(padj) & padj < padj_threshold) %>%
      ggplot(aes(x = log2FoldChange)) +
      geom_histogram(bins = 50, fill = 'salmon', color = 'white') +
      labs(
        title = paste0('Log2 Fold Change for DE Genes (padj < ', padj_threshold, ')'),
        x = 'Log2 Fold Change',
        y = 'Count'
      ) +
      theme_bw()
    
    return(p)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
    # Get the top N gene IDs sorted by ascending padj (most significant first)
    top_genes <- labeled_results %>%
      filter(!is.na(padj)) %>%
      arrange(padj) %>%
      head(num_genes) %>%
      pull(genes)
    
    # counts(dds, normalized=TRUE) divides raw counts by size factors
    # This makes samples comparable despite different sequencing depths
    norm_counts <- counts(dds_obj, normalized = TRUE)
    
    # Subset to top genes and reshape to long format for ggplot
    # pivot_longer converts from wide (one col per sample) to long (one row per sample per gene)
    plot_data <- norm_counts[top_genes, ] %>%
      as.data.frame() %>%
      rownames_to_column('genes') %>%
      pivot_longer(cols = -genes, names_to = 'samplename', values_to = 'norm_count')
    
    p <- plot_data %>%
      ggplot(aes(x = genes, y = log10(norm_count + 1), color = samplename)) +
      geom_jitter(width = 0.2, size = 2) +
      labs(
        title = paste('Normalized Counts: Top', num_genes, 'DE Genes by padj'),
        x = 'Gene',
        y = 'log10(Normalized Count + 1)',
        color = 'Sample'
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
    p <- labeled_results %>%
      filter(!is.na(padj)) %>%
      ggplot(aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
      geom_point(alpha = 0.6, size = 0.8) +
      scale_color_manual(
        values = c('UP' = 'firebrick', 'DOWN' = 'steelblue', 'NS' = 'grey60')
      ) +
      labs(
        title = 'Volcano Plot: Adult vs Postnatal Day 0',
        x = 'Log2 Fold Change',
        y = '-log10(Adjusted P-value)',
        color = 'Status'
      ) +
      theme_bw()
    
    return(p)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
    # Read the two-column ID->symbol lookup table
    id2gene <- read.table(id2gene_path, header = FALSE, sep = '\t',
                          col.names = c('genes', 'symbol'))
    
    ranked <- labeled_results %>%
      # Join MGI symbols onto results using the Ensembl gene ID column
      left_join(id2gene, by = 'genes') %>%
      # Remove genes with no matching symbol or missing log2FC
      filter(!is.na(symbol), !is.na(log2FoldChange)) %>%
      # For duplicate symbols (multiple Ensembl IDs -> same symbol),
      # keep the one with the largest absolute fold change
      group_by(symbol) %>%
      slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      # fgsea requires the vector sorted descending by the ranking metric
      arrange(desc(log2FoldChange))
    
    # setNames() builds the named vector: names=symbols, values=log2FC
    rnk_list <- setNames(ranked$log2FoldChange, ranked$symbol)
    
    return(rnk_list)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
    # gmtPathways() parses the GMT format into a named list of character vectors
    # Each element = one pathway, containing a vector of gene symbols
    pathways <- gmtPathways(gmt_file_path)
    
    # Set seed for reproducibility — fgsea uses random permutations internally
    set.seed(42)
    
    fgsea_res <- fgsea(
      pathways = pathways,
      stats    = rnk_list,
      minSize  = min_size,
      maxSize  = max_size
    )
    
    return(as_tibble(fgsea_res))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
    # Get top N by highest NES (enriched in upregulated genes)
    top_pos <- fgsea_results %>%
      arrange(desc(NES)) %>%
      head(num_paths)
    
    # Get top N by lowest NES (enriched in downregulated genes)
    top_neg <- fgsea_results %>%
      arrange(NES) %>%
      head(num_paths)
    
    plot_data <- bind_rows(top_pos, top_neg) %>%
      mutate(
        direction     = ifelse(NES > 0, 'Positive', 'Negative'),
        # Truncate long pathway names so they fit on the axis
        pathway_label = str_trunc(pathway, 60)
      )
    
    p <- plot_data %>%
      ggplot(aes(x = reorder(pathway_label, NES), y = NES, fill = direction)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = c('Positive' = 'firebrick', 'Negative' = 'steelblue')) +
      coord_flip() +
      labs(
        title = paste('Top', num_paths, 'Positive and Negative NES Pathways'),
        x     = 'Pathway',
        y     = 'Normalized Enrichment Score (NES)',
        fill  = 'Direction'
      ) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 7))
    
    return(p)
}


