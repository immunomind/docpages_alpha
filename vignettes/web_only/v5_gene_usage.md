---
title: Gene usage analysis
author: '<a href="https://immunomind.io">ImmunoMind</a>'
date: support@immunomind.io
output:
  html_document:
    fig_height: 8
    fig_width: 10
    theme: spacelab
    toc: 'yes'
  pdf_document:
    toc: 'yes'
  word_document:
    toc: 'yes'
---

# v5\_gene\_usage

\`\`\`{r setup, include=FALSE, echo=FALSE}

## knitr::knit\_hooks$set\(optipng = knitr::hook\_optipng\)

## knitr::opts\_chunk$set\(optipng = '-o7'\)

knitr::opts\_chunk$set\(echo = TRUE\) knitr::opts\_chunk$set\(fig.align = "center"\) knitr::opts\_chunk$set\(fig.width = 12\) knitr::opts\_chunk$set\(fig.height = 6\)

library\(immunarch\)

## source\("../R/testing.R"\)

## immdata = load\_test\_data\(\)

data\(immdata\)

```text
# Gene usage computation
`immunarch` comes with a gene segments data table containing known gene segments for several species following the [IMGT](http://www.imgt.org/IMGTrepertoire/LocusGenes/) nomenclature. In order to get the current statistics of genes, call the `gene_stats()` function:

```{r gene-usage}
gene_stats()
```

To compute the distributions of genes, `immunarch` includes the `geneUsage` function. It receives a repertoire or a list of repertoires as input and genes and species for which you want to get the statistics. E.g., if you plan to use `TRBV` genes of `Homo Sapiens`, you need to use the `hs.trbv` string in the function, where `hs` comes from the `alias` column and `trbv` is the gene name. Of if you plan to use `IGHJ` genes of `Mus Musculus`, you need to use `musmus.ighj`:

```text
# Next four function calls are equal. "hs" is from the "alias" column.
imm_gu <- geneUsage(immdata$data, "hs.trbv")
# imm_gu = geneUsage(immdata$data, "HomoSapiens.trbv")
# imm_gu = geneUsage(immdata$data, "hs.TRBV")
# imm_gu = geneUsage(immdata$data, "HomoSapiens.TRBV")

imm_gu
```

Gene distributions could be computed either using counts of individual clonotypes \(`.quant = "count"`\) or not using them \(`.quant = NA`\).

In order to compute allele-level or family-level distributions, change the `.type` parameter.

Parameter `.norm` controls whether `immunarch` will normalise the data to ensure the sum of all frequencies to be equal 1 or not.

You can visualise the histogram of gene usage in different ways:

\`\`\`{r, message=F, warning=FALSE, fig.width=12, fig.height=5}

## Compute the distribution of the first two samples

imm\_gu &lt;- geneUsage\(immdata$data\[c\(1, 2\)\], "hs.trbv", .norm = T\)

vis\(imm\_gu\)

```text
```{r, message=F, warning=FALSE, fig.width=12, fig.height=5}
imm_gu <- geneUsage(immdata$data, "hs.trbv", .norm = T)

vis(imm_gu, .by = "Status", .meta = immdata$meta)
```

\`\`\`{r, message=F, warning=F, fig.width=15, fig.height=8} vis\(imm\_gu, .grid = T\)

```text
Another practical approach to the visualisation of group distributions are box plots:

```{r, message=F, warning=FALSE, fig.width=12, fig.height=5}
vis(imm_gu, .by = "Status", .meta = immdata$meta, .plot = "box")
```

### Ambiguity of gene segment names

Due to the ambiguity of gene alignments for some clonotypes, `geneUsage` has the following options to deal with ambiguous data:

* `.ambig = "inc"` - includes all possible combinations of ambiguous gene alignments from the data. NOTE: ImmunoSEQ formats use non-standart gene segment names, so it is preferable to use this argument value with ImmunoSEQ formats. This argument is ON by default to ease the gene manipulation. Feel free to change it to `"exc"` in case of other data formats. It is ON by default, we recommend it to leave it that way.
* `.ambig = "exc"` - filters out all clonotypes with ambiguous gene alignments.
* `.ambig = "wei"` - introduces weighted approach \(divides by n \(`1/n`\) the frequency for each entry of the corresponding gene if there are `n` genes for a clonotype\).
* `.ambig = "maj"` - chooses only the first gene segment.

## Gene usage analysis

To analyse the gene usage `immunarch` introduces the `geneUsageAnalysis` function. The `.method` parameter controls how the data is going to be preprocessed and analysed. `geneUsageAnalysis` includes following methods for preprocessing:

* "js" - Jensen-Shannon Divergence.
* "cor" - correlation.
* "cosine" - cosine similarity.
* "pca" - principal component analysis.
* "mds" - multi-dimensional scaling.
* "tsne" - t-Distributed Stochastic Neighbor Embedding.

And a few methods for the actual analysis:

* "hclust" - clusters the data using hierarchical clustering.
* "kmeans" - clusters the data using K-means.
* "dbscan" - clusters the data using DBSCAN.
* "kruskall" - compute Kruskall for each gene separately on data splitted to groups \(without preprocessing\). Results could be used with Dunn test in order to detect significant differences between groups.

You can call several methods in a single line of code, which is probably the most powerful feature of the package. For instance, `"js+hclust"` first computes Jensen-Shannon divergence and then applies hierarchical clustering on the resulting distance matrix, whereas `"anova"` computes ANOVA on each gene separately after repertoires have been grouped:

\`\`\`{r, warning=FALSE, fig.width=12, fig.height=5} imm\_gu &lt;- geneUsage\(immdata$data, "hs.trbv", .norm = T\)

imm\_gu\_js &lt;- geneUsageAnalysis\(imm\_gu, .method = "js", .verbose = F\) imm\_gu\_cor &lt;- geneUsageAnalysis\(imm\_gu, .method = "cor", .verbose = F\)

p1 &lt;- vis\(imm\_gu\_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5\) p2 &lt;- vis\(imm\_gu\_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5\)

p1 + p2

```text
Now let us visualise the output after both preprocessing and analysis:

```{r, warning=FALSE, fig.width=10, fig.height=4}
imm_gu_js[is.na(imm_gu_js)] <- 0

vis(geneUsageAnalysis(imm_gu, "cosine+hclust", .verbose = F))

# vis(geneUsageAnalysis(imm_gu, "js+dbscan", .verbose = F))
```

On top of that you can add clustering:

\`\`\`{r, message=F, warning=F, fig.width=12, fig.height=4} imm\_cl\_pca &lt;- geneUsageAnalysis\(imm\_gu, "js+pca+kmeans", .verbose = F\) imm\_cl\_mds &lt;- geneUsageAnalysis\(imm\_gu, "js+mds+kmeans", .verbose = F\) imm\_cl\_tsne &lt;- geneUsageAnalysis\(imm\_gu, "js+tsne+kmeans", .perp = .01, .verbose = F\)

p1 &lt;- vis\(imm\_cl\_pca, .plot = "clust"\) p2 &lt;- vis\(imm\_cl\_mds, .plot = "clust"\) p3 &lt;- vis\(imm\_cl\_tsne, .plot = "clust"\) p1 + p2 + p3

```text
You can regulate the number of clusters as well:
```{r, message=F, warning=F, fig.width=8, fig.height=4}
imm_cl_pca2 <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .k = 3, .verbose = F)
vis(imm_cl_pca2)
```

## Spectratyping

Spectratype is a useful way to represent distributions of genes per sequence length. Parameter `.quant` controls the quantity that used to compute proportions of genes - either by clonotype \(`id`\) or by number of clones per clonotype \(`count`\). Parameter `.col` controls which column to choose, e.g., "nt" for lengths of CDR3 nucleotide sequences only \(without grouping by gene segments\), "aa+v" for lengths of CDR3 amino acid sequences \(grouped by V gene segments\).

\`\`\`{r spectr, fig.width=12, fig.height=4} p1 &lt;- vis\(spectratype\(immdata$data\[\[1\]\], .quant = "id", .col = "nt"\)\) p2 &lt;- vis\(spectratype\(immdata$data\[\[1\]\], .quant = "count", .col = "aa+v"\)\)

p1 + p2

\`\`\`

