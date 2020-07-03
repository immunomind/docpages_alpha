---
title: Repertoire overlap and public clonotypes
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

# v4\_overlap

\`\`\`{r setup, include=FALSE, echo=FALSE}

## knitr::knit\_hooks$set\(optipng = knitr::hook\_optipng\)

## knitr::opts\_chunk$set\(optipng = '-o7'\)

knitr::opts\_chunk$set\(echo = TRUE\) knitr::opts\_chunk$set\(fig.align = "center"\) knitr::opts\_chunk$set\(fig.width = 12\) knitr::opts\_chunk$set\(fig.height = 6\)

library\(immunarch\)

## source\("../R/testing.R"\)

## immdata = load\_test\_data\(\)

data\(immdata\)

```text
# Repertoire overlap
Repertoire overlap is the most common approach to measure repertoire similarity. It is achieved by computation of specific statistics on clonotypes shared between given repertoires, also called "public" clonotypes. `immunarch` provides several indices: 
- number of public clonotypes (`.method = "public"`) - a classic measure of overlap similarity.

- overlap coefficient (`.method = "overlap"`) - a normalised measure of overlap similarity. It is defined as the size of the intersection divided by the smaller of the size of the two sets.

- Jaccard index (`.method = "jaccard"`) - it measures similarity between finite sample sets, and is defined as the size of the intersection divided by the size of the union of the sample sets.

- Tversky index (`.method = "tversky"`) - an asymmetric similarity measure on sets that compares a variant to a prototype. If using default arguments, it's similar to Dice's coefficient.

- cosine similarity (`.method = "cosine"`) - a measure of similarity between two non-zero vectors

- Morisita's overlap index (`.method = "morisita"`) - a statistical measure of dispersion of individuals in a population. It is used to compare overlap among samples.

- incremental overlap - overlaps of the N most abundant clonotypes with incrementally growing N (`.method = "inc+METHOD"`, e.g., `"inc+public"` or `"inc+morisita"`).

The function that includes described methods is `repOverlap`. Again the output is easily visualised when passed to `vis()` function that does all the work:

```{r overlap, message=F, warning=FALSE, fig.width=12, fig.height=7}
imm_ov1 <- repOverlap(immdata$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(immdata$data, .method = "morisita", .verbose = F)

p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)

p1 + p2

vis(imm_ov1, "heatmap2")
```

You can easily change the number of significant digits: \`\`\`{r overlap-signif-digits, message=F, warning=FALSE, fig.width=12, fig.height=7} p1 &lt;- vis\(imm\_ov2, .text.size = 2.5, .signif.digits = 1\) p2 &lt;- vis\(imm\_ov2, .text.size = 2, .signif.digits = 2\)

p1 + p2

```text
<!---
Top-overlap [Work in Progress]
```{r, warning=TRUE, fig.width=14, fig.height=8}
warning("TODO")
```

--&gt;

To analyse the computed overlap measures function apply `repOverlapAnalysis`. \`\`\`{r overlap-1, warning=F, fig.width=8, fig.height=5}

## Apply different analysis algorithms to the matrix of public clonotypes:

## "mds" - Multi-dimensional Scaling

repOverlapAnalysis\(imm\_ov1, "mds"\)

## "tsne" - t-Stochastic Neighbor Embedding

repOverlapAnalysis\(imm\_ov1, "tsne"\)

## Visualise the results

repOverlapAnalysis\(imm\_ov1, "mds"\) %&gt;% vis\(\)

```text
```{r overlap-2, warning=F, fig.width=10, fig.height=5}
# Apply different analysis algorithms to the matrix of public clonotypes:
# "mds" - Multi-dimensional Scaling
repOverlapAnalysis(imm_ov1, "mds")
# "tsne" - t-Stochastic Neighbor Embedding
repOverlapAnalysis(imm_ov1, "tsne")

# Visualise the results
repOverlapAnalysis(imm_ov1, "mds") %>% vis()

# Clusterise the MDS resulting components using K-means
repOverlapAnalysis(imm_ov1, "mds+kmeans") %>% vis()
```

## Public repertoire

In order to build a massive table with all clonotypes from the list of repertoires use the `pubRep` function. \`\`\`{r, warning=TRUE, fig.width=14, fig.height=8}

## Pass "nt" as the second parameter to build the public repertoire table using CDR3 nucleotide sequences

pr.nt &lt;- pubRep\(immdata$data, "nt", .verbose = F\) pr.nt

```text
```{r, warning=TRUE, fig.width=14, fig.height=8}
# Pass "aa+v" as the second parameter to build the public repertoire table using CDR3 aminoacid sequences and V alleles
# In order to use only CDR3 aminoacid sequences, just pass "aa"
pr.aav <- pubRep(immdata$data, "aa+v", .verbose = F)
pr.aav
```

\`\`\`{r, eval=FALSE, warning=TRUE, fig.width=14, fig.height=8}

## You can also pass the ".coding" parameter to filter out all noncoding sequences first:

pr.aav.cod &lt;- pubRep\(immdata$data, "aa+v", .coding = T\)

```text
```{r, eval=FALSE, warning=TRUE, fig.width=14, fig.height=8}
# Create a public repertoire with coding-only sequences using both CDR3 amino acid sequences and V genes
pr <- pubRep(immdata$data, "aa+v", .coding = T, .verbose = F)

# Apply the filter subroutine to leave clonotypes presented only in healthy individuals
pr1 <- pubRepFilter(pr, immdata$meta, c(Status = "C"))

# Apply the filter subroutine to leave clonotypes presented only in diseased individuals
pr2 <- pubRepFilter(pr, immdata$meta, c(Status = "MS"))

# Divide one by another
pr3 <- pubRepApply(pr1, pr2)

# Plot it
p <- ggplot() +
  geom_jitter(aes(x = "Treatment", y = Result), data = pr3)
p
```
