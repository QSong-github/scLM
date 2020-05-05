# scLM
A R-based tool to do the automatic identification of co-expressed genes across mulitple single cell RNA-seq datasets simultaneously

# 1. Installation
```
devtools::github_install("scLM")
```

# 2. Input scRNA-seq Data File Format
scLM works with multiple single cell RNA-seq dataset as inputs. It also works with one single cell dataset. Bascially, the format looks like the following. Example data files can be found in the ```Data``` folder.

        | Cell1ID | Cell2ID | Cell3ID | Cell4ID | Cell5ID  | ... |
|----|--------|--------|--------|---------|-----|
| Gene1 | 12 | 0 | 0 | 0 | ... |
| Gene2 | 125 | 0 | 298 | 0  | ... |
| Gene3 | 0 | 0| 0 | 0  | ... |
|...    |...|...|...|...|...|

The gourd truth labels for cells in each dataset can also be input. The format is as following

| Cell1ID | Lable1 |
|----|--------|
| Cell2ID | Lable2 |
| Cell3ID | Lable3 |
| Cell4ID | Lable4 |
|...    |..


# 3. How to use

## 3.1 load scLM package
```
library(scLM)
```
recommend the linux system to run the codes

## 3.2 read in scRNA-seq datasets
```
#' load the example data
data("example1")
# or 
# data("example2")
```
In this example, we define 3 co-expression clusters for each dataset of the input list example1 and example2

## 3.3 read in ground truth cell labels (this is optional)
```
data(example1.member)
# data(example2.member)
```

## 3.4 run the function with specific lambda across multiple datasets
```
result1 <- Multi_NB(datalist=example1, N=3, K=nrow(example1[[1]]))
```
result1 contains the latent variables accompanied with other coefficients, the identified co-expression clusters

## 3.5 evaluate of clustering results using ground truth (this is optional)
```
# calculate the Adjusted Rand Index

library(clues)
adjustedRand(result1$clusters,example1.member)
```

# 4. Examples and reproducible results 
can be found using the example.R script

# 5. without prior knowledge, run the function with several lambda across multiple datasets 

## 5.1 How to use

```
for (lambda in 1:20)
{
    results <- Multi_NB(datalist=example1, N=lambda, K=nrow(example1[[1]]))
    save(results,file=paste0('path1',lambda,'results.RData')）
}
```
Output results in a designated path "path1"

## 5.2 Identify the optimal lambda and co-expression clusters

```
files <- list.files(path=path1, pattern='results.RData')

# load all the results in the resA with the list structure

count <- 0
resA <- vector('list')
for ( i in files){
    load(i)
    count <- count +1
    resA[[count]] <-  results
    names(resA[[count]]) <- paste0('res_',strsplit(strsplit(i,'_')[[1]][4],'results')[[1]][1])}


# identify the result with least BIC value

bicS <- lapply(1:length(resA),function(i){ Res <- resA[[i]][[1]]$BIC })
optimal.lambda <- grep(min(unlist(bicS)),unlist(bicS))

# optimal co-expression clusters

load(files[optimal.lambda])
opitmal.cluster <- results$clusters
```
