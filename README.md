### Test Run
```
git clone https://github.com/jpfeil/typhon.git
cd typhon
docker run --rm -v $(pwd):/data/ jpfeil/typhon:0.1.1 --RSEM typhon/test/rsem_genes.results --diagnosis MYCN-NA-Neuroblastoma
```

### Adding new models 
There are two analyses that are run for each diagnosis. The first analysis is the enrichment analysis that looks for multivariate expression of genes. The second anlaysis is a univariate analysis that looks for a set of druggable gene targets.

Multivariate hydra models should be placed under the `models/multivariate/<diagnosis>` directory. Each analysis should be uniquely named. The analysis directory needs the hydra model output directory and a features directory associated with the clusters in the model. If the directory structure matches the other diseases, then this new multivariate model should be incorporated into the typhon framework.

Similarly the gene-level analysis should be placed under the `models/univariate/<diagnosis>` directory. The hydra multimodal gene fits from the `filter` command should be copied. This directory also needs a druggable-genes directory with the training data from the hydra analysis and the DGIdb gene-drug interactions.
