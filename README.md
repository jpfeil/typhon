### Test Run
```
cd typhon/test
docker run --rm -v $(pwd):/data/ jpfeil/typhon:0.1.1 --RSEM rsem_genes.results --diagnosis MYCN-NA-Neuroblastoma
```
