### Test Run
```
cd typhon
docker run --rm -v $(pwd):/data/ jpfeil/typhon:0.1.1 --RSEM test/rsem_genes.results --diagnosis MYCN-NA-Neuroblastoma
```
