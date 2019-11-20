### Test Run
```
git clone https://github.com/jpfeil/typhon.git
cd typhon
docker run --rm -v $(pwd):/data/ jpfeil/typhon:0.1.1 --RSEM typhon/test/rsem_genes.results --diagnosis MYCN-NA-Neuroblastoma
```
