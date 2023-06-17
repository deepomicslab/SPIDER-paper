# SPIDER-paper

## Test small datasets
Run notebooks PDAD.ipynb and mouse_embryo.ipynb with adata provided in example_datasets. You only need to specify the R_path param, and the ourput files will be in example_datasets.

## Test large datasets
We provide one sample result for each datasets. The output files can be downloaded here: https://www.dropbox.com/scl/fi/4txttotme3lzgbmzqc02i/example_output.zip?dl=0&rlkey=ua8nokrbv7b8rjxpbhwdfnai4. The files contains pre-run results so you can run all notebooks without reproducing the results. Make sure you set the out_f parma correctly - if you extreact the downloaded zip file in this repo, then you should set:
``
out_f = f'../example_output/{ds}/{sample_name}/'
``

## Run all datasets
For interface construction and SVI identification, refer to scripts/01.run_all.py.
For interface feature extraction, clustering, and pseudotime analysis, refer to scripts/02.spatialpca.py.