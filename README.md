# SPIDER-paper

## Test small datasets
Run notebooks PDAD.ipynb and mouse_embryo.ipynb with adata provided in example_datasets. You only need to specify the R_path param, and the ourput files will be in example_datasets. Or you can run the .py file version in the scripts folders.

## Test large datasets
We provide one sample result for each datasets. The output files can be downloaded here: https://www.dropbox.com/scl/fi/4txttotme3lzgbmzqc02i/example_output.zip?dl=0&rlkey=ua8nokrbv7b8rjxpbhwdfnai4. The files contains pre-run results so you can run all notebooks without reproducing the results. Make sure you set the out_f parma correctly - if you extreact the downloaded zip file in this repo, then you should set:
``
out_f = f'../example_output/{ds}/{sample_name}/'
``
If you want to run the scripts for indicidual dataset sample and reproduce the results, you can download all the input adatas here: Download all the input adatas here: https://www.dropbox.com/scl/fi/viu7ahop16kaeovsymtx5/input_datasets.zip?dl=0&rlkey=jjomqnmddv0xiphov7yw4uw57. Make sure you set the out_f parma correctly - if you extreact the downloaded zip file in this repo, then you should set:
``
out_f = f'../input_datasets/{ds}/{sample_name}/'
``

## Run all datasets
Download all the input adatas here: https://www.dropbox.com/scl/fi/viu7ahop16kaeovsymtx5/input_datasets.zip?dl=0&rlkey=jjomqnmddv0xiphov7yw4uw57.
For interface construction and SVI identification, refer to scripts/01.run_all.py.
For example, to run SPIDER on the PDAC dataset:
``
python 01.run_all.py 1
``
For interface feature extraction, clustering, and pseudotime analysis, refer to scripts/02.spatialpca.py.
For example, to run SPIDER on the DLPFC sample 151673:
``
python 02.spatialpca.py 1
``