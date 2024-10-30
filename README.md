# SPIDER-paper

## Data availability
The user can download the example input from the following link and deposite them in the example_datasets folder.
- Example input data: stored in [input_datasets_small.zip](https://portland-my.sharepoint.com/:u:/g/personal/shiyingli7-c_my_cityu_edu_hk/EfGc0cv8EwlKjPNQmwi11ZwBEdoDQpCZdfBUn-WDNy_8zg?e=qbMLqy). 
- Full input data: stored in [input_datasets.zip](https://portland-my.sharepoint.com/:u:/g/personal/shiyingli7-c_my_cityu_edu_hk/EV3Xv7sSj75Fp6VBCHnOPAUBCdwTADbU70uqKv6fIq5rPw?e=yoGFM4).

## Example notebooks

### Case 1: Comparison between SPIDER and SpatialDM using the PDAC dataset
Run the notebook `PDAD.ipynb` with the provided adata from `./example_datasets/PDAC/PDAC_A/`. You only need to specify the `R_path` parameter. The default output folder will be `./results/PDAC/PDAC_A/`.

### Case 2: Analyses on SPIDER SVIs using the stereo-seq mouse brain dataset
The input adata is stored in the `/saw/mouse/` folder in [input_datasets_small.zip](https://portland-my.sharepoint.com/:u:/g/personal/shiyingli7-c_my_cityu_edu_hk/EfGc0cv8EwlKjPNQmwi11ZwBEdoDQpCZdfBUn-WDNy_8zg?e=qbMLqy). 
Run the notebook `SAW.ipynb` (need to specify the `R_path` parameter). The default output folder is `./results/saw/mouse/`. To skip the actual running time for this large dataset, the user can also download the output files in [saw.zip](https://portland-my.sharepoint.com/:u:/g/personal/shiyingli7-c_my_cityu_edu_hk/EWlpevc-GDtDj8q0RGZBPy4BRmw0rnwYIKf6wNBeOm6sfw?e=VKY5AF) and put the contents in `./results/saw/mouse/`.

### Case 3: Analyses on SPIDER SVI patterns using the slide-seq brain dataset
The input adata is stored in the `/slide-seq-v2/mousebrain/` folder in [input_datasets_small.zip](https://portland-my.sharepoint.com/:u:/g/personal/shiyingli7-c_my_cityu_edu_hk/EfGc0cv8EwlKjPNQmwi11ZwBEdoDQpCZdfBUn-WDNy_8zg?e=qbMLqy). 
Run notebooks `mouse_brain.ipynb` (need to specify the `R_path` parameter). The default output folder is `./results/slide-seq-v2/mousebrain/`. To skip the actual running time for this large dataset, the user can also download the output files in [slide-seq-v2.zip](https://portland-my.sharepoint.com/:u:/g/personal/shiyingli7-c_my_cityu_edu_hk/ESwmKf_lfBdJhGwWpWXmFvgB7e_vRy8QYuQgPZEl8iL9jA?e=PlCnbx) and put the contents in `./results/slide-seq-v2/mousebrain/`.