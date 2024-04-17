## This is adopted from StaG-mwc, a shotgun sequencing processing snakemake pipeline

### Add-on features:
1. Add kneaddata module in `workflow/rules/preproc/read_quality.smk`
2. Add a script `workflow/scripts/kneaddata_summary.py` to summarise pre-processing statistics for kneadata and modify `workflow/rules/preproc/preprocessing_summary.smk` to incoporate kneaddata summary.
3. Modify `workflow/rules/taxonomic_profiling/strainphlan.smk` to allow parallel run of StrainPhlAn on all available clades based on pre-specified prevalence in the sample batch.
4. Add `profiles/sge` and `profiles/pbs` to enable running on sge and pbs clusters. The sge cluster module was modified from https://github.com/Snakemake-Profiles/sge


### Check https://github.com/ctmrbio/stag-mwc for the latest git repo for stag-mwc
### Check https://stag-mwc.readthedocs.io for the full documentation of the lastest stag-mwc git repo

### Citing
If you find StaG-mwc useful in your research, please cite the Zenodo DOI:
https://zenodo.org/badge/latestdoi/125840716

