# Flatiron specific instructions
## Setting up conda environment
```
conda env create -n sfari -f hail_conda.yml
conda activate sfari
```

## Preparing input files
### Scripts require two input file to be present
`iWES.bam.list`  - each line has full path of bam file. 
`iWES_sample_ids.txt` - one sampleID per line
`
An example for chromosome 2
```
./create_key_by_sample.sh 2
./create_tsv_files.sh 2
./create_wdl_input.sh 2
```

## Running the WDL
Flatiron slurm was not made to run things like WDL. The way around this is to run the WDL in local mode on a node with a high number of CPUs. At the time the `rome` nodes with 128 CPUs served this purpose.

```
srun -N 1 --exclusive -p info --constraint=rome --job-name=intsesh --pty bash -i
```

* Edit workflow.options to change output directory (for each chromosome)
* Edit readviz_inputs_consolidated.json to change the input TSV file and the .chr to the chromosome (i.e. 20_1)
* `bash run.sh`









