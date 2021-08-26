# SFARI Readviz Pipeline
This is the calling pipeline that has been used to create de-identified bam files for variants in the SFARI browser to be visualized. It performs the follow tasks:
* Selects example heterozygous, homozygous or hemizygous variants from hail matrix table.
* Create per sample input files for GATK HaplotypeCaller to create bam files. This is the read data that HaplotypeCaller "sees" before calling a variant. 
* De-identify bam files by stripping all details from each read and then combine in individual bam files into a grouped bam file
* Create an sqlite file so that variants can be searched and mapped to corresponding grouped bam file 

All code is based the gnomAD readviz pipeline and scripts developed and tested for Google Cloud Platform and has been adapted to work on Institutional clusters. The gnomAD pipeline was implemented using a collection of scripts, while this implementation consolidates it into one WDL workflow.

## Requirements

* Java 1.8
* Python
* GATK (tested on 4.1.7.0)
* Cromwell (tested on v56)
* Sqlite3

### Python libraries required
```
hail
peewee
pysam
tqdm
```

### Running WDL using Cromwell
The Readviz WDL pipeline is run using cromwell using the following command line.
```
java -Dconfig.file=slurm.conf -jar \
/home/ml2529/shared/tools/jars/cromwell-56.jar run \
MultiSampleReadViz.wdl \
-i readviz_inputs.json  \
-o workflow.options

```

The `launch.sh` file is a batch script that can be used to launch in Slurm by the following command  
```
sbatch launch.sh
```

### Inputs
The main input file is specified in `readviz_inputs.json`. In summary this contains the location of the follow input files:
* Hail matrix table of variant calls that variants will be sampled
* sample IDs to include
* Location of bam/cram files for each sample
* hg38 reference sequence and associated index files
* location of required Python scripts

### Ouputs
There a two major outputs that can then be used by the SFARI browser (in the IGV window):
* Various group BAM files with read data for each variant. This can be served as a static file by the web server.
* Sqlite database file for each chromosome. This contains information on which file(s) a variant can be found, which can be queried by the GraphQL API.










