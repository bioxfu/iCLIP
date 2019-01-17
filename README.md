# iCLIP Workflow Tutorial
### 1. Install iCount (optional)
```
conda create -n icount icount python=3.6 snakemake=3.13.3 fastqc jellyfish=1.1.11 fastx_toolkit sra-tools igvtools
source activate icount
```

### 2. Preparing a genome index (optional)
```
# Get list of available ENSEMBL releases:
iCount releases
# Get list of available species for given ENSEMBL release:
iCount species -r 88

# Download the human/mouse genome sequence from release 88:
iCount genome homo_sapiens -r 88
iCount genome mus_musculus -r 88

# Download the annotation of human/mouse genome from release 88:
iCount annotation homo_sapiens -r 88
iCount annotation mus_musculus -r 88

# Generate human/mouse genome index that is used by STAR mapper:
mkdir hs88
iCount indexstar homo_sapiens.88.fa.gz hs88 --annotation homo_sapiens.88.gtf.gz --threads 20

mkdir mm88
iCount indexstar mus_musculus.88.fa.gz mm88 --annotation mus_musculus.88.gtf.gz --threads 20

# Generate a new annotation file with human/mouse genome segmentation:
iCount segment homo_sapiens.88.gtf.gz hs88seg.gtf.gz homo_sapiens.88.fa.gz.fai 
iCount segment mus_musculus.88.gtf.gz mm88seg.gtf.gz mus_musculus.88.fa.gz.fai 
```

### 3. Make project directory
```
# the project directory contains specific GROUP name and current DATE
GROUP=XXX
DATE=`date +"%Y%m%d"`
mkdir ~/Project/${GROUP}_${DATE}
```

### 4. Clone the repository
```
cd ~/Project/${GROUP}_${DATE}
git clone https://github.com/bioxfu/iCLIP
cd iCLIP
```

### 5. Create *config.yaml* and *Snakefile* based on the examples
```
cp example/example.config.yaml config.yaml
cp example/example.Snakefile Snakefile

# edit config.yaml 
```

### 6. Initiate the project
```
source init.sh
```

### 7. Dry run the workflow to check any mistakes
```
./dry_run.sh
```

### 8. Run the workflow
```
# if you are working on the HPC
./run_HPC.sh

# if you are working on the local machine (need large memory)
./run.sh

# check the workflow progress in nohup.out file
tail nohup.log 

# check the jobs on HPC
qstat

# if you get the error: Directory cannot be locked.
snakemake --unlock 
```

### 9. Remove the temporary files
```
./clean.sh
```

