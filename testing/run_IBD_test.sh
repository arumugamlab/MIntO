#!/usr/bin/env bash

LOCAL_DIR="/scratch/$USER/tmp2/"
TEST_DIR=`pwd`
MINTO_DIR="$TEST_DIR/MIntO"
CONDA_DIR="$TEST_DIR/conda_env"
SNAKE_PARAMS="--use-conda --restart-times 1 --keep-going --latency-wait 60 --conda-prefix $CONDA_DIR --jobs 16 --cores 160 --resources mem=1700"

# Get MIntO or pull the latest if it already exists

if [ -d "MIntO" ]; then
  cd MIntO
  git pull
else
  git clone https://github.com/arumugamlab/MIntO.git
fi

# Download dependencies

cat $MINTO_DIR/testing/dependencies.yaml.in | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__LOCAL_DIR__>@$LOCAL_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > dependencies.yaml
time (snakemake --snakefile $MINTO_DIR/smk/dependencies.smk --configfile dependencies.yaml $SNAKE_PARAMS >& dependencies.log && echo "OK")

# Download raw data

if [ ! -d "IBD_tutorial_raw" ]; then
  wget https://zenodo.org/record/6369313/files/IBD_tutorial_raw.tar.gz
  tar xfz IBD_tutorial_raw.tar.gz
fi

# Get data
mkdir -p IBD_tutorial
cd IBD_tutorial
cp $MINTO_DIR/tutorial/metadata/tutorial_metadata.txt .
cp $MINTO_DIR/tutorial/build_hg18_subset.fna .

# Run metaG steps

OMICS="metaG"
mkdir -p $OMICS
cd $OMICS
cat $MINTO_DIR/testing/QC_1.yaml.in | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__LOCAL_DIR__>@$LOCAL_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > QC_1.yaml
time (snakemake --snakefile $MINTO_DIR/smk/QC_1.smk --configfile QC_1.yaml $SNAKE_PARAMS >& QC_1.log && echo "OK")
if [ ! -f "QC_2.yaml.fixed" ]; then
  patch QC_2.yaml $MINTO_DIR/testing/QC_2.patch -o - | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__LOCAL_DIR__>@$LOCAL_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > QC_2.yaml.fixed
fi
time (snakemake --snakefile $MINTO_DIR/smk/QC_2.smk --configfile QC_2.yaml.fixed $SNAKE_PARAMS >& QC_2.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/assembly.smk --configfile assembly.yaml $SNAKE_PARAMS >& assembly.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/binning_preparation.smk --configfile assembly.yaml $SNAKE_PARAMS >& binning_prep.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/mags_generation.smk --configfile mags_generation.yaml $SNAKE_PARAMS >& mags.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/gene_annotation.smk --configfile mapping.yaml $SNAKE_PARAMS >& annotation.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/gene_abundance.smk --configfile mapping.yaml $SNAKE_PARAMS >& abundance.log && echo "OK")
cd ..

# Run metaT steps

OMICS="metaT"
mkdir -p $OMICS
cd $OMICS
cat ../metaG/QC_1.yaml | sed "s/metaG/metaT/g" > QC_1.yaml
time (snakemake --snakefile $MINTO_DIR/smk/QC_1.smk --configfile QC_1.yaml $SNAKE_PARAMS >& QC_1.log && echo "OK")
if [ ! -f "QC_2.yaml.fixed" ]; then
  patch QC_2.yaml $MINTO_DIR/testing/QC_2.patch -o - | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__LOCAL_DIR__>@$LOCAL_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > QC_2.yaml.fixed
fi
time (snakemake --snakefile $MINTO_DIR/smk/QC_2.smk --configfile QC_2.yaml.fixed $SNAKE_PARAMS >& QC_2.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/assembly.smk --configfile assembly.yaml $SNAKE_PARAMS >& assembly.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/binning_preparation.smk --configfile assembly.yaml $SNAKE_PARAMS >& binning_prep.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/mags_generation.smk --configfile mags_generation.yaml $SNAKE_PARAMS >& mags.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/gene_annotation.smk --configfile mapping.yaml $SNAKE_PARAMS >& annotation.log && echo "OK")
time (snakemake --snakefile $MINTO_DIR/smk/gene_abundance.smk --configfile mapping.yaml $SNAKE_PARAMS >& abundance.log && echo "OK")
cd ..

# Run integration

time (snakemake --snakefile $MINTO_DIR/smk/data_integration.smk --configfile data_integration.yaml $SNAKE_PARAMS >& integration.TPM.metaGT.log && echo "OK")
sed "s/abundance_normalization: MG/abundance_normalization: TPM/" data_integration.yaml > data_integration.yaml.TPM
time (snakemake --snakefile $MINTO_DIR/smk/data_integration.smk --configfile data_integration.yaml.TPM $SNAKE_PARAMS >& integration.MG.metaGT.log && echo "OK")
