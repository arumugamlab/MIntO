#!/usr/bin/env bash

# Which MIntO version are we using?
# Use specific tag by "tags/<TAG>" or "main"
# E.g.
# MINTO_STABLE_VERSION="tags/2.0.0-beta.2"

MINTO_STABLE_VERSION="main"

# Set MIntO and scratch locations

if [ ! -z "$COMPUTEROME_PROJ" ]; then
  # Danish computerome resource
  LOCAL_DIR="/home/projects/$COMPUTEROME_PROJ/scratch/$USER/MIntO/"
  MINTO_DIR="/home/projects/$COMPUTEROME_PROJ/apps/MIntO"
else
  LOCAL_DIR="/scratch/$USER/tmp/MIntO/"
  MINTO_DIR="$(pwd)/MIntO"
fi
CONDA_DIR="$MINTO_DIR/conda_env"

# Snakemake options

SNAKE_PARAMS="--use-conda --restart-times 1 --keep-going --latency-wait 60 --conda-prefix $CONDA_DIR --jobs 16 --cores 40 --resources mem=188"

# Where will the tutorial be tested?

TEST_DIR=$(pwd)

# Get MIntO or pull the latest if it already exists

if [ -d "$MINTO_DIR" ]; then
  cd $MINTO_DIR
  git checkout main
  git pull
  git checkout $MINTO_STABLE_VERSION
  cd $TEST_DIR
else
  cd $(dirname $MINTO_DIR)
  git clone https://github.com/arumugamlab/MIntO.git
  cd MIntO
  git checkout $MINTO_STABLE_VERSION
  cd $TEST_DIR
fi

# Download dependencies

echo -n "Downloading dependencies: "
cat $MINTO_DIR/testing/dependencies.yaml.in | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__LOCAL_DIR__>@$LOCAL_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > dependencies.yaml
time (snakemake --snakefile $MINTO_DIR/smk/dependencies.smk --configfile dependencies.yaml $SNAKE_PARAMS >& dependencies.log && echo "OK")

# Download raw data

if [ ! -d "IBD_tutorial_raw" ]; then
  echo -n "Downloading tutorial data: "
  wget https://zenodo.org/record/6369313/files/IBD_tutorial_raw.tar.gz
  tar xfz IBD_tutorial_raw.tar.gz
  echo "OK"
fi

# Get data
mkdir -p IBD_tutorial
cd IBD_tutorial
cp $MINTO_DIR/tutorial/metadata/tutorial_metadata.txt .
cp $MINTO_DIR/tutorial/build_hg18_subset.fna .

# Run metaG and metaT steps

OMICS="metaG"
for OMICS in metaG metaT; do
  echo ""
  echo "------------------"
  echo "Processing $OMICS:"
  echo "------------------"
  mkdir -p $OMICS
  cd $OMICS
  echo -n "QC_1: "
  cat $MINTO_DIR/testing/QC_1.yaml.in | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__LOCAL_DIR__>@$LOCAL_DIR@;s@<__TEST_DIR__>@$TEST_DIR@;s@<__OMICS__>@$OMICS@;" > QC_1.yaml
  time (snakemake --snakefile $MINTO_DIR/smk/QC_1.smk --configfile QC_1.yaml $SNAKE_PARAMS >& QC_1.log && echo "OK")
  if [ ! -f "QC_2.yaml.fixed" ]; then
    patch QC_2.yaml $MINTO_DIR/testing/QC_2.patch -o - | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__LOCAL_DIR__>@$LOCAL_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > QC_2.yaml.fixed
  fi
  echo -n "QC_2: "
  time (snakemake --snakefile $MINTO_DIR/smk/QC_2.smk --configfile QC_2.yaml.fixed $SNAKE_PARAMS >& QC_2.log && echo "OK")
  echo -n "ASSEMBLY: "
  time (snakemake --snakefile $MINTO_DIR/smk/assembly.smk --configfile assembly.yaml $SNAKE_PARAMS >& assembly.log && echo "OK")
  echo -n "BINNING_PREP: "
  time (snakemake --snakefile $MINTO_DIR/smk/binning_preparation.smk --configfile assembly.yaml $SNAKE_PARAMS >& binning_prep.log && echo "OK")
  echo -n "BINNING: "
  time (snakemake --snakefile $MINTO_DIR/smk/mags_generation.smk --configfile mags_generation.yaml $SNAKE_PARAMS >& mags.log && echo "OK")
  echo -n "GENE_ANNOTATION: "
  time (snakemake --snakefile $MINTO_DIR/smk/gene_annotation.smk --configfile mapping.yaml $SNAKE_PARAMS >& annotation.log && echo "OK")
  echo -n "GENE_ABUNDANCE: "
  time (snakemake --snakefile $MINTO_DIR/smk/gene_abundance.smk --configfile mapping.yaml $SNAKE_PARAMS >& abundance.log && echo "OK")
  cd ..
done

# Run integration

echo -n "DATA_INTEGRATION - TPM: "
time (snakemake --snakefile $MINTO_DIR/smk/data_integration.smk --configfile data_integration.yaml $SNAKE_PARAMS >& integration.TPM.metaGT.log && echo "OK")
sed "s/abundance_normalization: MG/abundance_normalization: TPM/" data_integration.yaml > data_integration.yaml.TPM
echo -n "DATA_INTEGRATION - MG: "
time (snakemake --snakefile $MINTO_DIR/smk/data_integration.smk --configfile data_integration.yaml.TPM $SNAKE_PARAMS >& integration.MG.metaGT.log && echo "OK")
