#!/usr/bin/env bash

# Which MIntO version are we using?
# Use specific tag by "tags/<TAG>" or "main"
# E.g.
# MINTO_STABLE_VERSION="tags/2.1.0"
# Developers use 'main' but users should stick to stable versions.

MINTO_STABLE_VERSION="tags/2.2.0"
TEST_SLURM="" # leave it empty to turn off slurm-testing

function profile_command() {
    local cmd=$1
    echo $cmd >> $COMMAND_LOG
    /usr/bin/time -f 'real\t%E\nuser\t%U\nsystem\t%S\navgmem\t%KKB\nmaxmem\t%MKB\ncpu\t%P' bash -c "$cmd && echo OK"
    echo ""
}

# Set MIntO and scratch locations

if [ ! -z "$COMPUTEROME_PROJ" ]; then
  # Danish computerome resource
  SHADOWDIR="/home/projects/$COMPUTEROME_PROJ/scratch/$USER/MIntO/"
  MINTO_DIR="/home/projects/$COMPUTEROME_PROJ/apps/MIntO"
else
  SHADOWDIR="/scratch/$USER/tmp/MIntO/"
  MINTO_DIR="$(pwd)/MIntO"
fi

# Where will the tutorial be tested?

TEST_DIR=$(pwd)
CONDA_DIR="$TEST_DIR/conda_env"

# Get MIntO or pull the latest if it already exists

echo "-------------"
echo "GETTING MINTO"
echo "-------------"

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

# Set up cluster definition

sed -i -e "s%^CLUSTER_WORKLOAD_MANAGER.*%CLUSTER_WORKLOAD_MANAGER = 'slurm'%" $MINTO_DIR/site/cluster_def.py
sed -i -e "s%^CLUSTER_LOCAL_DIR.*%CLUSTER_LOCAL_DIR = '/scratch/MIntO/mirror'%" $MINTO_DIR/site/cluster_def.py

# Record the snakemake commands that have been run

COMMAND_LOG='commands_dependencies.txt'

# Snakemake options
if [ ! -z "$COMPUTEROME_PROJ" ]; then
  SNAKE_PARAMS="--use-conda --restart-times 0 --keep-going --latency-wait 60 --conda-prefix $CONDA_DIR --shadow-prefix $SHADOWDIR --jobs 16 --default-resources gpu=0 mem=4 --cluster 'qsub -d $(pwd) -W group_list=$COMPUTEROME_PROJ -A $COMPUTEROME_PROJ -N {name} -l nodes=1:thinnode:ppn={threads},mem={resources.mem}gb,walltime=7200 -V -v TMPDIR=$SHADOWDIR' --local-cores 4"
elif [ ! -z "$TEST_SLURM" ]; then
  SNAKE_PARAMS="--use-conda --restart-times 0 --keep-going --latency-wait 60 --conda-prefix $CONDA_DIR --shadow-prefix $SHADOWDIR --local-cores 16 --jobs 16 --default-resources gpu=0 mem=4 \"qsub_args=''\" --cluster 'sbatch -J {name} --mem={resources.mem}G --gres=gpu:{resources.gpu} -c {threads} {resources.qsub_args} -e slurm-%x.e%A -o slurm-%x.o%A'"
else
  # Computerome thin nodes
  #SNAKE_PARAMS="--use-conda --restart-times 0 --keep-going --latency-wait 60 --conda-prefix $CONDA_DIR --shadow-prefix $SHADOWDIR --jobs 16 --cores 40 --resources mem=188"
  # Default
  SNAKE_PARAMS="--use-conda --restart-times 0 --keep-going --latency-wait 60 --conda-prefix $CONDA_DIR --shadow-prefix $SHADOWDIR --jobs 16 --cores 96 --resources mem=700"
fi

# Set code directory
CODE_DIR=$MINTO_DIR

# Download dependencies

echo ""
echo "------------"
echo "DEPENDENCIES"
echo "------------"

cat $MINTO_DIR/testing/dependencies.yaml.in | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > dependencies.yaml
echo "enable_GTDB: no" >> dependencies.yaml
cmd="snakemake --snakefile $CODE_DIR/smk/dependencies.smk --configfile dependencies.yaml $SNAKE_PARAMS >& dependencies.log"
profile_command "$cmd"

# Download raw data

echo ""
echo "-------------"
echo "TUTORIAL DATA"
echo "-------------"

if [ ! -d "IBD_tutorial_raw" ]; then
  echo -n "Downloading tutorial data: "
  wget --quiet https://zenodo.org/record/8320216/files/IBD_tutorial_raw_v2.0.0.tar.gz
  tar xfz IBD_tutorial_raw_v2.0.0.tar.gz
  echo OK
fi
echo ""

# Symlink samples for testing bulk mode

if [[ "$1" = "--bulk" ]]; then
  mkdir -p IBD_tutorial_bulk/metaG
  for SAMPLE in CD136 CD140 CD237 CD240;
  do
    find $PWD/IBD_tutorial_raw/metaG/${SAMPLE} -type f -exec ln -s {} IBD_tutorial_bulk/metaG/ \;
  done
fi

# Extract ref-genome

tar xfz $MINTO_DIR/tutorial/genomes.tar.gz

# Extract gene-catalog

tar xfz $MINTO_DIR/tutorial/gene_catalog.tar.gz

# Get data
mkdir -p IBD_tutorial
cd IBD_tutorial
cp $MINTO_DIR/tutorial/metadata/tutorial_metadata.txt .
cp $MINTO_DIR/tutorial/build_hg18_subset.fna .

# Run metaG and metaT steps

echo "---------------"
echo "DATA PROCESSING"
echo "---------------"

OMICS="metaG"
for OMICS in metaG metaT; do
  echo ""
  echo "------------------"
  echo "Processing $OMICS:"
  echo "------------------"
  mkdir -p $OMICS
  cd $OMICS

  COMMAND_LOG="commands_${OMICS}.txt"

  if [[ "$1" == "--bulk" && "$OMICS" == "metaG" ]]; then
      echo -n "QC_1 bulk: "
      cat $MINTO_DIR/testing/QC_1.bulk.yaml.in | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__TEST_DIR__>@$TEST_DIR@;s@<__OMICS__>@$OMICS@;" > QC_1.bulk.yaml
      cmd="snakemake --snakefile $CODE_DIR/smk/QC_1.smk --configfile QC_1.bulk.yaml $SNAKE_PARAMS >& QC_1_bulk.log"
      profile_command "$cmd"
  fi

  echo -n "QC_1: "
  cat $MINTO_DIR/testing/QC_1.yaml.in | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__TEST_DIR__>@$TEST_DIR@;s@<__OMICS__>@$OMICS@;" > QC_1.yaml
  cmd="snakemake --snakefile $CODE_DIR/smk/QC_1.smk --configfile QC_1.yaml $SNAKE_PARAMS >& QC_1.log"
  profile_command "$cmd"

  patch QC_2.yaml $CODE_DIR/testing/QC_2.patch -o - | sed "s@<__MINTO_DIR__>@$MINTO_DIR@;s@<__TEST_DIR__>@$TEST_DIR@" > QC_2.yaml.fixed
  echo -n "QC_2: "
  cmd="snakemake --snakefile $CODE_DIR/smk/QC_2.smk --configfile QC_2.yaml.fixed $SNAKE_PARAMS >& QC_2.log"
  profile_command "$cmd"

  echo -n "ASSEMBLY: "
  perl -pe "s/enable_COASSEMBLY: no/enable_COASSEMBLY: yes/; s/^# Contig-depth: bwa/EXCLUDE_ASSEMBLY_TYPES:\n - illumina_coas\n\n# Contig-depth: bwa/" < assembly.yaml > assembly.yaml.fixed
  cmd="snakemake --snakefile $CODE_DIR/smk/assembly.smk --configfile assembly.yaml.fixed $SNAKE_PARAMS >& assembly.log"
  profile_command "$cmd"

  echo -n "BINNING_PREP: "
  cmd="snakemake --snakefile $CODE_DIR/smk/binning_preparation.smk --configfile assembly.yaml.fixed $SNAKE_PARAMS >& binning_prep.log"
  profile_command "$cmd"

  echo -n "BINNING: "
  cmd="snakemake --snakefile $CODE_DIR/smk/mags_generation.smk --configfile mags_generation.yaml $SNAKE_PARAMS >& mags.log"
  profile_command "$cmd"

  echo -n "GENE_ANNOTATION - MAG: "
  sed "s/TAXONOMY_NAME: phylophlan,gtdb/TAXONOMY_NAME: phylophlan/" mapping.yaml > mapping.yaml.mag
  cmd="snakemake --snakefile $CODE_DIR/smk/gene_annotation.smk --configfile mapping.yaml.mag $SNAKE_PARAMS >& annotation.log"
  profile_command "$cmd"
  echo -n "GENE_ABUNDANCE - MAG: "
  cmd="snakemake --snakefile $CODE_DIR/smk/gene_abundance.smk --configfile mapping.yaml.mag $SNAKE_PARAMS >& abundance.log"
  profile_command "$cmd"

  sed "s@MINTO_MODE: MAG@MINTO_MODE: refgenome@; s@PATH_reference:@PATH_reference: $TEST_DIR/genomes@;" mapping.yaml.mag > mapping.yaml.refgenome
  echo -n "GENE_ANNOTATION - refgenome: "
  cmd="snakemake --snakefile $CODE_DIR/smk/gene_annotation.smk --configfile mapping.yaml.refgenome $SNAKE_PARAMS >& annotation.refgenome.log"
  profile_command "$cmd"
  echo -n "GENE_ABUNDANCE - refgenome: "
  cmd="snakemake --snakefile $CODE_DIR/smk/gene_abundance.smk --configfile mapping.yaml.refgenome $SNAKE_PARAMS >& abundance.refgenome.log"
  profile_command "$cmd"

  echo -n "GENE_ABUNDANCE - gene catalog: "
  sed "s@MINTO_MODE: MAG@MINTO_MODE: catalog@; s@PATH_reference:@PATH_reference: $TEST_DIR/gene_catalog@; s@NAME_reference:@NAME_reference: gene_catalog.fna@; s@abundance_normalization: TPM,MG@abundance_normalization: TPM@" mapping.yaml.mag > mapping.yaml.catalog
  cmd="snakemake --snakefile $CODE_DIR/smk/gene_abundance.smk --configfile mapping.yaml.catalog $SNAKE_PARAMS >& abundance.catalog.log"
  profile_command "$cmd"

  cd ..
done

# Run integration

echo ""
echo ""
echo "----------------"
echo "DATA INTEGRATION"
echo "----------------"

for OMICS in metaG_metaT metaG metaT; do

  echo ""
  echo "------------------"
  echo "Processing $OMICS:"
  echo "------------------"

  COMMAND_LOG="commands_integration_${OMICS}.txt"

  sed "s/omics: metaG_metaT/omics: $OMICS/; s/MINTO_MODE: .*/MINTO_MODE: MAG/" data_integration.yaml > data_integration.yaml.MG.$OMICS
  echo -n "MODE - MAG, MG: "
  mkdir -p $OMICS.MAG.MG
  cmd="cd $OMICS.MAG.MG && snakemake --snakefile $CODE_DIR/smk/data_integration.smk --configfile ../data_integration.yaml.MG.$OMICS $SNAKE_PARAMS >& integration.MAG.MG.$OMICS.log"
  profile_command "$cmd"
  sed "s/abundance_normalization: MG/abundance_normalization: TPM/" data_integration.yaml.MG.$OMICS > data_integration.yaml.TPM.$OMICS
  echo -n "MODE - MAG, TPM: "
  mkdir -p $OMICS.MAG.TPM
  cmd="cd $OMICS.MAG.TPM && snakemake --snakefile $CODE_DIR/smk/data_integration.smk --configfile ../data_integration.yaml.TPM.$OMICS $SNAKE_PARAMS >& integration.MAG.TPM.$OMICS.log"
  profile_command "$cmd"

  echo -n "MODE - refgenome, MG: "
  sed "s/MINTO_MODE: .*/MINTO_MODE: refgenome/" data_integration.yaml.MG.$OMICS > data_integration.yaml.refgenome.MG.$OMICS
  mkdir -p $OMICS.refgenome.MG
  cmd="cd $OMICS.refgenome.MG && snakemake --snakefile $CODE_DIR/smk/data_integration.smk --configfile ../data_integration.yaml.refgenome.MG.$OMICS $SNAKE_PARAMS >& integration.refgenome.MG.$OMICS.log"
  profile_command "$cmd"
  sed "s/abundance_normalization: MG/abundance_normalization: TPM/" data_integration.yaml.refgenome.MG.$OMICS > data_integration.yaml.refgenome.TPM.$OMICS
  echo -n "MODE - refgenome, TPM: "
  mkdir -p $OMICS.refgenome.TPM
  cmd="cd $OMICS.refgenome.TPM && snakemake --snakefile $CODE_DIR/smk/data_integration.smk --configfile ../data_integration.yaml.refgenome.TPM.$OMICS $SNAKE_PARAMS >& integration.refgenome.TPM.$OMICS.log"
  profile_command "$cmd"

  echo -n "MODE - gene-catalog, TPM: "
  sed "s/MINTO_MODE: .*/MINTO_MODE: catalog/; s@ANNOTATION_file:@ANNOTATION_file: $TEST_DIR/gene_catalog/gene_catalog.annotations.tsv@" data_integration.yaml.TPM.$OMICS > data_integration.yaml.catalog.TPM.$OMICS
  mkdir -p $OMICS.catalog.TPM
  cmd="cd $OMICS.catalog.TPM && snakemake --snakefile $CODE_DIR/smk/data_integration.smk --configfile ../data_integration.yaml.catalog.TPM.$OMICS $SNAKE_PARAMS >& integration.catalog.TPM.$OMICS.log"
  profile_command "$cmd"

done
