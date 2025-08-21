#!/bin/bash
#SBATCH --job-name=setup_tensorqtl
#SBATCH --output=logs/setup_tensorqtl.%j.out
#SBATCH --error=logs/setup_tensorqtl.%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Load conda if not already loaded
module load Anaconda3

# Create and activate a clean environment (Python 3.11 for tensorqtl)
ENV_NAME="tensorqtl_p3.11"
conda create -y -n $ENV_NAME python=3.11
source activate $ENV_NAME

# Install tensorqtl and its dependencies
pip install --upgrade pip

# Base installs
pip install tensorqtl pandas_plink
pip install 'rpy2<=3.5.12'

# Optional but recommended upgrades for conflicts Cindy saw
pip install scipy -U
pip install matplotlib -U
pip install multiqc -U
pip install requests -U
pip install importlib-metadata -U

# Ensure R/qvalue is available
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y r-essentials
conda install -c bioconda bioconductor-qvalue

# Fix for R shared libraries (specific to your HPC R build)
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH

# Test installation
echo "Testing tensorqtl..."
python -c "import tensorqtl; import pandas_plink; print('tensorqtl OK')"

# Export environment for reproducibility
conda env export > /wsu/home/ic/ic10/ic1032/tensorQTL/tensorqtl.environment.yaml
