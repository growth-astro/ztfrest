#!/bin/bash

# Runs ztfrest without slack processing.

# Run the pipeline (first half)
source /scr2/ztfrest/anaconda3/bin/conda activate ztfrest2
export PATH=/scr2/ztfrest/anaconda3/envs/ztfrest2/bin/:$PATH

cd /scr2/ztfrest/ZTF/sarah_test
python filter_kn.py --doWriteDb --doForcePhot --doAuxFp --doCLU

# trigger ForcePhotZTF which only works in an older python env
# writes this data to the db for later access
source /scr2/ztfrest/anaconda3/bin/conda activate ztfrest
export PATH=/scr2/ztfrest/anaconda3/envs/ztfrest/bin/:$PATH

python trigger_forcephotztf.py --doWriteDb

# run the pipeline (second half) - selection based on FPZTF / stacked FP
source /scr2/ztfrest/anaconda3/bin/conda activate ztfrest2
export PATH=/scr2/ztfrest/anaconda3/envs/ztfrest2/bin/:$PATH

python final_select_fpztf.py --doWriteDb

source /scr2/ztfrest/anaconda3/bin/conda deactivate

echo "Finished."
echo "Time:"
echo `date`
