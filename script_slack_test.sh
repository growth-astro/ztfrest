#!/bin/tcsh

# Run ztfrest with slack processing.

# set up python env
setenv LD_LIBRARY_PATH /data/ia/anaconda3/envs/ztfrest2025/lib
setenv LD_RUN_PATH /data/ia/anaconda3/envs/ztfrest2025/bin

# Run the pipeline
cd /data/ia/ztfrest_sarah/

echo "Starting pipeline..."
echo `date`

/data/ia/anaconda3/envs/ztfrest2025/bin/python /data/ia/ztfrest_sarah/filter_kn.py --doWriteDb --doAuxFp --doCLU

# Run ForcePhotZTF

echo "Working on ForcePhotZTF..."

cd /data/ia/ztfrest/

# set path for ForcePhotZTF package
setenv LD_LIBRARY_PATH /data/ia/anaconda3/lib
setenv LD_RUN_PATH /data/ia/anaconda3/bin

setenv PYTHONPATH /data/ia
setenv PYTHONPATH /data/ia/ForcePhotZTF:${PYTHONPATH}

/data/ia/anaconda3/bin/python /data/ia/ztfrest_sarah/trigger_forcephotztf.py --doWriteDb

# Select based on ForcePhotZTF and stacked forced photometry
setenv LD_LIBRARY_PATH /data/ia/anaconda3/envs/ztfrest2025/lib
setenv LD_RUN_PATH /data/ia/anaconda3/envs/ztfrest2025/bin

echo "Filtering ForcePhot..."
cd /data/ia/ztfrest_sarah

/data/ia/anaconda3/envs/ztfrest2025/bin/python /data/ia/ztfrest_sarah/final_select_fpztf.py --doWriteDb

# slack processing: TEST CHANNEL
echo "Sending to slack..."

/data/ia/anaconda3/envs/ztfrest2025/bin/python /data/ia/ztfrest_sarah/slack_bot_2025.py --channel test2025 -d --usenewdb

echo "Finished."
echo "Time:"
echo `date`
