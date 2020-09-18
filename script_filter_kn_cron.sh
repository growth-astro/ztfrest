#!/bin/tcsh

setenv LD_LIBRARY_PATH /data/ia/anaconda3/lib
setenv LD_RUN_PATH /data/ia/anaconda3/bin

# Set path for ZTF forced photometry
setenv PYTHONPATH /data/ia/
setenv PYTHONPATH /data/ia/ForcePhotZTF:${PYTHONPATH}

cd /data/ia/ztfrest/

echo "Time:"
echo `date`

# Run the pipeline
/data/ia/anaconda3/bin/python /data/ia/ztfrest/filter_kn.py --doWriteDb --doForcePhot --doCLU

# Send a slack message about the completion of the data processing
/data/ia/anaconda3/bin/python /data/ia/ztfrest/slack_bot_processing.py --channel test
/data/ia/anaconda3/bin/python /data/ia/ztfrest/slack_bot_processing.py --channel caltech
/data/ia/anaconda3/bin/python /data/ia/ztfrest/slack_bot_processing.py --channel partnership

echo "Finished cron job"
echo "Time:"
echo `date`
