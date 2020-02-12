# Download miniconda version
You can find it here - https://docs.conda.io/en/latest/miniconda.html#linux-installers

# Add the new miniconda file to your path
`export PATH="/home/ubuntu/miniconda3/bin/:$PATH"`

# Create r environment
`conda create -n r anaconda`

# Activate r environment
`conda activate r`

# Problems above
If you are having issues activating your R environment try the following
`source  /home/ubuntu/miniconda3/etc/profile.d/conda.sh`
Solution was found here - https://github.com/conda/conda/issues/7980

# R is installed
Once you activate your r environment do `which R` and you will see that you have R installed. But what if you want to change the R version?

# Check what R versions are available
`conda search r-base`
Then choose the verison you want:
`conda install -c r r=3.6.0`

