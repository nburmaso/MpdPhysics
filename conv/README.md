## Converting reconstruction into ROOT trees for analysis

### Prerequisites:

1) Install MpdRoot framework. It can be obtained from the official [website](http://mpd.jinr.ru/howto-install/) of the MPD collaboration.

2) Run event generation and reconstruction. See guides in `gen` and `reco`

### Step-by-step guide:

- Create a separate directory for your simulations
- Copy the following ROOT macros from `reco/urqmd/macros/`:
    * `tree.C`: collects events and tracks relevant to the analysis and filters them using kinematic cuts and
       cuts based on detector signals.
    * `conv.C`: wrapper for `tree.C` compilation
- If you are running simulation on a cluster with batch system, copy the following shell scripts from `gen/urqmd/macros/` to your simulation directory:
    * `submit.sh`: a script that can be used to submit a number of jobs at once
    * `srun.sh`: a script that submits each simulation job
    * `conv.sh`: a script that runs simulation for each subjob
- If you are running on a local PC, copy `conv_local.sh` to your simulation directory.
- Open `conv.sh` (or `conv_local.sh`) and set variables:
    * `MPDROOT_DIR`: MpdRoot installation path
    * `INDIR`: path to reconstruction results
    * `OUTDIR`: path to your results directory
- If running on a cluster: run `submit.sh`
- If running locally: run `conv_local.sh`

**Important**: 
- Check shell scripts for additional comments and modify scripts for your needs
- Macros used in this guide are very study-specific and were used for analysis of photon pairs. Use it as example for your analysis