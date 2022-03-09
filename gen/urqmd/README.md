## Running URQMD simulation

### Prerequisites:

Install URQMD event generator. It can be obtained from the official [website](https://urqmd.org/).

### Step-by-step guide:

- Create a separate directory for your simulations
- If you are running simulation on a cluster with batch system, copy the following shell scripts from `gen/urqmd/macros/` to your simulation directory:
    * `submit.sh`: a script that can be used to submit a number of jobs at once
    * `srun.sh`: a script that submits each simulation job
    * `urqmd.sh`: a script that runs simulation for each subjob
- If you are running on a local PC, copy `urqmd_local.sh` to your simulation directory.
- Create `inputfile` for URQMD and fill in desired parameters. Some examples are available in `gen/urqmd/inputfiles` (Note that generator seed that can be set using `rsd` parameter is set in `urqmd.sh` or `urqmd_local.sh` randomly in our scripts). In case of any uncertainty, check URQMD documentation. 
- Open `urqmd.sh` (or `urqmd.sh`) and set variables:
    * `URQMD_DIR`: URQMD installation path
    * `OUTDIR`: path to your simulation directory
    * `inputfile`: path to your inputfile
- If running on a cluster: run `submit.sh`
- If running locally: run `urqmd_local.sh`

**Important: check shell scripts for additional comments and modify scripts for your needs** 