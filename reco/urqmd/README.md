## Running full simulation and event reconstruction with MpdRoot 

### Prerequisites:

Install MpdRoot framework. It can be obtained from the official [website](http://mpd.jinr.ru/howto-install/) of the MPD collaboration.

### Step-by-step guide:

- Create a separate directory for your simulations
- Copy the following ROOT macros from `reco/urqmd/macros/`:
    * `geometry.C`: used to define geometry that will be used in simulation
    * `sim.C`: runs Geant simulation in MpdRoot
    * `rec.C`: runs track reconstruction using Kalman filter technique
- If you are running simulation on a cluster with batch system, copy the following shell scripts from `gen/urqmd/macros/` to your simulation directory:
    * `submit.sh`: a script that can be used to submit a number of jobs at once
    * `srun.sh`: a script that submits each simulation job
    * `run.sh`: a script that runs simulation for each subjob
- If you are running on a local PC, copy `run_local.sh` to your simulation directory.
- Open `run.sh` (or `run_local.sh`) and set variables:
    * `MPDROOT_DIR`: MpdRoot installation path
    * `URQMD_DIR`: path to URQMD-generated events
    * `OUTDIR`: path to your simulation directory
- If running on a cluster: run `submit.sh`
- If running locally: run `run_local.sh`

**Important: check shell scripts for additional comments and modify scripts for your needs**