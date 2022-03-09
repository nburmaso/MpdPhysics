## Photon pairs analysis with MPD

This is a collection of macros used for analysis of photon pairs in MPD

### List of analyses

- `photons_selection`:
    * `analysis.C`: A macro for selection of photon pairs using a variety of kinematic and detector cuts used. Requires reconstructed events stored in trees as input (check `reco` and `conv` for instructions).
    * Contains an example output for Bi-Bi simulation at ![](https://latex.codecogs.com/svg.image?\sqrt{s_{\mathrm{NN}}}=9.2&space;\mathrm{\~GeV}).

- `pions_yield`:
    * `draw_expected_yields.C`: A macro for estimating neutral pion yield. Requires results from `photons_selection/analysis.C` as input.
    * Contains an example of estimated neutral pion yields for Bi-Bi simulation at ![](https://latex.codecogs.com/svg.image?\sqrt{s_{\mathrm{NN}}}=9.2&space;\mathrm{\~GeV}).

- `pions_peak`:
    * `draw_pi_peak.C`: A macro for drawing diphoton mass spectrum and estimating parameters of neutral pion peak. Requires results from `photons_selection/analysis.C` as input.
    * Contains an example of estimated neutral pion peak for Bi-Bi simulation at ![](https://latex.codecogs.com/svg.image?\sqrt{s_{\mathrm{NN}}}=9.2&space;\mathrm{\~GeV}).
- `photon_reco_eff`:
    * `reconstruction_efficiency.C`: A macro for estimating photon reconstruction efficiency. Requires reconstructed events stored in trees as input (check `reco` and `conv` for instructions).
    * `draw_reconstruction_efficiency.C`: Visualization macro, takes resulting root-file from `reconstruction_efficiency.C` as input
    * Contains example estimates obtained using box simulations of photons with MPD (see `reco/box`).

- `tracks_par_resolution`:
    * `calc_resolution.C`: A macro that can be used to collect track parameters resolution for charged primary particles. Requires reconstructed events stored in trees as input (check `reco` and `conv` for instructions). An example macro for collecting events including additional primary particles can be found in `tree_prim.C`.
    * `draw_resolution.C`: Example macro for drawing parameter resolutions in dependence of ![](https://latex.codecogs.com/svg.image?p_{\mathrm&space;T}).
    * Contains example outputs for DCA and ![](https://latex.codecogs.com/svg.image?p_{\mathrm&space;T}).