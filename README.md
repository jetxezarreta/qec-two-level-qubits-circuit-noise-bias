# qec-two-level-qubits-circuit-noise-bias
<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
    <li><a href="#attribution">Attribution</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

This repository includes the Stim circuits used to simulate the performance of the rotated XZZX code under a hybrid biased-depolarizing (HBD) circuit-level noise model as described in "Leveraging biased noise for more efficient quantum error correction at the circuit-level with two-level qubits" (add arxiv number). It includes the Stim circuits containing the circuit-level error models considered, including a generic HBD noise model (with bias-preserving and non-bias-preserving CZ gates), a HBD model with CNOT gates exhibiting residual biased noise and a compilation of the XZZX code using only CZ gates that preserve bias. The files containing the obtained results for the plots presented in the paper are also accessible. Finally, a script used to determine pauli transfer matrices of a generic CNOT gate under biased noise is also accessible.

The repository is organized as follows:
* circuits
  * CZcompilation_XZZX_surface_code_HybridBiasCLN: is a script to generate the circuits for the XZZX code compiled only using CZ gates preserving the bias.
  * XZZX_surface_code_depolarizing_CZ: is a script to generate the circuits for the XZZX code with CZ (and CNOT) gates implemented in a non-bias-preserving way.
  * XZZX_surface_code_HybridBiasCLN_ResidualCNOT: is a script to generate the circuits for the XZZX code with CZ gates implemented in a bias-preserving way and CNOT gates exhibiting residual bias.
  * XZZX_surface_code_HybridBiasCLN: is a script to generate the circuits for the XZZX code with CZ gates implemented in a bias-preserving way and CNOT gates introducing depolarizing noise (generic HBD circuit-level noise).
* QTIP:
  *  CNOTbias: script to get the pauli transfer matrices of a noisy CNOT gate using a generic pulse for two-level qubits.
*  simulation/XZZX_code:
      * The folders inside include the STIM output files for the simulations presented in the paper. Those have comprehensive names.
      * plot_thresholds: jupyter notebook for plotting the thresholds obtained using the results folders.
      * simulation_XZZX_code_depolarizing_CZ: script to run simulations of the XZZX code with CZ gates (and CNOTs) implemented in a non-bias-preserving manner.
      * simulation_XZZX_code_HBD_CZbiased: script to run simulations of the all the considered HBD models with bias-preserving CZ gates.
      * teraquop_footprint: jupyter notebook to plot the footprints required by the XZZX code.
      * plot_utils: misc plotting functions.
* visualization:
  * visualize_XZZX_surface_code: jupyter notebook for visualizing the extraction circuits used for the XZZX codes.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- CONTACT -->
## Contact
This repository has been developed by:

* _Paul Schnabl_ - [@pjschna](https://x.com/pjschna) - pjschnabl@gmail.com
* _Josu Etxezarreta Martinez_ - [@katutxakur](https://x.com/katutxakur) - [@katutxakur.bsky.social](https://bsky.app/profile/katutxakur.bsky.social) - jetxezarreta@unav.es



Project Link: [https://github.com/jetxezarreta/qec-two-level-qubits-circuit-noise-bias](https://github.com/jetxezarreta/qec-two-level-qubits-circuit-noise-bias)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

Thanks to the following amazing projects and webs for the help, tools and information! Do not forget to visit and star/like their work also!

* [StimCircuits](https://github.com/oscarhiggott/StimCircuits) - Oscar Higgott
* [PyMatching](https://github.com/oscarhiggott/PyMatching) - Oscar Higgott
* [Stim](https://github.com/quantumlib/Stim) - Craig Gidney
* [Sinter](https://pypi.org/project/sinter/) - Craig Gidney
* [QuTiP](https://qutip.org/)
* [Best-README-Template](https://github.com/othneildrew/Best-README-Template) - Othneil Drew

The figures from the visualize_XZZX_surface_code file have been borrowed from:
* [Suppressing quantum errors by scaling a surface code logical qubit](https://www.nature.com/articles/s41586-022-05434-1)
* [Decoding algorithms for Surface Codes](https://quantum-journal.org/papers/q-2024-10-10-1498/)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ATTRIBUTION -->
## Attribution
Please cite the following article if you use our circuits or our findings:
```
@article{CLNbiasedNoiseTLSqubits,
    author = "{Etxezarreta Martinez}, Josu and {Schnabl}, Paul and {Oliva del Moral}, Javier and {Dastbasteh}, Reza and {Crespo} Pedro M. and {Otxoa}, Ruben M.",
    title = "{Leveraging biased noise for more efficient quantum error correction at the circuit-level with two-level qubits}",
    journal = {arXiv},
    pages = {2505.17718},
    archivePrefix = "arXiv",
    primaryClass = "quant-ph",
    month = may,
    year = {2025},
    url ={https://arxiv.org/abs/2505.17718}
}
