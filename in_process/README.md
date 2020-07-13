# RELAX
Relax implements a semi-analytic Fourier-domain solver and equivalent body forces to compute quasi-static relaxation of stress perturbation.

Note: Please download RELAX software from the website:
http://geodynamics.org/cig/software/relax/
OR downlaod RELAX software from github:
git clone: --recursive https://github.com/geodynamics/relax.git

Usage:
1. unzip the file
2. make one template for the studied case. In the template you should configure the basic parameters for you studied case.
3. use generate_script_RELAX_V2.py to generate templates batchly.
4. go to /Relax/examples/project_name/PBS_script
5. run ./run_RELAX to sub jobs to PBS/LSF or other cluster  

