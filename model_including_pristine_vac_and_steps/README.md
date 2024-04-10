`tpd-simulator`
--------------

# Installation

Create a new virtual environment following the instructions in the [Python documentation](https://docs.python.org/3/library/venv.html) or a more accessible guide such as [this one](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/).

Clone this repository from GitHub by clicking the green "Code" button and copying the URL. Then, in the terminal, navigate to the directory where you want to clone the repository and run the following command:

    git clone <link-from-url>

Activate the virtual environment and install the required packages:

```
    pip install -r requirements.txt
```

# Usage

Run the simulation from the model directory

```
    python run_model.py
```

You will have to set the following keys in the file `inputs.yaml`:

* `beta`: The heating rate
* `theta_sat`: The saturation coverage
* `grid_size`: The size of the grid for the numerical `mpmath.odefun` integration
* `Tmin`: The minimum temperature in the integration range
* `Tmax`: The maximum temperature in the integration range
* `mpmath_precision`: The precision of the `mpmath` integration
* `b`: Adsorbate-adsorbate interaction parameter 1--1
* `a`: Adsorbate-adsorbate interaction parameter 0--1
* `H0`: Enthalpy of desorption for site 0 in eV
* `H1`: Enthalpy of desorption for site 1 in eV
* `theta_change`: The coverage at which the adsorbate-adsorbate are "turned on" 
* `output_dir`: The directory where the output files will be saved

The enthalpy of desorption is computed as: 
    
        1.52  if  theta <= theta_change
        H0 - a * (theta) if  theta_change < theta <= 0.45
        H1 - b * (theta - 0.45) if theta > 0.45
        H0 otherwise

The free energy of desorption is computed as:

        Gx = Hx - T * S_config

where 

        S_config = k_B * ln(theta / (1 - theta))

For more information, look at the [paper here](https://pubs.rsc.org/en/content/articlelanding/2021/cp/d1cp01992a). The solver routine is not published yet, so please do not make the repository public.
