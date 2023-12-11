import os
from typing import List, Tuple, Dict, Union, Callable
import logging
import yaml

import numpy as np
import mpmath as mp

import matplotlib.pyplot as plt

from model import SecondOrderTPD

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

from plot_params import get_plot_params

get_plot_params()


class TPDOER(SecondOrderTPD):
    def __init__(
        self,
        beta: float,
        theta_sat: float,
        grid_size: int = 1000,
        Tmin: float = 100.0,
        Tmax: float = 1000.0,
        b: float = 0.5,
        theta_change: float = 0.5,
        H0: float = 0.0,
    ):
        super().__init__(beta, theta_sat, grid_size, Tmin, Tmax)

        # Make b and theta_change into mpmath floats
        self.b = mp.mpf(f"{b}")
        self.theta_change = mp.mpf(f"{theta_change}")
        self.H0 = mp.mpf(f"{H0}")

    def Ha(
        self,
        T: Union[float, List[float]],
        theta: Union[float, List[float]],
    ) -> Union[float, List]:
        """Activation energy for the _ADSORPTION_ reaction. Note that
        we would expect the enthalpy of the reaction, Ha to be positive
        for the desoprtion while the configurational entropy term
        depends on the value of theta."""

        if mp.mpf("0.0") <= theta <= self.theta_change:
            Ha = self.H0
        elif self.theta_change < theta < mp.mpmathify("1.0"):
            Ha = self.H0 - self.b * (theta - self.theta_change)
        else:
            logging.error(f"Theta is out of bounds: {theta}")
            raise ValueError("Theta must be between 0 and 1.")

        return Ha

    def Ga(
        self,
        T: Union[float, List[float]],
        theta: Union[float, List[float]],
    ) -> Union[float, List]:
        """Activation energy for the _DESORPTION_ reaction. Note that
        we would expect the enthalpy of the reaction, Ha to be positive
        for the desoprtion while the configurational entropy term
        depends on the value of theta."""

        # Get the enthalpy of desorption
        Ha = self.Ha(T=T, theta=theta)

        # The differential configurational entropy
        S_config = self.kB * mp.log(theta / (1 - theta))

        # The free energy of desorption
        Ga = Ha - T * S_config

        return Ga


if __name__ == "__main__":
    """Plot the rates as a function of coverage."""

    # Get the input parameters
    with open("inputs.yaml", "r") as f:
        inputs = yaml.safe_load(f)

    # Set mpmath precision
    mp.dps = inputs["mpmath_precision"]
    logger.info(f"mpmath precision set to {mp.dps}.")

    # Write the parameters to the log
    logger.info(f"beta = {inputs['beta']}")
    logger.info(f"theta_sat = {inputs['theta_sat']}")
    logger.info(f"grid_size = {inputs['grid_size']}")
    logger.info(f"Tmin = {inputs['Tmin']}")
    logger.info(f"Tmax = {inputs['Tmax']}")
    logger.info(f"b = {inputs['b']}")
    logger.info(f"theta_change = {inputs['theta_change']}")
    logger.info(f"H0 = {inputs['H0']}")
    logger.info(f"Model output will be saved to {inputs['output_dir']}.")

    # Create the output directory if it doesn't exist
    if not os.path.exists(inputs["output_dir"]):
        os.makedirs(inputs["output_dir"])

    # Set up the model
    model = TPDOER(
        beta=inputs["beta"],
        theta_sat=inputs["theta_sat"],
        grid_size=inputs["grid_size"],
        Tmin=inputs["Tmin"],
        Tmax=inputs["Tmax"],
        b=inputs["b"],
        theta_change=inputs["theta_change"],
        H0=inputs["H0"],
    )

    # Get the theta(T) values
    theta, Trange, additional_info = model.get_theta(get_additional_info=True)

    dtheta_dT_scaled = [i* (-1) for i in additional_info["dtheta_dT"]]

    # Plot the results
    fig, ax = plt.subplots(1, 3, figsize=(6, 2.5), constrained_layout=True)
    ax[0].plot(Trange, theta)
    ax[1].plot(theta, additional_info["Ga"], label="Ga")
    ax[1].plot(theta, additional_info["Ha"], label="Ha")
    ax[1].legend()
    ax[2].plot(Trange,dtheta_dT_scaled)

    ax[0].set_xlabel("Temperature (K)")
    ax[0].set_ylabel("Coverage")
    ax[1].set_xlabel("Coverage")
    ax[1].set_ylabel("Desorption energy (eV)")
    ax[2].set_xlabel("Temperature (K)")
    ax[2].set_ylabel(r"d$\theta$/dT")

    fig.savefig(f"{inputs['output_dir']}/model_output.png", dpi=300)
