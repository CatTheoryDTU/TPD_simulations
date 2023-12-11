import logging

import numpy as np
import mpmath as mp

from typing import List, Tuple, Dict, Union, Callable

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


class SecondOrderTPD:
    def __init__(
        self,
        beta: float,
        theta_sat: float,
        grid_size: int = 1000,
        Tmin: float = 100.0,
        Tmax: float = 1000.0,
    ):
        """Solve the ordinary differential equation
        for the second order rate equation to get
        theta(T). Once we have theta(T), we can determine
        dtheta/dt (which is the rate, and can be comapred to TPD.
        Meant to serve as a base class, subclass and implement
        a Ga function. Note that we will be using mpmath to perform
        all the arithmetic."""
        # Make the constants into mpmath floats
        self.beta = mp.mpf(f"{beta}")
        self.theta_sat = mp.mpf(f"{theta_sat}")
        self.grid_size = mp.mpf(f"{grid_size}")
        self.Tmin = mp.mpf(f"{Tmin}")
        self.Tmax = mp.mpf(f"{Tmax}")

        # Relevant constants
        self.kB = mp.mpf(8.617e-5)  # eV/K
        self.h = mp.mpf(4.135667662e-15)  # eV*s

        # Choose a fixed value of x_star = 1.
        self.x_star = mp.mpf("1.0")

    def get_dtheta_dT(
        self,
        T: Union[float, List],
        theta: Union[float, List],
    ) -> Union[float, List]:
        """Implement here the differential equation
        right hand side to get the dtheta/dT."""
        logger.debug(f"theta = {mp.nstr(theta, 5)}")
        logger.debug(f"T = {mp.nstr(T, 5)}")

        Ga = self.Ga(theta=theta, T=T)
        logger.debug(f"Ga = {mp.nstr(Ga, 5)}")

        dtheta_dT = -mp.mpf("2") * self.kB * T / self.h
        dtheta_dT *= mp.exp(-Ga / (self.kB * T))
        dtheta_dT *= theta**2
        dtheta_dT /= self.beta

        logger.debug(f"dtheta_dT = {mp.nstr(dtheta_dT, 5)}")
        
        return dtheta_dT

    def get_dx_dT(
        self, T: Union[float, List[float]], x: Union[float, List[float]]
    ) -> Union[float, List[float]]:
        """Implement here the differential equation
        right hand side to get the dx/dT."""

        # Define the coverage based on the number of sites.
        # This conversion implies that theta is always between 0 and 1.
        theta = x**2 / (self.x_star**2 + x**2)

        # Get dtheta_dT based on this coverage
        dtheta_dT = self.get_dtheta_dT(T=T, theta=theta)

        # We need to perform chain rule to get dx/dT
        dtheta_dx = mp.mpf("2.0") * x / (self.x_star**2 + x**2)
        dtheta_dx -= mp.mpf("2.0") * x**3 / (self.x_star**2 + x**2) ** 2

        dx_dT = dtheta_dT / dtheta_dx

        return dx_dT

    def Ha(
        self, theta: Union[float, List[float]], T: Union[float, List[float]]
    ) -> Union[float, List]:
        """Implement here the enthalpy as
        a function of coverage."""
        raise NotImplementedError("Implement this function in the child class.")

    def Ga(
        self, theta: Union[float, List[float]], T: Union[float, List[float]]
    ) -> Union[float, List]:
        """Implement here the activation energy as
        a function of coverage."""
        raise NotImplementedError("Implement this function in the child class.")

    def get_theta(
        self, get_additional_info: bool = True
    ) -> Union[Tuple[List, List], Tuple[List, List, Dict]]:
        """For a supplied Temperature range, find
        the corresponding theta(T) values by solving
        the ODE."""
        # Integrate the dtheta_dT function using mpmath integration
        # theta_T = mp.odefun(self.dtheta_dT, self.Tmin, self.theta_sat)
        # Saturation x value
        x_sat = mp.sqrt(self.theta_sat * self.x_star**2 / (1 - self.theta_sat))
        x_T = mp.odefun(self.get_dx_dT, self.Tmin, x_sat)

        Trange = mp.linspace(self.Tmin, self.Tmax, self.grid_size)

        # Store the coverage
        theta = []

        # Store some additional information for the purposes of debugging
        dtheta_dT = []
        Ga = []
        Ha = []

        for _T in Trange:
            # Get the coverage based on the number of sites.
            _theta = x_T(_T) ** 2 / (self.x_star**2 + x_T(_T) ** 2)
            theta.append(_theta)

            # Get the derivative of theta with respect to T
            _dtheta_dT = self.get_dtheta_dT(T=_T, theta=_theta)
            dtheta_dT.append(_dtheta_dT)

            # Get the desorption enthapy
            _Ha = self.Ha(theta=_theta, T=_T)
            Ha.append(_Ha)

            # Get the desorption free energy
            _Ga = self.Ga(theta=_theta, T=_T)
            Ga.append(_Ga)

        additional_info = {
            "dtheta_dT": dtheta_dT,
            "Ga": Ga,
            "Ha": Ha,
        }

        if get_additional_info:
            return theta, Trange, additional_info
        else:
            return theta, Trange
