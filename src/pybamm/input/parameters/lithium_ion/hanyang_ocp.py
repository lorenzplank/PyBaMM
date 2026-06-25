"""
Tabulated open-circuit-potential (OCP) curves digitized from
https://real.hanyang.ac.kr/materials

Each function returns a :class:`pybamm.Interpolant` of OCP [V] vs Li/Li+ as a
function of electrode stoichiometry, matching the tabulated-OCP convention used
by other PyBaMM parameter sets (e.g. ``Ai2020``). They are provided as a
*standalone resource* — the built-in parameter sets keep their own analytic
OCPs. To use one, override the relevant OCP entry, e.g.::

    import pybamm
    from pybamm.input.parameters.lithium_ion.hanyang_ocp import (
        nmc811_ocp_Hanyang, gr_ocp_Hanyang,
    )

    pv = pybamm.ParameterValues("Chen2020")
    pv.update({
        "Positive electrode OCP [V]": nmc811_ocp_Hanyang,
        "Negative electrode OCP [V]": gr_ocp_Hanyang,
    })

Available cathodes: LCO, LFP, NCA, NMC111, NMC532, NMC622, NMC811, LMO, LNMO.
Available anodes:   Gr, LTO, Si.

.. note::
    The interpolants are only valid over each material's digitized
    stoichiometry window (see the CSVs in ``data/``); PyBaMM will raise on
    extrapolation outside that range.
"""
import os

import pybamm

path, _ = os.path.split(os.path.abspath(__file__))


_lco_data = pybamm.parameters.process_1D_data("lco_ocp_Hanyang.csv", path=path)


def lco_ocp_Hanyang(sto):
    """Tabulated LCO (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lco_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_lfp_data = pybamm.parameters.process_1D_data("lfp_ocp_Hanyang.csv", path=path)


def lfp_ocp_Hanyang(sto):
    """Tabulated LFP (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lfp_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_nca_data = pybamm.parameters.process_1D_data("nca_ocp_Hanyang.csv", path=path)


def nca_ocp_Hanyang(sto):
    """Tabulated NCA (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nca_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_nmc111_data = pybamm.parameters.process_1D_data("nmc111_ocp_Hanyang.csv", path=path)


def nmc111_ocp_Hanyang(sto):
    """Tabulated NMC111 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc111_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_nmc532_data = pybamm.parameters.process_1D_data("nmc532_ocp_Hanyang.csv", path=path)


def nmc532_ocp_Hanyang(sto):
    """Tabulated NMC532 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc532_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_nmc622_data = pybamm.parameters.process_1D_data("nmc622_ocp_Hanyang.csv", path=path)


def nmc622_ocp_Hanyang(sto):
    """Tabulated NMC622 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc622_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_nmc811_data = pybamm.parameters.process_1D_data("nmc811_ocp_Hanyang.csv", path=path)


def nmc811_ocp_Hanyang(sto):
    """Tabulated NMC811 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc811_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_lmo_data = pybamm.parameters.process_1D_data("lmo_ocp_Hanyang.csv", path=path)


def lmo_ocp_Hanyang(sto):
    """Tabulated LMO (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lmo_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_lnmo_data = pybamm.parameters.process_1D_data("lnmo_ocp_Hanyang.csv", path=path)


def lnmo_ocp_Hanyang(sto):
    """Tabulated LNMO (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lnmo_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_gr_data = pybamm.parameters.process_1D_data("gr_ocp_Hanyang.csv", path=path)


def gr_ocp_Hanyang(sto):
    """Tabulated Gr (anode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _gr_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_lto_data = pybamm.parameters.process_1D_data("lto_ocp_Hanyang.csv", path=path)


def lto_ocp_Hanyang(sto):
    """Tabulated LTO (anode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lto_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


_si_data = pybamm.parameters.process_1D_data("si_ocp_Hanyang.csv", path=path)


def si_ocp_Hanyang(sto):
    """Tabulated Si (anode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _si_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")
