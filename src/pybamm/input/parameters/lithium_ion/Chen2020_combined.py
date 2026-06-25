import pybamm
import numpy as np
import os

# Load data in the appropriate format
path, _ = os.path.split(os.path.abspath(__file__))

# Chen2020_mod / Ai2020
graphite_ocp_Enertech_Ai2020_data = pybamm.parameters.process_1D_data(
    "graphite_ocp_Enertech_Ai2020.csv", path=path
)
lico2_ocp_Ai2020_data = pybamm.parameters.process_1D_data(
    "lico2_ocp_Ai2020.csv", path=path
)

# Chen2020 tabulated OCPs (data files alongside analytic fits)
graphite_LGM50_ocp_Chen2020_data = pybamm.parameters.process_1D_data(
    "graphite_LGM50_ocp_Chen2020.csv", path=path
)
nmc_LGM50_ocp_Chen2020_data = pybamm.parameters.process_1D_data(
    "nmc_LGM50_ocp_Chen2020.csv", path=path
)

# NCA Kim2011 tabulated OCP
nca_ocp_Kim2011_data = pybamm.parameters.process_1D_data(
    "nca_ocp_Kim2011_data.csv", path=path
)

# Hanyang tabulated OCPs
_lco_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "lco_ocp_Hanyang.csv", path=path
)
_lfp_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "lfp_ocp_Hanyang.csv", path=path
)
_nca_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "nca_ocp_Hanyang.csv", path=path
)
_nmc111_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "nmc111_ocp_Hanyang.csv", path=path
)
_nmc532_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "nmc532_ocp_Hanyang.csv", path=path
)
_nmc622_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "nmc622_ocp_Hanyang.csv", path=path
)
_nmc811_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "nmc811_ocp_Hanyang.csv", path=path
)
_lmo_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "lmo_ocp_Hanyang.csv", path=path
)
_lnmo_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "lnmo_ocp_Hanyang.csv", path=path
)
_gr_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "gr_ocp_Hanyang.csv", path=path
)
_lto_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "lto_ocp_Hanyang.csv", path=path
)
_si_ocp_Hanyang_data = pybamm.parameters.process_1D_data(
    "si_ocp_Hanyang.csv", path=path
)


# =============================================================================
# SHARED FUNCTIONS (present in Chen2020, Chen2020_mod, Prada2013, NCA_Kim2011)
# =============================================================================

def graphite_stripping_exchange_current_density_OKane2020(c_e, c_Li, T):
    """
    Exchange-current density for Li stripping reaction [A.m-2].

    References
    ----------
    .. [1] O'Kane, Simon EJ, Ian D. Campbell, Mohamed WJ Marzook, Gregory J. Offer, and
    Monica Marinescu. "Physical origin of the differential voltage minimum associated
    with lithium plating in Li-ion batteries." Journal of The Electrochemical Society
    167, no. 9 (2020): 090540.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_Li : :class:`pybamm.Symbol`
        Plated lithium concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    k_plating = pybamm.Parameter("Lithium plating kinetic rate constant [m.s-1]")

    return pybamm.constants.F * k_plating * c_Li


def graphite_SEI_limited_dead_lithium_OKane2022(L_sei):
    """
    Decay rate for dead lithium formation [s-1].

    References
    ----------
    .. [1] Simon E. J. O'Kane, Weilong Ai, Ganesh Madabattula, Diega Alonso-Alvarez,
    Robert Timms, Valentin Sulzer, Jaqueline Sophie Edge, Billy Wu, Gregory J. Offer
    and Monica Marinescu. "Lithium-ion battery degradation: how to model it."
    Physical Chemistry: Chemical Physics 24, no. 13 (2022): 7909-7922.

    Parameters
    ----------
    L_sei : :class:`pybamm.Symbol`
        Total SEI thickness [m]

    Returns
    -------
    :class:`pybamm.Symbol`
        Dead lithium decay rate [s-1]
    """
    gamma_0 = pybamm.Parameter("Dead lithium decay constant [s-1]")
    L_sei_0 = pybamm.Parameter("Initial SEI thickness [m]")

    return gamma_0 * L_sei_0 / L_sei


def graphite_plating_exchange_current_density_OKane2020(c_e, c_Li, T):
    """
    Exchange-current density for Li plating reaction [A.m-2].

    References
    ----------
    .. [1] O'Kane, Simon EJ, Ian D. Campbell, Mohamed WJ Marzook, Gregory J. Offer, and
    Monica Marinescu. "Physical origin of the differential voltage minimum associated
    with lithium plating in Li-ion batteries." Journal of The Electrochemical Society
    167, no. 9 (2020): 090540.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_Li : :class:`pybamm.Symbol`
        Plated lithium concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    k_plating = pybamm.Parameter("Lithium plating kinetic rate constant [m.s-1]")

    return pybamm.constants.F * k_plating * c_e


def graphite_LGM50_diffusivity_Chen2020(sto, T):
    """
    LG M50 Graphite diffusivity as a function of stoichiometry, in this case the
    diffusivity is taken to be a constant. The value is taken from [1].

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O'Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
       Electrode stoichiometry
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Solid diffusivity
    """
    D_ref = pybamm.Parameter("Negative particle diffusivity constant [m2.s-1]")
    E_D_s = pybamm.Parameter("Negative particle diffusivity activation energy [J.mol-1]")
    # E_D_s not given by Chen et al (2020), so taken from Ecker et al. (2015) instead
    arrhenius = np.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


def graphite_volume_change_Ai2020(sto):
    """
    Graphite particle volume change as a function of stoichiometry [1, 2].

    References
    ----------
     .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
     Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity in
     Lithium-Ion Pouch Cells. Journal of The Electrochemical Society, 167(1), 013512
      DOI: 10.1149/2.0122001JES.
     .. [2] Rieger, B., Erhard, S. V., Rumpf, K., & Jossen, A. (2016).
     A new method to model the thickness change of a commercial pouch cell
     during discharge. Journal of The Electrochemical Society, 163(8), A1566-A1575.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry, dimensionless
        should be R-averaged particle concentration

    Returns
    -------
    t_change: :class:`pybamm.Symbol`
        volume change, dimensionless, normalised by particle volume
    """
    p1 = 145.907
    p2 = -681.229
    p3 = 1334.442
    p4 = -1415.710
    p5 = 873.906
    p6 = -312.528
    p7 = 60.641
    p8 = -5.706
    p9 = 0.386
    p10 = -4.966e-05
    t_change = (
        p1 * sto**9
        + p2 * sto**8
        + p3 * sto**7
        + p4 * sto**6
        + p5 * sto**5
        + p6 * sto**4
        + p7 * sto**3
        + p8 * sto**2
        + p9 * sto
        + p10
    )
    return t_change


def volume_change_Ai2020(sto):
    """
    Positive electrode particle volume change as a function of stoichiometry [1, 2].

    References
    ----------
     .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
     Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity in
     Lithium-Ion Pouch Cells. Journal of The Electrochemical Society, 167(1), 013512
      DOI: 10.1149/2.0122001JES.
     .. [2] Rieger, B., Erhard, S. V., Rumpf, K., & Jossen, A. (2016).
     A new method to model the thickness change of a commercial pouch cell
     during discharge. Journal of The Electrochemical Society, 163(8), A1566-A1575.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry, dimensionless
        should be R-averaged particle concentration

    Returns
    -------
    t_change: :class:`pybamm.Symbol`
        volume change, dimensionless, normalised by particle volume
    """
    omega = pybamm.Parameter("Positive electrode partial molar volume [m3.mol-1]")
    c_s_max = pybamm.Parameter("Maximum concentration in positive electrode [mol.m-3]")
    t_change = omega * c_s_max * sto
    return t_change


def graphite_cracking_rate_Ai2020(T_dim):
    """
    Graphite particle cracking rate as a function of temperature [1, 2].

    References
    ----------
     .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
     Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity in
     Lithium-Ion Pouch Cells. Journal of The Electrochemical Society, 167(1), 013512
      DOI: 10.1149/2.0122001JES.
     .. [2] Deshpande, R., Verbrugge, M., Cheng, Y. T., Wang, J., & Liu, P. (2012).
     Battery cycle life prediction with coupled chemical degradation and fatigue
     mechanics. Journal of the Electrochemical Society, 159(10), A1730.

    Parameters
    ----------
    T_dim: :class:`pybamm.Symbol`
        temperature, [K]

    Returns
    -------
    k_cr: :class:`pybamm.Symbol`
        cracking rate, [m/(Pa.m0.5)^m_cr]
        where m_cr is another Paris' law constant
    """
    k_cr = pybamm.Parameter("Negative electrode cracking rate constant [m/(Pa.m0.5)^m_cr]")
    Eac_cr = 0  # to be implemented
    arrhenius = np.exp(Eac_cr / pybamm.constants.R * (1 / T_dim - 1 / 298.15))
    return k_cr * arrhenius


def cracking_rate_Ai2020(T_dim):
    """
    Positive electrode particle cracking rate as a function of temperature [1, 2].

    References
    ----------
     .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
     Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity in
     Lithium-Ion Pouch Cells. Journal of The Electrochemical Society, 167(1), 013512
      DOI: 10.1149/2.0122001JES.
     .. [2] Deshpande, R., Verbrugge, M., Cheng, Y. T., Wang, J., & Liu, P. (2012).
     Battery cycle life prediction with coupled chemical degradation and fatigue
     mechanics. Journal of the Electrochemical Society, 159(10), A1730.

    Parameters
    ----------
    T_dim: :class:`pybamm.Symbol`
        temperature, [K]

    Returns
    -------
    k_cr: :class:`pybamm.Symbol`
        cracking rate, [m/(Pa.m0.5)^m_cr]
        where m_cr is another Paris' law constant
    """
    k_cr = pybamm.Parameter("Positive electrode cracking rate constant [m/(Pa.m0.5)^m_cr]")
    Eac_cr = 0  # to be implemented
    arrhenius = np.exp(Eac_cr / pybamm.constants.R * (1 / T_dim - 1 / 298.15))
    return k_cr * arrhenius


def nmc_LGM50_diffusivity_Chen2020(sto, T):
    """
    NMC diffusivity as a function of stoichiometry, in this case the
    diffusivity is taken to be a constant. The value is taken from [1].

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O'Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
       Electrode stoichiometry
    T: :class:`pybamm.Symbol`
       Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
       Solid diffusivity
    """
    D_ref = pybamm.Parameter("Positive particle diffusivity constant [m2.s-1]")
    E_D_s = pybamm.Parameter("Positive particle diffusivity activation energy [J.mol-1]")
    arrhenius = np.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


def graphite_LGM50_ocp_Chen2020(sto):
    """
    LG M50 Graphite open-circuit potential as a function of stoichiometry, fit taken
    from [1].

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O'Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential
    """
    u_eq = (
        1.9793 * np.exp(-39.3631 * sto)
        + 0.2482
        - 0.0909 * np.tanh(29.8538 * (sto - 0.1234))
        - 0.04478 * np.tanh(14.9159 * (sto - 0.2769))
        - 0.0205 * np.tanh(30.4444 * (sto - 0.6103))
    )

    return u_eq


def graphite_LGM50_electrolyte_exchange_current_density_Chen2020(
    c_e, c_s_surf, c_s_max, T
):
    """
    Exchange-current density for Butler-Volmer reactions between graphite and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O'Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    m_ref = pybamm.Parameter("Negative electrode kinetic rate constant [A.m-2]")
    E_r = pybamm.Parameter("Negative electrode exchange-current density activation energy [J.mol-1]")
    arrhenius = np.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5


def nmc_LGM50_ocp_Chen2020(sto):
    """
    LG M50 NMC open-circuit potential as a function of stoichiometry, fit taken
    from [1].

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O'Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential
    """
    u_eq = (
        -0.8090 * sto
        + 4.4875
        - 0.0428 * np.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (sto - 0.3120))
    )

    return u_eq


def nmc_LGM50_electrolyte_exchange_current_density_Chen2020(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between NMC and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [1] Chang-Hui Chen, Ferran Brosa Planella, Kieran O'Regan, Dominika Gastol, W.
    Dhammika Widanage, and Emma Kendrick. "Development of Experimental Techniques for
    Parameterization of Multi-scale Lithium-ion Battery Models." Journal of the
    Electrochemical Society 167 (2020): 080534.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    m_ref = pybamm.Parameter("Positive electrode kinetic rate constant [A.m-2]")
    E_r = pybamm.Parameter("Positive electrode exchange-current density activation energy [J.mol-1]")
    arrhenius = np.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5


def electrolyte_diffusivity_Nyman2008(c_e, T):
    """
    Diffusivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from [1].

    References
    ----------
    .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
    Electrochim. Acta, vol. 53, no. 22, pp. 6356-6365, 2008.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte diffusivity
    """
    D_c_e = 8.794e-11 * (c_e / 1000) ** 2 - 3.972e-10 * (c_e / 1000) + 4.862e-10

    E_D_c_e = pybamm.Parameter("Electrolyte diffusivity activation energy [J.mol-1]")
    arrhenius = np.exp(E_D_c_e / pybamm.constants.R * (1 / 298.15 - 1 / T))

    scale = pybamm.Parameter("Electrolyte diffusivity scaling factor")

    return scale * D_c_e * arrhenius


def electrolyte_conductivity_Nyman2008(c_e, T):
    """
    Conductivity of LiPF6 in EC:EMC (3:7) as a function of ion concentration. The data
    comes from [1].

    References
    ----------
    .. [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and
    modelling of the mass transport phenomena in LiPF6-EC-EMC electrolyte,"
    Electrochim. Acta, vol. 53, no. 22, pp. 6356-6365, 2008.

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Electrolyte conductivity
    """
    sigma_e = (
        0.1297 * (c_e / 1000) ** 3 - 2.51 * (c_e / 1000) ** 1.5 + 3.329 * (c_e / 1000)
    )

    E_sigma_e = pybamm.Parameter("Electrolyte conductivity activation energy [J.mol-1]")
    arrhenius = np.exp(E_sigma_e / pybamm.constants.R * (1 / 298.15 - 1 / T))

    scale = pybamm.Parameter("Electrolyte conductivity scaling factor")

    return scale * sigma_e * arrhenius


# =============================================================================
# CHEN2020_MOD FUNCTIONS
# =============================================================================

def graphite_ocp_Enertech_Ai2020(sto):
    """
    Graphite open-circuit potential from Enertech/Ai2020 data, loaded as a cubic
    interpolant. Alternative to graphite_LGM50_ocp_Chen2020.

    References
    ----------
    .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
    Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity in
    Lithium-Ion Pouch Cells. Journal of The Electrochemical Society, 167(1), 013512.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential
    """
    name, (x, y) = graphite_ocp_Enertech_Ai2020_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


# =============================================================================
# AI2020 DATA OCP FUNCTIONS
# =============================================================================

def lico2_ocp_Ai2020(sto):
    """
    LiCoO2 open-circuit potential as a tabulated interpolant from Ai2020.

    References
    ----------
    .. [1] Ai, W., Kraft, L., Sturm, J., Jossen, A., & Wu, B. (2020).
    Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity in
    Lithium-Ion Pouch Cells. Journal of The Electrochemical Society, 167(1), 013512.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential [V]
    """
    name, (x, y) = lico2_ocp_Ai2020_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


# =============================================================================
# CHEN2020 TABULATED OCP FUNCTIONS (data-file versions of the analytic fits)
# =============================================================================

def graphite_LGM50_ocp_Chen2020_tabulated(sto):
    """
    LG M50 graphite OCP as a tabulated cubic interpolant loaded from the Chen2020
    data file. Covers the same stoichiometry range as the analytic
    ``graphite_LGM50_ocp_Chen2020`` but follows the measured data exactly.

    References
    ----------
    .. [1] Chang-Hui Chen et al., J. Electrochem. Soc. 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential [V]
    """
    name, (x, y) = graphite_LGM50_ocp_Chen2020_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def nmc_LGM50_ocp_Chen2020_tabulated(sto):
    """
    LG M50 NMC OCP as a tabulated cubic interpolant loaded from the Chen2020
    data file. Covers the same stoichiometry range as the analytic
    ``nmc_LGM50_ocp_Chen2020`` but follows the measured data exactly.

    References
    ----------
    .. [1] Chang-Hui Chen et al., J. Electrochem. Soc. 167 (2020): 080534.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential [V]
    """
    name, (x, y) = nmc_LGM50_ocp_Chen2020_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


# =============================================================================
# NCA KIM2011 TABULATED OCP
# =============================================================================

def nca_ocp_Kim2011_tabulated(sto):
    """
    NCA OCP as a tabulated cubic interpolant from digitized Kim2011 data.
    Alternative to the analytic polynomial ``nca_ocp_Kim2011``.

    References
    ----------
    .. [1] Kim, G. H. et al., J. Electrochem. Soc. 158 (2011): A955-A969.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential [V]
    """
    name, (x, y) = nca_ocp_Kim2011_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


# =============================================================================
# HANYANG TABULATED OCP FUNCTIONS
# Digitized from https://real.hanyang.ac.kr/materials
# Cathodes: LCO, LFP, NCA, NMC111, NMC532, NMC622, NMC811, LMO, LNMO
# Anodes:   Gr, LTO, Si
# =============================================================================

def lco_ocp_Hanyang(sto):
    """Tabulated LCO (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lco_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def lfp_ocp_Hanyang(sto):
    """Tabulated LFP (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lfp_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def nca_ocp_Hanyang(sto):
    """Tabulated NCA (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nca_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def nmc111_ocp_Hanyang(sto):
    """Tabulated NMC111 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc111_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def nmc532_ocp_Hanyang(sto):
    """Tabulated NMC532 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc532_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def nmc622_ocp_Hanyang(sto):
    """Tabulated NMC622 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc622_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def nmc811_ocp_Hanyang(sto):
    """Tabulated NMC811 (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _nmc811_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def lmo_ocp_Hanyang(sto):
    """Tabulated LMO (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lmo_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def lnmo_ocp_Hanyang(sto):
    """Tabulated LNMO (cathode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lnmo_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def gr_ocp_Hanyang(sto):
    """Tabulated Gr (anode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _gr_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def lto_ocp_Hanyang(sto):
    """Tabulated LTO (anode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _lto_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


def si_ocp_Hanyang(sto):
    """Tabulated Si (anode) OCP [V] vs Li/Li+ (digitized, Hanyang)."""
    name, (x, y) = _si_ocp_Hanyang_data
    return pybamm.Interpolant(x, y, sto, name=name, interpolator="cubic")


# =============================================================================
# PRADA2013 FUNCTIONS
# =============================================================================

def LFP_ocp_Afshar2017(sto):
    """
    Open-circuit potential for LFP as a function of stoichiometry.

    References
    ----------
    .. [1] Afshar, S., Morris, K., & Khajepour, A. (2017). Efficient electrochemical
    model for lithium-ion cells. arXiv preprint arXiv:1709.03970.

    Parameters
    ----------
    sto : :class:`pybamm.Symbol`
       stoichiometry of material (li-fraction)

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential [V]
    """
    c1 = -150 * sto
    c2 = -30 * (1 - sto)
    k = 3.4077 - 0.020269 * sto + 0.5 * np.exp(c1) - 0.9 * np.exp(c2)

    return k


def LFP_electrolyte_exchange_current_density_kashkooli2017(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between LFP and electrolyte.

    References
    ----------
    .. [1] Kashkooli, A. G., Amirfazli, A., Farhad, S., Lee, D. U., Felicelli, S., Park,
    H. W., ... & Chen, Z. (2017). Representative volume element model of lithium-ion
    battery electrodes based on X-ray nano-tomography. Journal of Applied
    Electrochemistry, 47(3), 281-293.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    m_ref = pybamm.Parameter("Positive electrode kinetic rate constant [A.m-2]")
    E_r = pybamm.Parameter("Positive electrode exchange-current density activation energy [J.mol-1]")
    arrhenius = np.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return m_ref * arrhenius * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5


# =============================================================================
# NCA_KIM2011 FUNCTIONS
# =============================================================================

def graphite_diffusivity_Kim2011(sto, T):
    """
    Graphite diffusivity as a function of stoichiometry and temperature [1].

    References
    ----------
    .. [1] Kim, G. H., Smith, K., Lee, K. J., Santhanagopalan, S., & Pesaran, A.
    (2011). Multi-domain modeling of lithium-ion batteries encompassing
    multi-physics in varied length scales. Journal of The Electrochemical
    Society, 158(8), A955-A969.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """
    D_ref = 9e-14
    E_D_s = 4e3
    arrhenius = np.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


def graphite_ocp_Kim2011(sto):
    """
    Graphite Open-circuit Potential (OCP) as a function of the stoichiometry [1].

    References
    ----------
    .. [1] Kim, G. H., Smith, K., Lee, K. J., Santhanagopalan, S., & Pesaran, A.
    (2011). Multi-domain modeling of lithium-ion batteries encompassing
    multi-physics in varied length scales. Journal of The Electrochemical
    Society, 158(8), A955-A969.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential [V]
    """
    u_eq = (
        0.124
        + 1.5 * np.exp(-70 * sto)
        - 0.0351 * np.tanh((sto - 0.286) / 0.083)
        - 0.0045 * np.tanh((sto - 0.9) / 0.119)
        - 0.035 * np.tanh((sto - 0.99) / 0.05)
        - 0.0147 * np.tanh((sto - 0.5) / 0.034)
        - 0.102 * np.tanh((sto - 0.194) / 0.142)
        - 0.022 * np.tanh((sto - 0.98) / 0.0164)
        - 0.011 * np.tanh((sto - 0.124) / 0.0226)
        + 0.0155 * np.tanh((sto - 0.105) / 0.029)
    )

    return u_eq


def graphite_electrolyte_exchange_current_density_Kim2011(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between graphite and LiPF6 in
    EC:DMC [1].

    References
    ----------
    .. [1] Kim, G. H., Smith, K., Lee, K. J., Santhanagopalan, S., & Pesaran, A.
    (2011). Multi-domain modeling of lithium-ion batteries encompassing
    multi-physics in varied length scales. Journal of The Electrochemical
    Society, 158(8), A955-A969.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    i0_ref = 36  # reference exchange current density at 100% SOC
    sto_ref = 0.36  # stoichiometry at 100% SOC
    c_s_n_ref = sto_ref * c_s_max  # reference electrode concentration
    c_e_ref = pybamm.Parameter("Initial concentration in electrolyte [mol.m-3]")
    alpha = 0.5  # charge transfer coefficient

    m_ref = i0_ref / (
        c_e_ref**alpha * (c_s_max - c_s_n_ref) ** alpha * c_s_n_ref**alpha
    )

    E_r = 3e4
    arrhenius = np.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return m_ref * arrhenius * c_e**alpha * c_s_surf**alpha * (c_s_max - c_s_surf) ** alpha


def nca_diffusivity_Kim2011(sto, T):
    """
    NCA diffusivity as a function of stoichiometry and temperature [1].

    References
    ----------
    .. [1] Kim, G. H., Smith, K., Lee, K. J., Santhanagopalan, S., & Pesaran, A.
    (2011). Multi-domain modeling of lithium-ion batteries encompassing
    multi-physics in varied length scales. Journal of The Electrochemical
    Society, 158(8), A955-A969.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """
    D_ref = 3e-15
    E_D_s = 2e4
    arrhenius = np.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


def nca_electrolyte_exchange_current_density_Kim2011(c_e, c_s_surf, c_s_max, T):
    """
    Exchange-current density for Butler-Volmer reactions between NCA and LiPF6 in
    EC:DMC [1].

    References
    ----------
    .. [1] Kim, G. H., Smith, K., Lee, K. J., Santhanagopalan, S., & Pesaran, A.
    (2011). Multi-domain modeling of lithium-ion batteries encompassing
    multi-physics in varied length scales. Journal of The Electrochemical
    Society, 158(8), A955-A969.

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    i0_ref = 4  # reference exchange current density at 100% SOC
    sto_ref = 0.41  # stoichiometry at 100% SOC
    c_s_ref = sto_ref * c_s_max  # reference electrode concentration
    c_e_ref = pybamm.Parameter("Initial concentration in electrolyte [mol.m-3]")
    alpha = 0.5  # charge transfer coefficient

    m_ref = i0_ref / (c_e_ref**alpha * (c_s_max - c_s_ref) ** alpha * c_s_ref**alpha)
    E_r = 3e4
    arrhenius = np.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return m_ref * arrhenius * c_e**alpha * c_s_surf**alpha * (c_s_max - c_s_surf) ** alpha


def nca_ocp_Kim2011(sto):
    """
    NCA Open Circuit Potential (OCP) as a function of the stoichiometry [1].

    References
    ----------
    .. [1] Kim, G. H., Smith, K., Lee, K. J., Santhanagopalan, S., & Pesaran, A.
    (2011). Multi-domain modeling of lithium-ion batteries encompassing
    multi-physics in varied length scales. Journal of The Electrochemical
    Society, 158(8), A955-A969.

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stoichiometry

    Returns
    -------
    :class:`pybamm.Symbol`
        Open-circuit potential [V]
    """
    U_posi = (
        1.638 * sto**10
        - 2.222 * sto**9
        + 15.056 * sto**8
        - 23.488 * sto**7
        + 81.246 * sto**6
        - 344.566 * sto**5
        + 621.3475 * sto**4
        - 554.774 * sto**3
        + 264.427 * sto**2
        - 66.3691 * sto
        + 11.8058
        - 0.61386 * np.exp(5.8201 * sto**136.4)
    )

    return U_posi


# =============================================================================
# PARAMETER VALUES
# =============================================================================

# Call dict via a function to avoid errors when editing in place
def get_parameter_values():
    """
    Combined parameter set based on Chen2020 (LG M50 NMC/graphite cell) with all
    functions from Chen2020, Chen2020_mod, Prada2013, and NCA_Kim2011 included.

    Base scalar values are taken from Chen2020. Additional parameters unique to
    Chen2020_mod (thermal LAM factors), Prada2013 (lithium interstitial molar volume),
    and NCA_Kim2011 (tab geometry, per-face heat transfer coefficients) are appended
    with their respective default values.

    References
    ----------
    .. [1] Chen2020 - Chang-Hui Chen et al., J. Electrochem. Soc. 167 (2020): 080534.
    .. [2] Prada2013 - E. Prada et al., J. Electrochem. Soc. 160 (2013): A616-A628.
    .. [3] NCA_Kim2011 - G. H. Kim et al., J. Electrochem. Soc. 158 (2011): A955-A969.
    .. [4] OKane2020 - S. E. J. O'Kane et al., J. Electrochem. Soc. 167 (2020): 090540.
    """

    return {
        # exchange current densities from OKane2020
        "Negative electrode exchange-current density activation energy [J.mol-1]": 35000.0,
        "Negative electrode kinetic rate constant [A.m-2]": 6.48e-7,
        "Positive electrode exchange-current density activation energy [J.mol-1]": 17800,
        "Positive electrode kinetic rate constant [A.m-2]": 3.42e-6,
        # lithium plating
        "Lithium plating potential sharpness": 100,
        "Lithium metal partial molar volume [m3.mol-1]": 1.3e-05,
        "Lithium plating kinetic rate constant [m.s-1]": 1e-09,
        "Exchange-current density for plating [A.m-2]": graphite_plating_exchange_current_density_OKane2020,
        "Exchange-current density for stripping [A.m-2]": graphite_stripping_exchange_current_density_OKane2020,
        "Initial plated lithium concentration [mol.m-3]": 0.0,
        "Typical plated lithium concentration [mol.m-3]": 1000.0,
        "Lithium plating transfer coefficient": 0.65,
        "Dead lithium decay constant [s-1]": 1e-06,
        "Dead lithium decay rate [s-1]": graphite_SEI_limited_dead_lithium_OKane2022,
        # negative electrode cracking (Ai2020)
        "Negative electrode volume change": graphite_volume_change_Ai2020,
        "Negative electrode initial crack length [m]": 2e-08,
        "Negative electrode initial crack width [m]": 1.5e-08,
        "Negative electrode number of cracks per unit area [m-2]": 3180000000000000.0,
        "Negative electrode Paris' law constant b": 1.12,
        "Negative electrode Paris' law constant m": 2.2,
        "Negative electrode cracking rate": graphite_cracking_rate_Ai2020,
        "Negative electrode cracking rate constant [m/(Pa.m0.5)^m_cr]": 3.9e-20,
        "Negative electrode LAM constant proportional term [s-1]": 2.7778e-07,
        "Negative electrode LAM constant exponential term": 2.0,
        "Negative electrode critical stress [Pa]": 60000000.0,
        "Negative electrode partial molar volume [m3.mol-1]": 3.1e-06,
        "Negative electrode Poisson's ratio": 0.3,
        "Negative electrode Young's modulus [Pa]": 15000000000.0,
        "Negative electrode reference concentration for free of deformation [mol.m-3]": 0,
        # negative electrode LAM thermal factors (from Chen2020_mod)
        "Negative electrode LAM thermal pre-exponential factor [s-1]": 1.0e-09,
        "Negative electrode LAM activation energy [J.mol-1]": 40000.0,
        "Negative electrode LAM time exponent": 0.0,
        # positive electrode cracking (Ai2020)
        "Positive electrode Poisson's ratio": 0.2,
        "Positive electrode Young's modulus [Pa]": 375000000000.0,
        "Positive electrode reference concentration for free of deformation [mol.m-3]": 0.0,
        "Positive electrode partial molar volume [m3.mol-1]": 1.25e-05,
        "Positive electrode volume change": volume_change_Ai2020,
        "Positive electrode initial crack length [m]": 2e-08,
        "Positive electrode initial crack width [m]": 1.5e-08,
        "Positive electrode number of cracks per unit area [m-2]": 3180000000000000.0,
        "Positive electrode Paris' law constant b": 1.12,
        "Positive electrode Paris' law constant m": 2.2,
        "Positive electrode cracking rate": cracking_rate_Ai2020,
        "Positive electrode cracking rate constant [m/(Pa.m0.5)^m_cr]": 3.9e-20,
        "Positive electrode LAM constant proportional term [s-1]": 2.7778e-07,
        "Positive electrode LAM constant exponential term": 2.0,
        "Positive electrode critical stress [Pa]": 375000000.0,
        # positive electrode LAM thermal factors (from Chen2020_mod)
        "Positive electrode LAM thermal pre-exponential factor [s-1]": 1.0e-09,
        "Positive electrode LAM activation energy [J.mol-1]": 40000.0,
        "Positive electrode LAM time exponent": 0.0,
        # general
        "chemistry": "lithium_ion",
        "Negative particle diffusivity constant [m2.s-1]": 3.3e-14,
        "Negative particle diffusivity activation energy [J.mol-1]": 30300.0,
        "Positive particle diffusivity constant [m2.s-1]": 4e-15,
        "Positive particle diffusivity activation energy [J.mol-1]": 25000.0,
        # SEI
        "Ratio of lithium moles to SEI moles": 2.0,
        "SEI partial molar volume [m3.mol-1]": 9.585e-05,
        "SEI reaction exchange current density [A.m-2]": 1.5e-07,
        "SEI resistivity [Ohm.m]": 200000.0,
        "SEI solvent diffusivity [m2.s-1]": 2.5e-22,
        "Bulk solvent concentration [mol.m-3]": 2636.0,
        "SEI open-circuit potential [V]": 0.4,
        "SEI electron conductivity [S.m-1]": 8.95e-14,
        "SEI lithium interstitial diffusivity [m2.s-1]": 1e-20,
        "Lithium interstitial reference concentration [mol.m-3]": 15.0,
        "Initial SEI thickness [m]": 5e-09,
        "EC initial concentration in electrolyte [mol.m-3]": 4541.0,
        "EC diffusivity [m2.s-1]": 2e-18,
        "SEI kinetic rate constant [m.s-1]": 1e-12,
        "SEI growth activation energy [J.mol-1]": 0.0,
        "Negative electrode reaction-driven LAM factor [m3.mol-1]": 0.0,
        "Positive electrode reaction-driven LAM factor [m3.mol-1]": 0.0,
        "Initial SEI on cracks thickness [m]": 5e-13,
        # lithium interstitial molar volume (from Prada2013)
        "Lithium interstitial molar volume [m3.mol-1]": 1.3e-05,
        # cell
        "Negative current collector thickness [m]": 1.2e-05,
        "Negative electrode thickness [m]": 8.52e-05,
        "Separator thickness [m]": 1.2e-05,
        "Positive electrode thickness [m]": 7.56e-05,
        "Positive current collector thickness [m]": 1.6e-05,
        "Electrode height [m]": 0.065,
        "Electrode width [m]": 1.58,
        "Cell cooling surface area [m2]": 0.00531,
        "Cell volume [m3]": 2.42e-05,
        "Cell thermal expansion coefficient [m.K-1]": 1.1e-06,
        "Negative current collector conductivity [S.m-1]": 58411000.0,
        "Positive current collector conductivity [S.m-1]": 36914000.0,
        "Negative current collector density [kg.m-3]": 8960.0,
        "Positive current collector density [kg.m-3]": 2700.0,
        "Negative current collector specific heat capacity [J.kg-1.K-1]": 385.0,
        "Positive current collector specific heat capacity [J.kg-1.K-1]": 897.0,
        "Negative current collector thermal conductivity [W.m-1.K-1]": 401.0,
        "Positive current collector thermal conductivity [W.m-1.K-1]": 237.0,
        "Nominal cell capacity [A.h]": 5.0,
        "Current function [A]": 5.0,
        "Contact resistance [Ohm]": 0,
        # tab geometry (from NCA_Kim2011, for 2D/3D thermal models)
        "Negative tab width [m]": 0.044,
        "Negative tab centre y-coordinate [m]": 0.013,
        "Negative tab centre z-coordinate [m]": 0.2,
        "Positive tab width [m]": 0.044,
        "Positive tab centre y-coordinate [m]": 0.137,
        "Positive tab centre z-coordinate [m]": 0.2,
        # negative electrode
        "Negative electrode conductivity [S.m-1]": 215.0,
        "Maximum concentration in negative electrode [mol.m-3]": 33133.0,
        "Negative particle diffusivity [m2.s-1]": graphite_LGM50_diffusivity_Chen2020,
        "Negative electrode OCP [V]": graphite_LGM50_ocp_Chen2020,
        "Negative electrode porosity": 0.25,
        "Negative electrode active material volume fraction": 0.75,
        "Negative particle radius [m]": 5.86e-06,
        "Negative electrode Bruggeman coefficient (electrolyte)": 1.5,
        "Negative electrode Bruggeman coefficient (electrode)": 0,
        "Negative electrode charge transfer coefficient": 0.5,
        "Negative electrode double-layer capacity [F.m-2]": 0.2,
        "Negative electrode exchange-current density [A.m-2]": graphite_LGM50_electrolyte_exchange_current_density_Chen2020,
        "Negative electrode density [kg.m-3]": 1657.0,
        "Negative electrode specific heat capacity [J.kg-1.K-1]": 700.0,
        "Negative electrode thermal conductivity [W.m-1.K-1]": 1.7,
        "Negative electrode OCP entropic change [V.K-1]": 0.0,
        # positive electrode
        "Positive electrode conductivity [S.m-1]": 0.18,
        "Maximum concentration in positive electrode [mol.m-3]": 63104.0,
        "Positive particle diffusivity [m2.s-1]": nmc_LGM50_diffusivity_Chen2020,
        "Positive electrode OCP [V]": nmc_LGM50_ocp_Chen2020,
        "Positive electrode porosity": 0.335,
        "Positive electrode active material volume fraction": 0.665,
        "Positive particle radius [m]": 5.22e-06,
        "Positive electrode Bruggeman coefficient (electrolyte)": 1.5,
        "Positive electrode Bruggeman coefficient (electrode)": 0,
        "Positive electrode charge transfer coefficient": 0.5,
        "Positive electrode double-layer capacity [F.m-2]": 0.2,
        "Positive electrode exchange-current density [A.m-2]": nmc_LGM50_electrolyte_exchange_current_density_Chen2020,
        "Positive electrode density [kg.m-3]": 3262.0,
        "Positive electrode specific heat capacity [J.kg-1.K-1]": 700.0,
        "Positive electrode thermal conductivity [W.m-1.K-1]": 2.1,
        "Positive electrode OCP entropic change [V.K-1]": 0.0,
        # separator
        "Separator porosity": 0.47,
        "Separator Bruggeman coefficient (electrolyte)": 1.5,
        "Separator density [kg.m-3]": 397.0,
        "Separator specific heat capacity [J.kg-1.K-1]": 700.0,
        "Separator thermal conductivity [W.m-1.K-1]": 0.16,
        # electrolyte
        "Initial concentration in electrolyte [mol.m-3]": 1000.0,
        "Cation transference number": 0.2594,
        "Thermodynamic factor": 1.0,
        "Electrolyte diffusivity [m2.s-1]": electrolyte_diffusivity_Nyman2008,
        "Electrolyte diffusivity scaling factor": 1.0,
        "Electrolyte diffusivity activation energy [J.mol-1]": 17000.0,
        "Electrolyte conductivity [S.m-1]": electrolyte_conductivity_Nyman2008,
        "Electrolyte conductivity scaling factor": 1.0,
        "Electrolyte conductivity activation energy [J.mol-1]": 17000.0,
        # experiment
        "Reference temperature [K]": 298.15,
        "Total heat transfer coefficient [W.m-2.K-1]": 10.0,
        "Ambient temperature [K]": 298.15,
        "Number of electrodes connected in parallel to make a cell": 1.0,
        "Number of cells connected in series to make a battery": 1.0,
        "Lower voltage cut-off [V]": 2.5,
        "Upper voltage cut-off [V]": 4.2,
        "Open-circuit voltage at 0% SOC [V]": 2.5,
        "Open-circuit voltage at 100% SOC [V]": 4.2,
        "Initial concentration in negative electrode [mol.m-3]": 29866.0,
        "Initial concentration in positive electrode [mol.m-3]": 17038.0,
        "Initial temperature [K]": 298.15,
        # per-face heat transfer coefficients (from NCA_Kim2011, for 2D/3D thermal models)
        "Negative current collector surface heat transfer coefficient [W.m-2.K-1]": 0.0,
        "Positive current collector surface heat transfer coefficient [W.m-2.K-1]": 0.0,
        "Negative tab heat transfer coefficient [W.m-2.K-1]": 25.0,
        "Positive tab heat transfer coefficient [W.m-2.K-1]": 25.0,
        "Edge heat transfer coefficient [W.m-2.K-1]": 0.3,
        # citations
        "citations": ["Chen2020", "Prada2013", "Kim2011"],
    }
