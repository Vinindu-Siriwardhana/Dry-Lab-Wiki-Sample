"""
params.py -- the single source of truth for every parameter used in the figures.

No value is hard-coded inside a solver; everything enters through this file, so
changing a number here changes every figure that depends on it.  Each entry
carries its provenance TIER, following the team constants register
(constants_register.xlsx, updated 2026-10-09):

    T1  measured in house
    T2  mechanistic correlation
    T3  named literature analogue, matching criterion stated
    T4  order-of-magnitude prior (>= 1 decade)
    D   derived from other entries
    L   locked design assumption

Units: cm, s, mol cm^-3.  1 uM = 1e-9 mol cm^-3.
"""
import math

UM = 1e-9
DAY = 86400.0
HOUR = 3600.0
LN2 = math.log(2.0)

# --- physical constants (register sheet 0) ---------------------------------
T_BODY       = 310.0        # K    37 C                                    D
K_B          = 1.380649e-23 # J/K                                          D
ETA_WATER    = 0.00069      # Pa s water at 37 C                           T3

# --- geometry ---------------------------------------------------------------
L_GEL        = 0.015    # cm   (150 um) as-manufactured 5% GelMA          T1
L_GEL_ALT    = 0.05     # cm   (500 um) V3-document baseline, comparison  design
L_TIS        = 0.15     # cm   (1.5 mm) dermal depth to the systemic sink T3
L_TIS_RANGE  = (0.10, 0.20)   # cm, band carried as sensitivity           T3
DELTA        = 0.05     # cm   (500 um) fibroblast-active / target depth  T4 (sweep 300-1000 um)
R_W_DEFAULT  = 0.2      # cm   (2 mm) wound radius (grid 1-5 mm)          T4

# --- species ----------------------------------------------------------------
MW_V14       = 1355.59  # Da                                               T1
MW_FGF       = 16487.87 # Da                                               T1
A_V14_NM     = 0.73     # nm   hydrodynamic radius (0.87 = sensitivity)    T2
A_FGF_NM     = 1.68     # nm   (2.34 = sensitivity bound)                  T2


def stokes_einstein(a_nm):
    """Free diffusivity in cm2/s from the Stokes-Einstein relation at 37 C."""
    return K_B * T_BODY / (6 * math.pi * ETA_WATER * a_nm * 1e-9) * 1e4


D0_V14       = 4.5e-6   # cm2/s Stokes-Einstein from A_V14_NM, ETA     D  (doc band 2.5-4.5e-6)
D0_FGF       = 2.0e-6   # cm2/s Stokes-Einstein from A_FGF_NM, ETA     D  (doc 1.3-1.7e-6)
T_HALF_V14   = 6.1 * HOUR   # s  serum half-life proxy                     T4
T_HALF_FGF   = 168 * HOUR   # s  functional half-life at 37 C              T1

# --- hydrogel (Model 1) -----------------------------------------------------
A_FIBRE_NM   = 0.6      # nm   gelatin chain radius                        T3
XI_NM        = 10.0     # nm   mesh size (range 5-50 swept)                T4
EPS_GEL      = 0.618    # --   5% GelMA porosity                           T3
TAU_GEL      = EPS_GEL ** -0.5   # 1.27, Bruggeman                         T2
K_DEG_GEL    = LN2 / T_HALF_V14  # 1/s  in-gel degradation V14 (3.15e-5)   D
K_DEG_GEL_F  = LN2 / T_HALF_FGF  # 1/s  in-gel degradation FGF (1.15e-6)   D
S_MAX        = 0.0      # mol/cm3 heparin/HS site density -- LOCKED: no heparin   L
R_V          = 1.0      # --   V14 sorption retardation -- LOCKED            L

# --- tissue (Model 2) -------------------------------------------------------
EPS_TIS      = 0.80     # --   tissue porosity                             T3
K_P          = 1.0      # --   gel/tissue partition coefficient            T4
D_TIS_V14    = 1.2e-6   # cm2/s tissue diffusivity (range 0.7-1.8e-6)      T3
D_TIS_V14_LO = 0.7e-6
D_TIS_V14_HI = 1.8e-6
D_TIS_FGF    = 5e-7     # cm2/s                                            T4
K_PROT_V14   = LN2 / T_HALF_V14  # 1/s (3.15e-5) -- HIGHEST-PRIORITY GAP   T4
K_PROT_FGF   = LN2 / T_HALF_FGF  # 1/s (1.15e-6)                           T4
K_PROT_RANGE = (1e-3, 1e-4, K_PROT_V14, 1e-6)   # register sweep 1e-6..1e-3
K_CL_FGF     = 0.000115 # 1/s  NOT in the PDE (Dirichlet sink replaces it) T3
K_CL_V14     = 0.0035   # 1/s  NOT in the PDE; unverified                  T4

# --- dosing (Model 7) -------------------------------------------------------
C_TARGET_V14 = 25 * UM     # mol/cm3  efficacious tissue concentration    T4
C_SOLUBILITY = 738 * UM    # mol/cm3  1 mg/mL saturation ceiling          T4 -- MEASURE THIS
F_REL        = 0.5         # --  fraction of loaded dose released         T4
REAPPLY_D    = 1.0         # d   dressing-change interval                 T4

# --- receptor / signalling (Models 3-4) -------------------------------------
KD_L         = 6.5e-11  # mol/cm3 (65 nM) LPS-MD2                          T3
L_LPS        = 1.0 * UM # mol/cm3 (1 uM) wound endotoxin load              T4
A_LPS        = L_LPS / KD_L   # -- scaled endotoxin load, 15.4             D
THETA_UNTREATED = A_LPS / (1 + A_LPS)   # 0.94, closure calibration anchor D
TLR4_TOT     = 1.0      # nM                                               T4
GAMMA        = 10.0     # --  scaled receptor density                      T4
HILL_H       = 2.0      # --  NF-kB p65 Hill coefficient                   T3
K_NF         = 10.0     # nM                                               T4
P65_BAS, P65_MAX = 0.5, 50.0   # nM                                        T4
PSI_RANGE    = (1e-2, 1e4)

# --- nitric oxide (Learning Cycle 1 and Model 5) ---------------------------
D0_NO        = 3.0e-5   # cm2/s NO in water                                T3
NO_HALFLIFE  = (0.09, 0.5, 2.0)   # s, normoxic tissue, report as a RANGE:
                                  # autoxidation is second order, so t_1/2 is
                                  # concentration dependent                T3
INV_TAU      = (1.0, 0.5, 0.3)    # swept group 1/tau                      T2
STACK_MM     = (0.5, 3.0)         # V2 cellulose + hydrogel stack, mm      design (V2)
K_SCAV       = 0.01     # 1/s  NO scavenging                               T3
Q_CHRONIC    = 5.0e-19  # mol/cell/s sustained iNOS rate (3.5-8.1e-19)     T3
Q_MAX_ACUTE  = 8.17e-17 # mol/cell/s                                       T3
N_MAC        = 262      # cells/mm2 chronic baseline M1 density            T3
MAC_BAND     = 0.005    # cm                                               T4
K_TURN       = LN2 / (4 * HOUR)   # 1/s iNOS turnover                      T4
K_INOS_REL   = 0.3      # x p65_max                                        T4
P_INOS       = 2.0      # --                                               T4
PHI_NO2      = 0.9      # --  NO -> nitrite yield                          T4

# --- closure (Model 6) -- every entry here is Tier 4 ------------------------
D_N          = 1e-9     # cm2/s fibroblast random motility                 T4
D_N_RANGE    = (1e-10, 1e-8)
K_CARRY      = 1.0      # normalised carrying capacity (5e6 cells/cm3)
R0           = 0.8 / DAY   # 1/s basal proliferation, chronic              T4
D0_DEATH     = 0.24 / DAY  # 1/s basal death, chronic (calibrated 2026-10-09, was 1.0)
DELTA_D      = 0.3 / DAY   # 1/s max death increment (calibrated, was 0.7) T4
E_MAX        = 2.0         # --  FGF max efficacy                          T4
EC50_F       = 0.5         # ng/mL                                         T4
BETA_R       = 0.5         # --  inflammation penalty on proliferation     T4
K_I_R        = 0.3         # --                                            T4
K_I_D        = 0.3         # --                                            T4
S_I          = 2.0         # --  death Hill exponent                       T4
KE0          = 1.0 / (18 * HOUR)  # 1/s effect-compartment rate            T4
W_SMOOTH     = 75e-4       # cm  IC smoothing width (50-100 um)            T4


def rp_of(C_e, I_e):
    """Proliferation r_p(C_e, I_e) in 1/s, register closures."""
    return R0 * (1 + E_MAX * C_e / (EC50_F + C_e)) * (1 - BETA_R * I_e / (K_I_R + I_e))


def d_of(I_e):
    """Death d(I_e) in 1/s, register closures."""
    return D0_DEATH + DELTA_D * I_e ** S_I / (K_I_D ** S_I + I_e ** S_I)


# chronic baseline: untreated wound, no FGF, theta_LPS = 0.94
I_CHRONIC    = THETA_UNTREATED
R_CHRONIC    = rp_of(0.0, I_CHRONIC)      # 0.497 /d
D_CHRONIC    = d_of(I_CHRONIC)            # 0.512 /d  -> net -0.015 /d: stall
# therapeutic state (scenario): FGF2-G3 at EC50, V14 at psi_50 (half the
# active receptor suppressed, so I_e = theta_0 / 2)
I_TREATED    = THETA_UNTREATED / 2
R_PULSE      = rp_of(EC50_F, I_TREATED)   # 1.11 /d
D_PULSE      = d_of(I_TREATED)            # 0.45 /d
TAU_RELAX    = 1.0 / KE0                  # s, drug-effect relaxation (18 h)

PROVENANCE = {
    # module: counts of T1-T4 entries in the constants register, and how many
    # of them have no planned experiment.  D (derived) and L (locked) rows are
    # not counted.
    "M1 release":      dict(T1=4, T2=3, T3=3, T4=3, unplanned=1),
    "M2 tissue":       dict(T1=0, T2=0, T3=4, T4=6, unplanned=5),
    "M3 receptor":     dict(T1=0, T2=0, T3=2, T4=3, unplanned=2),
    "M4 signalling":   dict(T1=0, T2=0, T3=1, T4=4, unplanned=4),
    "M5 NO readout":   dict(T1=0, T2=0, T3=9, T4=5, unplanned=5),
    "M6 closure":      dict(T1=0, T2=0, T3=0, T4=15, unplanned=15),
    "M7 design":       dict(T1=0, T2=0, T3=0, T4=3, unplanned=2),
}
