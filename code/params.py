"""
params.py -- the single source of truth for every parameter used in the figures.

No value is hard-coded inside a solver; everything enters through this file, so
changing a number here changes every figure that depends on it.  Each entry
carries its provenance TIER, following the team constants register
(constants_register.xlsx, updated 2026-10-16; gel changed to 1% Ca-alginate
on 2026-10-09, and to the 2.5 mm / 3 mg/mL V14 / 200 ng/mL FGF2-G3 design with
f_rel = 0.5 / 0.15 on 2026-10-16; f_rel,FGF = 0.31 from the BCA fit on
2026-10-16):

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
L_GEL        = 0.25     # cm   (2.5 mm) as-manufactured 1% Ca-alginate,
                        #      team design 2026-10-16 (was 500 um)          T1
L_GEL_ALT    = 0.015    # cm   (150 um) earlier GelMA build, comparison   design
L_TIS        = 0.15     # cm   (1.5 mm) dermal depth to the systemic sink T3
L_TIS_RANGE  = (0.10, 0.20)   # cm, band carried as sensitivity           T3
DELTA        = 0.05     # cm   (500 um) fibroblast-active / target depth  T4 (sweep 300-1000 um)
R_W_DEFAULT  = 0.2      # cm   (2 mm) wound radius (grid 1-5 mm)          T4

# --- species ----------------------------------------------------------------
MW_V14       = 1355.59  # Da                                               T1
MW_FGF       = 16487.87 # Da                                               T1
A_V14_NM     = 0.87     # nm   hydrodynamic radius, fitted correlation     T2
                        #      a_h = 0.0515 M^0.392 (doc App. A); band 0.8-1.0.
                        #      0.73 = anhydrous-sphere lower bound, not used
A_FGF_NM     = 2.34     # nm   same correlation; band 1.9-2.4 (Wilkins ~2.0) T2
                        #      1.68 = anhydrous-sphere lower bound, not used


def stokes_einstein(a_nm):
    """Free diffusivity in cm2/s from the Stokes-Einstein relation at 37 C."""
    return K_B * T_BODY / (6 * math.pi * ETA_WATER * a_nm * 1e-9) * 1e4


D0_V14       = 3.6e-6   # cm2/s consensus of 3 correlations, doc App. A T2 (band 2.5-4.5e-6)
D0_FGF       = 1.5e-6   # cm2/s same                                    T2 (band 1.3-1.7e-6)
                        # stokes_einstein() at the anhydrous radii gives 4.5 / 2.0e-6:
                        # upper bounds only
T_HALF_V14   = 6.1 * HOUR   # s  serum half-life proxy                     T4
T_HALF_FGF   = 168 * HOUR   # s  functional half-life at 37 C              T1

# --- hydrogel (Model 1) -----------------------------------------------------
# 1% w/v sodium alginate, ionically cross-linked in 100 mM CaCl2, in HBS
# (10 mM HEPES, 150 mM NaCl, pH 7.4).  Team recipe.                         T1
A_FIBRE_NM   = 0.6      # nm   alginate chain radius; not measured for this gel.
                        #      at A_FGF_NM = 2.34 and XI_NM = 14, raising it
                        #      0.6 -> 2 nm lowers H(FGF) 0.89 -> 0.83        T4
XI_NM        = 14.0     # nm   1% Ca-alginate mesh, Grassi et al. 2009,
                        #      Molecules 14:3003 (14 mechanical, 26 release) T3
EPS_GEL      = 0.99     # --   1% w/v polymer -> ~0.6% v/v solid           D
TAU_GEL      = EPS_GEL ** -0.5   # 1.005, Bruggeman                        T2
K_DEG_GEL    = LN2 / T_HALF_V14  # 1/s  in-gel degradation V14 (3.15e-5)   D
K_DEG_GEL_F  = LN2 / T_HALF_FGF  # 1/s  in-gel degradation FGF (1.15e-6)   D
S_MAX        = 0.0      # mol/cm3 binding sites -- LOCKED: plain alginate does not
                        #      bind FGFs (Freeman et al. 2008, SPR)          L
R_V          = 1.0      # --   V14 sorption -- LOCKED: alginate is hydrophilic L

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
K_CL_FGF     = 0.000115 # 1/s  register value only: NOT used anywhere (Dirichlet sink) T3
K_CL_V14     = 0.0035   # 1/s  register value only: NOT used anywhere; unverified T4

# --- dosing (Model 7) -------------------------------------------------------
# Target from the 14 Oct 2026 flow panel (replaces the earlier 25 uM T4 guess):
#   gate-frequency detection limit 27.8 % (t-test, n = 3)
#     -> occupancy drop 13.1 % via v3model.occupancy_from_gate at SIGMA_FLOW_LOG
#     -> K_D,P >= 10.08 uM (one-sided LOWER BOUND)
#     -> C_target = K_D,P (a - 1) = 145 uM: the dose at which 50 % of the MD2 pool
#        is LPS-bound (psi_50 = 1 + a, which halves theta_LPS, is 165 uM).
SIGMA_FLOW_LOG = 0.7       # --  log-SD of per-cell FITC intensity       T4 -- needs FCS files
KD_P_LOWER   = 10.08 * UM  # mol/cm3  lower bound on V14 K_D,P            D
C_TARGET_V14 = 145 * UM    # mol/cm3  efficacious tissue concentration   D
# Patch design, 2026-10-16.  Loadings are TOTAL per gel volume; only the
# released (mobile) fraction f_rel diffuses, the rest stays in the gel for the
# life of the patch, so the free initial concentration is f_rel C_load / R_V
# (twodomain.free_from_total).
C_LOAD_V14   = 3.0e-3 / MW_V14   # mol/cm3  3 mg/mL = 2213 uM                T1 (design)
C_LOAD_FGF_NGML = 200.0          # ng/mL    FGF2-G3 loading (12.1 nM)        T1 (design)
C_LOAD_FGF   = C_LOAD_FGF_NGML * 1e-9 / MW_FGF   # mol/cm3                    D
C_SOLUBILITY = C_LOAD_V14  # mol/cm3  LOWER BOUND: 3 mg/mL dissolves (was a
                           #          1 mg/mL = 738 uM T4 guess). Saturation
                           #          solubility still unmeasured.           T1 (bound)
F_REL_V14    = 0.5         # --  released fraction of V14 (= f_mob)         T4
F_REL_FGF    = 0.31        # --  released fraction of FGF2-G3 (= f_mob)     T1
                           #     least-squares fit of f_rel x Crank to the 30 Sep
                           #     BCA 30-min read (G03): 0.31 +/- 0.07 (SE, n = 12).
                           #     Was a 0.15 T4 guess. Assay is near its noise floor.
F_REL_FGF_SE = 0.07        # --  standard error of F_REL_FGF                T1
F_REL        = F_REL_V14   # --  alias kept for the shipped code/figures scripts
REAPPLY_D    = 1.0         # d   dressing-change interval                 T4

# --- receptor / signalling (Models 3-4) -------------------------------------
KD_L         = 6.5e-11  # mol/cm3 (65 nM) LPS-MD2                          T3
L_LPS        = 1.0 * UM # mol/cm3 (1 uM) wound endotoxin load              T4
A_LPS        = L_LPS / KD_L   # -- scaled endotoxin load, 15.4             D
THETA_UNTREATED = A_LPS / (1 + A_LPS)   # 0.94, closure calibration anchor D
TLR4_TOT     = 1.0      # nM                                               T4
HILL_H       = 2.0      # --  NF-kB p65 Hill coefficient                   T3
K_NF         = 10.0     # nM                                               T4
GAMMA        = TLR4_TOT / K_NF   # -- scaled receptor density, 0.1         D
P65_BAS, P65_MAX = 0.5, 50.0   # nM                                        T4
PSI_RANGE    = (1e-2, 1e4)

# --- nitric oxide (Learning Cycle 1 and Model 5) ---------------------------
# TWO DIFFERENT RATES, and conflating them was an error in earlier drafts.
# NO_HALFLIFE is NO consumption in TISSUE (cell respiration + haem), which is
# what sets the penetration depth L_d in the legacy sensing analysis.
# K_SCAV is autoxidation in CULTURE MEDIUM, which is ~2 decades slower and is
# what Eq. (no-well) uses, because Model 5 is a well, not tissue.
D0_NO        = 3.0e-5   # cm2/s NO in water                                T3
NO_HALFLIFE  = (0.09, 0.5, 2.0)   # s, normoxic TISSUE; report as a RANGE:
                                  # consumption is second order, so t_1/2 is
                                  # concentration dependent.  Thomas 2001,
                                  # PNAS 98:355 (PMID 11134509)            T3
INV_TAU      = (1.0, 0.5, 0.3)    # swept group 1/tau                      T2
STACK_MM     = (0.5, 3.0)         # V2 cellulose + hydrogel stack, mm      design (V2)
K_SCAV       = 0.01     # 1/s  NO autoxidation in CULTURE MEDIUM (t_1/2 =
                        #      ln2/K_SCAV = 69 s).  NOT the tissue rate
                        #      above; still fast vs iNOS turnover (4 h), so
                        #      the quasi-steady free-NO form holds.         T3
Q_CHRONIC    = 5.0e-19  # mol/cell/s sustained iNOS rate (3.5-8.1e-19)     T3
Q_MAX_ACUTE  = 8.17e-17 # mol/cell/s                                       T3
N_MAC        = 262      # cells/mm2 chronic baseline M1 density            T3
MAC_BAND     = 0.0364   # cm   M1 band depth, CALIBRATED so that the
                        #      quasi-steady balance reproduces the team's
                        #      36 nM wound-NO figure:
                        #        [NO] = (N_MAC/MAC_BAND) Q_CHRONIC / K_SCAV
                        #      Was 0.005 cm, which gives 262 nM and so
                        #      contradicted the 36 nM it was meant to match.
                        #      Not used by any current figure.             T4
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
    "M7 design":       dict(T1=1, T2=0, T3=0, T4=2, unplanned=2),   # f_rel,FGF T4 -> T1 (2026-10-16)
}
