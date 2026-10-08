"""
params.py -- the single source of truth for every parameter used in the figures.

No value is hard-coded inside a solver; everything enters through this file, so
changing a number here changes every figure that depends on it.  Each entry
carries its provenance TIER, following Appendix "Parameter provenance":

    T1  measured in house
    T2  mechanistic correlation
    T3  named literature analogue, matching criterion stated
    T4  order-of-magnitude prior (>= 1 decade)

Units: cm, s, mol cm^-3.  1 uM = 1e-9 mol cm^-3.
"""

UM = 1e-9
DAY = 86400.0

# --- geometry ---------------------------------------------------------------
L_GEL        = 0.1      # cm   (1 mm) gel thickness                       design
L_GEL_THIN   = 500e-4   # cm   (500 um) thin-gel case                     design
L_TIS        = 0.3      # cm   (3 mm) dermal depth to the systemic sink   T3
DELTA        = 0.1      # cm   (1 mm) fibroblast-active / target depth    T4 (choice)
R_W_DEFAULT  = 0.5      # cm   (5 mm) wound radius                        design

# --- transport --------------------------------------------------------------
D0_V14       = 3.6e-6   # cm2/s free-solution diffusivity, 37 C           T2
D0_FGF       = 1.5e-6   # cm2/s                                           T2
A_V14_NM     = 0.87     # nm   hydrodynamic radius                        T2
A_FGF_NM     = 2.34     # nm                                              T2
A_FIBRE_NM   = 0.6      # nm   gelatin chain radius                       T3
XI_NM        = 10.0     # nm   mesh size (range 5-50 swept)               T3
EPS_GEL      = 0.90     # --   10% w/v GelMA at equilibrium hydration     T3
EPS_TIS      = 0.20     # --   tissue interstitial volume fraction        T3
K_P          = 1.0      # --   gel/tissue partition coefficient           T4
D_TIS_V14    = 1.2e-6   # cm2/s tissue diffusivity (range 0.7-1.8e-6)     T2
D_TIS_V14_LO = 0.7e-6
D_TIS_V14_HI = 1.8e-6

# --- loss rates -------------------------------------------------------------
K_DEG_GEL    = 1e-5     # 1/s  in-gel degradation (1e-5 to 1e-4)          T4
K_PROT_V14   = 1e-5     # 1/s  tissue proteolysis -- HIGHEST-PRIORITY GAP T4
K_PROT_RANGE = (1e-4, 1e-5, 1e-6)

# --- dosing -----------------------------------------------------------------
C_TARGET_V14 = 25 * UM     # mol/cm3  efficacious tissue concentration    T4
C_SOLUBILITY = 738 * UM    # mol/cm3  1 mg/mL saturation ceiling          T4 -- MEASURE THIS
MW_V14       = 1355.59     # Da
F_REL        = 0.5         # --  fraction of loaded dose released         T4

# --- nitric oxide (Learning Cycle 1) ---------------------------------------
D0_NO        = 3.3e-5   # cm2/s                                           T3
NO_HALFLIFE  = (3.0, 5.0, 10.0)   # s, report as a RANGE: autoxidation is
                                  # second order, so t_1/2 is concentration
                                  # dependent                             T3
INV_TAU      = (1.0, 0.5, 0.3)    # swept group 1/tau                     T2

# --- closure (Model 6) -- every entry here is Tier 4 ------------------------
D_N          = 1e-9     # cm2/s fibroblast random motility                T4
D_N_RANGE    = (1e-10, 1e-8)
R0           = 1.0 / DAY   # 1/s basal proliferation                      T4
D0_DEATH     = 0.3 / DAY   # 1/s basal death                              T4
K_CARRY      = 1.0         # normalised carrying capacity
W_SMOOTH     = 75e-4       # cm  IC smoothing width (50-100 um)

# chronic baseline and therapeutic pulse used in the durability section
R_CHRONIC    = 0.8 / DAY
D_CHRONIC    = 1.0 / DAY
R_PULSE      = 2.0 / DAY
D_PULSE      = 0.4 / DAY
TAU_RELAX    = 0.5 * DAY   # s, drug-effect relaxation time constant

# --- receptor / signalling --------------------------------------------------
A_LPS        = 10.0     # --  scaled endotoxin load [L]/KD_L              T4
GAMMA        = 10.0     # --  scaled receptor density                     T4
HILL_H       = 1.5      # --                                              T4
PSI_RANGE    = (1e-2, 1e4)

PROVENANCE = {
    # module: (T1, T2, T3, T4, no_planned_experiment)
    "M1 release":      dict(T1=2, T2=4, T3=3, T4=2, unplanned=1),
    "M2 tissue":       dict(T1=0, T2=3, T3=2, T4=4, unplanned=3),
    "M3 receptor":     dict(T1=1, T2=0, T3=0, T4=3, unplanned=2),
    "M4 signalling":   dict(T1=0, T2=0, T3=0, T4=5, unplanned=4),
    "M5 NO readout":   dict(T1=1, T2=0, T3=2, T4=3, unplanned=2),
    "M6 closure":      dict(T1=0, T2=0, T3=0, T4=11, unplanned=11),
}
