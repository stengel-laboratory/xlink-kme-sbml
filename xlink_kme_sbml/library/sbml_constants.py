# REACTIONS
S_REACT_XL = "XL"
S_REACT_XL_TRANS = "XLTrans"
S_REACT_MONO_HYDRO = "MonoHydro"
S_REACT_MONO = "Mono"
S_REACT_MONO_TRANS = "MonoTrans"
S_REACT_MONO_TRANS_HYDRO = "MonoTransHydro"

# KINETIC CONSTANTS
S_K_OFF = "koff"
S_K_ON = "kon"
S_K_HYDROLYSIS = "kh"
S_K_LYS = "klys"
S_K_ON_XL = 'kon_xl'
S_K_GENERIC = 'rate_const'

S_K_MAX_K_LYS = 'rate_const_max_klys'
S_K_MAX_K_LYS_P1 = S_K_MAX_K_LYS + '_p1'
S_K_MAX_K_LYS_P2 = S_K_MAX_K_LYS + '_p2'
S_K_MAX_K_ON_XL = 'rate_const_max_kon_xl'
S_K_MAX_K_ON_XL_P1 = S_K_MAX_K_ON_XL + 'p1'
S_K_MAX_K_ON_XL_P2 = S_K_MAX_K_ON_XL + 'p2'
S_K_RATE_GT_K_LYS_MAX = 'rate_const_greater_than_klys_max'
S_K_RATE_GT_K_ON_XL_MAX = 'rate_const_greater_than_kon_xl_max'

# UNITS
S_UNIT_LITRE_PER_MOLE_PER_SECOND = "litre_per_mole_per_second"

# SPECIES
S_CROSSLINKER_BI_HYDROLIZED = "CrosslinkerBiHydrolized"
S_CROSSLINKER_MONO_HYDROLIZED = "CrosslinkerMonoHydrolized"
S_CROSSLINKER = "Crosslinker"
S_LYS = "LYS"

# MISC
S_NONE = "None"
S_EXP = 'exp'

# MODEL STRINGS
D_NAME = "name"
D_TYPE = "type"
D_ID = "id"
D_LOCATION_ID = "location_id"
D_POSITION = "position"
D_PROTEIN = "protein"
D_PRECURSOR_LIST = "precursor_list"
D_LOCATION_LIST = "location_list"
D_POSITION_LYSINE = "pos_lys"
D_CONCENTRATION = 'conc'

# PLOTS
S_EXP_NAME = 'UniProtID'
S_SIM_NAME = 'SimName'
S_UXID = "uxID"
S_SUFFIX_EXP = "_exp"
S_SUFFIX_SIM = "_sim"
S_VALUE = "value"
S_LOG2RATIO = "log2ratio"
S_UID_REV = "uID_rev"
S_UID = "uID"
S_UID_SHORT = "uID_short"
S_LINK_TYPE = "link_type"
S_VAR = "variable"
S_NXL1 = 'nxl1'
S_NXL2 = 'nxl2'
S_NXL_SUM = 'nxl_sum'

# list of explorable variables
var_explore_list = [S_K_OFF, S_K_ON, S_K_HYDROLYSIS, S_CROSSLINKER]
