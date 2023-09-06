"""
2023-09 Kai-Michael Kammer
Constants used for kinetic model creation
"""
# REACTIONS
S_REACT_XL = "XL"
S_REACT_XL_TRANS = "XLTrans"
S_REACT_MONO_HYDRO = "MonoHydro"
S_REACT_MONO = "Mono"
S_REACT_MONO_SUM = "MonoSum"
S_REACT_MONO_TRANS = "MonoTrans"
S_REACT_MONO_TRANS_HYDRO = "MonoTransHydro"

S_REACT_DISPLAY_MONO = 'Mono-link'
S_REACT_DISPLAY_XL = 'Cross-link'

# KINETIC CONSTANTS
S_K_OFF = "koff"
S_K_ON = "kon"
S_K_HYDROLYSIS = "kh"
S_K_LYS = "klys"
S_K_ON_XL = 'kon_xl'
S_K_GENERIC = 'rate_const'
S_K_ZERO = 'k_zero'

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
S_MONOLNKER_HYDROLIZED = "MonolinkerHydrolized"
S_MONOLINKER = "Monolinker"
S_LYS = "LYS"
S_RATIO_CROSSLINKER_LYSINE = "linker/lys ratio"
S_MOLECULAR_WEIGHT = 'MolecularWeight'
S_NUM_LYS = 'NumLys'

# MISC
S_NONE = "None"
S_EXP = 'exp'
S_PROT = 'protein'
S_IMPUTED = 'isImputed'

# MODEL STRINGS
D_NAME = "name"
D_TYPE = "type"
D_ID = "id"
D_LOCATION_ID = "location_id"
D_PROTEIN_LIST = "protein_list"
D_PROTEIN = "prot"
D_CHAIN = "chain"
D_PRECURSOR_LIST = "precursor_list"
D_LOCATION_LIST = "location_list"
D_POSITION_LIST = "pos_list"
D_CONCENTRATION = 'conc'
D_MIN_DIST = 'minimum_xl_dist'
D_MAX_DIST = 'maximum_xl_dist'
D_REACTIVITY_DATA_MONO = 'reactivity_dict_mono'
D_REACTIVITY_DATA_XL = 'reactivity_dict_xl'

# PLOTS
S_EXP_NAME = 'UniProtID'
S_SIM_NAME = 'SimName'
S_UXID = "uxID"
S_SUFFIX_EXP = "_exp"
S_SUFFIX_SIM = "_sim"
S_SUFFIX_SCALED = "_scaled"
S_VALUE = "value"
S_VALUE_SCALED = S_VALUE + S_SUFFIX_SCALED
S_VALUE_GL = "value_g_per_l"
S_LOG2RATIO = "log2ratio"
S_REACTION_RATIO = "Reaction Ratio"
S_UID_REV = "uID_rev"
S_UID = "uID"
S_UID_SHORT = "uID_short"
S_LINK_TYPE = "link_type"
S_VAR = "variable"
S_NXL1 = 'nxl1'
S_NXL2 = 'nxl2'
S_NXL_SUM = 'nxl_sum'
S_COND = 'cond'

# SUPP_EXP
S_SUFFIX_REF = "_ref"
S_SUFFIX_CHANGE = "_change"
S_SUFFIX_INC_ABS = "_inc_abs"
S_VALUE_REF = S_VALUE + S_SUFFIX_REF
S_VALUE_SCALED_REF = S_VALUE_SCALED + S_SUFFIX_REF
S_K_GENERIC_REF = S_K_GENERIC + S_SUFFIX_REF
S_COND_REF = S_COND + S_SUFFIX_REF
S_VALUE_CHANGE = S_VALUE + S_SUFFIX_CHANGE
S_VALUE_INC_ABS = S_VALUE + S_SUFFIX_INC_ABS
S_VALUE_SCALED_CHANGE = S_VALUE_SCALED + S_SUFFIX_CHANGE
S_VALUE_SCALED_INC_ABS = S_VALUE_SCALED + S_SUFFIX_INC_ABS
S_LOG2RATIO_SCALED = S_LOG2RATIO + S_SUFFIX_SCALED
S_REACTION_RATIO_SCALED = S_REACTION_RATIO + S_SUFFIX_SCALED

var_supp_ref_list = [S_VALUE_REF, S_VALUE_SCALED_REF, S_K_GENERIC_REF, S_COND_REF]



# list of explorable variables
var_explore_list = [S_K_OFF, S_K_ON, S_K_HYDROLYSIS, S_CROSSLINKER]
