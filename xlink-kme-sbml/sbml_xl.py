# %%
import simplesbml as sbml
import tesbml
import tesbml.libsbml
from typing import Callable, Iterator, Union, Optional, List, Tuple, Dict

d_name = "name"
d_type = "type"
d_id = "id"
d_location_id = "location_id"
d_position = "position"
d_protein = "protein"
d_precursor_list = "precursor_list"
d_location_list = "location_list"

s_str_lys = "LYS"


def get_user_data_dict(
    s_type, s_pos=None, s_prot=None, s_precursor_list=None, sort_precursors=True
):
    user_data_dict = {}
    user_data_dict[d_precursor_list] = s_precursor_list
    user_data_dict[d_type] = s_type
    if s_pos and not s_precursor_list:
        user_data_dict[d_location_list] = [(str(s_prot), s_pos)]
    elif s_precursor_list and not s_pos:
        for precursor in s_precursor_list:
            user_data_dict.setdefault(d_location_list, []).extend(
                precursor.UserData[d_location_list]
            )
    else:
        print("WARNING: Invalid combination of arguments")
    if sort_precursors:
        user_data_dict[d_location_list] = sorted(user_data_dict[d_location_list])
    s_loc_id = ""
    for location in user_data_dict[d_location_list]:
        if location[0] != "None":
            s_loc_id += f"{location[0]}_"  # protein
        s_loc_id += f"{location[1]}_"  # sequence postion
    s_loc_id = s_loc_id[:-1]
    user_data_dict[d_location_id] = s_loc_id
    user_data_dict[d_id] = f"{s_type}_{s_loc_id}"
    return user_data_dict


# %%
class XLReaction(object):
    """ 
    Template Class for Crosslinker Reactions
    Possible Reactions:
    1) 
        Multiple Species + Single Crosslinker <-> Multiple Transistion Monolinks
        Single kon <-> Single koff
        2 Reactants
    2) 
        Multiple Transistion Monolinks -> Multiple Monolinks
        Mutiple k, no backwards k
        1 Reactant
    3)
        Multiple Monolinks -> Multiple Hydrolyzed Monolinks
        Single kh, no backwards k
        1 Reactant
    4) 
        Multiple Monolinks x Multiple Species <-> Multiple Transistion Crosslinks
        Multiple kon <-> Single koff
        2 Reactants
    5) 
        Multiple Transition Crosslinks -> Multiple Crosslinks
        Multiple k, no backwards k
        1 Reactant
    """

    def __init__(
        self,
        reactants: List[tesbml.libsbml.Species],
        param_forward,
        reaction_string: str,
        sbml_model: sbml.sbmlModel,
    ):
        super().__init__()
        self.sbml_model = sbml_model
        self.reactants = reactants
        self.param_forward = param_forward
        self.reaction_string = reaction_string
        self.products = []  # type: List[tesbml.libsbml.Species]
        self.reactions = []  # type: List[tesbml.libsbml.Reaction
        self.__create_reactions()

    def __create_reactions(self):
        for reactants in self._get_reactants():
            prod = self._create_product(reactants)
            param_forward = self._get_param_forward(reactants)
            param_backward = self._get_param_backward(reactants, prod)
            expr_forward = self._get_expr_forward(reactants, param_forward)
            expr_backward = self._get_expr_backward(prod, param_backward)
            expr = expr_forward + expr_backward
            r = self.sbml_model.addReaction(
                reactants=self._get_reactants_for_model(reactants),
                products=[prod.getId()],
                expression=expr,
            )  # type: tesbml.libsbml.Reaction
            self.reactions.append(r.getKineticLaw())

    def _create_product(self, reactants) -> tesbml.libsbml.Species:
        product_type = self.reaction_string
        if isinstance(reactants, list):
            user_data_dict = get_user_data_dict(s_type=self.reaction_string, s_precursor_list=reactants)
        else:
            user_data_dict = get_user_data_dict(s_type=self.reaction_string, s_precursor_list=[reactants])
        product_loc_id = user_data_dict[d_location_id]
        product_name = f"{self.reaction_string}_{product_loc_id}"
        s = self.sbml_model.addSpecies(product_name, 0)
        s.UserData = user_data_dict
        s.setSpeciesType(product_type)
        self.products.append(s)
        return s

    def _create_reactions(self):
        pass

    def _get_reactants(self):
        pass

    def _get_reactants_for_model(self):
        pass

    def _get_param_forward(self, reactants):
        pass

    def _get_param_backward(self, reactants, product):
        pass

    def _get_expr_forward(self, reactants, param_forward):
        pass

    def _get_expr_backward(self, product, param_backward):
        return ""


class XLReactionMonoTrans(XLReaction):
    """  
    Class representing the reaction of
    Multiple Lysines -> Multiple Transition Monolinks; single kon, single koff
    Lys_1 + LinkerFree -> Mono_trans_1; kon*Lys_1*LinkerFree - koff*Mono_trans_1
    """

    def __init__(
        self,
        reactants: List[tesbml.libsbml.Species],
        param_forward: tesbml.libsbml.Parameter,
        reaction_string: str,
        sbml_model: sbml.sbmlModel,
        crosslinker: tesbml.libsbml.Species,
        param_backward: tesbml.libsbml.Parameter,
    ):
        print("INFO: Creating Mono Trans Reactions")
        self.crosslinker = crosslinker
        self.param_backward = param_backward
        super().__init__(reactants, param_forward, reaction_string, sbml_model)

    def _get_reactants(self):
        return self.reactants

    def _get_expr_forward(self, reactants, param_forward):
        return f"{param_forward.getId()}*{self.crosslinker.getId()}*{reactants.getId()}"

    def _get_expr_backward(self, product, param_backward):
        return f"-{param_backward.getId()}*{product.getId()}"

    def _get_param_forward(self, reactants):
        return self.param_forward

    def _get_param_backward(self, reactants, product):
        return self.param_backward

    def _get_reactants_for_model(self, reactants):
        return [reactants.getId(), self.crosslinker.getId()]


class XLReactionMono(XLReaction):
    """  
    Class representing the reaction of
    Multiple Transistion Monolinks -> Multiple Monolinks; multiple klys, no backwards reaction
    Mono_trans_1 -> Mono_1; klys_1*Mono_trans_1
    """

    def __init__(
        self,
        reactants: List[tesbml.libsbml.Species],
        param_forward: Dict[str, tesbml.libsbml.Parameter],
        reaction_string: str,
        sbml_model: sbml.sbmlModel,
    ):
        print("INFO: Creating Mono Reactions")
        self.param_lys_list = []
        super().__init__(reactants, param_forward, reaction_string, sbml_model)

    def _get_reactants(self):
        return self.reactants

    def _get_reactants_for_model(self, reactants):
        return [reactants.getId()]

    def _get_param_forward(self, reactants):
        return self.param_forward[reactants.getId()]

    def _get_expr_forward(self, reactants, param_forward):
        return f"{param_forward.getId()}*{reactants.getId()}"


class XLReactionMonoHydrolized(XLReaction):
    """  
    Class representing the reaction of
    Multiple Monolinks -> Multiple Hydrolized Monolinks; single kh, no backwards reaction
    Mono_1 -> Mono_Hydrolized_1; kh*Mono_1
    """

    def __init__(
        self,
        reactants: List[tesbml.libsbml.Species],
        param_forward: tesbml.libsbml.Parameter,
        reaction_string: str,
        sbml_model: sbml.sbmlModel,
    ):
        print("INFO: Creating Mono Hydrolized Reactions")
        super().__init__(reactants, param_forward, reaction_string, sbml_model)

    def _get_reactants(self):
        return self.reactants

    def _get_reactants_for_model(self, reactants):
        return [reactants.getId()]

    def _get_param_forward(self, reactants):
        return self.param_forward

    def _get_expr_forward(self, reactants, param_forward):
        return f"{param_forward.getId()}*{reactants.getId()}"


class XLReactionXLTrans(XLReaction):
    """  
    Class representing the reaction of
    Mutiple Monolinks x Multiple Lysines -> Multiple Transition Crosslinks; multiple kon_xl, single koff
    Lys_1 + Mono_3 -> XL_trans_1_3; kon_xl_1_3*Lys_1*Mono_3 - koff*XL_trans_1_3
    """

    def __init__(
        self,
        reactants: List[tesbml.libsbml.Species],
        param_forward: Dict[str, tesbml.libsbml.Parameter],
        reaction_string: str,
        sbml_model: sbml.sbmlModel,
        param_backward: tesbml.libsbml.Parameter,
        reactants_2: List[tesbml.libsbml.Species],
    ):
        print("INFO: Creating XL Trans Reactions")
        self.param_backward = param_backward
        self.reactants_2 = reactants_2
        self.param_zero = sbml_model.addParameter(
            "k_zero", 0
        )  # used for non-viable reactions
        super().__init__(reactants, param_forward, reaction_string, sbml_model)

    def _create_product(self, reactants) -> tesbml.libsbml.Species:
        user_data_dict = get_user_data_dict(
            self.reaction_string, s_precursor_list=reactants, sort_precursors=False
        )
        product_type = self.reaction_string
        product_loc_id = user_data_dict[d_location_id]
        product_name = f"{product_type}_{product_loc_id}"
        s = self.sbml_model.addSpecies(product_name, 0)
        s.UserData = user_data_dict
        s.setSpeciesType(product_type)
        self.products.append(s)
        return s

    def _get_expr_forward(self, reactants, param_forward):
        return f"{param_forward.getId()}*{reactants[0].getId()}*{reactants[1].getId()}"

    def _get_expr_backward(self, product, param_backward):
        return f"-{param_backward.getId()}*{product.getId()}"

    def _get_param_forward(self, reactants):
        user_data_dict = get_user_data_dict(self.reaction_string, s_precursor_list=reactants)
        location_id = user_data_dict[d_location_id]
        if location_id in self.param_forward:
            return self.param_forward[location_id]
        else:
            print(f"WARNING: No kon_xl found for {location_id}. Set to 0 instead")
            return self.param_zero

    def _get_param_backward(self, reactants, product):
        return self.param_backward

    def _get_reactants(self):
        tuple_pair_list = []  # type: List[Tuple[tesbml.libsbml.Species]]
        for react_1 in self.reactants:
            id_1 = react_1.UserData[d_location_id]
            for react_2 in self.reactants_2:
                id_2 = react_2.UserData[d_location_id]
                if id_1 != id_2:
                    if self._get_param_forward([react_1, react_2]) != self.param_zero:
                        tuple_pair_list.append((react_1, react_2))
        return tuple_pair_list

    def _get_reactants_for_model(self, reactants):
        return [reactants[0].getId(), reactants[1].getId()]


class XLReactionXL(XLReaction):
    """  
    Class representing the reaction of
    Multiple Transistion Crosslinks -> Multiple Crosslinks; multiple klys, no backwards reaction
    XL_trans_1_2 -> XL_1_2; klys_1*XL_trans_1_2
    XL_trans_2_1 -> XL_1_2; klys_2*XL_trans_2_1
    """

    def __init__(
        self,
        reactants: List[tesbml.libsbml.Species],
        param_forward: List[tesbml.libsbml.Parameter],
        reaction_string: str,
        sbml_model: sbml.sbmlModel,
    ):
        print("INFO: Creating XL Reactions")
        self.unique_id_prod_dict = {}  # type: Dict[str, tesbml.libsbml.Parameter]
        super().__init__(reactants, param_forward, reaction_string, sbml_model)

    def _create_product(self, reactants) -> tesbml.libsbml.Species:
        product_type = self.reaction_string
        user_data_dict = get_user_data_dict(self.reaction_string, s_precursor_list=[reactants])
        location_id = user_data_dict[d_location_id]
        if location_id in self.unique_id_prod_dict:
            return self.unique_id_prod_dict[location_id]
        product_name = f"{product_type}_{location_id}"
        s = self.sbml_model.addSpecies(product_name, 0)
        s.UserData = user_data_dict
        s.setSpeciesType(product_type)
        self.products.append(s)
        self.unique_id_prod_dict[location_id] = s
        return s

    def _get_reactants(self):
        return self.reactants

    def _get_reactants_for_model(self, reactants):
        return [reactants.getId()]

    def _get_param_forward(self, reactants):
        id_xl_trans = reactants.UserData[d_location_id]
        for p_lys in self.param_forward:
            location_id_klys = p_lys.UserData[d_location_id]
            for precursor in reactants.UserData[d_precursor_list]:
                if precursor.getSpeciesType() == s_str_lys:
                    if location_id_klys == precursor.UserData[d_location_id]:
                        return p_lys
        print("Warning: No matching klys found")
        return None

    def _get_expr_forward(self, reactants, param_forward):
        return f"{param_forward.getId()}*{reactants.getId()}"


# %%
class MinimalModel(object):
    def __init__(self, c_xl=1, kh=1e-6, kon=10e-4, koff=10e-1):
        self.sbml_model = sbml.sbmlModel()
        self.species_xl = []
        self.xl_xx = self.sbml_model.addSpecies("Crosslinker", c_xl)
        self.xl_xx.setSpeciesType("XL_FREE")

        self.xl_xh = self.sbml_model.addSpecies("CrosslinkerMonoHydrolized", 0)
        self.xl_xh.setSpeciesType("XL_MONO_HYDRO")

        self.xl_hh = self.sbml_model.addSpecies("CrosslinkerBiHydrolized", 0)
        self.xl_hh.setSpeciesType("XL_BI_HYDRO")

        self.p_kh = self.sbml_model.addParameter("kh", kh)

        self.sbml_model.addReaction(
            reactants=[self.xl_xx.getId()],
            products=[self.xl_xh.getId()],
            expression=f"{self.p_kh.getId()}*{self.xl_xx.getId()}",
        )
        self.sbml_model.addReaction(
            reactants=[self.xl_xh.getId()],
            products=[self.xl_hh.getId()],
            expression=f"{self.p_kh.getId()}*{self.xl_xh.getId()}",
        )  # type: tesbml.libsbml.Reaction

        # create custom unit: litre per mole per second: L/(mole*sec)
        # needed for 2nd order kinetics
        ud = self.sbml_model.getModel().createUnitDefinition()
        ud.setIdAttribute("litre_per_mole_per_second")
        u_m = ud.createUnit()
        u_m.setKind(tesbml.libsbml.UNIT_KIND_MOLE)
        u_m.setExponent(-1)
        u_m.setScale(0)
        u_m.setMultiplier(1)
        u_l = ud.createUnit()
        u_l.setKind(tesbml.libsbml.UNIT_KIND_LITRE)
        u_l.setExponent(1)
        u_l.setScale(0)
        u_l.setMultiplier(1)
        u_s = ud.createUnit()
        u_s.setKind(tesbml.libsbml.UNIT_KIND_SECOND)
        u_s.setExponent(-1)
        u_s.setScale(0)
        u_s.setMultiplier(1)
        self.p_kon = self.sbml_model.addParameter(
            "kon", kon, units="litre_per_mole_per_second"
        )  # type: tesbml.libsbml.Parameter
        self.p_koff = self.sbml_model.addParameter("koff", koff)


# %%
class AllXLReactions(object):
    def __init__(self, min_model: MinimalModel, params: Dict):
        self.min_model = min_model
        self.params = params
        self.species_lys = self.add_lys()
        self.react_mono_trans = self.create_mono_trans()
        self.mono_trans_id_to_param_dict = self.add_lys_params()
        self.react_mono = self.create_mono()
        self.react_mono_hydro = self.create_mono_hydro()
        self.xl_trans_dict = self.add_xl_trans_params()
        self.react_xl_trans = self.create_xl_trans()
        self.react_xl = self.create_xl()

    def create_mono_trans(self):
        return XLReactionMonoTrans(
            reactants=self.species_lys,
            param_forward=self.min_model.p_kon,
            reaction_string="MonoTrans",
            sbml_model=self.min_model.sbml_model,
            crosslinker=self.min_model.xl_xx,
            param_backward=self.min_model.p_koff,
        )
    def create_mono(self):
        return XLReactionMono(
            reactants=self.react_mono_trans.products,
            param_forward=self.mono_trans_id_to_param_dict,
            reaction_string="Mono",
            sbml_model=self.min_model.sbml_model,
        )
    def create_mono_hydro(self):
        return XLReactionMonoHydrolized(
            reactants=self.react_mono.products,
            param_forward=self.min_model.p_kh,
            reaction_string="MonoHydro",
            sbml_model=self.min_model.sbml_model,
        )
    def create_xl_trans(self):
        return XLReactionXLTrans(
            reactants=self.species_lys,
            param_forward=self.xl_trans_dict,
            reaction_string="XLTrans",
            sbml_model=self.min_model.sbml_model,
            param_backward=self.min_model.p_koff,
            reactants_2=self.react_mono.products,
        )
    def create_xl(self):
        return XLReactionXL(
            reactants=self.react_xl_trans.products,
            param_forward=self.react_mono.param_forward.values(),
            reaction_string="XL",
            sbml_model=self.min_model.sbml_model,
        )

    def add_lys(self) -> List[tesbml.libsbml.Species]:
        """
        Must return a list of species. The species should be created here.
        """
        lys_list = []
        return lys_list

    def add_lys_params(self, params=None) -> Dict[str, tesbml.libsbml.Parameter]:
        """
        Must return a dict mapping mono_trans.getId() to klys. The klys parameters should be created here.
        """
        mono_trans_to_param_dict = {}
        return mono_trans_to_param_dict

    def add_xl_trans_params(self, params=None) -> Dict[str, tesbml.libsbml.Parameter]:
        """
        Must return a dict mapping the location_id of xl_trans to kon_xl. The kon_xl parameters should be created here.
        """
        unique_id_kon_dict = {}
        return unique_id_kon_dict
