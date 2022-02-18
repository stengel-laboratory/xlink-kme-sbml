import simplesbml.simplesbml as sbml
import tesbml
import tesbml.libsbml
import xlink_kme_sbml.library.sbml_constants as const
from typing import List, Tuple, Dict


def get_user_data_dict(
        s_type, s_pos=None, s_prot=None, s_precursor_list=None, sort_precursors=True
):
    user_data_dict = {const.D_PRECURSOR_LIST: s_precursor_list, const.D_TYPE: s_type}
    if s_pos and not s_precursor_list:
        user_data_dict[const.D_LOCATION_LIST] = [(str(s_prot), s_pos)]
    elif s_precursor_list and not s_pos:
        for precursor in s_precursor_list:
            user_data_dict.setdefault(const.D_LOCATION_LIST, []).extend(
                precursor.UserData[const.D_LOCATION_LIST]
            )
    else:
        print("WARNING: Invalid combination of arguments")
    if sort_precursors:
        user_data_dict[const.D_LOCATION_LIST] = sorted(user_data_dict[const.D_LOCATION_LIST])
    s_loc_id = ""
    for location in user_data_dict[const.D_LOCATION_LIST]:
        if location[0] != const.S_NONE:
            s_loc_id += f"{location[0]}_"  # protein
        s_loc_id += f"{location[1]}_"  # sequence postion
    s_loc_id = s_loc_id[:-1]
    user_data_dict[const.D_LOCATION_ID] = s_loc_id
    user_data_dict[const.D_ID] = f"{s_type}_{s_loc_id}"
    return user_data_dict


def get_location_id_dict(s_list: List[tesbml.libsbml.Species]):
    loc_id_dict = {}
    for s in s_list:
        loc_id_dict[s.UserData[const.D_LOCATION_ID]] = s
    return loc_id_dict


class XLReaction(object):
    """ 
    Template Class for Crosslinker Reactions
    Correct usage is to overwrite the methods provided by this class
    Currently defined reactions:
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
            sbml_model: sbml.SbmlModel,
            reactants_loc_id_to_products_dict=None,
    ):
        super().__init__()
        self.sbml_model = sbml_model
        self.reactants = reactants
        self.param_forward = param_forward
        self.reaction_string = reaction_string
        # should map a single or multiple (in a tuple) location ids to a product
        self.reactants_loc_id_to_products_dict = reactants_loc_id_to_products_dict  # type: Dict[str ,tesbml.libsbml.Species]
        self.products = []  # type: List[tesbml.libsbml.Species]
        self.reactions = []  # type: List[tesbml.libsbml.Reaction]
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
        if self.reactants_loc_id_to_products_dict:
            return self.reactants_loc_id_to_products_dict[reactants.UserData[const.D_LOCATION_ID]]
        else:
            product_type = self.reaction_string
            if isinstance(reactants, list):
                user_data_dict = get_user_data_dict(s_type=self.reaction_string, s_precursor_list=reactants)
            else:
                user_data_dict = get_user_data_dict(s_type=self.reaction_string, s_precursor_list=[reactants])
            product_loc_id = user_data_dict[const.D_LOCATION_ID]
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

    def _get_reactants_for_model(self, reactants):
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
            sbml_model: sbml.SbmlModel,
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
            sbml_model: sbml.SbmlModel,
            products=None,
    ):
        print("INFO: Creating Mono Reactions")
        super().__init__(reactants, param_forward, reaction_string, sbml_model, products)

    def _get_reactants(self):
        return self.reactants

    def _get_reactants_for_model(self, reactants):
        return [reactants.getId()]

    def _get_param_forward(self, reactants):
        # the location id of a lysine is the key in the param_forward dict
        return self.param_forward[reactants.UserData[const.D_LOCATION_ID]]

    def _get_expr_forward(self, reactants, param_forward):
        return f"{param_forward.getId()}*{reactants.getId()}"


class XLReactionMonoEff(XLReaction):
    """
    Class representing the reaction of
    Multiple Lysines -> Multiple Monolinks; multiple kon, no backwards reaction
    Lys_1 + LinkerFree -> Mono_1; keff1*Lys_1*LinkerFree
    This class includes diffusion: keff = (kon*klys1)/(koff+klys1)
    """

    def __init__(
            self,
            reactants: List[tesbml.libsbml.Species],
            param_forward: Dict[str, tesbml.libsbml.Parameter],
            reaction_string: str,
            sbml_model: sbml.SbmlModel,
            crosslinker: tesbml.libsbml.Species,
            param_kon: tesbml.libsbml.Parameter,
            param_koff: tesbml.libsbml.Parameter,
            products=None,
    ):
        print("INFO: Creating Mono Eff Reactions")
        self.crosslinker = crosslinker
        self.param_kon = param_kon
        self.param_koff = param_koff
        super().__init__(reactants, param_forward, reaction_string, sbml_model, products)

    def _get_reactants(self):
        return self.reactants

    def _get_param_forward(self, reactants):
        # the location id of a lysine is the key in the param_forward dict
        return self.param_forward[reactants.UserData[const.D_LOCATION_ID]]

    def __get_k_eff(self, param_forward):
        return f"(({param_forward.getId()}*{self.param_kon.getId()})/({param_forward.getId()}+{self.param_koff.getId()}))"

    def _get_expr_forward(self, reactants, param_forward):
        return f"{self.__get_k_eff(param_forward)}*{self.crosslinker.getId()}*{reactants.getId()}"

    def _get_reactants_for_model(self, reactants):
        return [reactants.getId(), self.crosslinker.getId()]


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
            sbml_model: sbml.SbmlModel,
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
            sbml_model: sbml.SbmlModel,
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
        product_loc_id = user_data_dict[const.D_LOCATION_ID]
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
        location_id = user_data_dict[const.D_LOCATION_ID]
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
            id_1 = react_1.UserData[const.D_LOCATION_ID]
            for react_2 in self.reactants_2:
                id_2 = react_2.UserData[const.D_LOCATION_ID]
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
            sbml_model: sbml.SbmlModel,
    ):
        print("INFO: Creating XL Reactions")
        self.unique_id_prod_dict = {}  # type: Dict[str, tesbml.libsbml.Parameter]
        super().__init__(reactants, param_forward, reaction_string, sbml_model)

    def _create_product(self, reactants) -> tesbml.libsbml.Species:
        product_type = self.reaction_string
        user_data_dict = get_user_data_dict(self.reaction_string, s_precursor_list=[reactants])
        location_id = user_data_dict[const.D_LOCATION_ID]
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
        for p_lys in self.param_forward:
            location_id_klys = p_lys.UserData[const.D_LOCATION_ID]
            for precursor in reactants.UserData[const.D_PRECURSOR_LIST]:
                if precursor.getSpeciesType() == const.S_LYS:
                    if location_id_klys == precursor.UserData[const.D_LOCATION_ID]:
                        return p_lys
        print("Warning: No matching klys found")
        return None

    def _get_expr_forward(self, reactants, param_forward):
        return f"{param_forward.getId()}*{reactants.getId()}"


class XLReactionXLEff(XLReactionXL):
    """
    Class representing the reaction of
    Mutiple Monolinks x Multiple Lysines -> Multiple Crosslinks; multiple kon_xl, no back backwards reaction
    Lys_1 + Mono_3 -> XL_1_3; keff_xl_1_3*Lys_1*Mono_3
    This class includes diffusion: keff_xl = (kon_xl*klys)/(koff+klys)
    """

    def __init__(
            self,
            reactants: List[tesbml.libsbml.Species],
            param_forward: Dict[str, tesbml.libsbml.Parameter], # kon_xl dict
            reaction_string: str,
            sbml_model: sbml.SbmlModel,
            param_klys: List[tesbml.libsbml.Parameter], # list of klys params
            param_koff: tesbml.libsbml.Parameter,
            reactants_2: List[tesbml.libsbml.Species],
    ):
        print("INFO: Creating XL Eff Reactions")
        self.reactants_2 = reactants_2
        self.param_klys = param_klys
        self.param_koff = param_koff
        self.param_zero = sbml_model.addParameter(
            "k_zero", 0
        )  # used for non-viable reactions
        self.unique_id_prod_dict = {}
        super().__init__(reactants, param_forward, reaction_string, sbml_model)

    def _create_product(self, reactants) -> tesbml.libsbml.Species:
        product_type = self.reaction_string
        user_data_dict = get_user_data_dict(self.reaction_string, s_precursor_list=reactants)
        location_id = user_data_dict[const.D_LOCATION_ID]
        if location_id in self.unique_id_prod_dict:
            return self.unique_id_prod_dict[location_id]
        product_name = f"{product_type}_{location_id}"
        s = self.sbml_model.addSpecies(product_name, 0)
        s.UserData = user_data_dict
        s.setSpeciesType(product_type)
        self.products.append(s)
        self.unique_id_prod_dict[location_id] = s
        return s
    # def _create_product(self, reactants) -> tesbml.libsbml.Species:
    #     product_type = self.reaction_string
    #     user_data_dict = get_user_data_dict(product_type, s_precursor_list=reactants)
    #     product_loc_id = user_data_dict[const.D_LOCATION_ID]
    #     product_name = f"{product_type}_{product_loc_id}"
    #     s = self.sbml_model.addSpecies(product_name, 0)
    #     s.UserData = user_data_dict
    #     s.setSpeciesType(product_type)
    #     self.products.append(s)
    #     return s

    def _get_param_forward(self, reactants):
        user_data_dict = get_user_data_dict(self.reaction_string, s_precursor_list=reactants)
        location_id = user_data_dict[const.D_LOCATION_ID]
        if location_id in self.param_forward:
            return self.param_forward[location_id]

        else:
            print(f"WARNING: No kon_xl found for {location_id}. Set to 0 instead")
            return self.param_zero

    def _get_reactants(self):
        tuple_pair_list = []  # type: List[Tuple[tesbml.libsbml.Species]]
        for react_1 in self.reactants:
            id_1 = react_1.UserData[const.D_LOCATION_ID]
            for react_2 in self.reactants_2:
                id_2 = react_2.UserData[const.D_LOCATION_ID]
                if id_1 != id_2:
                    if self._get_param_forward([react_1, react_2]) != self.param_zero:
                        tuple_pair_list.append((react_1, react_2))
        return tuple_pair_list

    def _get_expr_forward(self, reactants, param_forward):
        p_klys =  self.__get_k_lys(reactants)
        return f"{self.__get_k_eff(param_forward, p_klys)}*{reactants[0].getId()}*{reactants[1].getId()}"

    def _get_reactants_for_model(self, reactants):
        return [reactants[0].getId(), reactants[1].getId()]

    def __get_k_eff(self, param_kon_xl, param_klys):
        # param_forward is kon_xl
        return f"(({param_kon_xl.getId()}*{param_klys.getId()})/({param_klys.getId()}+{self.param_koff.getId()}))"

    def __get_k_lys(self, reactants):
        for p_lys in self.param_klys:
            location_id_klys = p_lys.UserData[const.D_LOCATION_ID]
            for species in reactants:
                if species.getSpeciesType() == const.S_LYS:
                    if location_id_klys == species.UserData[const.D_LOCATION_ID]:
                        return p_lys
        print("Warning: No matching klys found")
        return None


class MinimalModel(object):
    """
    Defines a minimal SBML model for a crosslinking reaction.
    Defined are the
    - crosslinker (and its initial concentration), its mono- and bihydrolized form
    - the hydrolysis constant
    - the kinetic constants kon and koff
    - the unit for kon in the form of L/(mole*sec)
    """
    def __init__(
            self,
            kinetic_params: Dict[str, float] = None,# params dict
            c_linker: float = 1):
        self.sbml_model = sbml.SbmlModel()
        self.c_linker = c_linker
        if kinetic_params is None:
            kinetic_params = self._get_default_params()
        self.kinetic_params = kinetic_params
        self.kh = kinetic_params[const.S_K_HYDROLYSIS]
        self._add_linker()
        self._create_second_order_unit_def()
        self._add_params()
        # self.sbml_model.getCompartment('c1').setUnits('dimensionless')

    def _get_default_params(self):
        return {const.S_K_HYDROLYSIS: 1e-6, const.S_K_ON: 10e-4, const.S_K_OFF: 10e-1}

    def _add_linker(self):
        self.xl_xx = self.sbml_model.addSpecies(const.S_CROSSLINKER, self.c_linker, )
        # self.xl_xx.setSpeciesType(S_CROSSLINKER)

        self.xl_xh = self.sbml_model.addSpecies(const.S_CROSSLINKER_MONO_HYDROLIZED, 0, )
        # self.xl_xh.setSpeciesType("XL_MONO_HYDRO")

        self.xl_hh = self.sbml_model.addSpecies(const.S_CROSSLINKER_BI_HYDROLIZED, 0, )
        # self.xl_hh.setSpeciesType("XL_BI_HYDRO")

        self.p_kh = self.sbml_model.addParameter(const.S_K_HYDROLYSIS, self.kinetic_params[const.S_K_HYDROLYSIS])

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

    def _create_second_order_unit_def(self):
        # create custom unit: litre per mole per second: L/(mole*sec)
        # needed for 2nd order kinetics
        ud = self.sbml_model.getModel().createUnitDefinition()
        ud.setIdAttribute(const.S_UNIT_LITRE_PER_MOLE_PER_SECOND)
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

    def _add_params(self):
        self.p_kon = self.sbml_model.addParameter(
            const.S_K_ON, self.kinetic_params[const.S_K_ON], units=const.S_UNIT_LITRE_PER_MOLE_PER_SECOND
        )  # type: tesbml.libsbml.Parameter
        self.p_koff = self.sbml_model.addParameter(const.S_K_OFF, self.kinetic_params[const.S_K_OFF])

    def add_compartments(self, n_comp):
        def _set_comp_id(obj, comp_name):
            obj.setId(obj.getName() + '_' + comp_name)

        def _update_obj_comp(obj, comp_name):
            if comp_name == 'c1':
                obj.setName(obj.getId())
                _set_comp_id(obj, comp_name)
                return None
            else:
                obj_clone = obj.clone()  # type: tesbml.libsbml.Species
                _set_comp_id(obj_clone, comp_name)
                return obj_clone

        def _update_reaction_comp(r, comp_name):
            def _updt(s_ref, f):
                s_name = s_ref.getSpecies()
                s_name_new = s_name + '_' + comp_name
                s_ref.setSpecies(s_name_new)
                if s_name in f:
                    f = f.replace(s_name, s_name_new)
                return f
            formula = r.getKineticLaw().getFormula()
            for r_ref in r.getListOfReactants():
                formula = _updt(r_ref, formula)
            for p_ref in r.getListOfProducts():
                formula = _updt(p_ref, formula)
            math = tesbml.libsbml.parseL3Formula(formula)
            r.createKineticLaw()
            r.getKineticLaw().setMath(math)

        new_species = []
        new_reactions = []
        for i in range(n_comp):
            self.sbml_model.addCompartment()
        for comp in self.sbml_model.getListOfCompartments():
            comp = comp.getId()
            for s in self.sbml_model.getListOfSpecies():
                s_new = _update_obj_comp(s, comp)
                if s_new:
                    s_new.setCompartment(comp)
                    new_species.append(s_new)
            for r in self.sbml_model.getListOfReactions():
                r_new = _update_obj_comp(r, comp)
                if r_new:
                    _update_reaction_comp(r_new, comp)
                    new_reactions.append(r_new)
        for s_new in new_species:
            self.sbml_model.model.addSpecies(s_new)
        for r in self.sbml_model.getListOfReactions():
            _update_reaction_comp(r, 'c1')
        for r_new in new_reactions:
            self.sbml_model.model.addReaction(r_new)


class MinimalModelMonolinker(MinimalModel):
    def __init__(
            self,
            kinetic_params: Dict[str, float] = None,# params dict
            c_linker: float = 1):
        super().__init__(kinetic_params, c_linker)

    def _add_linker(self):
        self.xl_xx = self.sbml_model.addSpecies(const.S_MONOLINKER, self.c_linker, )

        self.xl_xh = self.sbml_model.addSpecies(const.S_MONOLNKER_HYDROLIZED, 0, )

        self.p_kh = self.sbml_model.addParameter(const.S_K_HYDROLYSIS, self.kinetic_params[const.S_K_HYDROLYSIS])

        self.sbml_model.addReaction(
            reactants=[self.xl_xx.getId()],
            products=[self.xl_xh.getId()],
            expression=f"{self.p_kh.getId()}*{self.xl_xx.getId()}",
        )

class AllXLReactions(object):
    """
    This class needs to be overwritten to define the following functions:
    - add_lys(): define number and position of lysines
    - add_lys_params(): define the reactivity of each lysine
    - add_xl_trans_params(): define the reactivity of each crosslink
    This allows the above parameters to be supplied externally.
    It is also possible to overwrite the create_* methods with a pass state to avoid the formation of a species.
    I.E.: def create_xl_trans(self): pass means to no crosslinks will be formed
    """
    def __init__(self, min_model: MinimalModel, params: Dict):
        self.min_model = min_model
        self.params = params
        self.species_lys = self.add_lys()
        self.react_mono_trans = self.create_mono_trans()
        self.react_mono_hydro_trans = self.create_mono_hydro_trans()
        self.lys_id_to_param_dict = self.add_lys_params()
        self.react_mono = self.create_mono()
        self.react_mono_hydro = self.create_mono_hydro()
        self.mono_hydro_trans_prod_dict = {}
        for prod_mono_hydro in self.react_mono_hydro.products:
            loc_id = prod_mono_hydro.UserData[const.D_LOCATION_ID]
            self.mono_hydro_trans_prod_dict[loc_id] = prod_mono_hydro
        self.react_mono_hydro_alt = self.create_mono_hydro_alt()
        self.xl_trans_dict = self.add_xl_trans_params()
        self.react_xl_trans = self.create_xl_trans()
        self.react_xl = self.create_xl()

    def create_mono_trans(self):
        return XLReactionMonoTrans(
            reactants=self.species_lys,
            param_forward=self.min_model.p_kon,
            reaction_string=const.S_REACT_MONO_TRANS,
            sbml_model=self.min_model.sbml_model,
            crosslinker=self.min_model.xl_xx,
            param_backward=self.min_model.p_koff,
        )

    def create_mono_hydro_trans(self):
        return XLReactionMonoTrans(
            reactants=self.species_lys,
            param_forward=self.min_model.p_kon,
            reaction_string=const.S_REACT_MONO_TRANS_HYDRO,
            sbml_model=self.min_model.sbml_model,
            crosslinker=self.min_model.xl_xh,
            param_backward=self.min_model.p_koff,
        )

    def create_mono_hydro_alt(self):
        return XLReactionMono(
            reactants=self.react_mono_hydro_trans.products,
            param_forward=self.lys_id_to_param_dict,
            reaction_string=const.S_REACT_MONO_HYDRO,
            sbml_model=self.min_model.sbml_model,
            products=self.mono_hydro_trans_prod_dict,
        )

    def create_mono(self):
        return XLReactionMono(
            reactants=self.react_mono_trans.products,
            param_forward=self.lys_id_to_param_dict,
            reaction_string=const.S_REACT_MONO,
            sbml_model=self.min_model.sbml_model,
        )

    def create_mono_hydro(self):
        return XLReactionMonoHydrolized(
            reactants=self.react_mono.products,
            param_forward=self.min_model.p_kh,
            reaction_string=const.S_REACT_MONO_HYDRO,
            sbml_model=self.min_model.sbml_model,
        )

    def create_xl_trans(self):
        return XLReactionXLTrans(
            reactants=self.species_lys,
            param_forward=self.xl_trans_dict,
            reaction_string=const.S_REACT_XL_TRANS,
            sbml_model=self.min_model.sbml_model,
            param_backward=self.min_model.p_koff,
            reactants_2=self.react_mono.products,
        )

    def create_xl(self):
        return XLReactionXL(
            reactants=self.react_xl_trans.products,
            param_forward=self.react_mono.param_forward.values(),
            reaction_string=const.S_REACT_XL,
            sbml_model=self.min_model.sbml_model,
        )

    def add_lys(self) -> List[tesbml.libsbml.Species]:
        """
        Must return a list of species. The species should be created here.
        """
        lys_list = []
        return lys_list

    def add_lys_params(self) -> Dict[str, tesbml.libsbml.Parameter]:
        """
        Must return a dict mapping the location_id of a lys to klys. The klys parameters should be created here.
        """
        lys_to_param_dict = {}
        return lys_to_param_dict

    def add_xl_trans_params(self) -> Dict[str, tesbml.libsbml.Parameter]:
        """
        Must return a dict mapping the location_id of xl_trans (two lysine positions) to kon_xl.
        The kon_xl parameters should be created here.
        """
        unique_id_kon_dict = {}
        return unique_id_kon_dict


class AllXLReactionsNoDiff(AllXLReactions):
    """
    This class needs to be overwritten to define the following functions:
    - add_lys(): define number and position of lysines
    - add_lys_params(): define the reactivity of each lysine
    - add_xl_trans_params(): define the reactivity of each crosslink
    This allows the above parameters to be supplied externally.
    It is also possible to overwrite the create_* methods with a pass statement to avoid the formation of a species.
    I.E.: def create_xl_trans(self): pass means to no crosslinks will be formed
    """
    def __init__(self, min_model: MinimalModel, params: Dict):
        self.min_model = min_model
        self.params = params
        self.species_lys = self.add_lys()
        self.lys_id_to_param_dict = self.add_lys_params()
        self.react_mono = self.create_mono()
        self.react_mono_hydro = self.create_mono_hydro()
        self.mono_hydro_trans_prod_dict = {}
        for prod_mono_hydro in self.react_mono_hydro.products:
            loc_id = prod_mono_hydro.UserData[const.D_LOCATION_ID]
            self.mono_hydro_trans_prod_dict[loc_id] = prod_mono_hydro
        self.react_mono_hydro_alt = self.create_mono_hydro_alt()
        self.xl_trans_dict = self.add_xl_trans_params()
        self.react_xl = self.create_xl()

    def create_mono(self):
        return XLReactionMonoEff(
            reactants=self.species_lys,
            param_forward=self.lys_id_to_param_dict,
            reaction_string=const.S_REACT_MONO,
            sbml_model=self.min_model.sbml_model,
            crosslinker=self.min_model.xl_xx,
            param_kon=self.min_model.p_kon,
            param_koff=self.min_model.p_koff,
        )

    def create_mono_hydro(self):
        return XLReactionMonoHydrolized(
            reactants=self.react_mono.products,
            param_forward=self.min_model.p_kh,
            reaction_string=const.S_REACT_MONO_HYDRO,
            sbml_model=self.min_model.sbml_model,
        )

    def create_mono_hydro_alt(self):
        return XLReactionMonoEff(
            reactants=self.species_lys,
            param_forward=self.lys_id_to_param_dict,
            reaction_string=const.S_REACT_MONO_HYDRO,
            sbml_model=self.min_model.sbml_model,
            crosslinker=self.min_model.xl_xh,
            param_kon=self.min_model.p_kon,
            param_koff=self.min_model.p_koff,
            products=self.mono_hydro_trans_prod_dict,
        )

    def create_xl(self):
        return XLReactionXLEff(
            reactants=self.species_lys,
            param_forward=self.xl_trans_dict,
            reaction_string=const.S_REACT_XL,
            sbml_model=self.min_model.sbml_model,
            param_klys=self.react_mono.param_forward.values(),
            param_koff=self.min_model.p_koff,
            reactants_2=self.react_mono.products,
        )


class AllMonoReactionsNoDiff(AllXLReactionsNoDiff):
    """
    Class for creating monolinks only without explicit diffusion and only a combined effective reaction rate keff
    This class needs to be overwritten to define the following functions:
    - add_lys(): define number and position of lysines
    - add_lys_params(): define the reactivity of each lysine
    This allows the above parameters to be supplied externally.
    It is also possible to overwrite the create_* methods with a pass statement to avoid the formation of a species.
    """
    def __init__(self, min_model: MinimalModel, params: Dict):
        self.min_model = min_model
        self.params = params
        self.species_lys = self.add_lys()
        self.lys_id_to_param_dict = self.add_lys_params()
        self.react_mono = self.create_mono()
        self.mono_hydro_trans_prod_dict = {}

    def create_mono(self):
        return XLReactionMonoEff(
            reactants=self.species_lys,
            param_forward=self.lys_id_to_param_dict,
            reaction_string=const.S_REACT_MONO_HYDRO,
            sbml_model=self.min_model.sbml_model,
            crosslinker=self.min_model.xl_xx,
            param_kon=self.min_model.p_kon,
            param_koff=self.min_model.p_koff,
        )

    def create_mono_hydro(self):
        pass

    def create_mono_hydro_alt(self):
        pass

    def create_xl(self):
        pass