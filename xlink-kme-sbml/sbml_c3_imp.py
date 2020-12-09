# %%
import IMP
import IMP.core
import IMP.atom
import numpy as np

import sbml_constants
import sbml_xl

# TODO: Read fasta with PMI, build coarse grained model

# %%
def get_random_reactivity():
    return 0.1#np.random.exponential(scale=0.1)

# normal distribution centered around 23 A
def get_kon(dist, mu=23, sigma=5):
    if dist > 35 or dist < 5:
        return 0
    return (1/(sigma*np.sqrt(2*np.pi)) * np.exp(-0.5*((dist-mu)/sigma)**2))

model = IMP.Model()
hroot_c3 = IMP.atom.read_pdb(
    "/home/kai/Coding/xlink_kme/xlink_kme/pdb/2A73.pdb",
    model,
    IMP.atom.CAlphaPDBSelector(),
)
hroot_c3b = IMP.atom.read_pdb(
    "/home/kai/Coding/xlink_kme/xlink_kme/pdb/2I07.pdb",
    model,
    IMP.atom.CAlphaPDBSelector(),
)

# IMP.atom.show_molecular_hierarchy(hroot)
# print(IMP.atom.get_chain_id(hroot))
# IMP.atom.Selection()

ps_c3 = IMP.atom.Selection(hroot_c3).get_selected_particles()
ps_c3b = IMP.atom.Selection(hroot_c3b).get_selected_particles()

# cs_c3 = IMP.core.XYZs(ps_c3)
# cs_c3b = IMP.core.XYZs(ps_c3b)

rb_c3 = IMP.core.RigidBody.setup_particle(IMP.Particle(model), ps_c3)
rb_c3b = IMP.core.RigidBody.setup_particle(IMP.Particle(model), ps_c3b)
print(rb_c3)

ps_lys_c3 = IMP.atom.Selection(
    hroot_c3, residue_type=IMP.atom.LYS
).get_selected_particles()
ps_lys_c3b = IMP.atom.Selection(
    hroot_c3b, residue_type=IMP.atom.LYS
).get_selected_particles()


# (p1,p2)->kon
calc_dist = set()
dist_list_c3 = []

for pA in ps_lys_c3:
    for pB in ps_lys_c3:
        if pA != pB:
            parts = sorted([pA, pB])
            parts_t = tuple(parts)
            if parts_t not in calc_dist:
                dist = IMP.core.get_distance(IMP.core.XYZ(pA), IMP.core.XYZ(pB))
                dist_list_c3.append(dist)
                calc_dist.add(parts_t)

dist_list_c3 = np.array(dist_list_c3)
# optionally plot dist_list with sns.distplot
# normally distributed but right-skewed
# print(kon_xl)
# %%
d_imp_ps = 'imp_ps'
def get_imp_ps_user_data_dict(ps):
    offset = 22
    res = IMP.atom.Residue(IMP.atom.Hierarchy(ps).get_parent())
    chain_id = IMP.atom.Chain(res.get_parent()).get_id()
    name_atom = str(res.get_residue_type()).replace('"', "")
    pos = res.get_index() + offset 
    user_data_dict = sbml_xl.get_user_data_dict(s_type=name_atom,s_pos=pos, s_prot=chain_id) 
    user_data_dict[d_imp_ps] = ps
    return user_data_dict
# %%

res = IMP.atom.get_residue(hroot_c3, 1)
IMP.atom.get_residue_type(res)


# %%
# IMP.atom.show_molecular_hierarchy(hroot_c3)

# %%
class AllXLReactionsIMP(sbml_xl.AllXLReactions):
    
    # def create_xl_trans(self):
        # return None
    
    def create_xl(self):
        return None
   
    def add_lys(self):
        lys_list = []
        offset = 22 # pdb to fasta offset
        for ps in self.params["ps_lys"]:
            user_data_dict = get_imp_ps_user_data_dict(ps)
            s = self.min_model.sbml_model.addSpecies(user_data_dict[sbml_constants.D_ID], 1, )
            s.setSpeciesType(user_data_dict[sbml_constants.D_TYPE])
            s.UserData = user_data_dict
            lys_list.append(s)
        return lys_list

    def add_lys_params(self):
        mono_trans_to_param_dict = {}
        for mono_trans in self.react_mono_trans.products:
            user_data_dict = sbml_xl.get_user_data_dict(s_type='klys', s_precursor_list=[mono_trans])
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[sbml_constants.D_ID], get_random_reactivity()
            ) # type: tesbml.libsbml.Parameter
            p_klys.UserData = user_data_dict 
            mono_trans_to_param_dict[mono_trans.getId()] = p_klys
        return mono_trans_to_param_dict

    def add_xl_trans_params(self):
        unique_id_kon_dict = {}
        for lys in self.species_lys:
            location_id_lys = lys.UserData[sbml_constants.D_LOCATION_ID]
            for mono in self.react_mono.products:
                location_id_mono = mono.UserData[sbml_constants.D_LOCATION_ID]
                if location_id_lys == location_id_mono:
                    continue
                user_data_dict = sbml_xl.get_user_data_dict(s_type="kon_xl", s_precursor_list=[lys, mono])
                location_id = user_data_dict[sbml_constants.D_LOCATION_ID]
                if location_id in unique_id_kon_dict:
                    continue
                ps_lys = lys.UserData[d_imp_ps]
                # super hacky -> make it nice
                ps_mono = mono.UserData[sbml_constants.D_PRECURSOR_LIST][0].UserData[sbml_constants.D_PRECURSOR_LIST][0].UserData[d_imp_ps]
                dist = IMP.core.get_distance(IMP.core.XYZ(ps_lys), IMP.core.XYZ(ps_mono))
                if dist > 5 and dist < 35:
                    p_kon_xl = self.min_model.sbml_model.addParameter(
                        user_data_dict[sbml_constants.D_ID],
                        get_kon(dist),
                        units="litre_per_mole_per_second",
                    )
                    p_kon_xl.UserData = user_data_dict
                    unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict

# %%
min_model = sbml_xl.MinimalModel(c_xl=111)
# %%
params_c3 = {
    "ps_lys": ps_lys_c3,
}
all_reactions = AllXLReactionsIMP(min_model, params_c3)
# %%
with open(
    "/home/kai/Coding/xlink_kme/xlink_kme/model_c3.xml", "w"
) as f:
    f.write(min_model.sbml_model.toSBML())


# %%
