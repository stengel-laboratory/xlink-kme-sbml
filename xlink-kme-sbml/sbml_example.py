# %%
import sbml_xl
import tellurium as te
import pandas as p
import numpy as np
import pandas as pd

# %%
class AllXLReactionsExample(sbml_xl.AllXLReactions):
    def add_lys(self):
        lys_list = []
        for i in range(self.params["n_lys"]):
            pos = i + 1
            prot = "A"
            user_data_dict = sbml_xl.get_user_data_dict(s_type=sbml_xl.s_str_lys, s_pos=pos)
            s = self.min_model.sbml_model.addSpecies(user_data_dict[sbml_xl.d_id], 1)
            s.setSpeciesType(sbml_xl.s_str_lys)
            s.UserData = user_data_dict
            lys_list.append(s)
        return lys_list

    def add_lys_params(self, params=None):
        mono_trans_to_param_dict = {}
        for mono_trans in self.react_mono_trans.products:
            user_data_dict = sbml_xl.get_user_data_dict(s_type='klys', s_precursor_list=[mono_trans])
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[sbml_xl.d_id], np.random.exponential(scale=0.1)
            ) # type: tesbml.libsbml.Parameter
            p_klys.UserData = user_data_dict 
            mono_trans_to_param_dict[mono_trans.getId()] = p_klys
        return mono_trans_to_param_dict

    def add_xl_trans_params(self, params=None):
        unique_id_kon_dict = {}
        for lys in self.species_lys:
            location_id_lys = lys.UserData[sbml_xl.d_location_id]
            for mono in self.react_mono.products:
                location_id_mono = mono.UserData[sbml_xl.d_location_id]
                if location_id_lys == location_id_mono:
                    continue
                user_data_dict = sbml_xl.get_user_data_dict(s_type='kon_xl', s_precursor_list=[lys, mono])
                location_id = user_data_dict[sbml_xl.d_location_id]
                if location_id in unique_id_kon_dict:
                    continue
                p_kon_xl = self.min_model.sbml_model.addParameter(
                    user_data_dict[sbml_xl.d_id],
                    np.random.exponential(scale=0.1),
                    units="litre_per_mole_per_second",
                )
                p_kon_xl.UserData = user_data_dict
                unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict



# %%
min_model = sbml_xl.MinimalModel(c_xl=5)
# %%
all_reactions = AllXLReactionsExample(min_model, params={"n_lys":5})
# %%
with open(
    "/home/kai/Coding/xlink_kme/xlink_kme/model.xml", "w"
) as f:
    f.write(min_model.sbml_model.toSBML())


# %%
rr = te.loadSBMLModel(min_model.sbml_model.toSBML())

# %%
rr.draw(layout="dot")

# %%
results = rr.simulate(0, 50000000)
df_res = pd.DataFrame(results, columns=results.colnames)
print("convergence", max(abs(rr.dv())))


# %%
df_res.tail()

# %%
with open(
    "/home/kai/Coding/xlink_kme/xlink_kme/model_antimony.txt", "w"
) as f:
    f.write(rr.getAntimony())

# %%