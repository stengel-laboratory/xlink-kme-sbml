#!/usr/bin/env python
from xlink_kme_sbml.library import sbml_constants as const, sbml_xl
import tellurium as te
import numpy as np


# %%
class AllXLReactionsExample(sbml_xl.AllXLReactionsImplicitDiffusionSimple):
#    def create_xl_trans(self):
#        pass
#
#    def create_xl(self):
#        pass

    def add_lys(self):
        lys_list = []
        no_lys = self.params["n_lys"]
        cnt = 0
        amount = self.params[const.D_CONCENTRATION]
        for i in range(no_lys):
            pos = i + 1
            cnt += 1
            prot = "Prot"
            chain = 'A'
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_LYS, s_pos=pos, s_prot=prot, s_chain=chain)
            s = self.min_model.sbml_model.addSpecies(user_data_dict[const.D_ID], amount, )
            s.setSpeciesType(const.S_LYS)
            s.UserData = user_data_dict
            lys_list.append(s)
            # if cnt == 1:
            #     amount *= 10
            #     cnt = 0
        return lys_list
    def add_lys_params(self, params=None):
        lys_to_params_dict = {}
        for lys in self.species_lys:
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_K_LYS, s_precursor_list=[lys])
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[const.D_ID], 0.1, #np.random.exponential(scale=0.1)
            ) # type: tesbml.libsbml.Parameter
            p_klys.UserData = user_data_dict 
            lys_to_params_dict[user_data_dict[const.D_LOCATION_ID]] = p_klys
        return lys_to_params_dict

    def add_xl_trans_params(self, params=None):
        unique_id_kon_dict = {}
        skip_percentage = (self.params["n_lys"]*2)/((self.params["n_lys"])**1.89) # use 2 instead of 1.9 for linear scale; 1.9 means that more monolinks lead to exponential more crosslinks
        for lys in self.species_lys:
            location_id_lys = lys.UserData[const.D_LOCATION_ID]
            for mono in self.react_mono.products:
                location_id_mono = mono.UserData[const.D_LOCATION_ID]
                if location_id_lys == location_id_mono:
                    continue
                if skip_percentage <= np.random.uniform():
                    continue
                user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_K_ON_XL, s_precursor_list=[lys, mono])
                location_id = user_data_dict[const.D_LOCATION_ID]
                if location_id in unique_id_kon_dict:
                    continue
                p_kon_xl = self.min_model.sbml_model.addParameter(
                    user_data_dict[const.D_ID],
                    5e8,#np.random.exponential(scale=5e8),
                    units=const.S_UNIT_LITRE_PER_MOLE_PER_SECOND,
                )
                p_kon_xl.UserData = user_data_dict
                unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict

    # def add_xl_trans_params(self, params=None):
    #     return {}


class MonoReactionsNoDiffExample(sbml_xl.AllMonoReactionsNoDiff):
    def add_lys(self):
        lys_list = []
        cnt = 0
        amount = 1
        for i in range(self.params["n_lys"]):
            pos = i + 1
            cnt += 1
            # prot = f"A{i + 1}"
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_LYS, s_pos=pos)#, s_prot=prot)
            s = self.min_model.sbml_model.addSpecies(user_data_dict[const.D_ID], amount, )
            s.setSpeciesType(const.S_LYS)
            s.UserData = user_data_dict
            lys_list.append(s)
            # if cnt == 1:
            #     amount *= 10
            #     cnt = 0
        return lys_list

    def add_lys_params(self, params=None):
        lys_to_params_dict = {}
        for lys in self.species_lys:
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_K_LYS, s_precursor_list=[lys])
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[const.D_ID], np.random.exponential(scale=0.01)
            ) # type: tesbml.libsbml.Parameter
            p_klys.UserData = user_data_dict
            lys_to_params_dict[user_data_dict[const.D_LOCATION_ID]] = p_klys
        return lys_to_params_dict

# %%

for no_lys in [20, 40, 60, 80, 100]:
    mw_kda = 1.65 * no_lys * 1e3  # the number of lysines times 1.65 roughly gives a molecular weight in kDa
    lys_conc = 2/mw_kda # 2 g/L
    xl_conc = no_lys/mw_kda # 1 g/L * no_lys
    min_model = sbml_xl.MinimalModel(kinetic_params={'kh': 2.5e-5, 'koff': 1e9, 'kon': 1e7}, c_linker=xl_conc)
    params = {
        const.D_CONCENTRATION: lys_conc,
        "n_lys": no_lys,
    }
    #min_model = sbml_xl.MinimalModelMonolinker(kinetic_params={'kh': 1e-5, 'koff': 1e-1, 'kon': 1e-3}, c_linker=25)
    # %%
    all_reactions = AllXLReactionsExample(min_model, params=params)
    #all_reactions = MonoReactionsNoDiffExample(min_model, params={"n_lys": 100})
    #all_reactions = MonoReactionsNoDiffExample(min_model, params={"n_lys": 25})

    # %%
    with open(
            f"model_random_same_react_{no_lys}.xml", "w"
    ) as f:
        f.write(min_model.sbml_model.toSBML())

    # rr = te.loadSBMLModel(min_model.sbml_model.toSBML())
    # with open(
    #         f"model_random_{no_lys}.txt", "w"
    # ) as f:
    #     f.write(rr.getAntimony())
# %%

# with open(
#         "/home/kai/Nextcloud-unikn/Documents/latex/kinetics/kinetic_models/model_120lys.txt", "w"
# ) as f:
#     f.write(rr.getAntimony())
#
# # %%
# rr.draw(layout="dot")
#
# # %%
# results = rr.simulate(0, 50000000)
# df_res = pd.DataFrame(results, columns=results.colnames)
# print("convergence", max(abs(rr.dv())))
#
#
# # %%
# df_res.tail()
#
# # %%

# %%