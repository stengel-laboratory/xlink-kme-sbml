import pandas as pd
import altair as alt
import numpy as np
#import seaborn as sns
import xlink_kme_sbml.library.sbml_constants as const
import xlink_kme_sbml.library.sbml_sim_helper as helper
import link_library.plot_library as plib
alt.data_transformers.disable_max_rows()

"""Kai Kammer - 2020-08. 
Library of classes to visualize kinetic crosslink simulations
"""



class PlotMasterVariableExplorer(object):
    def __init__(self, df, rr, out_dir='plots'):
        self.df = df
        self.rr = rr
        self.var_name = helper.get_explored_var_name(df)
        if not self.var_name:
            print("Warning: No variable to explore found in input. Exiting")
            exit(1)
        self.out_dir = out_dir
        self.df = helper.add_reactions_constants(self.df, self.rr)
        self.df = helper.add_initial_concentration(self.df, self.rr)
        # print(self.df[const.S_LINK_TYPE])
        if const.S_REACT_XL in self.df[const.S_LINK_TYPE].values:
            print("Add nxl and connected rate constants")
            self.df = helper.add_nxl(self.df)
            self.df = helper.add_max_connected_rate_constant(self.df)

    def plot_all(self):
        def _save(plot, name):
            print(f"Saving {name} plot")
            plib.save_g(plot, f"kme_explore_{self.df[const.S_EXP][0]}_{self.var_name}_{name}", [self.df],
                        out_dir=self.out_dir)
        plot_var = self.get_var_plot()
        _save(plot_var, 'var')
        plot_corr_rate = self.get_corr_plot(const.S_K_GENERIC)
        _save(plot_corr_rate, 'corr_rate')
        plot_corr_nxl = self.get_corr_plot(const.S_NXL_SUM)
        _save(plot_corr_nxl, 'corr_nxl')
        plot_corr_yield_mono = self.get_explorer_plot(const.S_REACT_MONO_HYDRO)
        _save(plot_corr_yield_mono, 'yield_mono')
        plot_corr_yield_xl = self.get_explorer_plot(const.S_REACT_XL)
        _save(plot_corr_yield_xl, 'yield_xl')

    def get_corr_df(self, x, agg='mean'):
        df_reg = self.df.groupby([x, const.S_LINK_TYPE, self.var_name])[[const.S_VALUE]].agg(agg).reset_index()
        df_reg = (df_reg.groupby([const.S_LINK_TYPE, self.var_name]).corr()).reset_index().drop(
            columns=['level_2', x]).drop_duplicates(subset=[const.S_LINK_TYPE, self.var_name]).rename(
            columns={const.S_VALUE: 'corr'})
        df = pd.merge(self.df, df_reg, on=[const.S_LINK_TYPE, self.var_name])
        return df

    def get_corr_plot(self, x, agg='mean'):
        y = const.S_VALUE
        df = self.get_corr_df(x, agg)
        base = alt.Chart(df)
        c = base.mark_point().encode(
            x=alt.X(x),
            y=f'{agg}({y})',
            tooltip=[x, f'{agg}({y})', self.var_name]
        )
        regression = c.transform_regression(x, y).mark_line()
        params = base.mark_text(align="left", size=13).encode(
            x=alt.value(2),  # pixels from left
            y=alt.value(10),  # pixels from top
            text=alt.Text("corr:N"),
            detail='count()'
        ).transform_calculate(
            corr=f'"r=" + format(datum.corr,".2f")',
        )
        if x == const.S_NXL_SUM:
            c += c.mark_errorbar(extent='ci', opacity=0.5)
        c += regression + params
        c = c.facet(
            row=const.S_LINK_TYPE,
            column=self.var_name
        ).resolve_scale(y='independent', x='independent')
        return c

    def get_var_plot(self, agg='mean'):
        c = alt.Chart(self.df).mark_point().encode(
            x=alt.X(self.var_name, scale=alt.Scale(base=10, type='log')),
            y=f'{agg}({const.S_VALUE})',
        )
        #c_error_bars =
        c += c.mark_line() + c.mark_errorbar(extent='ci', opacity=0.5)

        c = c.facet(column=const.S_LINK_TYPE).resolve_scale(y='independent')
        return c

    def get_explorer_plot(self, link_type, sort_by=const.S_K_GENERIC, log_y=False, independent_y=False):
        sort_by = alt.EncodingSortField(field=sort_by, order='descending')
        base = alt.Chart(self.df[self.df[const.S_LINK_TYPE] == link_type])
        bars = base.mark_bar().encode(
            x=f'{self.var_name}:N',
            y=alt.Y(const.S_VALUE),
            color=alt.Color(f'{self.var_name}:N', legend=alt.Legend(orient='top')),
        )
        if log_y:
            bars = bars.encode(
                y=alt.Y(const.S_VALUE, scale=alt.Scale(type='log', base=10)),
            )
        text_reaction_const = base.mark_text(dx=-0, y=10, color='black').encode(
            text=alt.Text(f'{const.S_K_GENERIC}:N'),
            detail='count()'  # const.S_UID_SHORT,
        ).transform_calculate(
            **{const.S_K_GENERIC: f'"k=" + format(datum.{const.S_K_GENERIC},".2e")'}
        )
        text_nxl1 = base.mark_text(dx=-0, y=25, color='black').encode(
            text=alt.Text(f'{const.S_NXL1}:N'),
            detail='count()'  # const.S_UID_SHORT,
        ).transform_calculate(
            **{const.S_NXL1: f'"{const.S_NXL1}=" + datum.{const.S_NXL1}'}
        )
        text_nxl2 = base.mark_text(dx=-0, y=40, color='black').encode(
            text=alt.Text(f'{const.S_NXL2}:N'),
            detail='count()'  # const.S_UID_SHORT,
        ).transform_calculate(
            **{const.S_NXL2: f'"{const.S_NXL2}=" + datum.{const.S_NXL2}'}
        )
        text_max_klys = base.mark_text(dx=-0, y=55, color='black').encode(
            # x=const.S_UID_SHORT,
            # const.S_VALUE,#const.S_UID_SHORT,
            text=alt.Text(f'{const.S_K_MAX_K_LYS}:N'),
            # detail encoding conflicts with condition
            # detail='count()',
            color=alt.condition(
                alt.FieldEqualPredicate(field=const.S_K_RATE_GT_K_LYS_MAX, equal=True),
                alt.value("green"),  # The positive color
                alt.value("red")  # The negative color
            )
        ).transform_calculate(
            # rate_const=alt.value('abc ') + alt.datum.rate_const,
            **{const.S_K_MAX_K_LYS: f'"max klys=" + format(datum.{const.S_K_MAX_K_LYS},".2e")'},
        )
        text_max_kon_xl = base.mark_text(dx=-0, y=70, color='black').encode(
            # x=const.S_UID_SHORT,
            # const.S_VALUE,#const.S_UID_SHORT,
            text=alt.Text(f'{const.S_K_MAX_K_ON_XL}:N'),
            # detail encoding conflicts with condition
            # detail=f'mean({const.S_UID_SHORT})',
            color=alt.condition(
                alt.FieldEqualPredicate(field=const.S_K_RATE_GT_K_ON_XL_MAX, equal=True),
                alt.value("green"),  # The positive color
                alt.value("red")  # The negative color
            )
        ).transform_calculate(
            # rate_const=alt.value('abc ') + alt.datum.rate_const,
            **{const.S_K_MAX_K_ON_XL: f'"max kon_xl=" + format(datum.{const.S_K_MAX_K_ON_XL},".2e")'}
        )
        if link_type == const.S_REACT_XL:
            c = alt.layer(bars, text_reaction_const, text_nxl1, text_nxl2, text_max_klys, text_max_kon_xl)
        else:
            c = alt.layer(bars, text_reaction_const, text_nxl1, text_max_kon_xl)
        c = c.facet(alt.Column(const.S_UID_SHORT, sort=sort_by), columns=14)
        if independent_y:
            c = c.resolve_scale(y='independent')
        return c


class PlotMasterTwoStates(object):
    def __init__(self, df_s1: pd.DataFrame, df_s2: pd.DataFrame, name_s1: str, name_s2: str, mapper_dict: dict,
                 df_xtract=None, df_dist=None, out_dir='plots'):
        self.name_s1 = name_s1
        self.name_s2 = name_s2
        self.df_s1_melt = helper.prepare_df(df_s1, self.name_s1, mapper_dict)
        self.df_s2_melt = helper.prepare_df(df_s2, self.name_s2, mapper_dict)
        self.df_log2 = self.get_log2_df()
        self.df_dist = df_dist
        self.df_xtract = df_xtract
        self.df_concat = pd.concat([self.df_s1_melt, self.df_s2_melt])
        self.out_dir = out_dir

        if self.df_xtract is not None:
            dfm_f = pd.merge(
                self.df_log2[[const.S_UID, const.S_LOG2RATIO, const.S_LINK_TYPE]],
                self.df_xtract[[const.S_UID, const.S_LOG2RATIO]],
                on=[const.S_UID],
                suffixes=[const.S_SUFFIX_SIM, const.S_SUFFIX_EXP],
            )
            dfm_rev = pd.merge(
                self.df_log2[[const.S_UID_REV, const.S_LOG2RATIO, const.S_LINK_TYPE]].rename(columns={const.S_UID_REV: const.S_UID}),
                self.df_xtract[[const.S_UID, const.S_LOG2RATIO]],
                on=[const.S_UID],
                suffixes=[const.S_SUFFIX_SIM, const.S_SUFFIX_EXP],
            )
            self.df_xtract_sim_merge = pd.concat([dfm_f, dfm_rev]).drop_duplicates().reset_index(drop=True)

    @staticmethod
    def plot_regression(df, x, y, col=None):
        c = (
            alt.Chart(df)
                .mark_point()
                .encode(x=alt.X(x), y=alt.Y(y))
        )
        regression = c.transform_regression(x, y).mark_line()
        params = (
            c.transform_regression(x, y, params=True)
                .mark_text(align="left")
                .encode(
                x=alt.value(20),  # pixels from left
                y=alt.value(20),  # pixels from top
                text=alt.Text("rSquared:N", format=".2e"),
            )
        )
        c += regression + params
        if col:
            c = c.facet(row=const.S_LINK_TYPE, column=col)
        else:
            c = c.facet(row=const.S_LINK_TYPE)
        return c.resolve_scale(x="independent", y="independent")

    def plot_yield_abs(self):
        if hasattr(self, 'var_name'):
            color = self.var_name
        else:
            color = const.S_EXP
        chart = alt.Chart(self.df_concat).mark_point(size=50).encode(
            x=const.S_VAR, y=alt.Y(const.S_VALUE), row=const.S_LINK_TYPE,
            color=alt.Color(color, legend=alt.Legend(orient='top'))
        ).resolve_scale(x="independent", y="independent")
        return chart

    def plot_log2_exp_vs_sim_regression(self, log2_cutoff_sim=0.5):
        fil_mono = (self.df_xtract_sim_merge[const.S_LINK_TYPE] == "MonoHydro") & (
                (self.df_xtract_sim_merge[const.S_LOG2RATIO + const.S_SUFFIX_SIM] > log2_cutoff_sim) | (
                self.df_xtract_sim_merge[const.S_LOG2RATIO + const.S_SUFFIX_SIM] < -log2_cutoff_sim)
        )
        fil_xl = self.df_xtract_sim_merge[const.S_LINK_TYPE] == "XL"
        df_sim_xtract_merge = self.df_xtract_sim_merge[fil_mono | fil_xl]
        chart = self.plot_regression(df_sim_xtract_merge, x=const.S_LOG2RATIO + const.S_SUFFIX_SIM,
                                     y=const.S_LOG2RATIO + const.S_SUFFIX_EXP)
        return chart

    def plot_log2ratio_by_link_type(self, log2_cutoff=1):
        chart = alt.Chart(
            self.df_log2[
                (self.df_log2[const.S_LOG2RATIO] > log2_cutoff)
                | (self.df_log2[const.S_LOG2RATIO] < -log2_cutoff)
                ]
        ).mark_bar().encode(
            x=const.S_VAR, y=alt.Y(const.S_LOG2RATIO, stack="zero"), column=const.S_LINK_TYPE,
        ).resolve_scale(
            x="independent", y="independent"
        )
        return chart

    def plot_log2ratio_vs_dist_delta(self):
        df_delta_dist = self.df_dist.groupby(const.S_UXID).apply(self.get_delta_dist).dropna()
        df_delta_dist = pd.DataFrame(df_delta_dist).reset_index()
        df_delta_dist = df_delta_dist.rename(columns={const.S_UXID: const.S_UID, 0: "delta_dist"})
        dfdd = pd.merge(
            df_delta_dist, self.df_xtract[[const.S_UID, const.S_LOG2RATIO, const.S_LINK_TYPE]], on=[const.S_UID]
        )
        chart = self.plot_regression(dfdd, x="delta_dist", y=const.S_LOG2RATIO)
        return chart

    def plot_yield_abs_sim_vs_exp(self):
        dfs = pd.merge(
            self.df_concat[[const.S_UID, const.S_VALUE, const.S_EXP, const.S_LINK_TYPE]],
            self.df_xtract[[const.S_UID, "ms1_area_sum", "ms1_area_sum_ref"]],
            on=[const.S_UID],
        )
        dfs_rev = pd.merge(
            self.df_concat[[const.S_UID_REV, const.S_VALUE, const.S_EXP, const.S_LINK_TYPE]].rename(
                columns={const.S_UID_REV: const.S_UID}),
            self.df_xtract[[const.S_UID, "ms1_area_sum", "ms1_area_sum_ref"]],
            on=[const.S_UID],
        )
        dfs = pd.concat([dfs, dfs_rev]).drop_duplicates().reset_index(drop=True)
        # dfs = dfs[~(dfs["value"] >= 0.99)].reset_index(drop=True)
        # assign correct ms1 areas to the corresponding simulation values
        mask_c3 = dfs[const.S_EXP] == self.name_s1
        dfs.loc[mask_c3, 'ms1_area_sum_exp'] = dfs['ms1_area_sum_ref']
        dfs.loc[~mask_c3, 'ms1_area_sum_exp'] = dfs['ms1_area_sum']
        chart = self.plot_regression(dfs, x=const.S_VALUE, y="ms1_area_sum_exp")
        return chart

    def get_log2_df(self):
        df_merge = pd.merge(
            self.df_s1_melt,
            self.df_s2_melt,
            on=[const.S_VAR, const.S_LINK_TYPE, const.S_UID, const.S_UID_REV],
            suffixes=["_" + self.name_s1, "_" + self.name_s2],
        )
        df_merge[const.S_LOG2RATIO] = np.log2(
            df_merge[f'{const.S_VALUE}_{self.name_s2}'] / df_merge[f'{const.S_VALUE}_{self.name_s1}'])
        df_merge = (
            df_merge.replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
        )
        return df_merge

    def get_delta_dist(self, x, metric="SASD"):
        exp_name = "exp_name"
        if len(x) == 2:
            d_ref_exp = x[x[exp_name] == self.name_s1][metric].values[0]
            d_exp = x[x[exp_name] == self.name_s2][metric].values[0]
            return d_exp - d_ref_exp
        else:
            return None

    def plot_log2_exp_vs_sim_abs(self):
        dfmm = pd.melt(
            self.df_xtract_sim_merge,
            value_vars=[const.S_LOG2RATIO + const.S_SUFFIX_EXP, const.S_LOG2RATIO + const.S_SUFFIX_SIM],
            id_vars=[const.S_UID, const.S_LINK_TYPE]
        )
        c = alt.Chart(dfmm).mark_point(size=50).encode(
            x=const.S_UID, y=alt.Y(const.S_VALUE, stack="zero"), color=const.S_VAR, row=const.S_LINK_TYPE
        ).resolve_scale(x="independent", y="independent")
        return c

    def plot_all(self):
        plot_sim_abs = self.plot_yield_abs()
        plib.save_g(plot_sim_abs, "kme_sim_yield_abs", [self.df_s1_melt, self.df_s2_melt], out_dir=self.out_dir)

        plot_sim_log2 = self.plot_log2ratio_by_link_type()
        plib.save_g(plot_sim_log2, "kme_sim_log2ratio", [self.df_log2], out_dir=self.out_dir)

        if self.df_xtract is not None:
            plot_sim_vs_exp_log2_abs = self.plot_log2_exp_vs_sim_abs()
            plib.save_g(plot_sim_vs_exp_log2_abs, "kme_sim_vs_exp_log2_abs", [self.df_xtract_sim_merge],
                        out_dir=self.out_dir)
            plot_sim_vs_exp_log2_regr = self.plot_log2_exp_vs_sim_regression()
            plib.save_g(plot_sim_vs_exp_log2_regr, "kme_sim_vs_exp_log2_regression", [self.df_xtract_sim_merge],
                        out_dir=self.out_dir)
            if "ms1_area_sum" in self.df_xtract:
                plot_sim_vs_exp_abs_yield = self.plot_yield_abs_sim_vs_exp()
                plib.save_g(plot_sim_vs_exp_abs_yield, "kme_sim_vs_exp_abs_yield", [self.df_concat, self.df_xtract],
                            out_dir=self.out_dir)
            if self.df_dist is not None:
                plot_exp_log2_vs_delta_dist = self.plot_log2ratio_vs_dist_delta()
                plib.save_g(plot_exp_log2_vs_delta_dist, "kme_exp_log2_vs_delta_dist_regression",
                            [self.df_dist, self.df_xtract],
                            out_dir=self.out_dir)


