import pandas as pd
import numpy as np
from xlink_kme_sbml.library import sbml_sim_helper as helper
from xlink_kme_sbml.library import sbml_constants as const
import altair as alt
import scipy.stats

"""
2023-09 Kai-Michael Kammer
Library of classes to visualize kinetic crosslink simulations
"""


def get_signi_plot_rr(df, min_val_abs=0, min_inc_abs=0, min_log2ratio_abs=0, exclude_imputed=False, y_label=None):
    if not y_label:
        y_label = const.S_VAR
    df = helper.get_rr_significance_columns(df, min_val_abs=min_val_abs, min_inc_abs=min_inc_abs,
                                            min_log2ratio_abs=min_log2ratio_abs)
    df = df[df[const.S_IS_SIGNIFICANT]]
    if exclude_imputed:
        df = df[~df[const.S_IMPUTED]]
    base = alt.Chart(df).mark_bar().encode(
        y=alt.Y(f'count({const.S_VAR})').title(y_label),
        x=alt.X(const.S_COND,
                sort=alt.EncodingSortField(field=const.S_IS_SIGNIFICANT, op="sum", order='descending')).title(""),
    )
    text = base.mark_text(dx=0, dy=-5, color='black', size=10).encode(
        text=alt.Text(f'sum({const.S_IS_SIGNIFICANT})'),
    )
    c = base + text
    c = c.properties(title="Suppression Experiment")
    c = c.configure_title(orient='top', anchor='middle', fontSize=14)
    c = prepare_for_print(c)
    # c = c.resolve_scale(
    # x='independent'
    # )
    return c


def get_xl_fraction_plot(df_xl_fraction, species, species_title, domain=None):
    df_filtered = df_xl_fraction[df_xl_fraction[const.S_TYPE].str.fullmatch(species)]
    y = alt.Y(const.S_VALUE).title(species_title)
    if domain:
        y = y.scale(domain=domain)
    c = alt.Chart(df_filtered).mark_line().encode(
        x=alt.X(const.S_EXP).title("Protein"),
        y=y,
        color=alt.Color(const.S_COND).title("Condition"),
    )
    c += c.mark_point()
    c = prepare_for_print(c)
    c = c.properties(width=200)
    return c


def get_area_plot_time(df, species="", val_col=const.S_VALUE, min_val=0, fullmatch=True, column=None,
                       color=const.S_LINK_TYPE, metric='sum'):
    df = df[df[const.S_VALUE] >= min_val]
    groupby_list = [const.S_TIME, color]
    if column:
        col = alt.Column(column)
        groupby_list.append(column)
    else:
        col = alt.Column()
    if fullmatch:
        df_filtered = df[df[const.S_TYPE].str.fullmatch(species)]
    else:
        df_filtered = df[df[const.S_TYPE].str.contains(species)]
    df_mean = df_filtered.groupby(groupby_list)[val_col].agg(
        ['mean', 'sum', 'median', 'std', 'count', ci95]).reset_index()
    c = alt.Chart(df_mean).mark_area(opacity=0.8).encode(
        y=alt.Y(metric).stack("normalize").title("Quantity Percentage"),
        color=color + ':N',
        x=const.S_TIME,
    )
    c = c.facet(
        column=col,
    ).resolve_scale(
        y='independent',
    )
    return c


def get_val_plot_time(df, species="", val_col=const.S_VALUE, min_val=0, fullmatch=True, column=None, color=const.S_EXP,
                      metric='mean', show_variation=True):
    df = df[df[const.S_VALUE] >= min_val]
    groupby_list = [const.S_TIME, color]
    if column:
        col = alt.Column(column)
        groupby_list.append(column)
    else:
        col = alt.Column()
    if fullmatch:
        df_filtered = df[df[const.S_TYPE].str.fullmatch(species)]
    else:
        df_filtered = df[df[const.S_TYPE].str.contains(species)]
    df_mean = df_filtered.groupby(groupby_list)[val_col].agg(
        ['mean', 'sum', 'median', 'std', 'count', ci95]).reset_index()
    df_mean['mean_p_ci'] = df_mean['mean'] + df_mean['ci95']
    df_mean['mean_m_ci'] = df_mean['mean'] - df_mean['ci95']
    c = alt.Chart(df_mean).mark_line().encode(
        y=metric,
        color=color + ':N',
        x=const.S_TIME,
    )
    if (metric == 'mean') and show_variation:
        c += c.mark_area(opacity=0.5).encode(y='mean_p_ci', y2='mean_m_ci')
    c = c.facet(
        column=col,
    ).resolve_scale(
        y='independent',
    )
    return c


def get_signi_plot_line(df, min_val=0, val_col=const.S_VALUE, metric='count'):
    df = df[df[val_col] >= min_val]
    species = [const.S_REACT_MONO_SUM, const.S_REACT_XL]
    species = '|'.join(species)
    df = df[df[const.S_TYPE].str.fullmatch(species)]

    c = alt.Chart(df).mark_point().encode(
        y=alt.Y(f'{metric}({val_col})').title(metric),
        x=alt.X(f'{const.S_EXP}').title("Protein"),
        color=const.S_TYPE,
    )
    return c


def get_signi_plot(df, min_val_abs=0, val_col=const.S_VALUE, show_signi_only=True, metric='count'):
    species = [const.S_REACT_MONO_SUM, const.S_REACT_XL]
    species = '|'.join(species)
    df = df[df[const.S_TYPE].str.fullmatch(species)]
    df = df.copy()

    df[const.S_IS_SIGNIFICANT] = False
    df.loc[df[val_col].abs() >= min_val_abs, const.S_IS_SIGNIFICANT] = True

    base = alt.Chart(df).mark_bar().encode(
        y=alt.Y(f'{metric}({val_col})').title(metric),
        x=const.S_TYPE,
    )
    bars = base.encode(
        color=const.S_IS_SIGNIFICANT,
        # tooltip=[str_var]
    )
    text_label = alt.Text(f'{metric}({val_col})').format('.0f')
    if show_signi_only:
        text = base.mark_text(dx=0, dy=-5, color='black', size=10).encode(
            text=text_label,
            detail=alt.Detail(const.S_IS_SIGNIFICANT),
        ).transform_filter(
            alt.FieldEqualPredicate(field=const.S_IS_SIGNIFICANT, equal=True)
        )
    else:
        text = base.mark_text(dx=0, dy=-5, color='black', size=10).encode(
            text=text_label,
        )
    c = bars + text
    c = c.facet(
        column=alt.Column(const.S_EXP).title('Protein'),
    )
    return c


def get_corr_plot(df, x, y, group_by, species=None, exp=None, facet=None, min_val_y=0, fit=True, name_x='', name_y='',
                  axis_format='.1e'):
    if not name_x:
        name_x = x
    if not name_y:
        name_y = y
    df = df[df[y] >= min_val_y]
    if species:
        species = '|'.join(species)
        df = df[df[const.S_VAR].str.contains(species, regex=True)]
    if exp:
        exp = '|'.join(exp)
        df = df[df[const.S_EXP].str.contains(exp, regex=True)]
    df = helper.get_corr_df(df, x, y, group_by=group_by)

    base = alt.Chart(df)
    c = base.mark_point().encode(
        x=alt.X(x, axis=alt.Axis(format=axis_format)).title(name_x),
        y=alt.Y(y).title(name_y).scale(domainMax=1.0),
        tooltip=[x, y],
    )
    if df[const.S_EXP].nunique() > 1:
        c = c.encode(
            color=alt.Color(f'{const.S_EXP}:N').sort(helper.get_sorted_cats_list(df[const.S_EXP])).title("Protein"),
        )
    regression = c.transform_regression(x, y, method='linear').mark_line()
    # regression_params = c.transform_regression(x, y, method='quad', params=True)
    # text='params:N'
    # ).transform_calculate(
    # params=f'"r²=" + format({alt.datum.rSquared},".2f") + "    y=" + round({alt.datum.coef}, 2)',
    # )
    params = base.mark_text(align="left", size=13).encode(
        x=alt.value(2),  # pixels from left
        y=alt.value(10),  # pixels from top
        text=alt.Text("corr:N"),
        detail='count()'
    ).transform_calculate(
        corr=f'"r²=" + format({alt.datum.corr ** 2},".2f")',
    )
    if fit:
        c += regression + params
    if x == const.S_NXL_SUM:
        c += c.mark_errorbar(extent='ci', opacity=0.5)
    if df[const.S_TYPE].nunique() > 1:
        c = c.facet(
            column=alt.Column(const.S_DISPLAY_NAME + ':N').title("Species").header(labelFontSize=12, titleFontSize=14).sort("descending")
        )
    c = c.resolve_scale(
        y='independent',
        x='independent'
    )
    c = update_display_name_species(c)
    c = prepare_for_print(c)
    return c


def update_display_name_species(c, target_str=const.S_TYPE):
    d = {const.S_REACT_MONO_SUM: const.S_REACT_DISPLAY_MONO,
         const.S_REACT_MONO_HYDRO: const.S_REACT_DISPLAY_MONO, const.S_REACT_XL: const.S_REACT_DISPLAY_XL}
    c = c.transform_calculate(
        display_name=f"{d}[datum.{target_str}]",
    )
    return c


def prepare_for_print(c):
    c = c.configure_axis(
        labelFontSize=12,
        titleFontSize=14,
    ).configure_legend(
        labelFontSize=12,
        titleFontSize=14,
    )
    return c


def get_suppress_plot(df, var_col=const.S_VAR, val_col=const.S_VALUE, color="", col="", min_val_abs=0, min_inc_abs=0,
                      min_log2ratio_abs=0, exclude_imputed=False):
    df = helper.get_rr_significance_columns(df, min_val_abs=min_val_abs, min_inc_abs=min_inc_abs,
                                            min_log2ratio_abs=min_log2ratio_abs)
    df = df[df[const.S_IS_SIGNIFICANT]]
    if exclude_imputed:
        df = df[~df[const.S_IMPUTED]]
    c = alt.Chart(df).mark_bar().encode(
        y=val_col,
        x=var_col,
    )
    if color:
        c = c.encode(
            color=color
        )
    if col:
        c = c.facet(
            column=alt.Column(col, sort=alt.EncodingSortField(field=var_col, op="count", order='descending')),
        )
        c = c.resolve_scale(
            x='independent'
        )
    return c


def get_connectivity_plot(df, min_val=0, val_col=const.S_VALUE, show_chains=True):
    df = helper.get_xl_connectivity_count_df(df, min_val=min_val, val_col=val_col)
    df = df[df[const.S_COUNT] > 0]
    groupby_list = [const.S_EXP]
    if show_chains:
        x_val = const.S_CHAIN_ID
        groupby_list.append(const.S_CHAIN_ID)
    else:
        x_val = const.S_EXP
    df[const.S_COUNT] = df.groupby(groupby_list)[const.S_COUNT].transform('mean')
    df = df.drop(labels=[const.S_POS], axis=1)
    df = df.drop_duplicates()
    # print(df)
    # return
    c = alt.Chart(df).mark_bar().encode(
        y=f'mean({const.S_COUNT})',
        x=x_val,
        color=const.S_IS_SIGNIFICANT
    )
    if show_chains:
        c = c.facet(
            column=const.S_EXP,
        ).resolve_scale(
            x='independent'
        )
    # c += c.mark_errorbar(extent='ci', opacity=1)
    return c


def get_val_sum_plot(df, species="", val_col=const.S_VALUE, min_val_abs=0, fullmatch=True):
    df = df.copy()
    df[const.S_IS_SIGNIFICANT] = False
    df.loc[df[val_col].abs() >= min_val_abs, const.S_IS_SIGNIFICANT] = True
    if fullmatch:
        df_filtered = df[df[const.S_TYPE].str.fullmatch(species)]
    else:
        df_filtered = df[df[const.S_TYPE].str.contains(species)]
    base = alt.Chart(df_filtered).mark_bar().encode(
        y=f'sum({val_col})',
        x=const.S_TYPE,
    )
    text = base.mark_text(dx=25, dy=0, color='white', size=10, angle=90).encode(
        # y=alt.Y(f'sum({str_is_significant})'),
        text=alt.Text(f'sum({val_col})', format='.4e'),
    )
    bars = base.encode(
        color=const.S_IS_SIGNIFICANT
    )
    c = bars + text
    c = c.facet(
        column=const.S_EXP,
    )
    return c


def get_val_percentage_plot(df, species="", val_col=const.S_VALUE, fullmatch=True, min_val=0):
    df = df.copy()
    df = df[df[val_col] >= min_val]
    if fullmatch:
        df_filtered = df[df[const.S_TYPE].str.fullmatch(species)]
    else:
        df_filtered = df[df[const.S_TYPE].str.contains(species)]
    base = alt.Chart(df_filtered).mark_bar().encode(
        y=alt.Y(f'sum({val_col})').stack('normalize').title('Quantity Percentage'),
    )
    text = base.transform_joinaggregate(
        yield_sum=f'sum({val_col})',
        groupby=[const.S_EXP]
    ).transform_joinaggregate(
        yield_sum_species=f'sum({val_col})',
        groupby=[const.S_EXP, const.S_TYPE]
    ).transform_calculate(
        yield_norm=f'datum.yield_sum_species / datum.yield_sum'
    ).mark_text(dx=0, dy=-5, color='black', size=10, angle=0).encode(
        # y=alt.Y(f'sum({str_is_significant})'),
        # text=alt.Text(f'sum({val_col})', format='.4e'),
        text=alt.Text('yield_norm:Q').format('.1%'),
        # detail=str_type,
    )
    bars = base.encode(
        color=const.S_TYPE,
    )
    c = bars + text
    c = c.facet(
        column=alt.Column(const.S_EXP).title('Protein'),
    )
    return c


def get_val_plot(df, species="", val_col=const.S_VALUE, min_val=0, fullmatch=True):
    df = df[df[val_col] >= min_val]
    if fullmatch:
        df_filtered = df[df[const.S_TYPE].str.fullmatch(species)]
    else:
        df_filtered = df[df[const.S_TYPE].str.contains(species)]
        y = val_col
        x = const.S_VAR
    c = alt.Chart(df_filtered).mark_bar().encode(
        y=val_col,
        x=const.S_VAR,
        color=const.S_TYPE
    )
    c = c.facet(
        column=const.S_EXP
    )
    return c


def ci95(x):
    return scipy.stats.sem(x, ddof=1) * 1.96  # use mean +(-) ci to get 95% conf interval


def get_explore_chart(df, var, ref_val=None, share_y_axis=False, metric='mean', min_val=0, x_axis_format="",
                      ref_val_label=None, val_col=const.S_VALUE, label_var="", label_val=""):
    if not label_var:
        label_var = var
    if not label_val:
        label_val = f'{metric}({val_col})'
    df = df[df[val_col] >= min_val]
    group_by = [const.S_EXP, const.S_LINK_TYPE, var]
    if 'cond' in df:
        group_by.append('cond')
    df_mean = df.groupby(group_by)[val_col].agg(['mean', 'sum', 'median', 'std', 'count', ci95]).reset_index()
    df_mean['mean_p_ci'] = df_mean['mean'] + df_mean['ci95']
    df_mean['mean_m_ci'] = df_mean['mean'] - df_mean['ci95']
    if ref_val_label:
        x_axis_title = [label_var, f'({ref_val_label:.2e})']
    else:
        x_axis_title = label_var
    c = alt.Chart(df_mean).mark_circle().encode(
        x=alt.X(var, scale=alt.Scale(base=10, type='log'), axis=alt.Axis(format=f"{x_axis_format}"),
                title=x_axis_title),
        # title=f'{str_species_crosslinker} Concentration/M'),
        y=alt.Y(f'{metric}', axis=alt.Axis(title=label_val), scale=alt.Scale(domainMin=0)),
        # y=alt.Y(f'{metric}({str_value})', axis=alt.Axis(title=f'{metric}')),
        color=alt.Color(f'{const.S_EXP}:N').title("Protein").sort(helper.get_sorted_cats_list(df_mean[const.S_EXP])),
    )
    # c_error_bars =
    if 'cond' in df:
        # row = alt.Row('cond', title=None, header=alt.Header(labelFontSize=14), sort=sorted_cats_list)
        c = c.mark_circle().encode(
            color=const.S_COND + ':N',
        )
        row = alt.Row()

    else:
        row = alt.Row()
    if share_y_axis:
        y_axis = 'shared'
    else:
        y_axis = 'independent'
    if metric == 'mean':
        c += c.mark_line() + c.mark_errorband(opacity=0.35, borders=True).encode(
            y=alt.Y('mean_p_ci').title(label_val),
            y2='mean_m_ci',
        )
    else:
        c += c.mark_line()
    if ref_val is not None:
        c_rule = alt.Chart(pd.DataFrame({var: [ref_val]})).mark_rule(strokeDash=[12, 6], size=2, color='grey').encode(
            x=var)
        c += c_rule

    c = c.facet(
        column=alt.Column(const.S_DISPLAY_NAME + ":N", title=None, header=alt.Header(labelFontSize=14)).sort("descending"),
        row=row
    ).resolve_scale(y=y_axis)
    c = prepare_for_print(c)
    c = update_display_name_species(c, target_str=const.S_LINK_TYPE)
    return c


def get_explore_chart_2d(df, x, y, min_val=0, z=None, metric='mean', log_x=True, log_y=True, format_x_scientific=True,
                         format_y_scientific=True, ref_val_x=None, ref_val_y=None, ref_val_label_x=None,
                         ref_val_label_y=None, val_col=const.S_VALUE, axis_title_x=None, axis_title_y=None,
                         legend_title=None):
    def _get_rule(var, val):
        c_rule = alt.Chart(pd.DataFrame({var: [val]})).mark_rule(strokeDash=[12, 6], size=2, color='grey', opacity=0.8)
        return c_rule

    df = df[df[val_col] >= min_val]
    group_by = [const.S_LINK_TYPE, x, y]
    if z is not None:
        group_by.append(z)
    df_mean = df.groupby(group_by)[val_col].agg(['mean', 'sum', 'median', 'std', 'count', ci95]).reset_index()
    df_mean['value_spread'] = 100 * df_mean['ci95'] / df_mean['mean']

    if not axis_title_x:
        axis_title_x = x
    axis_x = alt.Axis(title=axis_title_x)
    if ref_val_label_x:
        axis_title_x = [axis_title_x, f'({ref_val_label_x:.2e})']
    axis_science_x = alt.Axis(format='.1e', title=axis_title_x)

    if not axis_title_y:
        axis_title_y = y
    axis_y = alt.Axis(title=axis_title_y)
    if ref_val_label_y:
        axis_title_y = [axis_title_y, f'({ref_val_label_y:.2e})']
    axis_science_y = alt.Axis(format='.1e', title=axis_title_y)

    scale_log = alt.Scale(base=10, type='log')

    if not log_x and not format_x_scientific:
        ax_x = alt.X(x, axis=axis_x)
    if not log_y and not format_y_scientific:
        ax_y = alt.Y(y, axis=axis_y)

    if log_x and not format_x_scientific:
        ax_x = alt.X(x, scale=scale_log, axis=axis_x)
    if log_y and not format_y_scientific:
        ax_y = alt.Y(y, scale=scale_log, axis=axis_y)

    if not log_x and format_x_scientific:
        ax_x = alt.X(x, axis=axis_science_x)
    if not log_y and format_y_scientific:
        ax_y = alt.Y(y, axis=axis_science_y)

    if log_x and format_x_scientific:
        ax_x = alt.X(x, axis=axis_science_x, scale=scale_log)
    if log_y and format_y_scientific:
        ax_y = alt.Y(y, axis=axis_science_y, scale=scale_log)

    if not legend_title:
        legend_title = metric

    c = alt.Chart(df_mean).mark_square(opacity=1.0, size=75).encode(
        x=ax_x,
        y=ax_y,
        color=alt.Color(metric).legend(title=legend_title).scale(scheme='spectral'),
        tooltip=[x, y, metric]
    )
    c += c.mark_circle()  # stupid hack to combine the rule_chart (i.e. y or x ref) with the original chart and then be able to facet

    if ref_val_x is not None:
        c_rule = _get_rule(x, ref_val_x)
        c += c_rule.encode(x=x)
    if ref_val_y is not None:
        c_rule = _get_rule(y, ref_val_y)
        c += c_rule.encode(y=y)

    if z is not None:
        c = c.facet(
            column=alt.Column(const.S_DISPLAY_NAME + ':N').title(None).header(labelFontSize=14).sort("descending"),
            row=f'{z}:O'
        )
    else:
        c = c.facet(
            column=alt.Column(const.S_DISPLAY_NAME + ':N').title(None).header(labelFontSize=14).sort("descending"),
        )
    c = c.resolve_scale(
        color='independent'
    )
    c = update_display_name_species(c, target_str=const.S_LINK_TYPE)
    c = prepare_for_print(c)
    # title=(f'{x} ({ref_val_label_x})')

    return c


def get_single_supp_plot(df, var_col, val_col):
    c = alt.Chart(df).mark_bar().encode(
        y=alt.Y(val_col),
        x=alt.X(var_col).title(""),
        color=alt.Color(const.S_DISPLAY_NAME + ':N').title("Species").sort("descending"),
    )
    c = c.facet(
        column=alt.Column(val_col + '_sign').title("").header(labelFontSize=12, titleFontSize=14)
    ).resolve_scale(
        x='independent',
        y='independent'
    )
    c = c.properties(title=df['cond'].iloc[0].replace("S", "Suppression at"))
    c = c.configure_title(orient='top', anchor='middle', fontSize=14)
    c = prepare_for_print(c)
    c = update_display_name_species(c, target_str=const.S_TYPE)
    return c
