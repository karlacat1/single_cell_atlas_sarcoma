"""
Kaplan-Meier event-free survival by group.

Input: CSV with columns time_months, event (1 observed, 0 censored) and group.
Output: <prefix>_km.pdf and <prefix>_medians.csv in --output-dir.

Usage:
  python survival_km.py --data episodes.csv --output-dir results/ \
      --groups C2 C5 C6 C9 C10 C11 --prefix archetype
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts

DEFAULT_COLOR = '#bfbfbf'


def fit_groups(df, groups, time_col, event_col, group_col, colors):
    kmfs = []
    for group_name, group_df in df.groupby(group_col, observed=True):
        if len(group_df) == 0:
            continue
        if groups is not None and group_name not in groups:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(durations=group_df[time_col], event_observed=group_df[event_col],
                label=f'{group_name}')
        kmf.plot_color = colors.get(group_name, DEFAULT_COLOR)
        kmf.group_name = group_name
        kmfs.append(kmf)

    if not kmfs:
        raise ValueError('No groups left to plot.')
    return kmfs


def plot_km(kmfs, output_path, xlim, title):
    fig, ax = plt.subplots(figsize=(8, 6))

    for kmf in kmfs:
        kmf.plot_survival_function(ax=ax, ci_show=False, color=kmf.plot_color)

    ax.axhline(0.5, color='grey', ls=':', lw=1, alpha=0.7)

    # median survival drop-lines
    for kmf in kmfs:
        if kmf.plot_color == DEFAULT_COLOR:
            continue
        med = kmf.median_survival_time_
        if np.isfinite(med):
            ax.vlines(med, 0, 0.5, color=kmf.plot_color, ls='--', lw=1.5, alpha=0.9)
            ax.plot(med, 0.5, marker='o', color=kmf.plot_color, ms=6, zorder=5)
            ax.annotate(f'{med:.0f}', xy=(med, 0.5), xytext=(0, 8),
                        textcoords='offset points', ha='center', va='bottom',
                        fontsize=10, color=kmf.plot_color, fontweight='bold')

    # set ticks before adding the at-risk table
    ax.set_xticks(list(np.arange(0, xlim + 1, 12)))
    ax.set_xlim(0, xlim)

    add_at_risk_counts(*kmfs, ax=ax, fontsize=10, rows_to_show=['At risk'])

    ax.set_xlabel('Time (months)')
    ax.set_ylabel('Event-free survival probability')
    ax.set_title(title)
    ax.set_ylim(0, 1.05)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Kaplan-Meier event-free survival by group.')
    parser.add_argument('--data', required=True, help='CSV with time_months, event and group')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--prefix', default='survival')
    parser.add_argument('--groups', nargs='+', default=None,
                        help='Groups to include. Omit to include all.')
    parser.add_argument('--time-col', default='time_months')
    parser.add_argument('--event-col', default='event')
    parser.add_argument('--group-col', default='group')
    parser.add_argument('--colors', default=None,
                        help='Optional CSV with columns group and color')
    parser.add_argument('--xlim', type=float, default=60)
    parser.add_argument('--title', default='Kaplan-Meier — event-free survival')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.data)
    for col in (args.time_col, args.event_col, args.group_col):
        if col not in df.columns:
            raise ValueError(f'Column not found in {args.data}: {col}')

    colors = {}
    if args.colors is not None:
        color_df = pd.read_csv(args.colors)
        colors = dict(zip(color_df['group'].astype(str), color_df['color']))

    kmfs = fit_groups(df, args.groups, args.time_col, args.event_col,
                      args.group_col, colors)

    plot_km(kmfs, os.path.join(args.output_dir, f'{args.prefix}_km.pdf'),
            args.xlim, args.title)

    medians = pd.DataFrame([{
        'group': k.group_name,
        'n': int(k.event_observed.shape[0]),
        'n_events': int(k.event_observed.sum()),
        'median_survival_months': k.median_survival_time_ if np.isfinite(k.median_survival_time_) else np.nan,
    } for k in kmfs])
    medians.to_csv(os.path.join(args.output_dir, f'{args.prefix}_medians.csv'), index=False)
    print(medians.to_string(index=False))


if __name__ == '__main__':
    main()
