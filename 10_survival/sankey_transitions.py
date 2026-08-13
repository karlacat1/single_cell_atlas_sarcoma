"""
Sankey plot of group transitions across longitudinal timepoints.

Input: CSV with one row per sample, holding a sample id, a group label, an
episode date and an event indicator. Patient identity is taken from the first
--patient-id-parts underscore-separated fields of the sample id.

Output: <prefix>_sankey.pdf and <prefix>_transitions.csv in --output-dir.

Usage:
  python sankey_transitions.py --data metadata.csv --output-dir results/ \
      --group-col archetype_cluster_12 --prefix archetype
"""

import argparse
import os

import pandas as pd
import plotly.graph_objects as go

FALLBACK_COLOR = '#4c78a8'


def build_timepoints(df, sample_col, group_col, date_col, event_col,
                     missing_event_code, patient_id_parts, max_timepoints):
    """Assign consecutive timepoints per patient and pivot to one row per patient."""
    d = df[[sample_col, group_col, date_col, event_col]].copy()
    d.columns = ['sample_id', 'group', 'date_current', 'event']
    d = d.drop_duplicates(subset='sample_id').dropna(subset=['group'])

    d['patient_id'] = d['sample_id'].str.split('_').str[:patient_id_parts].str.join('_')

    d = d[d['event'] != missing_event_code]
    d['date_current'] = pd.to_datetime(d['date_current'], errors='coerce')
    d = d.dropna(subset=['date_current'])

    d = d.sort_values(['patient_id', 'date_current'])
    d['timepoint'] = d.groupby('patient_id').cumcount() + 1

    max_timepoints = min(max_timepoints, int(d['timepoint'].max()))
    d = d[d['timepoint'] <= max_timepoints]

    # only patients that can contribute a transition
    tp_per_patient = d.groupby('patient_id')['timepoint'].nunique()
    d = d[d['patient_id'].isin(tp_per_patient[tp_per_patient >= 2].index)]

    if d.empty:
        raise ValueError('No patients with at least two usable timepoints.')

    wide = d.pivot_table(index='patient_id', columns='timepoint',
                         values='group', aggfunc='first')
    print(f'{wide.shape[0]} patients with longitudinal samples, {len(wide.columns)} timepoints')
    return wide, d


def hex_to_rgba(h, alpha=0.4):
    h = h.lstrip('#')
    if len(h) == 3:
        h = ''.join(c * 2 for c in h)
    return f'rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{alpha})'


def main():
    parser = argparse.ArgumentParser(description='Sankey plot of group transitions.')
    parser.add_argument('--data', required=True, help='CSV of sample metadata')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--group-col', required=True)
    parser.add_argument('--prefix', default='transitions')
    parser.add_argument('--sample-col', default='Project #ID')
    parser.add_argument('--date-col', default='Date of current episode')
    parser.add_argument('--event-col',
                        default='event (1=progression/death, 0=zensiert, 9=data not available)')
    parser.add_argument('--missing-event-code', type=int, default=9)
    parser.add_argument('--patient-id-parts', type=int, default=3)
    parser.add_argument('--max-timepoints', type=int, default=2)
    parser.add_argument('--colors', default=None,
                        help='Optional CSV with columns group and color')
    parser.add_argument('--title', default='Group transitions across timepoints')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.data)
    for col in (args.sample_col, args.group_col, args.date_col, args.event_col):
        if col not in df.columns:
            raise ValueError(f'Column not found in {args.data}: {col}')

    wide, d = build_timepoints(df, args.sample_col, args.group_col, args.date_col,
                               args.event_col, args.missing_event_code,
                               args.patient_id_parts, args.max_timepoints)

    colors = {}
    if args.colors is not None:
        color_df = pd.read_csv(args.colors)
        colors = dict(zip(color_df['group'].astype(str), color_df['color']))

    timepoints = sorted(wide.columns)
    max_tp = max(timepoints)

    def tp_x(tp):
        return 0.001 if max_tp == 1 else 0.05 + 0.9 * (tp - 1) / (max_tp - 1)

    # nodes: heights proportional to patient count, each column centred
    groups_present = sorted(d['group'].unique())
    global_max = max(int(wide[tp].notna().sum()) for tp in timepoints)
    scale, pad = 0.9, 0.015

    node_map, labels, node_colors, xs, ys = {}, [], [], [], []
    for tp in timepoints:
        present = [g for g in groups_present if g in set(wide[tp].dropna())]
        counts = {g: int((wide[tp] == g).sum()) for g in present}
        heights = {g: counts[g] / global_max * scale for g in present}
        cursor = (1 - (sum(heights.values()) + pad * (len(present) - 1))) / 2
        for g in present:
            node_map[(tp, g)] = len(labels)
            labels.append(f'{g} (n={counts[g]})')
            node_colors.append(colors.get(g, FALLBACK_COLOR))
            xs.append(tp_x(tp))
            ys.append(cursor + heights[g] / 2)
            cursor += heights[g] + pad

    # links: patient flow between consecutive timepoints
    src, tgt, val, link_colors, records = [], [], [], [], []
    for k in timepoints:
        if (k + 1) not in wide.columns:
            continue
        pair = wide[[k, k + 1]].dropna()
        for (g0, g1), n in pair.groupby([k, k + 1]).size().items():
            src.append(node_map[(k, g0)])
            tgt.append(node_map[(k + 1, g1)])
            val.append(int(n))
            link_colors.append(hex_to_rgba(colors.get(g0, FALLBACK_COLOR), 0.4))
            records.append({'from_timepoint': k, 'to_timepoint': k + 1,
                            'from_group': g0, 'to_group': g1, 'n_patients': int(n)})

    fig = go.Figure(go.Sankey(
        arrangement='snap',
        node=dict(label=labels, color=node_colors, x=xs,
                  pad=15, thickness=18, line=dict(color='white', width=0.5)),
        link=dict(source=src, target=tgt, value=val, color=link_colors),
    ))

    tp_counts = {tp: int(wide[tp].notna().sum()) for tp in timepoints}
    fig.update_layout(
        title=args.title, height=700, width=950, font=dict(size=12),
        annotations=[dict(x=tp_x(tp), y=1.08, xref='paper', yref='paper',
                          showarrow=False, font=dict(size=13),
                          text=f'<b>T{tp}</b>  (n={tp_counts[tp]})')
                     for tp in timepoints],
    )
    fig.write_image(os.path.join(args.output_dir, f'{args.prefix}_sankey.pdf'), scale=3)

    transitions = pd.DataFrame(records)
    transitions.to_csv(os.path.join(args.output_dir, f'{args.prefix}_transitions.csv'), index=False)
    print(transitions.to_string(index=False))


if __name__ == '__main__':
    main()
