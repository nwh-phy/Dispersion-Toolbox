"""Recompute descriptive checks from the pasted CSV excerpts only.
Does NOT load raw EELS, run MATLAB, verify source SHA, or assess identifiability.
Usage: python recompute_review_tables.py (requires numpy and pandas).
"""
from pathlib import Path
import json, re
import numpy as np
import pandas as pd

def main() -> None:
    root = Path(__file__).resolve().parent
    ev, out = root / 'evidence', root / 'recomputed'
    out.mkdir(exist_ok=True)
    parent = pd.read_csv(ev/'boundary_summary.csv')
    assert len(parent) == 108 and not parent.duplicated(['key','N']).any()
    assert parent.groupby('key')['N'].apply(lambda s: sorted(s)==[1,2]).all()
    grouped = parent.groupby(['peak_model','N','status']).size().rename('count').reset_index()
    grouped = grouped.rename(columns={'N':'n_components'})
    grouped.to_csv(out/'parent_status_counts.csv',index=False)
    pp = parent.pivot(index='key',columns='N',values='normalized_sse')
    pp.columns = ['Q_n1','Q_n2']; pp['delta_Q']=pp.Q_n1-pp.Q_n2
    pp['Q2_le_Q1']=pp.Q_n2 <= pp.Q_n1+1e-12
    pp.to_csv(out/'parent_pairwise_objectives.csv')
    new = pd.read_csv(ev/'component_parameters.csv')
    assert len(new)==12 and not new.duplicated(['target_q','peak_model','N_components']).any()
    pairs=new.pivot(index=['target_q','peak_model'],columns='N_components',values='normalized_sse')
    pairs.columns=['Q_n1','Q_n2']; pairs['delta_Q']=pairs.Q_n1-pairs.Q_n2
    pairs['relative_reduction_percent']=100*pairs.delta_Q/pairs.Q_n1
    pairs.to_csv(out/'new_pairwise_objectives.csv')
    # Compare only the low-q N3 pair, whose center/membership is unchanged in the source reports.
    comparisons=[]
    for row in new[new.target_q == -.0025].itertuples():
        old=parent[(parent.key==f'N3_R1_W1800_{row.peak_model}') & (parent.N==row.N_components)].iloc[0]
        comparisons.append({'peak_model':row.peak_model,'n_components':row.N_components,
            'parent_Q':old.normalized_sse,'new_Q':row.normalized_sse,
            'absolute_difference':abs(old.normalized_sse-row.normalized_sse),
            'parameter_stability':'not_assessed_parameters_missing'})
    pd.DataFrame(comparisons).to_csv(out/'unchanged_lowq_objective_check.csv',index=False)
    bins=pd.read_csv(ev/'centered_bins.csv'); dq=.0005; qzero=236
    br=[]
    for row in bins.itertuples():
        members=np.array([int(x) for x in re.findall(r'\d+',str(row.source_channel))])
        qs=(members-qzero)*dq
        assert len(members)==row.N and np.all(np.diff(members)==1)
        assert abs(qs.mean()-row.actual_q)<1e-12
        assert abs(row.actual_q-row.target_q)<1e-12
        assert np.all(np.sign(qs)==np.sign(row.target_q))
        br.append({'target_q':row.target_q,'N_q':row.N,'center_native':members[len(members)//2],
            'q_left':qs[0]-dq/2,'q_right':qs[-1]+dq/2,'width':len(members)*dq,
            'center_check':True,'basis':'conditional on dq=0.0005 and native qzero=236; not instrumental response'})
    pd.DataFrame(br).to_csv(out/'centered_bin_geometry.csv',index=False)
    shifts=pd.read_csv(ev/'alignment_shifts.csv'); blocks=pd.read_csv(ev/'frame_block_summary.csv')
    assert len(shifts)==300 and np.array_equal(shifts.frame,np.arange(1,301))
    sh=shifts.A1_common_shift_pixels.to_numpy(); blockchecks=[]
    # Infer a common reference from each block median, then check all reported extrema.
    refs=[]
    for b in blocks.itertuples():
        v=sh[b.first_frame-1:b.last_frame]
        refs.append(b.zlp_median_pixel-float(np.median(v)))
    assert np.ptp(refs)==0
    reference=refs[0]
    assert np.median(sh)==0
    for b in blocks.itertuples():
        z=sh[b.first_frame-1:b.last_frame]+reference
        check=(np.median(z)==b.zlp_median_pixel and min(z)==b.zlp_min_pixel and max(z)==b.zlp_max_pixel)
        assert check
        blockchecks.append({'block':b.block,'n':len(z),'median_z':float(np.median(z)),
            'min_z':float(min(z)),'max_z':float(max(z)),'reported_values_match':bool(check)})
    pd.DataFrame(blockchecks).to_csv(out/'frame_shift_internal_consistency.csv',index=False)
    freq=pd.Series(sh).value_counts().sort_index().rename_axis('observed_offset_pixels').rename('count')
    freq.to_csv(out/'shift_histogram_counts.csv')
    summary={
        'scope':'Arithmetic of transcribed user tables; not independent raw/MAT or MATLAB validation',
        'parent_model_rows':len(parent), 'parent_pair_keys':len(pp),
        'parent_boundary_rows':int((parent.status=='boundary').sum()),
        'parent_nonboundary_converged_rows':int((parent.status=='converged').sum()),
        'parent_all_pairs_Q2_le_Q1':bool(pp.Q2_le_Q1.all()),
        'candidate_rows_expected_if_all_24_exist':54*2*24,
        'parameter_rows_expected_if_all_candidates_have_valid_shapes':54*24*(5+8),
        'selected_parameter_rows_expected':54*(5+8),
        'candidate_completeness_independently_verified':False,
        'new_model_rows':len(new),'new_boundary_rows':int(new.selected_boundary.sum()),
        'new_double_boundary_rows':int(new.loc[new.N_components==2,'selected_boundary'].sum()),
        'new_all_pairs_Q2_le_Q1':bool((pairs.delta_Q>=-1e-12).all()),
        'new_optimized_starts_expected_if_all_attempts_executed':12*24,
        'new_parameter_rows_expected_for_full_export':3*2*(1+2),
        'new_full_fit_arrays_provided':False,
        'same_center_geometry_rows_checked':len(br),
        'sequence_indices_checked':len(sh),'inferred_zlp_reference_pixel':reference,
        'observed_offset_min_pixels':int(sh.min()),'observed_offset_max_pixels':int(sh.max()),
        'observed_range_meV_conditional_on_parent_dE4':int((sh.max()-sh.min())*4),
        'lag1_pearson_sample_correlation_of_offsets':float(np.corrcoef(sh[:-1],sh[1:])[0,1]),
        'lag1_interpretation':'descriptive serial association only; may contain trend, not a stationary noise model or ESS estimate',
        'first_to_last_block_integrated_signal_drop_percent':float(100*(1-blocks.integrated_signal_mean.iloc[-1]/blocks.integrated_signal_mean.iloc[0])),
        'block_means_strictly_decreasing':bool(np.all(np.diff(blocks.integrated_signal_mean)<0)),
        'A1_transform_executed':False,
        'tests_reported_passed':int(pd.read_csv(ev/'p4p5_test_status.csv').passed.sum()),
        'matlab_tests_run_by_reviewer':False,
    }
    (out/'review_numbers.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    print(pairs.to_string())

if __name__=='__main__':
    main()
