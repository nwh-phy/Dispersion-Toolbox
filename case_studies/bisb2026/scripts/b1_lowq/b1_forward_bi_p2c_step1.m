function b1_forward_bi_p2c_step1(run)
%B1_FORWARD_BI_P2C_STEP1  Record the detector-model choice before low-window fits.
root=bisb_find_project_root(fileparts(mfilename('fullpath')));
run=fullfile(root,'paper_results','b1_forward_bi_20261119','p2c',run);
ids={'Bi1107_10w','Bi1107_20w','Bi1117_10w','Bi1117_20w'};
paths=cell(size(ids));
models=cell(size(ids));
rows={};
messages={};
for j=1:numel(ids)
    id=ids{j};
    paths{j}=fullfile(run,'asymmetry',id);
    A=readtable(fullfile(paths{j},'asymmetry_fits.csv'),'TextType','string');
    D=readtable(fullfile(paths{j},'dispersion.csv'),'TextType','string');
    O=readtable(fullfile(run,'odd',id,'odd_bands.csv'));
    S=readtable(fullfile(run,'odd',id,'source_predictors.csv'),'TextType','string');
    f=@(n) A.held_chi2n(A.model==n);
    base=f("none"); ge=f("gainE"); sh=f("shift");
    % Freeze the empirical response below the B1 window only when its two
    % held-block scores improve by at least 0.01. Retain symmetric and shift
    % equivalents as detector-model sensitivity wherever feasible.
    if all(base-ge>=.01) && mean(ge)<=mean(sh)+.01
        models{j}='gainE';
    elseif all(base-sh>=.01)
        models{j}='shift';
    else
        models{j}='none';
    end
    q=D.qcut==2; e=D(q,:);
    v=A(A.model=="shift" & A.held_block==6,:);
    msg=sprintf(['%s: P2b held 6/1 %.4f/%.4f, odd gain %.4f/%.4f, shift proxy %.4f/%.4f; ' ...
        'delta0=%+.6g A^-1 delta1=%+.5f effective theta=%+.2f deg; ' ...
        'E0(0.001) %.1f--%.1f meV, E0(0.0025) %.1f--%.1f meV, ' ...
        'Delta envelope %.1f--%.1f meV (Q>=.001), rho %.1f--%.1f A; low-window primary %s.'], ...
        id,base,ge,sh,v.delta0_Ainv,v.delta1,v.effective_angle_deg,min(e.E0_at_001),max(e.E0_at_001), ...
        min(e.E0_at_0025),max(e.E0_at_0025),min(e.Delta_lo),max(e.Delta_hi),min(e.rho_lower),max(e.rho_lower),models{j});
    messages{end+1}=msg; %#ok<AGROW>
    for k=1:max(O.band)
        sub=O(O.band==k,:);
        g1=sub(sub.abs_q<=.001+eps,:); g2=sub(sub.abs_q>=.003-eps,:);
        qmid=groupsummary(sub,'block','median','raw_ratio');
        rows(end+1,:)={id,k,sub.lo_meV(1),sub.hi_meV(1),median(sub.raw_ratio), ...
            median(g1.raw_ratio),median(g2.raw_ratio),min(qmid.median_raw_ratio),max(qmid.median_raw_ratio), ...
            median(sub.ratio_over_elastic_area_ratio),median(sub.odd_sigma,'omitnan'), ...
            median(sub.far_R_1),median(sub.far_R_2)}; %#ok<AGROW>
    end
    m=S.held_block==6;
    messages{end+1}=sprintf('%s predictor hold-6 odd chi2/pair [%s] = [%s].',id, ...
        strjoin(S.model(m),', '),strjoin(string(round(S.held_odd_chi2_per_pair(m),4)),', ')); %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'id','band','lo','hi','median_R','median_R_q_le_001','median_R_q_ge_003', ...
    'block_median_R_min','block_median_R_max','R_over_elastic_area','median_residual_odd_sigma','far_R_055_060','far_R_065_085'});
out={fullfile(run,'step1_bands.csv'),fullfile(run,'step1_choice.csv'),fullfile(run,'step1_conclusion.log')};
assert(~any(cellfun(@isfile,out)),'Step 1 conclusion already exists.');
writetable(T,out{1});
writetable(table(string(ids(:)),string(models(:)),string(paths(:)), ...
    'VariableNames',{'id','model','fit_folder'}),out{2});
fid=fopen(out{3},'w');
assert(fid>=0);
nowUTC=datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z''');
fprintf(fid,'Step 1 completed %s.\n',char(nowUTC));
fprintf(fid,'%s\n',messages{:});
fprintf(fid,['Across acquisitions, a constant q offset alone and a static angular gain leave systematic odd residual. ' ...
    'An energy-dependent shift proxy and a phenomenological angular response both predict held blocks better. ' ...
    'The shift uses a derivative of the even P2b exit spectrum and yields an effective angle. ' ...
    'These data do not determine a physical tilt or detector gain across the measurement system. ' ...
    'Elastic wing and far-axis controls are tabulated in step1_bands.csv. Odd correction is fitted to the B1 material counts; ' ...
    'low-window response is frozen at 250 meV and the symmetric response is a separate sensitivity run. ' ...
    'P2b E0 nodes are conditional on geometry and q cut, with the region-A zigzag retained.\n']);
fclose(fid);
type(out{3});
end
