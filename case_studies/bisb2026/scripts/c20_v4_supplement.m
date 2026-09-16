function c20_v4_supplement(out)
% Additional diagnostics derived only from already persisted v4 arrays.
dest=fullfile(out,'590_PL2_10w'); loaded=load(fullfile(dest,'sequence_block_spectra.mat')); f=loaded.frames;
rows={}; qc=readtable(fullfile(dest,'frame_qc.csv'));
for t=1:height(qc)
 [w0,pk0]=width0(f.native_E,f.reference_profiles_A0(:,t));
 [w1,pk1]=width0(f.E,f.reference_profiles_A1(:,t));
 rows(end+1,:)={t,w0,w1,pk0,pk1}; %#ok<AGROW>
end
writetable(cell2table(rows,'VariableNames',{'sequence','reference_zero_FWHM_A0_meV','reference_zero_FWHM_A1_meV','apex_A0_meV','apex_A1_meV'}),fullfile(dest,'ZLP_width_diagnostics.csv'));
v=qc.measured_offset_pixels; signal=qc.integrated_signal;
rawcorr=corrcoef(v(1:end-1),v(2:end)); dv=detrend(v); dcorr=corrcoef(dv(1:end-1),dv(2:end));
signal_corr=corrcoef(signal(1:end-1),signal(2:end)); ds=detrend(signal); signal_dcorr=corrcoef(ds(1:end-1),ds(2:end));
[w0,p0]=width0(f.native_E,sum(f.reference_profiles_A0,2)); [w1,p1]=width0(f.E,sum(f.reference_profiles_A1,2));
summary=struct('A0_L1_max_error',f.A0_identity_error,'valid_shift_frames',nnz(f.alignment.valid), ...
 'observed_offset_range_pixels',[min(v) max(v)],'offset_lag1',rawcorr(1,2), ...
 'linear_detrended_offset_lag1',dcorr(1,2),'signal_lag1',signal_corr(1,2), ...
 'linear_detrended_signal_lag1',signal_dcorr(1,2),'reference_ZLP_FWHM_A0_meV',w0, ...
 'reference_ZLP_FWHM_A1_meV',w1,'reference_ZLP_apex_A0_meV',p0,'reference_ZLP_apex_A1_meV',p1, ...
 'A1_absolute_reference','rounded median ZLP peak; parent energy axis retained, no global recentering', ...
 'interpretation','quantized sequence diagnostics, not resolution deconvolution, iid test or effective sample size');
c20_v4_io('json',fullfile(dest,'A0_A1_diagnostics.json'),summary);
brows={};
for j=1:3
 mask=f.E>=300&f.E<=1800; mean_shape=mean(f.block_per_frame_mean_A0(mask,:,j),2);
 for b=1:6
  y=f.block_per_frame_mean_A0(mask,b,j);
  brows(end+1,:)={j,b,f.targets(j),f.block_count(b),trapz(f.E(mask),y),sqrt(mean((y-mean_shape).^2))}; %#ok<AGROW>
 end
end
writetable(cell2table(brows,'VariableNames',{'target','block','q','valid_frames','B1_area_per_frame','RMS_difference_from_mean_block'}),fullfile(dest,'block_spectral_diagnostics.csv'));
all=load(fullfile(dest,'solver_background_comparison.mat')); sets={all.legacy,all.independent,all.background}; solver={}; compare={};
for mode=1:3
 ds=sets{mode};
 for j=1:numel(ds)
  for fit=ds{j}.fits
   for k=1:numel(fit.candidates)
    c=fit.candidates(k);
    solver(end+1,:)={mode,string(ds{j}.key),fit.n_components,k,c.exitflag,c.objective,c.firstorderopt,c.iterations,c.funcCount, ...
     string(c.candidate_type),c.jacobian_condition,string(mat2str(c.singular_values,8)), ...
     string(mat2str(c.column_scaled_singular_values,8)),fit.success&&c.exitflag<=0&&c.objective<fit.normalized_sse}; %#ok<AGROW>
   end
   old=all.legacy{j}.fits(fit.n_components);
   compare(end+1,:)={mode,string(ds{j}.key),fit.n_components,fit.normalized_sse-old.normalized_sse, ...
    max(abs(fit.prediction-old.prediction))/max(1,max(abs(old.observed))),string(fit.baseline_mode)}; %#ok<AGROW>
  end
 end
end
writetable(cell2table(solver,'VariableNames',{'mode','key','n_components','start','exitflag','objective','firstorderopt','iterations','funcCount', ...
 'candidate_type','raw_J_condition','raw_J_singular_values','column_scaled_J_singular_values','unconverged_below_selected'}),fullfile(dest,'solver_diagnostics.csv'));
writetable(cell2table(compare,'VariableNames',{'mode','key','n_components','delta_Q_vs_legacy','relative_max_curve_change_vs_legacy','baseline_mode'}),fullfile(dest,'solver_background_differences.csv'));
p=jsondecode(fileread(fullfile(out,'provenance','parents.json')));
rec=jsondecode(fileread(fullfile(p.v2,'input_manifest.resolved.yaml'))); rec=rec.sessions(strcmp({rec.sessions.session_id},'590_PL2_10w'));
raw=rec.files(find(endsWith(string({rec.files.path}),'.npy'),1));
fid=fopen(raw.path,'rb'); cl=onCleanup(@()fclose(fid)); header=char(fread(fid,256,'*uint8').');
line_end=find(header==char(10),1); assert(~isempty(line_end));
c20_v4_io('text',fullfile(out,'provenance','NPY_header.txt'),header(11:line_end));
diffrows={}; old=dir(fullfile(p.v2,'source_snapshot','**','*.m'));
for k=1:numel(old)
 path=fullfile(old(k).folder,old(k).name); rel=extractAfter(path,[fullfile(p.v2,'source_snapshot') filesep]);
 now=fullfile(pwd,rel);
 if isfile(now), diffrows(end+1,:)={string(rel),string(c20_v4_io('hash',path)),string(c20_v4_io('hash',now))}; end %#ok<AGROW>
end
writetable(cell2table(diffrows,'VariableNames',{'path','parent_snapshot_sha256','current_sha256'}),fullfile(out,'provenance','parent_code_comparison.csv'));
end
function [width,apex]=width0(E,Y)
mask=E>=-100&E<=100; E=E(mask); Y=Y(mask); width=NaN; apex=NaN;
if any(~isfinite(Y))||isempty(Y)||max(Y)<=0, return; end
[pk,ip]=max(Y); apex=E(ip); l=find(Y(1:ip)<=pk/2,1,'last'); r=ip-1+find(Y(ip:end)<=pk/2,1,'first');
if isempty(l)||isempty(r)||l==ip||r==ip, return; end
width=interp1(Y(r-1:r),E(r-1:r),pk/2)-interp1(Y(l:l+1),E(l:l+1),pk/2);
end
