function b1_fullq_stage0(out, labels)
% B1 full-q task, stage 0: does the energy-dependent kinematic prefactor fix
% the |q| = 0.0025 misfit, and which ZLP component count / perpendicular
% acceptance to carry forward? v7 A1 spectra of the six mirror-q targets,
% DL n=2 B1 peaks. Writes one CSV per target into out. Parent outputs stay read-only.
root=bisb_find_project_root(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
p2=fullfile(root,'paper_results','b1_components_v2','20260912T021813081Z_31fcc543_18af79b6');
p7=fullfile(root,'paper_results','b1_mirror_q_v7','20261007T074250Z');
pp=fullfile(root,'paper_results','b1_physics_preview','20260912T064410245Z_5756ece8_93247099');
cfg=struct('model','lorentz','n_peaks',2,'J_starts',12,'seed',20260912,'window',[300 1800],'beam_kV',30);
all_labels={'R1','R2','R3','R1m','R2m','R3m'}; seed_offset=[1 2 3 1 2 3];
if nargin<2, labels=all_labels; end
if ~isfolder(out), mkdir(out); end
l=load(fullfile(p2,'590_PL2_10w','L1_minimal.mat')); qe=l.d.qe;
m=load(fullfile(pp,'590_PL2_10w','A1_full_q_map.mat')); sp=probe_sigma(m.map);
v7=load(fullfile(p7,'mirror_fits.mat')); keys=cellfun(@(a)a.key,v7.A,'UniformOutput',false);
K={'none','',0;'2d_h0','2d',0;'2d_h05','2d',0.0005;'2d_h10','2d',0.001;'3d_h0','3d',0};
aux={'aux1',[30 300];'aux2f',[30 80;80 300]};
for li=1:numel(labels)
 r=find(strcmp(all_labels,labels{li})); d7=v7.A{strcmp(keys,[labels{li} '_' cfg.model])};
 E=d7.E1(:); Y=d7.Y1(:); mq=qe.q_Ainv(d7.n3_channels); seed=cfg.seed+seed_offset(r);
 rows={};
 for nz=[2 3]
  for a=1:size(aux,1)
   for k=1:size(K,1)
    pre=[];
    if ~isempty(K{k,2})
     pre=qe_kinematic_prefactor(mq,form=K{k,2},h_perp=K{k,3},sigma_probe=sp,beam_kV=cfg.beam_kV);
    end
    t0=tic; f=qe_zlp_joint_fit(E,Y,'fit_window',[-Inf cfg.window(2)],'signal_window',cfg.window, ...
     'peak_model',cfg.model,'n_peaks',cfg.n_peaks,'n_starts',cfg.J_starts,'seed',seed, ...
     'aux_windows',aux{a,2},'n_zlp',nz,'prefactor',pre); secs=toc(t0);
    rows(end+1,:)=summary_row(labels{li},mean(mq),nz,aux{a,1},K{k,1},f,secs); %#ok<AGROW>
    fprintf('%s nz=%d %s %s: chi2 gain %.2f low %.2f sig %.2f | E0 %.0f/%.0f W %.0f/%.0f | wlow %.3f | %.0fs\n', ...
     labels{li},nz,aux{a,1},K{k,1},rows{end,6},rows{end,7},rows{end,8},rows{end,9},rows{end,11},rows{end,10},rows{end,12},rows{end,13},secs);
   end
  end
 end
 T=cell2table(rows,'VariableNames',{'region','q_mean','n_zlp','aux','kinematic','chi2_gain','chi2_low', ...
  'chi2_signal','E0_1','W_1','E0_2','W_2','weight_low','status','seconds'});
 writetable(T,fullfile(out,sprintf('stage0_%s.csv',labels{li})));
end
fprintf('sigma_probe = %.2e A^-1\n',sp);
end

function row=summary_row(label,q,nz,aux,kin,f,secs)
p=nan(2,2); cg=NaN; cl=NaN; cs=NaN; w=NaN;
if f.success
 p=f.parameters(:,1:2); cg=f.chi2_red_gain; cs=f.chi2_red_signal;
 x=f.energy_meV-f.zlp_parameters.center_meV; lo=x>=30&x<=300&isfinite(f.residual);
 cl=mean((f.residual(lo)./f.sigma(lo)).^2);
 Ep=f.energy_meV>0; a=trapz(f.energy_meV(Ep),f.loss_function(Ep,:)); w=a(1)/sum(a);
end
row={label,q,nz,aux,kin,cg,cl,cs,p(1,1),p(1,2),p(2,1),p(2,2),w,f.numerical_status,secs};
end

function s=probe_sigma(map)
% Probe angular spread from the elastic q-profile at q = 0 (|E| <= 12 meV):
% second moment of the central five channels minus the channel box variance.
el=sum(map.A1(abs(map.E)<=12,:),1); el(map.invalid_q)=0; [~,i0]=max(el);
k=-2:2; w=el(i0+k); w=w/sum(w); mu=sum(w.*k);
s=sqrt(max(sum(w.*(k-mu).^2)-1/12,0))*median(diff(map.q));
end
