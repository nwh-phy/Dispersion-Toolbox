function archive=c20_v6_finish(out)
a=load(fullfile(out,'appendix','products.mat')); products=a.products; cfg=a.cfg;
checks={}; candidate_count=0; reused_count=0;
for si=1:numel(products)
 p=products{si}; if isempty(p), continue; end
 for di=1:numel(p.details)
  d=p.details{di}; [ok,ix]=ismember(d.unit.source_channel,p.map.source_channel); assert(all(ok));
  Y=mean(p.map.A1(:,ix),2);
  for f=d.fits
   candidate_count=candidate_count+numel(f.candidates); reused_count=reused_count+d.reused_v5;
   mask=p.map.E>=min(f.energy_meV)&p.map.E<=max(f.energy_meV); err=max(abs(Y(mask)-f.observed)); assert(err<1e-8);
   model=peak_models(f.peak_model);
   for c=f.candidates
    if ~all(isfinite(c.p)), continue; end
    pred=c.p(1)*f.scale*(f.energy_meV/1000).^(-c.p(2));
    for k=1:f.n_components
     ii=3+3*(k-1); pred=pred+model.model_fn(c.p(ii)*1000,c.p(ii+1)*1000,c.p(ii+2)*f.ampunit*f.scale,f.energy_meV);
    end
    q=sum((pred-f.observed).^2)/f.scale^2; assert(abs(q-c.objective)<1e-7*max(1,q));
   end
   if f.success
    curves=zeros(size(f.components));
    for k=1:f.n_components, v=f.parameters(k,:); curves(:,k)=model.model_fn(v(1),v(2),v(3),f.energy_meV); end
    assert(max(abs(curves-f.components),[],'all')<1e-8*max(f.observed));
    assert(max(abs(f.observed-f.prediction-f.residual))<1e-8);
   end
   checks(end+1,:)={string(p.session),string(d.key),f.n_components,d.reused_v5,err,f.success}; %#ok<AGROW>
  end
 end
end
assert(size(checks,1)==56&&candidate_count==1680&&reused_count==12);
writetable(cell2table(checks,'VariableNames',{'session','key','n','reused','input_error','success'}),fullfile(out,'appendix','actual_array_checks.csv'));
b=load(fullfile(out,'appendix','parent_hashes_before.mat'));
assert(isequal(b.before{1},c20_v4_io('inventory',b.p2))&&isequal(b.before{2},c20_v4_io('inventory',b.p5)));
r=runtests('tests/C20PhysicsPreviewTest.m'); assert(all([r.Passed])); writetable(table(r),fullfile(out,'appendix','final_tests.csv'));
sources={which('c20_v6_physics_figures'),which('qe_spectral_centroids'),which('c20_v6_finish'),which('c20_packet_readback')}; rows={};
for k=1:numel(sources)
 file=sources{k}; h=c20_v4_io('hash',file); [~,name,ext]=fileparts(file); copyfile(file,fullfile(out,'appendix',[name '_final' ext]));
 rows(end+1,:)={string(file),string(h)}; %#ok<AGROW>
end
writetable(cell2table(rows,'VariableNames',{'path','sha256'}),fullfile(out,'appendix','final_export_source_hashes.csv'));
c20_v4_io('json',fullfile(out,'stage_status.json'),struct('session590_sparse_points',10,'repeat_points',4,'reused_models',12, ...
 'new_models',44,'saved_models',56,'candidates_checked',candidate_count,'centroid_identity','passed', ...
 'frozen_shapes','descriptive only; reversed apex order or missing reference skipped', ...
 'literature','six DOI-verified originals; five full texts, one author abstract only', ...
 'physical_origin','candidate hypotheses, not identified quasiparticles','CI','not_computed','old_runs_unchanged',true));
packet=fullfile(out,'physics_packet'); assert(~isfolder(packet)); mkdir(packet);
for folder={'figures','590_PL2_10w','n0_PL2_10w_repeat'}, copyfile(fullfile(out,folder{1}),fullfile(packet,folder{1})); end
for file={'physics_readout.md','literature_comparison.md','experiment_calculation_brief.md','physics_candidates.csv','spectral_weight_readout.csv','frozen_shape_counterfactuals.mat','config_resolved.json','stage_status.json'}
 copyfile(fullfile(out,file{1}),packet);
end
mkdir(fullfile(packet,'appendix'));
files=dir(fullfile(out,'appendix','*')); files=files(~[files.isdir]);
for k=1:numel(files)
 name=files(k).name;
 if endsWith(name,'.pdf')||ismember(name,{'Timrov2017.txt','Torbatian2020.txt','products.mat'}), continue; end
 copyfile(fullfile(files(k).folder,name),fullfile(packet,'appendix',name));
end
manifest=c20_v4_io('inventory',packet); writetable(manifest,fullfile(packet,'FILE_MANIFEST.csv'));
archive=fullfile(out,'physics_packet.zip'); assert(~isfile(archive)); zip(archive,{'*'},packet);
c20_packet_readback(archive,fullfile(out,'packet_readback_storage.json'),@(decoded)persistReadback(decoded,out));
disp(['PHYSICS_PACKET=' archive]);
end

function result=persistReadback(decoded,out)
count=0;
for sid={'590_PL2_10w','n0_PL2_10w_repeat'}
 a=load(fullfile(decoded,sid{1},'fit_details.mat'));
 for d=a.details
  d=d{1};
  for f=d.fits
   count=count+1; assert(max(abs(f.background+sum(f.components,2)-f.prediction))<1e-8*max(f.observed));
   assert(abs(sum(f.residual.^2)-f.sse)<1e-7*max(1,f.sse));
  end
 end
end
manifest=readtable(fullfile(decoded,'FILE_MANIFEST.csv'),TextType='string');
assert(count==56); result=struct('payload_hashes',height(manifest),'fit_arrays',count,'passed',true);
c20_v4_io('json',fullfile(out,'packet_readback.json'),result);
end
