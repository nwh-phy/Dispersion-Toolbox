function tests=test_b1_component_pilot_v2
tests=functiontests(localfunctions);
end
function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(fullfile(root,'src')));
end
function q=fixture(dq)
q=struct('q_Ainv',(-.015:dq:.015),'dq_Ainv',dq,'energy_meV',(300:4:1800)');
q.intensity=repmat(1:numel(q.q_Ainv),numel(q.energy_meV),1);
end
function testConservationRealDq(t)
for dq=[.0005 .00025]
 q=fixture(dq);
 for N=[1 3 5]
  b=qe_prepare_count_bins(q,N,variance=2*ones(size(q.intensity)),variance_source="synthetic independent");
  verifyEqual(t,sum(b.sum,2),sum(q.intensity(:,b.valid_q),2),'AbsTol',1e-10);
  members=[b.units.q_indices]; verifyEqual(t,sort(members),find(b.valid_q));
  verifyEqual(t,numel(unique(members)),numel(members));
  for k=1:numel(b.units)
   n=b.units(k).source_q_count;
   verifyEqual(t,b.mean(:,k)*n,b.sum(:,k),'AbsTol',1e-10);
   verifyEqual(t,b.variance_sum(:,k),2*n*ones(numel(q.energy_meV),1));
   verifyEqual(t,b.variance_mean(:,k),2/n*ones(numel(q.energy_meV),1),'AbsTol',1e-12);
   verifyEqual(t,b.units(k).q_right-b.units(k).q_left,n*dq,'AbsTol',1e-12);
   verifyEqual(t,numel(unique(sign(q.q_Ainv(b.units(k).q_indices)))),1);
  end
  if N>1, verifyTrue(t,any([b.units.source_q_count]==N)); end
 end
end
end
function testGapsInvalidAndPartial(t)
q=fixture(.0005); q.source_channel=1:size(q.intensity,2);
q.intensity(10,40)=NaN; q.source_channel(45:end)=q.source_channel(45:end)+1;
invalid=false(size(q.q_Ainv)); invalid(47)=true;
b=qe_prepare_count_bins(q,5,invalid_q=invalid);
verifyFalse(t,any(ismember([b.units.q_indices],[40 47])));
verifyTrue(t,any([b.units.partial]));
for u=b.units
 verifyTrue(t,all(diff(u.q_indices)==1)); verifyTrue(t,all(diff(u.source_channel)==1));
end
verifyTrue(t,all(isnan(b.variance_sum),'all'));
end
function testMatchedWidths(t)
a=qe_prepare_count_bins(fixture(.0005),3); b=qe_prepare_count_bins(fixture(.00025),6);
verifyEqual(t,a.units(1).q_width,b.units(1).q_width,'AbsTol',1e-12);
a=qe_prepare_count_bins(fixture(.0005),5); b=qe_prepare_count_bins(fixture(.00025),10);
verifyEqual(t,a.units(1).q_width,b.units(1).q_width,'AbsTol',1e-12);
end
function testExistingExtractorFixedBinCompatibility(t)
E=(300:8:1800)'; m=peak_models('lorentz_symmetric');
Y=.1*(E/1000).^-1.2+m.model_fn(750,180,300,E)+m.model_fn(1250,300,420,E);
q=struct('energy_meV',E,'q_Ainv',.001:.0005:.0035,'dq_Ainv',.0005,'intensity',repmat(Y,1,6));
opts=struct('binning_policy','fixed_nonoverlapping','bin_size',3,'q_range_Ainv',[-.015 .015], ...
 'q_skip_Ainv',.0005,'energy_window_meV',[300 1800],'peak_model','lorentz_symmetric', ...
 'min_peak_amplitude_fraction',0,'bootstrap_ci_samples',0);
b=qe_prepare_count_bins(q,3); a=b1_double_peak_binning_extract(q,q,opts);
verifyEqual(t,a.binning_map.source_q_count,[b.units.source_q_count]');
verifyEqual(t,a.binning_map.q_Ainv,[b.units.q_Ainv]','AbsTol',1e-12);
end
function testModelsAndContracts(t)
E=(300:8:1800)'; m=peak_models('lorentz_symmetric'); old=peak_models('lorentz');
verifyEqual(t,old.model_fn(900,200,1e6,E),1e6*E*200./((E.^2-900^2).^2+E.^2*200^2));
verifyEqual(t,m.model_fn(900,200,1,900),1/(pi*100),'AbsTol',1e-14);
Y=.1*(E/1000).^-1.3+m.model_fn(750,180,300,E)+m.model_fn(1250,300,420,E);
a=qe_compare_component_models(E,Y,n_starts=12);
b=qe_compare_component_models(E,3*Y,n_starts=12);
verifyTrue(t,all([a.success])); verifyEqual(t,a(2).parameters(:,1),[750;1250],'AbsTol',2);
verifyEqual(t,a(2).parameters(:,2),[180;300],'AbsTol',3);
verifyEqual(t,a(2).parameters(:,3),[300;420],'AbsTol',3);
for j=1:2
 verifyEqual(t,a(j).prediction,a(j).background+sum(a(j).components,2),'AbsTol',1e-12);
 verifyEqual(t,a(j).residual,Y-a(j).prediction,'AbsTol',1e-12);
 verifyEqual(t,a(j).sse,sum(a(j).residual.^2),'AbsTol',1e-12);
 verifyEqual(t,b(j).parameters(:,1:2),a(j).parameters(:,1:2),'AbsTol',.1);
 verifyEqual(t,b(j).prediction,3*a(j).prediction,'AbsTol',1e-5);
 for c=a(j).candidates
  verifyTrue(t,all(c.p0>=a(j).lb)&all(c.p0<=a(j).ub));
 end
end
single=qe_compare_component_models(E,.1*(E/1000).^-1.3+m.model_fn(950,320,400,E),n_starts=12);
verifyLessThan(t,single(1).sse,1e-8);
verifyEqual(t,single(2).scientific_status,'unresolved_or_invalid');
end
function testNoCacheReadOrWriteAndExplicitIdentity(t)
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s'));
raw=ones(16,512,20); raw(:,10,:)=10; energy=(0:511)'*4;
p=fullfile(folder,'synthetic_10w_raw.mat'); save(p,'raw','energy');
a3=999*ones(512,20); e=energy; cache=fullfile(folder,'eq3D_processed.mat'); save(cache,'a3','e');
d=load_qe_dataset(p,.0005,q_crop=[1 20],use_cache=false,write_cache=false);
verifyEqual(t,d.source_path,p); verifyLessThan(t,max(d.qe.intensity,[],'all'),999);
a=load(cache); verifyEqual(t,a.a3,a3);
explicit=load_qe_dataset(cache,.0005); verifyEqual(t,explicit.qe.intensity,a3);
% Same bytes and restored mtime cannot affect an explicitly no-cache read.
stamp=java.io.File(p).lastModified(); raw=2*raw; save(p,'raw','energy'); java.io.File(p).setLastModified(stamp);
d2=load_qe_dataset(p,.0005,q_crop=[1 20],use_cache=false,write_cache=false);
verifyEqual(t,d2.qe.intensity,2*d.qe.intensity);
end
function testRawSentinelCannotSetCalibration(t)
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s'));
p=fullfile(folder,'sentinel_10w.npy'); raw=ones(3,20,64,'single');
raw(:,:,10)=20; raw(:,8,10)=100; raw(1,19,50)=single(2^32);
header="{'descr': '<f4', 'fortran_order': True, 'shape': (3, 20, 64), }";
header=char(header); header=[header repmat(' ',1,mod(-(10+numel(header)+1),16)) newline];
fid=fopen(p,'wb','ieee-le'); fwrite(fid,[147 double('NUMPY') 1 0],'uint8');
fwrite(fid,numel(header),'uint16'); fwrite(fid,header,'char'); fwrite(fid,raw,'single'); fclose(fid);
meta=struct('is_sequence',true,'spatial_calibrations',struct('units','eV','scale',.004,'offset',0));
fid=fopen(replace(p,'.npy','.json'),'w'); fprintf(fid,'%s',jsonencode(meta)); fclose(fid);
d=load_raw_session(p,q_crop=[1 20],dq_Ainv=.0005,write_cache=false,align_zlp=false,show_progress=false,invalid_count_threshold=double(intmax('uint32')));
verifyEqual(t,d.qe.q_zero_index,8); verifyEqual(t,d.qe.energy_meV(10),0);
verifyEqual(t,d.invalid_native_q,19); verifyTrue(t,all(isnan(d.qe.intensity(:,19))));
verifyFalse(t,isfile(fullfile(folder,'eq3D_processed.mat')));
verifyEqual(t,d.frame_diagnostics.zlp_energy_pixel,10*ones(3,1));
end
