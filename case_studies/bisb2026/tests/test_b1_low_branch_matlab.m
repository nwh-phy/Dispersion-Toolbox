function tests=test_b1_low_branch_matlab
tests=functiontests(localfunctions);
end

function setupOnce(tc)
case_root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(case_root,'scripts'));
root=bisb_find_project_root(case_root);run(fullfile(root,'startup.m'));
tc.TestData.root=root;
end

function testExistingJointFitterRecoversTwoKnownResponses(tc)
E=(-180:8:1800).';K=qe_kinematic_prefactor(.005,dq=.0005,h_perp=.0066,sigma_probe=.0004);
z=5e4*(1+(E/12).^2/2).^(-2)+2000*(1+(E/35).^2).^(-1);
a=abs(E);L=zeros(numel(E),2);
for j=1:2
 e0=[650 1350];g=[250 350];
 L(:,j)=a*g(j)./((a.^2-e0(j)^2).^2+(a*g(j)).^2).*K(max(a,250));
 L(E<0,j)=L(E<0,j).*exp(-a(E<0)/25.852);
 L(:,j)=L(:,j)/max(L(:,j));
end
Y=z+L*[400;300];
f=qe_zlp_joint_fit(E,Y,signal_window=[300 1800],peak_model='lorentz', ...
 n_peaks=2,n_zlp=2,n_starts=4,kernel_mode='none',prefactor=K, ...
 prefactor_floor_meV=250,noise_sigma=sqrt(max(Y,1)));
verifyTrue(tc,f.success);
verifyEqual(tc,f.parameters(:,1),[650;1350],'AbsTol',15);
verifyEqual(tc,f.parameters(:,2),[250;350],'AbsTol',30);
verifyLessThan(tc,max(abs(f.prediction-f.background-sum(f.peaks,2))),1e-8*max(f.prediction));
end

function testPreparedBinsKeepNativeLowQAndOnlyThreeFiveHighQ(tc)
D=prepared(tc);q=[D.units.q_Ainv];N=[D.units.N];
verifyEqual(tc,unique(N(abs(q)<=.006+1e-12)),1);
verifyEqual(tc,unique(N(abs(q)>.006+1e-12)),[3 5]);
verifyLessThan(tc,D.max_weight_error,1e-8);
for n=[1 3 5]
 v=sort(q(N==n));verifyEqual(tc,v,-fliplr(v),'AbsTol',1e-12);
end
end

function testPreparedCountsMatchIndependentNativeFrameRead(tc)
D=prepared(tc);raw=read_npy(D.signature.source);
meta=jsondecode(fileread(strrep(D.signature.source,'.npy','.json')));
ax=meta.spatial_calibrations;ie=find(contains(string({ax.units}),'eV'),1,'last');
E=(ax(ie).offset+(0:size(raw,3)-1)*ax(ie).scale)'*1000-D.energy_origin_meV;
[ok,ei]=ismembertol(D.E,E,1e-8,'DataScale',1);verifyTrue(tc,all(ok));
q=((1:size(raw,2))-D.signature.session.zero)*D.dq;
for N=[1 3 5]
 ix=find([D.units.N]==N&[D.units.q_Ainv]>0);b=ix(ceil(numel(ix)/2));u=D.units(b);
 for t=[1 100 250]
  center=u.q_Ainv+D.drift_Ainv(t);
  w=max(0,min(q+D.dq/2,center+u.q_width/2)-max(q-D.dq/2,center-u.q_width/2))/D.dq;
  used=find(w>1e-12);eraw=ei+D.alignment.measured_offset_pixels(t);
  measured=reshape(double(raw(t,used,eraw)),numel(used),numel(ei)).'*w(used).';
  verifyEqual(tc,double(D.frames(:,t,b)),measured,'RelTol',2e-6,'AbsTol',1e-4);
 end
end
end

function D=prepared(tc)
file=fullfile(tc.TestData.root,'paper_results','b1_low_branch_matlab_20261024','590','prepared.mat');
assumeTrue(tc,isfile(file),'Local raw-data preparation absent; run MATLAB prepare stage for real-data tests');
L=load(file,'D');D=L.D;
end
