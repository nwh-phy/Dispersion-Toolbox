function s = c20_v4_simulate(E,q_members,scenario,seed)
% Named hypothetical means, independent Gaussian engineering noise only.
stream=RandStream('mt19937ar','Seed',seed); m=peak_models('lorentz_symmetric');
bg=.12*(E/1000).^-1.3; truth=[900 600 900];
members=repmat(bg,1,numel(q_members));
switch scenario
 case 'static_single'
  members=members+repmat(m.model_fn(900,600,900,E),1,numel(q_members));
 case 'q_mixed_single'
  shift=160*(q_members-mean(q_members))/.0005;
  for j=1:numel(q_members), members(:,j)=members(:,j)+m.model_fn(900+shift(j),280,900,E); end
  truth=[900 280 900];
 case 'overlapping_double'
  truth=[620 800 800;1140 1260 900];
  members=members+repmat(m.model_fn(620,800,800,E)+m.model_fn(1140,1260,900,E),1,numel(q_members));
 case 'background_mismatch'
  members=members+repmat(m.model_fn(900,600,900,E)+.16,1,numel(q_members));
 otherwise, error('Unknown simulation scenario');
end
mu=mean(members,2); sigma=.01*max(mu); noise=sigma*randn(stream,size(mu));
s=struct('scenario',scenario,'seed',seed,'E',E,'q_members',q_members,'member_truth',members, ...
 'truth_parameters',truth,'mean',mu,'noise',noise,'observed',mu+noise,'sigma',sigma, ...
 'noise_provenance','hypothetical independent Gaussian, 1 percent max ordinate; not measured noise', ...
 'role','engineering_smoke_or_pilot_not_experimental_false_positive_rate');
end
