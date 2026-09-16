function bins = qe_centered_bins(qe,targets,Ns)
% Locked centers, reusing fixed-bin validation and sum contract.
bins=struct([]);
for t=targets
 for N=Ns
  [distance,c]=min(abs(qe.q_Ainv-t)); idx=c-(N-1)/2:c+(N-1)/2;
  u=struct('target_q',t,'N',N,'valid',false,'reason','invalid_members', ...
   'q_Ainv',NaN,'q_left',NaN,'q_right',NaN,'q_width',NaN,'indices',[], ...
   'source_channel',[],'q_members',[],'sum',[],'mean',[],'variance_mean',[], ...
   'valid_mask',[],'partial',false,'processing_level','L1','independent_of_other_N',false);
  if distance<qe.dq_Ainv*1e-6 && mod(N,2)==1 && all(idx>=1 & idx<=numel(qe.q_Ainv))
   local=qe; local.intensity=qe.intensity(:,idx); local.q_Ainv=qe.q_Ainv(idx);
   local.source_channel=qe.source_channel(idx); b=qe_prepare_count_bins(local,N);
   if numel(b.units)==1 && b.units.source_q_count==N
    v=b.units; u.valid=true; u.reason='valid'; u.q_Ainv=v.q_Ainv;
    u.q_left=v.q_left; u.q_right=v.q_right; u.q_width=v.q_width;
    u.indices=idx; u.source_channel=local.source_channel; u.q_members=local.q_Ainv;
    u.sum=b.sum; u.mean=b.mean; u.variance_mean=b.variance_mean; u.valid_mask=b.valid_q;
   end
  end
  if isempty(bins), bins=u; else, bins(end+1)=u; end %#ok<AGROW>
 end
end
end
