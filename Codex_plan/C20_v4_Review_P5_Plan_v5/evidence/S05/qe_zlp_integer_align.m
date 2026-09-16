function a = qe_zlp_integer_align(E,reference_profiles,X,window)
% X is E x sequence x member. Estimate only in fixed elastic window.
assert(size(reference_profiles,1)==numel(E)&&size(X,1)==numel(E));
idx=find(E>=window(1)&E<=window(2)); T=size(X,2);
peak=nan(1,T); valid=false(1,T);
for t=1:T
 y=reference_profiles(idx,t);
 if all(isfinite(y)) && max(y)>0 && max(y)>3*max(median(y),eps)
  [~,k]=max(y); valid(t)=k>1&&k<numel(idx);
  if valid(t), peak(t)=idx(k); end
 end
end
assert(any(valid),'qe_zlp_integer_align:NoZLP','No valid elastic reference');
ref=round(median(peak(valid))); offset=peak-ref;
support=(1:numel(E))'; keep=all(support+offset(valid)>=1 & support+offset(valid)<=numel(E),2);
support=support(keep); aligned=nan(numel(support),T,size(X,3));
for t=find(valid), aligned(:,t,:)=X(support+offset(t),t,:); end
a=struct('E',E(support),'support',support,'valid',valid,'reference_pixel',ref, ...
 'measured_offset_pixels',offset,'correction_pixels',-offset,'aligned',aligned, ...
 'window_meV',window,'interpolation','none; integer noncircular; common support','peak_pixel',peak);
end
