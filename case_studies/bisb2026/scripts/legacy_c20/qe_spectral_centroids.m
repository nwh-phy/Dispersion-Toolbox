function r=qe_spectral_centroids(E,components,reference_components)
% Fixed-window identity and frozen-reference shapes, without extra fitting.
area=trapz(E,components,1); mu=trapz(E,E.*components,1)./area;
fraction=area/sum(area); combined=sum(components,2);
direct=trapz(E,E.*combined)/trapz(E,combined); identity=sum(fraction.*mu);
r=struct('area',area,'fraction',fraction,'component_centroid',mu, ...
 'combined_centroid',direct,'weighted_centroid',identity,'identity_error',abs(direct-identity), ...
 'frozen_centroid',NaN,'frozen_curve',nan(size(E)),'reference_applicable',false);
if nargin>=3 && ~isempty(reference_components) && isequal(size(reference_components),size(components)) && all(isfinite(reference_components),'all')
 ar=trapz(E,reference_components,1);
 if all(ar>0)
  frozen=sum((reference_components./ar).*fraction,2);
  r.frozen_curve=frozen; r.frozen_centroid=trapz(E,E.*frozen)/trapz(E,frozen); r.reference_applicable=true;
 end
end
end
