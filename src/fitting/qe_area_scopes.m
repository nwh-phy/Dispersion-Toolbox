function a = qe_area_scopes(E,Y,reference)
% Integrate sampled curves only on observed support; no extrapolation.
if nargin<3, reference=[300 1800]; end
E=E(:); assert(size(Y,1)==numel(E)&&all(diff(E)>0));
mask=E>=reference(1)&E<=reference(2); full=min(E)<=reference(1)&&max(E)>=reference(2);
area=nan(1,size(Y,2)); support=[NaN NaN];
if nnz(mask)>=2, area=trapz(E(mask),Y(mask,:),1); support=[min(E(mask)) max(E(mask))]; end
a=struct('area_fit_window',trapz(E,Y,1),'area_reference_window',area, ...
 'fit_window_meV',[min(E) max(E)],'reference_window_meV',reference, ...
 'reference_actual_support_meV',support,'reference_window_fully_observed',full, ...
 'integration_method','trapezoid on observed samples only','area_unit','corrected-count ordinate * meV');
end
