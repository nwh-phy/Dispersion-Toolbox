function denoised=BM3D_QRS(raw,BMfactor,undersamplingRatio)
% Block-matching and 3D filtering algorithm. Code downloaded from
% http://www.cs.tut.fi/~foi/GCF-BM3D and then modified by Ruishi Qi, Peking
% University. 
% Usage: 
% denoisedImage=BM3D_QRS(rawImage,BMfactor)
% Input arguments: 
% rawImage: raw image you want to denoise
% BMfactor: [optional] noise level overestimation factor. Default: 1 (no over-estimation)
if nargin<2
    BMfactor=1;
end
if nargin<3
    undersamplingRatio =1;
end
% raw=double(raw);
raw(isnan(raw)|isinf(raw))=0;
rawsize=size(raw);
% if min(rawsize)<8
%     disp('Minimum resolution of 8x8 needed!')
%     denoised =raw;
%     return
% end
if undersamplingRatio<1
    raw=imresize(raw,undersamplingRatio);
end

shift=min(raw(:));
raw=raw-shift;
scalefactor=255/max(raw(:));
raw=raw*scalefactor;
FT=fftshift(fft2(raw));
[x, y]=ndgrid(1:size(FT,1),1:size(FT,2));
x=x-mean(x(:));y=y-mean(y(:));
ok=x.^2/size(FT,1)^2+y.^2/size(FT,2)^2<0.05;
smth=abs(ifft2(FT.*ok));
% img(raw,[],smth,[],raw-smth,[])
% pause
diff=raw(:)-smth(:);
diff=sort(diff);
diff=diff(round(0.1*length(diff)):round(0.9*length(diff)));
sigma=rms(diff)*sqrt(2);
%
[~,denoised]=BM3D(1, raw, sigma*BMfactor);
denoised=denoised/scalefactor*255;
denoised=denoised+shift;
if undersamplingRatio<1
    denoised=imresize(denoised,rawsize);
end
%%
% for ii=1:size(data,3)
%     subplot 121,imagesc(data(:,:,ii)),setfig
%         set(gca,'clim',[0 max(data(:,:,ii),[],'all')])
%     subplot 122,imagesc(DenoisedData(:,:,ii)),setfig
%     set(gca,'clim',[0 max(data(:,:,ii),[],'all')])
% pause
% end