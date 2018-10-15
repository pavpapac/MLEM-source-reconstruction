function [nrec_std,FWHMstd,FWTMstd]=srcrec_errors(field,pro,x,DOSEav,DOSEstd,PSFav,PSFstd,N)
% FUNCTION SCRREC_ERRORS: Reconstructes the source with uncertainties due
% to jaw positioning, PSF and dose measurements. The number of iterations 
% (N) refers to the number of times the reconstructiong is repeated in
% order to achieve a robust estimate of the total uncertainty. 


nrec=zeros(1601,N,'double');
FWHM=zeros(1,N,'double');
FWTM=zeros(1,N,'double');

for j=1:N
    
    jaw_std=0.2; %assume a 0.2 mm std for jaw re-positioning
    jaw=10*(field+jaw_std*randn); % multiply field by 10 to make it an integer so you can use the round
    % After that divide by 10 again to go back to mm.
    jaw=round(jaw)./10;
    d=DOSEav+DOSEstd.*randn; %caluclate the dose profile with random variation based on measured profile std
    psf=PSFav+PSFstd.*randn; %caluclate the dose profile with random variation based on measured profile std
    [A,nrec(:,j),~,~,~,~,~,FWHM(j),FWTM(j),~,~]=RecSource(jaw,x,d,psf,pro,1);
end

%fianlly calculate the std of all reconstructed sources and FWHM,FWTM
%metrics.
nrec_std=std(nrec,0,2);
FWHMstd=std(FWHM);
FWTMstd=std(FWTM);