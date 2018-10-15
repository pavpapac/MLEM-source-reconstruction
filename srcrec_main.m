function [A,field_opt,snrec,snrec_std,FWHMrec,FWHMstd,FWTMrec,FWTMstd,pro_err90_10,src_RMSE100_10]=srcrec_main(x,dose,dose_std)
%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION srcrec_main %%%%%%%%%%%%%%%%%%%%%%%
% This is the main code to run for reconstructing the source. 
% INPUT: x=off-axis positions (mm), dose = measured dose profile 
% with Pb foil, pro='cro' or 'in', N: number of reconstructions to estimate
% uncertainties due to jaw repositioning (set at 30 or more), PSF and film measurements.
% OUTPUT: field_opt (mm) = optimum field size selected, 
% FWHMrec (mm) = FWHM of reconstructed source, TWHMrec (mm)=
% TWHM of rec source, pro_err90_10 = mean absolute error between input and 
% rec profile in the 90-10 % region, src_RMSE100_10= rmse between the rec 
% source and a gaussian fit in the 100 - 10 % region.

 
%% INITIALIZATION %%

% First initialize the field array. This contains all possible jaw
% positions to be tested for the 5 mm field (4-6 mm). 
field=4.0:0.1:6.0;
N=length(field);

%Now initialize the FWHM, k50 and 90-10% error matrix

k50=zeros(1,N,'double');
err90_10=zeros(1,N,'double');
FWHM=zeros(1,N,'double');

% Use the lower jaws (crossplane -> 'cro') for all measurements to be consistent in both
% orientations. However, one may want to switch to inplane ('in') if the upper jaws
% are used instead. 

pro = 'cro';

% Define the number of reconstructions to be performed for uncertainty
% estimation. Recommended value N = 30 or larger

Niter=50;

%% LOAD KERNELS %%

load('PSF','PSFx_av','PSFy_av','PSFx_std','PSFy_av')

if strcmp(pro,'cro')
    psf=PSFx_av;
    psf_std=PSFx_std;
else
    psf=PSFy_av;
    psf_std=PSFy_std;
end  

%% ITERATIVE RECONSTRUCTION %%

% Now reconstruct the source for each posible jaw position. The jaw
% position that minimizes the 90-10% error matrix will be considered the
% closest to the true. 

for i=1:N
    
   [~,~,~,~,~,~,k50_temp,FWHM_temp,~,err90_10_temp,~]=RecSource(field(i),x,dose,psf,pro,1);
   k50(i)=k50_temp(1);
   FWHM(i)=FWHM_temp(1);
   err90_10(i)=err90_10_temp(1);
       
end

% Find the minimum of the 90-10% error matrix and
% extract the source distribution for the optimum jaw position.

min90_10=find(err90_10==min(err90_10),1);
field_opt=field(min90_10);
[A,nrec,Xn,q,m,Xq,~,FWHMrec,FWTMrec,pro_err90_10,~]=RecSource(field_opt,x,dose,psf,pro,1);

%% FIT GAUSSIAN FUNCTION ON RECONSTRUCTED SOURCE %%

gausseqn='exp(-x.^2./(2.*(sigma.^2)))';
fit_rec=fit(Xn',nrec,gausseqn);
FWHMfit=2*sqrt(2*log(2)).*fit_rec.sigma;
sigmafit=FWHMfit./(2*sqrt(2*log(2)));
sigmarec=FWHMrec./(2*sqrt(2*log(2)));
fgauss=@(x,sigmafit) exp(-(x.^2)./(2.*(sigmafit).^2));
fgauss_rec=@(x,sigmarec) exp(-(x.^2)./(2.*(sigmarec).^2));
nfit=fgauss(Xn,sigmafit)';
ngauss=fgauss_rec(Xn,sigmarec)';

%% COMPARISON WITH MC SIMULATIONS (OPTIONAL -> UNCOMMENT) %%

%MC-tuned model electron and photon source sizes

% switch(pro)
%     case 'cro'
%         FWHMemc=2.0;
%         FWHMpmc=1.53;
%     case 'in'
%         FWHMemc=2.0;
%         FWHMpmc=1.53;
% end

% sigmaemc=FWHMemc./(2*sqrt(2*log(2)));
% sigmapmc=FWHMpmc./(2*sqrt(2*log(2)));
% nemc=fgauss(Xn,sigmaemc);
% npmc=fgauss(Xn,sigmapmc);

%% Now calculate the RMSE between the gaussian fit and the reconstructed source. 
% Smooth and re-normalize to the middle

c=(nrec-nfit).^2;
elem100=find(nrec==1,1);
elem10=find(nrec>0.10,1);
src_RMSE100_10=sqrt(mean(c(elem10:elem100)))./mean(nrec(elem10:elem100));
snrec=smooth(nrec,15);
snrec=snrec./snrec(round(length(snrec)/2));

%% Now calculate errors due to jaw re-positioning. Assuming a random variation
% on jaw position of 0.2 mm (1 sigma) following a normal distribution. 

[snrec_std,FWHMstd,FWTMstd]=srcrec_errors(field_opt,pro,x,dose,dose_std,psf,psf_std,Niter);

%% PLOTS %%

%Reconstructed source and Gaussian source
figure(1);
subplot(1,3,1);
plot(Xn,snrec,'k-');
hold on;
plot(Xn,nfit,'r-');
%plot(Xn,ngauss,'r-');
plot(Xn,snrec+snrec_std,'k--');
plot(Xn,snrec-snrec_std,'k--');
title('Reconstructed source','fontsize',15);
xlabel('off-axis (mm)','fontsize',15);
ylabel('Relative fluence','fontsize',15);
legend('MLEM rec source','Gaussian fit','\pm 1 \sigma');
xlim([0 3]);
grid on;
% MC input electrona nd photon sources (optional)
%plot(Xn,nemc,'r--');
%plot(Xn,npmc,'b-.');

% Measured and reconstructed dose profiles
figure(1);
subplot(1,3,2);
plot(Xq,q,'k-');
hold on;
plot(Xq,m,'r--');
title('Dose profiles','fontsize',15);
xlabel('off-axis (mm)','fontsize',15);
ylabel('Relative dose','fontsize',15);
xlim([0 7])
grid on;
legend('MLEM rec profile','Input dose profile')

% Field optimization: jaw position vs 90-10% error matrix 
figure(1);
subplot(1,3,3);
plot(field,err90_10,'k-');
hold on;
title('Field optimization','fontsize',15);
xlabel('field size (mm)','fontsize',15);
ylabel('Error','fontsize',15);
grid on;

