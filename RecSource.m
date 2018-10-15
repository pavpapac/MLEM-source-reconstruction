
function [A,nrec,Xn,q,m,Xq,k50,FWHM,TWHM,err90_10,i]=RecSource(field,x,dose,kernel,pro,withKernel)
%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION RECSOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstructs the source for a selected field size and system matrix. 
% INPUT: field=field size (mm) as projected at SSD=100 cm, x = off-axis 
% positions (mm), dose = dose profile measured on Pb foil, kernel = PSF 
% used for de-convolving the dose profile to fluence profile, pro='cro' 
% or 'in' for crossplane and inplane orientations respectively, withKernel 
% = 1, then use convolve the system matrix with PSF, otherwise skip.

% First derive the system matrix 
[A,n,Xn,Xq,Zn,Zq,res,range]=ExtrSystemMat(field,pro);

% Use a spline intrpolation on the image to sample data to the needed
% resolution. Validate that the choices you make give you a smooth system
% matrix. 

dose=spline(x,dose,-range/2:res:range/2);
kernel=spline(x,kernel,-range/2:res:range/2);

%Convolve each column with the kernel while keeping the same dimensions 
%of matrix A.
switch(withKernel)
    case 1
        for j=1:size(A,2)
            a_temp=conv(A(:,j),kernel,'same'); 
            a_temp=a_temp./max(a_temp);
            a_temp(isnan(a_temp))=0;
            A(:,j)=a_temp;
        end
end

% Now perform a first reconstruction assuming a uniform distribution 
% on the source plane and using the system matrix A

[nrec,m,FWHM,TWHM,q,k50,err90_10]=MLEM(A,n,dose,res);

% iterative reconstruction. Repeat until the FWHM converges to a value
% of absolute diff < 2 % compared to last 30 previous iterations or untill a maximum
% number of iterations is reached (max=1000)

for i=1:1000
     FWHM_old=FWHM;
     [nrec,m,FWHM,TWHM,q,k50,err90_10]=MLEM(A,nrec,dose,res);
     k(i)=abs(1-(FWHM./FWHM_old)); %% calculate the correction to the previous estimation
     if i>30
          %if in the last 30 measurements there was no correction larger
          %than 2 % then consider that the FWHM has converged and stop the
          %iterations.
          last30=k(length(k):-1:length(k)-29);
          lim=find(last30>0.02, 1);
          if isempty(lim)
              break;
          end
     end
end

for j=1:i
    [nrec,m,FWHM,TWHM,q,k50,err90_10]=MLEM(A,nrec,dose,res);
end

% (OPTIONAL) UNCOMMENT TO FIT EITH GAUSSIAN FUNCTION, COMPARE WITH MC
% INPUT SOURCES AND PLOT RESULTS. 

% gausseqn='exp(-x.^2./(2.*(sigma.^2)))';
% fit_rec=fit(Xn',nrec,gausseqn);
% FWHMfit=2*sqrt(2*log(2)).*fit_rec.sigma;
% sigmafit=FWHMfit./(2*sqrt(2*log(2)));
% fgauss=@(Xn,sigmafit) exp(-(Xn'.^2)./(2.*(sigmafit).^2));
% nfit=fgauss(Xn,sigmafit)';

% % Now calculate the RMSE between the gaussian fit and the reconstructed source. 
% 
% se=(nrec-fgauss(Xn,sigmafit)).^2;
% elem100=find(nrec==1,1);
% elem10=find(nrec>0.10,1);
% src_mse100_10=mean(se(elem10:elem100))./mean(nrec(elem10:elem100));


%MC-tuned model source size

% switch(pro)
%     case 'cro'
%         FWHMmc=0.75;
%     case 'in'
%         FWHMmc=0.75;
% end
% % 
% sigmamc=FWHMmc./(2*sqrt(2*log(2)));

% %Now plot the reconstructed source and the gaussian fit
% figure(1);
% plot(Xn,nrec,'k-');
% hold
% plot(Xn,fgauss(Xn,sigmafit),'r--');
% %plot(Xn,fgauss(Xn,sigmamc),'b--');
% legend('MLEM reconstructed','Gaussian fit','Monte Carlo');

% figure(3);
% plot(src_RMSE100_10,'bo')




% %Now plot the measured and reconstructed profiles
% figure(2);
% plot(Xq,m,'r--');
% hold
% plot(Xq,q,'k-');
% legend('Measured profile', 'MLEM reconstructed profile')
% %hold
