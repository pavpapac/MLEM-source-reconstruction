function [nrec,m,FWHMrec,TWHMrec,q,k50,err90_10]=MLEM(A,n,dose,res)
%%%%%%%%%%%% FUNCTION MLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the 'core' of the algorithm implementing 1 iteration of the 
% MLEM optimization for a given accelerator geometry. INPUT: A: system matrix, n= initial source
% distribution estimation, dose= measured dose profile using the Pb foil, 
% res= image and source resolution (mm). OUTPUT:nrec= new estimation of
% source distribution, m=measured profile, FWHMrec=FWHM of rec source,
% TWHMrec=TWHM of rec source, q=expected profile, k50=mean error in the
% 50-50 % region, err90_10= mean error between rec profile and
% initial profile in the 90-10 % region

%% MLEM iteration %%

%First, derive the expected image by applying the system matrix
% on the initial source distribution n. The first guess is usually a
% uniform distribution.

q=A*n;

% Measurement vector on the image plane

m=dose';

% Normalize to the mid point
q=q./q(round(length(q)/2));
m=m./m(round(length(m)/2));

% Now derive image pixel corrections as the ratios of the true measurement
% pixel values to the calculated pixel values. 

r=m./q;

% Remove any Inf or NaN points

r(r==Inf)=0;
r(r==-Inf)=0;
r(isnan(r))=0;

% Now backward map the projection bin corrections to the source plane. This
% give us a correction matrix for the source distribution.

c=A'*r;

% Normalize to the sensitivity matrix s=sum(A). This step basically
% normalizes the correction to the number of projection bin
% contributions. Essentially we are calculating the average correction for
% each source pixel from all the image pixel corrections.

s=sum(A)';
c=c./s;

% Now extract the new source estimation by applying the pixel
% corrections to the previous estimate.

nrec=n.*c;

% Finally re-normalize to the mid point
nrec=nrec./nrec(round(length(nrec)/2));
nrec(isnan(nrec))=0;

%% SOURCE AND PROFILE METRICS

%Find the FWHM and TWHM of the reconstructed source and FWHM of 
% and initial  and reconstructed profile. Calculate also the mean error 
% in the dose range 50-50% and 90-10%, which can be used as metrics to tune
% the jaw positions.

elem50=find(q>0.50,1);
elemMAX=find(q==1);
FWHMq=2*(abs(elemMAX-elem50).*res);

elem50=find(m>0.50,1);
elemMAX=find(m==1,1);
FWHMm=2*(abs(elemMAX-elem50).*res);

elem50=find(nrec>0.50,1);
elemMAX=find(nrec==1,1);
FWHMrec=2*(abs(elemMAX-elem50).*res);

elem10=find(nrec>0.10,1);
elemMAX=find(nrec==1,1);
TWHMrec=2*(abs(elemMAX-elem10).*res);

k50=FWHMm./FWHMq;
k50=abs(1-k50);
elem90=find(m>0.9,1);
elem10=find(m>0.1,1);
err90_10=mean(abs(r(elem10:elem90)-1));


