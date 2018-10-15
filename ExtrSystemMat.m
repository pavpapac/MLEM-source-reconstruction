function [A,n,Xn,Xq,Zn,Zq,res,range]=ExtrSystemMat(field,pro)

%% FUNCTION EXTRSYSTEMMAT: Derives the system matrix for a field size 
% as defined at SSD=100 cm and a profile orientation ('cro' or 'inp'). 
% OUTPUT: system matrix A, initial source estimate n, source coordinates
% Xn, image coordinates Xq, source z position Zn, image z position Zq,
% res=image and source resolution (mm), range = source and image range (mm) 
% To be used as an input to the MLEM algorithm. 

%% INITIALIZATION

Zn=0.0; % Source z position (mm)
Zq=1050; % image z position (mm)
res=0.01; % choose the image and source resolution. Validate that the choices
%you make give you a smooth system matrix. 
range=16; % choose the image and source range. 
Xn=-range/2:res:range/2; % source range and resolution
n=ones(length(Xn),1,'double');% initial estimate of source: uniform distribution
Xq=-range/2:res:range/2; % image range and coordinates (mm)
q=zeros(1,length(Xq),'double'); % projection image matrix
Ln=length(n); 
Lq=length(q);
Xq_left=zeros(Ln,1,'double');
Xq_right=zeros(Ln,1,'double');
hits_left_bot=zeros(Ln,1,'uint8'); % hits arrays
hits_left_top=zeros(Ln,1,'uint8');
hits_right_bot=zeros(Ln,1,'uint8');
hits_right_top=zeros(Ln,1,'uint8');
A=zeros(Lq,Ln,'double'); % System matrix A

%% JAWS

%Z-coordinates

switch(pro)
    case 'cro'
        Zjtop=366.1; % Novalis: 367 TrueBeam : 366.1
        Zjbot=444.1; % Novalis: 443.485 TrueBeam: 444.1 
        
    case 'in'
        Zjtop=278.9; % Novalis: 280 TrueBeam: 278.9
        Zjbot=356.6; % Novalis: 356.485 TruBeam: 356.6
end

% MC tuned field openings (width=0.47 cm at SSD=100 cm)
%Xjtop=0.08624;
%Xjbot=0.10422;

Xjtop=(field/2).*(Zjtop./1000);
Xjbot=(field/2).*(Zjbot./1000);

% lateral distance to top field opening
Dxright_top=Xjtop-Xn;
Dxleft_top=-Xjtop-Xn;

% lateral distance to bottom field opening
Dxright_bot=Xjbot-Xn;
Dxleft_bot=-Xjbot-Xn;

% tangents of rays tracing the edges of the field on the top and bottom 
% collimation level
TanThetaL_top=Dxleft_top./(Zjtop-Zn);
TanThetaR_top=Dxright_top./(Zjtop-Zn);

TanThetaL_bot=Dxleft_bot./(Zjbot-Zn);
TanThetaR_bot=Dxright_bot./(Zjbot-Zn);

%% RAY-TRACING AND JAW COLLIMATION

% First, ray-trace to the bottom collimation level (Xjbot) using the tangent
% of the top collimation level and check if there is any collision. For all the
% rays that pass through the opening, calculate the final position on the
% image plane. 

Xjbot_left_top=Xn+TanThetaL_top.*(Zjbot-Zn);
hits_left_bot(abs(Xjbot_left_top)>Xjbot)=1;
Xq_left(~hits_left_bot)=Xn(~hits_left_bot)+TanThetaL_top(~hits_left_bot).*(Zq-Zn);

% Now ray-trace to the top collimation level (Xjtop) using the tangent of 
% the bottom collimation level and check if there is any collision. 
% For all the rays that pass the top openings and did not experience a hit
% calculate the final position on the image plane. 

Xjtop_left_bot=Xn+TanThetaL_bot.*(Zjtop-Zn);
hits_left_top(abs(Xjtop_left_bot)>Xjtop)=1;
Xq_left(~hits_left_top)=Xn(~hits_left_top)+TanThetaL_bot(~hits_left_top).*(Zq-Zn);

% Ok, now repeat for the right collimation

Xjbot_right_top=Xn+TanThetaR_top.*(Zjbot-Zn);
hits_right_bot(abs(Xjbot_right_top)>Xjbot)=1;
Xq_right(~hits_right_bot)=Xn(~hits_right_bot)+TanThetaR_top(~hits_right_bot).*(Zq-Zn);

Xjtop_right_bot=Xn+TanThetaR_bot.*(Zjtop-Zn);
hits_right_top(abs(Xjtop_right_bot)>Xjtop)=1;
Xq_right(~hits_right_top)=Xn(~hits_right_top)+TanThetaR_bot(~hits_right_top).*(Zq-Zn);

% Now, find which pixels of the image matrix are within the field and assign
% a value of 1 to all of them. The resulting vector shows us which image  
% pixels each source pixel contributes. Thus, this is a COLUMN of the system
% matrix A. Repeat for all source pixels to create the system matrix
% A=Ln*Lq.

for i=1:length(Xq_left)
    
       Xq_pixels=find(Xq>Xq_left(i) & Xq<Xq_right(i));
       q(Xq_pixels)=1;
       A(:,i)=q;
       q=zeros(1,length(Xq),'double');
       
end
