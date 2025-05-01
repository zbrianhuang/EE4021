clc;close all;clear;
%map size&time
max_time=550 ;
max_space=201;
%parameters
mu = pi*4e-7;
elp = 8.85e-12;
c = sqrt(1/mu/elp);
wavelength = 523e-9;   %500nm
freq = c/wavelength;
%Unit
dx = wavelength/20;
dy = dx;
dt = 1/c*((1/dx)^2+(1/dy)^2)^-0.5; %Stability Conditions
%E H Conditions
Ez=zeros(max_space+1,1);
Hy=zeros(max_space+1,1);
Eeta=dt/dx/elp;     %updating variable  
Heta=dt/dx/mu;      %updating variable
%Gaussian Source
t0=max_space/2;                                                   %Gaussian  start point
spread=15;                                                %Gaussian  width
Ez(2:max_space,1)= exp(-(t0-(2:max_space)).^2/spread^2);  %Gaussian  function
Hy  = -Ez/sqrt(mu/elp)
%graph init
h_fig = figure;
im = plot(Ez)
for n=1:max_time
     Hy(1:max_space-1)=Hy(1:max_space-1)+Heta*(-Ez(1:max_space-1)+Ez(2:max_space)); %t=0.5          %equation
     Ez(2:max_space)=Ez(2:max_space)+Eeta*(Hy(2:max_space)-Hy(1:max_space-1));     %t=1              %equation
    
     Ez(1) =0;
     Hy(max_space) = Hy(max_space)+Heta*(0-Ez(max_space));
   
     %exportgraphics(h_fig, 'e1.gif', 'Append', true);
    
     plot(Ez);
     axis([1 max_space -1 1])
     title(['1D-FDTD Right & Left   ','t = ',num2str(n)],'FontSize',18);
     pause(0.001)
    
  
end   
%}

