clc;close all;clear;
%map size&time
max_time=300 ;
max_space=150;

%parameters
mu = pi*4e-7;
elp = 8.85e-12;
c = sqrt(1/mu/elp); 
wavelength = 523e-9;   %500nm
freq = c/wavelength;

%Unit 
dx = wavelength/20;

dt = dx/(2*c) %stability

%E H Conditions
Ez=zeros(max_space+1,1);
Hy=zeros(max_space+1,1);
 
Eeta=dt/dx/elp;     %updating variable   
Heta=dt/dx/mu;      %updating variable

%Gaussian Source
t0=max_space/2;                                                   %Gaussian  start point
spread=7;                                                %Gaussian  width
Ez(2:max_space,1)= exp(-(t0-(2:max_space)).^2/spread^2);  %Gaussian  function
%Hy = -Ez/sqrt(mu/elp)


 H3 = 0;H2 = 0;H1 = 0;
 E3 = 0;E2 = 0; E1 = 0;


h_fig = figure;
ax = axes(h_fig);
im = plot(ax, Ez);



for n=1:max_time
      Hy(1:max_space-1) = Hy(1:max_space-1)+Heta*(Ez(2:max_space)-Ez(1:max_space-1));
      
      Hy(max_space) = Hy(max_space)+ Heta*(E3-Ez(max_space));
      H3 = H2; H2 = H1; H1 = Hy(1);

      %Hy(max_space) = Hy(max_space)+Heta*(0-Ez(max_space)); 
    %mHx = Heta
    %mEy = Eeta
      %Ez(1) = Ez(1)+Eeta*(Hy(1)-0); 
     
      
      Ez(2:max_space)=Ez(2:max_space)+Eeta*(Hy(2:max_space)-Hy(1:max_space-1));     %t=1              %equation 
      

      Ez(1) = Ez(1)+ Eeta*(Hy(1)-H3);
     
      E3 = E2; E2 = E1; E1 = Ez(max_space);
      
      

    title(ax, ['Time Step = ', num2str(n)]);
    
    % Append frame to GIF
   % exportgraphics(h_fig, 'wave_simulation3.gif', 'Append', true);
    %

      plot(Ez);
      axis([1 max_space -1 1]);
      %title(['1D-FDTD Right & Left   ','t = ',num2str(n)],'FontSize',18);
      pause(0.001)
      
      
    
end    
%}
