clc;close all;clear;
t = 350
Nx = 200; Ny = 200;
%param

wavelength = 4e-9;  
da = wavelength/20;

mu0 = pi*4e-7;
e0 = 8.85e-12; 
c =sqrt(1/mu0/e0)

dt = 1/c*((1/da)^2+(1/da)^2)^-0.5; %Stability Conditions
%dt = da/(2*c) %new stability condition
Hx = zeros(Nx, Ny); 
Hy = zeros(Nx, Ny); 

E0 = 1.0;

x0 = Nx-Nx/3; y0 = Ny/2; 
sigma_x = 3; sigma_y = 3; 
[X, Y] = meshgrid(1:Ny, 1:Nx);
Ez = E0 * exp(-((X - x0).^2 / (2 * sigma_x^2) + (Y - y0).^2 / (2 * sigma_y^2)));

Z0 = sqrt(mu0 / e0);

theta =pi/8;  %propagation angle


Hx = (sin(theta+pi/2) / Z0) * Ez; 
Hy = (-cos(theta+pi/2) / Z0) * Ez; 


Heta = dt/mu0/da;
Eeta = dt/e0/da;

% --- Prepare Figure ---


%graph init
h_fig = figure;
ax = axes(h_fig);
im = imagesc(ax, Ez);
caxis(ax, [-0.1 0.1]);
colorbar;
axis equal tight;





for n=1:t
  
   Hx(1:Nx-1, 1:Ny-1) = Hx(1:Nx-1, 1:Ny-1) - Heta * (Ez(1:Nx-1, 2:Ny) - Ez(1:Nx-1, 1:Ny-1));
   Hy(1:Nx-1, 1:Ny-1) = Hy(1:Nx-1, 1:Ny-1) + Heta * (Ez(2:Nx, 1:Ny-1) - Ez(1:Nx-1, 1:Ny-1));
   Ez(2:Nx-1, 2:Ny-1) = Ez(2:Nx-1, 2:Ny-1) + Eeta * ((Hy(2:Nx-1, 2:Ny-1) - Hy(1:Nx-2, 2:Ny-1)) - (Hx(2:Nx-1, 2:Ny-1) - Hx(2:Nx-1, 1:Ny-2)));


   Ez(1, :) = 0; Ez(Nx, :) = 0;
   Ez(:, 1) = 0; Ez(:, Ny) = 0;

  
      
 set(im, 'CData', Ez);
    title(ax, ['Time Step = ', num2str(n)]);
    

   %exportgraphics(h_fig, 'e2_2d60.gif', 'Append', true);
    
shading flat;
 c = colorbar;
 vis_limit =0.1; 
c.Limits = [-vis_limit vis_limit]; 
caxis([-vis_limit vis_limit]); 
 pause(0.02);

   % --- End of Visualization Modification ---

end