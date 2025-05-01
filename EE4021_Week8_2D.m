clc;
close all;
clear;
t = 250



%grid parameters
Nx = 75; Ny = 75; 


%const
wavelength = 4e-9;   %500nm
da = wavelength/20; 
mu = pi*4e-7;
e0 = 8.85e-12;

E0 = 1.0;
x0 = 70; y0 = 70;
c0 = 1 / sqrt(e0 * mu); % Speed of light in vacuum
c =sqrt(1/mu/e0);

dt = 0.1/c*((1/da)^2+(1/da)^2)^-0.5; %Stability Conditions
dt = (da / (c0 * sqrt(2)))*0.9 ;


%double the size because the yee grid is staggered
Nx2 = 2*Nx;
Ny2 = 2*Ny;

NPML = [0 20 20 20];

sigx = zeros(Nx2,Ny2);
for nx = 1 : 2*NPML(1)
nx1 = 2*NPML(1) - nx + 1;
sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
nx1 = Nx2 - 2*NPML(2) + nx;
sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end
sigy = zeros(Nx2,Ny2);
for ny = 1 : 2*NPML(3)
ny1 = 2*NPML(3) - ny + 1;
sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1 : 2*NPML(4)
ny1 = Ny2 - 2*NPML(4) + ny;
sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end




%% Calculating Update Coefficients

% Interpolate conductivities onto the 1x Yee grid locations


% Assume free space: URxx=1, URyy=1, EPSrxx=1, EPSryy=1
URxx = 1.0; URyy = 1.0;  EPSrzz=1.0; % Relative Permeability/Permittivity

sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = - c0./URxx./mHx0;
mHx3 = - (c0*dt/e0) * sigHx./URxx ./ mHx0;


sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./URyy./mHy0;
mHy3 = - (c0*dt/e0) * sigHy./URyy ./ mHy0;


sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c0./mDz0;
mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;


mEz1 = 1/EPSrzz;
% initialize fields
Hx = zeros(Nx, Ny);
Hy = zeros(Nx, Ny);
Dz = zeros(Nx, Ny);
Ez = zeros(Nx,Ny);


%initalize curl arrays
CEx = zeros(Nx, Ny); % Stores (Ez(i,j+1)-Ez(i,j))/da
CEy = zeros(Nx, Ny); % Stores (Ez(i+1,j)-Ez(i,j))/da
CHz = zeros(Nx,Ny); % Stores (dHy/dx - dHx/dy)

% initalize integration arrays (Auxiliary fields for PML)
ICEx = zeros(Nx, Ny);
ICEy = zeros(Nx, Ny);
IDz = zeros(Nx, Ny); % Initialize auxiliary field for Dz update if using mDz3 term

dHy_dx = zeros(Nx,Ny);
dHx_dy = zeros(Nx,Ny);


% Source Parameters (for injection)

source_x=fix(Nx/2); 
source_y=fix(Ny/2); 
source_t0=40.0;
source_spread=7.0; 



%graph init
h_fig = figure;
ax = axes(h_fig);
im = imagesc(ax, Ez);
caxis(ax, [-0.05 0.05]);
colorbar;
axis equal tight;


%% Main loop



for n = 1:t
   % --- Update H field (using E and I fields from previous step) ---

   %compute curl of E
    CEx(:, 1:Ny-1) = (Ez(:, 2:Ny) - Ez(:, 1:Ny-1)) / da;
  % CEx(:, 1:Ny-1) = (Ez(:, 2:Ny) - Ez(:, 1:Ny-1)) / da;
    CEy(1:Nx-1, :) = -(Ez(2:Nx, :) - Ez(1:Nx-1, :)) / da;

   %update h integration
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;

    %Update H Field
    Hx = mHx1 .* Hx + mHx2 .* CEx + mHx3 .* ICEx;
    Hy = mHy1 .* Hy + mHy2 .* CEy + mHy3 .* ICEy; % Note the '+' sign for mHy2 term

    % --- Update E field (using H and I fields from current step) ---
    
   
%compute Curl of H
    dHy_dx(2:Nx, :) = (Hy(2:Nx, :) - Hy(1:Nx-1, :))/da;
    dHx_dy(:, 2:Ny) = (Hx(:, 2:Ny) - Hx(:, 1:Nx-1))/da;
    CHz = dHy_dx - dHx_dy;
            

source_pulse=-2.0*((source_t0-n)./source_spread).*exp(-1.*((source_t0-n)./source_spread)^2);
    Dz = (mDz1 .* Dz + mDz2 .* CHz + mDz4 .* IDz); % Applying formula from slide 4 with integrated IDz
    Dz(source_x,source_y) = Dz(source_x,source_y)+source_pulse;
    % --- Inject Source ---

    % Update Ez Field from Dz

    Ez = mEz1.*Dz;
   

    % --- Visualization ---
 set(im, 'CData', Ez);
    title(ax, ['Time Step = ', num2str(n)]);
    
    % Append frame to GIF
   % exportgraphics(h_fig, 'wave_simulation2.gif', 'Append', true);
    

 c = colorbar;
 vis_limit =0.05; 
c.Limits = [-vis_limit vis_limit]; % Use renamed vis_limit
caxis([-vis_limit vis_limit]); % Use renamed vis_limit
 pause(0.02);

end