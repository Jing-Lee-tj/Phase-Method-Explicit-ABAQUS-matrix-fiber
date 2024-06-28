E1 = 26.5e9; %Longitudinal Young%s modulus
E2 = 2.6e9;  %Transverse Young’s modulus
G12 = 1.3e9; %Shear Modulus
v12 = 0.35;  %Poisson's ration
rho = 1.19e-6*1e9; %density
L_min = 0.2/1000; % minimum element size
eta_L = 3e-5 *1000; %artificial viscosity parameters for fiber
eta_m = eta_L;
l0_L  = 0.4 * 1e-3; %internal length scale
l0_m  = l0_L;
amp_displacement = 0.3 / 1000; %total displacement
%%%%%%%%%%%%%%%%% other variables %%%%%%%%%%%%%%%%%
c0 = integral(@(d) 4* sqrt(2*d-d.*d), 0, 1); %normalization constant
%%%%%%%%%%%%%%%%% other modelling setting
simulation_time = 0.001;
dtime = calculate_dtime(E1,v12,rho,L_min,c0,eta_L,l0_L,eta_m,l0_m);
t_d_curve = displacement_curve(amp_displacement,simulation_time,true);


%%%%%%%%%%%%%%% output %%%%%%%%%%%%%
% disp(['dtime=' num2str(dtime)]);