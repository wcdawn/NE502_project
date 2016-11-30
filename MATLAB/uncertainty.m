clear

% Take data from files containing points for correlations in F & S
fname_f_upper = 'f_upper.csv';
fname_f_lower = 'f_lower.csv';
fname_s_upper = 's_upper.csv';
fname_s_lower = 's_lower.csv';

fid_f_upper = fopen(fname_f_upper,'r');
fid_f_lower = fopen(fname_f_lower,'r');
fid_s_upper = fopen(fname_s_upper,'r');
fid_s_lower = fopen(fname_s_lower,'r');

f_upper = fscanf(fid_f_upper,'%e,%e',[2 Inf]);
f_lower = fscanf(fid_f_lower,'%e,%e',[2 Inf]);
s_upper = fscanf(fid_s_upper,'%e,%e',[2 Inf]);
s_lower = fscanf(fid_s_lower,'%e,%e',[2 Inf]);

fclose(fid_f_upper);
fclose(fid_f_lower);
fclose(fid_s_upper);
fclose(fid_s_lower);

% Data from Doster
P = 1040; % psia
G = 1.21e6; % lbm/(hr*ft^2)
hin = 527.9; % Btu/lbm
Tin = 533.063; % deg_F
qpp_bar = 144032; % Btu/(hr*ft^2)
pitch = 0.640 / 12; % ft
d_o = 0.493 / 12; % ft
r_o = d_o / 2;
th = 0.034 / 12; % ft
d_i = d_o - 2 * th;
r_i = d_i / 2;
k_clad = 9.6; % Btu/(hr*ft*R)
d_pellet = 0.416 / 12; % ft
h_gap = 1000; % Btu/(hr*ft^2*R)
H = 148 / 12; % ft
gamma = 0.97;
Fz = 1.4;

% Data from SteamTab
mu_f = 0.219226; % lbm/(ft*hr)
mu_g = 0.0460478; % lbm/(ft*hr)
cp_f = 1.29859 ; % Btu/(lbm*R)
k_f  = 0.328877; % Btu/(hr*ft*R)
gc = 32.2 * 3600^2;
J = 778;
sigma = 0.0017; % lbf/ft (approx)
h_f = 548.746; % Btu/hr
h_g = 1191.05; % Btu/hr
h_fg = h_g - h_f;
nu_f = 0.0217438; % ft^3/lbm
nu_g = 0.426878 ; % ft^3/lbm
nu_fg = nu_g - nu_f;
rho_f = 1 / nu_f;
rho_g = 1 / nu_g;
Tsat = 549.432; % deg_F
Tsat_abs = Tsat + 459.67; % R 

length_z = 1e6;
z_mesh = linspace(0,H,length_z);

% shape_fun = zeros(length_z,1);
% for i = 1:length(z_mesh)
	% shape_fun(i) = ((pi * (H - z_mesh(i))) / (H)) * sin((pi * (H - z_mesh(i))) / (H));
% end
% diff_shape_fun = diff(shape_fun);
% max_loc = 0;
% diff_min = 10;
% shape_max = 0;
% for i = 2:((length(z_mesh) - 2))
	% if (abs(diff_shape_fun(i))) < diff_min 
		% diff_min = diff_shape_fun(i);
		% max_loc = i;
		% shape_max = shape_fun(i);
	% end
% end
% zmax = z_mesh(max_loc);

% Normalization, et cetera
syms z
eqn = 0 == diff(((pi * (H - z)) / (H)) * sin((pi * (H - z)) / (H)),z);
zmax = double(vpasolve(eqn,z));
shape_max = ((pi * (H - zmax)) / (H)) * sin((pi * (H - zmax)) / (H));

syms lambda z
eqn = Fz == H * shape_max / int(((pi * (H + lambda - z)) / (H + 2 * lambda)) * sin((pi * (H + lambda - z)) / (H + 2 * lambda)),z,0,H);
lambda = vpasolve(eqn,lambda);
lambda = double(lambda);
He = H + 2 * lambda;
q0pp = qpp_bar * Fz / shape_max;

% Based on Weisman form for Liquid Only
xC = 0.042 * (pitch / d_o) - 0.024;
Ax = pitch^2 - (pi * d_o^2) / 4;
Pw = pi * d_o;
De = (4 * Ax) / Pw;
mdot = G * Ax;

% From now on, I'm only working at the location of maximum heat flux
qpp = q0pp * ((pi * (H + lambda - zmax)) / (He)) * sin((pi * (H + lambda - zmax)) / (He));
syms z
int_temp = double(int(q0pp * ((pi * (H + lambda - z)) / (He)) * sin((pi * (H + lambda - z)) / (He)),z,0,zmax));
hinf = hin + (1 / (mdot * gamma)) * pi * d_o * int_temp;
Tinf = Tin + (1 / (mdot * cp_f * gamma)) * pi * d_o * int_temp;
if Tinf > Tsat
	Tinf = Tsat;
end
x = (hinf - h_f) / h_fg;

num_w = 2e4;
Tw_dist = zeros(num_w,1);
T0_dist = zeros(num_w,1);
kint = @(T) (3978.1 * log((692.6 + T)/692.6) + (6.02366e-12 / 4) * ((T + 460)^4 - 460^4));
for w = 1:num_w
	ooxtt = (x / (1 - x))^0.9 * (rho_f / rho_g)^0.5 * (mu_g / mu_f)^0.1;
	f_upper_point = interpolate(f_upper,ooxtt);
	f_lower_point = interpolate(f_lower,ooxtt);
	weight_f = rand;
	f_distr = (f_upper_point - f_lower_point) * weight_f + f_lower_point;
	F = f_distr;

	hlo = 0.023 * ((G * (1 - x) * De) / mu_f)^0.8 * (cp_f * mu_f / k_f)^0.4 * (k_f/De) * F;

	Re_LP = ((G * (1 - x) * De) / mu_f);
	Re_2f = Re_LP * F^1.25;
	if (Re_2f > s_upper(1,size(s_upper,2))) || (Re_2f < s_upper(1,1)) || (Re_2f > s_lower(1,size(s_lower,2))) || (Re_2f < s_lower(1,1))
		S = 0.9622 - 0.5822 * atan(Re_2f / 6.18e4);
		% continue
	else
		s_upper_point = interpolate(s_upper,Re_2f);
		s_lower_point = interpolate(s_lower,Re_2f);
		weight_s = rand;
		s_distr = (s_upper_point - s_lower_point) * weight_s + s_lower_point;
		S = s_distr;
	end

	syms Tco
	h2f = 0.00122 * ((k_f^0.79 * cp_f^0.45 * rho_f^0.49 * gc^0.25) / (sigma^0.5 * mu_f^0.29 * h_fg^0.24 * rho_g^0.24)) * ((h_fg * J) / (Tsat_abs * nu_fg))^0.75 * (Tco - Tsat)^0.99 * S;
	eqn = qpp == hlo * (Tco - Tinf) + h2f * (Tco - Tsat);
	Tw_dist(w) = double(vpasolve(eqn,Tco,Tinf));
	
	Tci = Tw_dist(w) + ((qpp * r_i) / k_clad) * log(r_o / r_i);
	Ts  = Tci + qpp * (r_o / r_i) * (1 / h_gap);
	surface_condint = kint(Ts);
	center_condint = surface_condint + qpp * (d_o / 4);
	syms Tcl
	eqn = center_condint == kint(Tcl);
	T0_dist(w) = double(vpasolve(eqn,Tcl,Ts));
	
	disp([w])
end
%%

close all
figure(1)
histogram(Tw_dist,100)
title('T_W Probability Distribution')
xlabel('Temperature [°F]')
ylabel('Frequency')
figure(2)
histogram(T0_dist,100)
title('T_{CL} Probability Distribution')
xlabel('Temperature [°F]')
ylabel('Frequency')

%%
fprintf('Tw_dist \n');
fprintf('mean   = %.4f \n',mean(Tw_dist));
fprintf('stdev  = %.4f \n',std(Tw_dist));
fprintf('stdev%% = %.4f%%\n',(std(Tw_dist) * 100 / mean(Tw_dist)));
fprintf('median = %.4f \n',median(Tw_dist));
fprintf('min    = %.4f \n',min(Tw_dist));
fprintf('max    = %.4f \n',max(Tw_dist));
fprintf('spread = %.4f \n',(max(Tw_dist) - min(Tw_dist)));
fprintf('\n');
fprintf('T0_dist \n');
fprintf('mean   = %.4f \n',mean(T0_dist));
fprintf('stdev  = %.4f \n',std(T0_dist));
fprintf('stdev%% = %.4f%%\n',(std(T0_dist) * 100 / mean(T0_dist)));
fprintf('median = %.4f \n',median(T0_dist));
fprintf('min    = %.4f \n',min(T0_dist));
fprintf('max    = %.4f \n',max(T0_dist));
fprintf('spread = %.4f \n',(max(T0_dist) - min(T0_dist)));