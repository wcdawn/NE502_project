program uncertainty
use tableio
use calculations
IMPLICIT NONE

integer :: i, k, ios

! CONSTANT
real(8),parameter :: pi = 3.1415926535987d0

! CORRELATION INFORMATION
integer :: length_f_upper, length_f_lower, length_s_upper, length_s_lower
real(8),dimension(:,:),allocatable :: f_upper, f_lower, s_upper, s_lower
real(8) :: f_upper_point, f_lower_point, f_distr, s_upper_point, s_lower_point, s_distr
real(8) :: weight_f, weight_s

! PROPERTIES
real(8) :: P, G, hin, Tin, qpp_bar, pitch, d_o, d_i, th, k_clad, d_pellet, h_gap, H, gamma, Fz
real(8) :: mu_f, mu_g, cp_f, k_f, gc, J, sigma, h_f, h_g, h_fg, nu_f, nu_g, nu_fg, rho_f, rho_g, Tsat, Tsat_abs

! WEISMANN & MORE PROPERTIES/DIMENSIONS
real(8) :: C, Ax, Pw, De, mdot

! CALCULATIONS
integer :: length_z, max_loc, top, bot
real(8) :: q0pp, shape_max, diff_min
real(8) :: dz
real(8) :: search_up, search_dn, search, tol, sol_up, sol_dn
real(8) :: lambda, He, lambda_up, lambda_dn, lambda_a_up, lambda_a_dn, lambda_b_up, lambda_b_dn
real(8) :: int_temp
real(8) :: Tw_lo, Tw_hi, Tw_up, Tw_dn, RHS_up, RHS_dn, LHS, max_Tw
real(8) :: ooxtt, hlo, h2f_up, h2f_dn, h2f_top, h2f_bot, F, S, Re_LP, Re_2f
real(8),dimension(:),allocatable :: z_mesh, shape_fun, diff_shape_fun, lambda_fun_up, lambda_fun_dn, int_lambda_fun
real(8),dimension(:),allocatable :: qpp, hinf, Tinf, x, Tw, T0

character(80) :: fname_f_upper, fname_f_lower, fname_s_upper, fname_s_lower, fname_Tw, fname_Tinf, fname_qpp

101 format(a)
102 format(e12.6)
103 format(i3)

! Take data from files containing points for correlations in F & S
fname_f_upper = 'f_upper.csv'
fname_f_lower = 'f_lower.csv'
fname_s_upper = 's_upper.csv'
fname_s_lower = 's_lower.csv'

call datainput(fname_f_upper,f_upper,length_f_upper)
call datainput(fname_f_lower,f_lower,length_f_lower)
call datainput(fname_s_upper,s_upper,length_s_upper)
call datainput(fname_s_lower,s_lower,length_s_lower)

! Data from Doster
P = 1040.0d0 ! psia
G = 1.21d6 ! lbm/(hr*ft^2)
hin = 527.9d0 ! Btu/lbm
Tin = 533.063d0 ! deg_F
qpp_bar = 144032.0d0 ! Btu/(hr*ft^2)
pitch = 0.640d0 / 12.0d0 ! ft
d_o = 0.493d0 / 12.0d0 ! ft
th = 0.034d0 / 12.0d0 ! ft
d_i = d_o - 2 * th
k_clad = 9.6d0 ! Btu/(hr*ft*R)
d_pellet = 0.416d0 / 12.0d0 ! ft
h_gap = 1000.0d0 ! Btu/(hr*ft^2*R)
H = 148.0d0 / 12.0d0 ! ft
gamma = 0.97d0
Fz = 1.4d0

! Data from SteamTab
mu_f = 0.219226d0 ! lbm/(ft*hr)
mu_g = 0.0460478d0 ! lbm/(ft*hr)
cp_f = 1.29859d0  ! Btu/(lbm*R)
k_f  = 0.328877d0 ! Btu/(hr*ft*R)
gc = 32.2d0 * 3600.0d0 ** 2.0d0
J = 778.0d0
sigma = 0.0017 ! lbf/ft (approx)
h_f = 548.746d0 ! Btu/hr
h_g = 1191.05d0 ! Btu/hr
h_fg = h_g - h_f
nu_f = 0.0217438d0 ! ft^3/lbm
nu_g = 0.426878d0  ! ft^3/lbm
nu_fg = nu_g - nu_f
rho_f = 1.0d0 / nu_f
rho_g = 1.0d0 / nu_g
Tsat = 549.432d0 ! deg_F
Tsat_abs = Tsat + 459.67 ! R

! Build basic functions & meshes, etc.
length_z = 100000
allocate(z_mesh(length_z))
allocate(shape_fun(length_z))

dz = (H - 0.0d0) / (real(length_z,8) - 1.0d0)
z_mesh(1) = 0.0d0
do i = 2,length_z
	z_mesh(i) = z_mesh(i - 1) + dz
	shape_fun(i) = ((pi * (H - z_mesh(i))) / (H)) * sin((pi * (H - z_mesh(i))) / (H))
enddo
call forwarddiff(shape_fun,length_z,dz,diff_shape_fun)
max_loc = 0
diff_min = 10.0d0
shape_max = 0.0d0
do i = 2,length_z-2
	if (abs(diff_shape_fun(i)) .lt. diff_min) then
		diff_min = diff_shape_fun(i)
		max_loc = i
		shape_max = shape_fun(i)
	endif
enddo

tol = 1.0d-5
search = 1.0d1
k = length_z / 2
allocate(lambda_fun_up(length_z))
allocate(lambda_fun_dn(length_z))
top = length_z
bot = 1

do while (search .gt. tol)
	lambda_up = z_mesh((top + k) / 2)
	lambda_dn = z_mesh((bot + k) / 2)
	do i = 1,length_z
		lambda_a_up = ((pi * (H + lambda_up - z_mesh(i))) / (H + 2.0d0 * lambda_up))
		lambda_b_up = sin(((pi * (H + lambda_up - z_mesh(i))) / (H + 2.0d0 * lambda_up)))
		lambda_fun_up(i) = lambda_a_up * lambda_b_up
		lambda_a_dn = ((pi * (H + lambda_dn - z_mesh(i))) / (H + 2.0d0 * lambda_dn))
		lambda_b_dn = sin(((pi * (H + lambda_dn - z_mesh(i))) / (H + 2.0d0 * lambda_dn)))
		lambda_fun_dn(i) = lambda_a_dn * lambda_b_dn
	enddo
	call trapz(lambda_fun_up,length_z,z_mesh,0.0d0,(H - 1.0d-10),sol_up)
	call trapz(lambda_fun_dn,length_z,z_mesh,0.0d0,(H - 1.0d-10),sol_dn)
	sol_up = (H * shape_max) / sol_up
	sol_dn = (H * shape_max) / sol_dn
	search_up = abs(Fz - sol_up) / Fz
	search_dn = abs(Fz - sol_dn) / Fz
	if (search_dn .lt. search_up) then
		search = search_dn
		bot = bot
		top = k		
		k = (top + bot) / 2
		lambda = lambda_dn
	else
		search = search_up
		bot = k
		top = top
		k = (top + bot) / 2
		lambda = lambda_up
	endif
	write(*,102) lambda
enddo
He = H + 2.0d0 * lambda

! Based on Weisman form for Liquid Only
C = 0.042d0 * (pitch / d_o) - 0.024
Ax = pitch ** 2.0d0 - (pi * d_o ** 2.0d0) / 4.0d0
Pw = pi * d_o
De = (4.0d0 * Ax) / Pw
mdot = G * Ax

allocate(qpp(length_z))
allocate(hinf(length_z))
allocate(Tinf(length_z))
allocate(x(length_z))
q0pp = qpp_bar * Fz / shape_max
do i = 2,length_z
	qpp(i) = q0pp * ((pi * (H + lambda - z_mesh(i))) / (He)) * sin((pi * (H + lambda - z_mesh(i))) / (He))
	call trapz(qpp,length_z,z_mesh,0.0d0,z_mesh(i),int_temp)
	hinf(i) = hin + (1.0d0 / (mdot * gamma)) * pi * d_o * int_temp
	Tinf(i) = Tin + (1.0d0 / (mdot * cp_f * gamma)) * pi * d_o * int_temp
	if (Tinf(i) .gt. Tsat) then
		Tinf(i) = Tsat
	endif
	x(i) = (hinf(i) - h_f) / (h_fg)
enddo






write(*,101)






allocate(Tw(length_z))

do i = 2,length_z
	if (x(i) .lt. 0.0d0) then
		cycle
	endif
	ooxtt = (x(i) / (1 - x(i))) ** 0.9d0 * (rho_f / rho_g) ** 0.5d0 * (mu_g / mu_f) ** 0.1d0
	if (ooxtt .le. 0.10d0) then
		F = 1.0d0
	else
		F = 2.35d0 * (ooxtt + 0.213d0) ** 0.736d0
	endif
	
	hlo = 0.023d0 * ((G * (1 - x(i)) * De) / mu_f) ** 0.8d0 * (cp_f * mu_f / k_f) ** 0.4d0 * (k_f / De) * F
	
	Re_LP = ((G * (1 - x(i))) * De) / mu_f
	Re_2f = Re_LP * F ** 1.25d0
	
	S = 0.9622d0 - 0.5822d0 * atan(Re_2f / 6.18d4)
	
	! guess
	Tw_lo = Tinf(i)
	Tw_hi = 6.5d2
	Tw(i) = (Tw_lo + Tw_hi) / 2.0d0
	search = 1.0d1
	tol = 1.0d-6
	k = 0
	! if (i .eq. 1386) then 
		! stop
	! endif
	! write(*,'(i12)') i
	do while (search .gt. tol)
		k = k + 1
		if (k .ge. 1001) then
			write(*,'(i12)') i
			exit
		endif
		Tw_up = (Tw(i) + Tw_hi) / 2.0d0
		Tw_dn = (Tw(i) + Tw_lo) / 2.0d0
		h2f_top = k_f ** 0.79d0 * cp_f ** 0.45d0 * rho_f ** 0.49d0 * gc ** 0.25d0
		h2f_bot = sigma ** 0.5d0 * mu_f ** 0.29d0 * h_fg ** 0.24d0 * rho_g ** 0.24d0
		h2f_up = 0.00122d0 * ((h2f_top) / (h2f_bot)) * ((h_fg * J) / (Tsat * nu_fg)) ** 0.75d0 * (Tw_up - Tsat) ** 0.99 * S
		h2f_dn = 0.00122d0 * ((h2f_top) / (h2f_bot)) * ((h_fg * J) / (Tsat * nu_fg)) ** 0.75d0 * (Tw_dn - Tsat) ** 0.99 * S
		RHS_up = hlo * (Tw_up - Tinf(i)) + h2f_up * (Tw_up - Tsat)
		RHS_dn = hlo * (Tw_dn - Tinf(i)) + h2f_dn * (Tw_dn - Tsat)
		LHS = qpp(i)
		search_up = abs(RHS_up - LHS) / LHS
		search_dn = abs(RHS_dn - LHS) / LHS
		if (search_dn .lt. search_up) then
			search = search_dn
			Tw_lo = Tw_lo
			Tw_hi = Tw(i)
			Tw(i) = Tw_dn
		else
			search = search_up
			Tw_lo = Tw(i)
			Tw_hi = Tw_hi
			Tw(i) = Tw_up
		endif
		! write(*,102) Tw(i)
	enddo
	! write(*,103) k
	
enddo

fname_qpp = 'qpp.csv'
call csvout(fname_qpp,z_mesh,qpp,length_z)
fname_Tinf = 'Tinf.csv'
call csvout(fname_Tinf,z_mesh,Tinf,length_z)
fname_tw = 'tw.csv'
call csvout(fname_tw,z_mesh,Tw,length_z)

max_Tw = 0.0d0
do i = 1,length_z
	if (Tw(i) .gt. max_Tw) then
		max_Tw = Tw(i)
	endif
enddo
	
write(*,102) max_Tw
write(*,102) Tw(max_loc)
write(*,'(i12)') max_loc

endprogram uncertainty