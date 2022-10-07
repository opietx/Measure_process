%% Constants
clear all

pi = 3.14159265359 ;
ub = 5.788381e-5 ;
e = 1.60217662e-19 ;
g=2 ;
m0 = 9.10938356e-31 ;
kb = 8.6173303e-5 ;
ht = 6.582119e-16 ;
omegaDebye = 13.81e13 ; 
rho = 2200 ; 
vl = 3750 ;
vt = 5900 ;

PhiU = 8.77 ;

PhiD = 5 ;

It = 16/105 * PhiU^2 ;
Il = 4*(PhiU^2/35 + 2*PhiD*PhiU/15 + PhiD^2/3) ;
mt = 0.19*m0 ;

bose = @(omega, T) 1./(exp(ht.*omega./(kb.*T/1000))-1) ;

%% Input Parameters
%Johnson noise parameters
l=100e-9 ; 
R=10e3 ;

Evs = 275e-6  ; %valley splitting
r = 69e-9 ; %dipole moment

Delta = 0.40e-9 ; %splitting at the anticrossing point
Rk = 25.813e3 ; %quantum resistance

GammaValley_ph = @(omega) (ht.^2./(2.*pi)).*(ht.*omega).^5 .* r.^2 ./ (4.*pi.*rho.*ht.^6) .* (It./vt.^7 + 2.*Il./vl.^7) .* e ;
GammaValley_jh = @(omega) 2.* R / Rk * omega * ht^2 * 3*r.^2 ./ l^2 ;

%% Simulation for fixed magnetic field
B=1.5 ;
omega = g*ub*B/ht ;

% Matrix of coefficients between states
a1 = -(Evs-ht*omega) ./ sqrt((Evs-ht*omega).^2+Delta^2) ;
a2 = -(Evs+ht*omega) ./ sqrt((Evs+ht*omega).^2+Delta^2) ;

SinGamma = @(x) sqrt((1-x)/2) ;
CosGamma = @(x) sqrt((1+x)/2) ;

valley_coupling = zeros(4,4) ;
valley_coupling(1,2) = -CosGamma(a1)*SinGamma(a2) - SinGamma(a1)*CosGamma(a2) ;
valley_coupling(1,3) = SinGamma(a1)*SinGamma(a2) - CosGamma(a1)*CosGamma(a2) ;
valley_coupling(1,4) = SinGamma(a2)*CosGamma(a2) ;
valley_coupling(2,3) = SinGamma(a1)*CosGamma(a1) ;
valley_coupling(2,4) = -CosGamma(a1)*CosGamma(a2) + SinGamma(a1)*SinGamma(a2) ;
valley_coupling(3,4) = SinGamma(a1)*CosGamma(a2) + CosGamma(a1)*SinGamma(a2) ;

valley_coupling = valley_coupling + valley_coupling' ;

% Energy of states
E = zeros(5) ;
E(1) = -(Evs + ht*omega)./2 ;
E(2) = -sqrt((Evs - ht*omega).^2+Delta^2) / 2 ;
E(3) = +sqrt((Evs - ht*omega)^2+Delta^2) / 2 ;
E(4) = (Evs + ht*omega)./2 ;

% Getting Rates
index=[1,2,3,4] ;
temp_res = 100 ;
Gamma = zeros(4,4) ;
rate = zeros(4,4, temp_res) ; % ij means transition from i to j
rate_first_order = zeros(4,4, temp_res) ;
rate_second_order = zeros(4,4, temp_res) ;
rate_second_order_jh = zeros(4,4, temp_res) ;
rate_second_order_twophonon = zeros(4,4, temp_res) ;
rate_orbit= zeros(1, temp_res) ;
Gammac=zeros(4, temp_res) ;
T = linspace(10, 2000, temp_res) ;
T_list = T ;

%% Calculate relaxation rates.
% Calculate first order relaxation rate, it is needed to evaluate the
% second order

h = waitbar(0, '1', 'Name', 'Getting T1 ...') ;

for ii=1:length(T_list)
    %fprintf('\nPercentage: %.2f', ii)
    progress = ii / length(T_list) ;
    waitbar(progress, h, sprintf('Progress ... %.1f %%', 100*progress))
    T=T_list(ii) ;
    
    for i=1:4
        for j=1:4
            if j~=i

                if j>i
                    rate_phonon = 2.*pi./ht^2 .* valley_coupling(i, j)^2 .* GammaValley_ph(abs(E(i)-E(j))/ht) .* bose(abs(E(i)-E(j))/ht, T) ;
                    rate_johnson = 2.*pi./ht^2 .* valley_coupling(i, j)^2 .* GammaValley_jh(abs(E(i)-E(j))/ht) .* bose(abs(E(i)-E(j))/ht, T) ;
                else
                    rate_phonon = 2.*pi./ht^2 .* valley_coupling(i, j)^2 .* GammaValley_ph(abs(E(i)-E(j))/ht) .* (1+bose(abs(E(i)-E(j))/ht, T)) ;
                    rate_johnson =  2.*pi./ht^2 .* valley_coupling(i, j)^2 .* GammaValley_jh(abs(E(i)-E(j))/ht) .* (1+bose(abs(E(i)-E(j))/ht, T)) ;
                end

                rate_first_order(i, j, ii) = rate_phonon  + rate_johnson;

            end
        end         
    end

    % Get state lifetime
    tmp = sum(rate_first_order, 2) ;
    Gammac(:, ii) = tmp(:,1,ii) + 1000 ;  %% 1 ms is the maximum lifetime limited by readout

    omega_min = 100  ;

    % Calculate Raman, Orbach process
    for i=1:4
        for j=1:4
            if j~=i
                %fprintf('Indeces %d, %d', i, j)
                m=index(index~=i) ;
                m=m(m~=j) ;
                map = map_Raman([int2str(i), int2str(j)] ) ;

                if length(map(:,1))==1

                    coeff1 = @(omega1) valley_coupling(j,m(1)).*valley_coupling(m(1),i)./(E(m(1))-E(i)-sign(E(m(1))-E(i)).*ht.*omega1 + 0.5.*ht.*1i.*Gammac(m(1), ii)) ;
                    coeff2 = @(omega1) valley_coupling(j,m(2)).*valley_coupling(m(2),i)./(E(m(2))-E(i)-sign(E(m(2))-E(i)).*ht.*omega1 + 0.5.*ht.*1i.*Gammac(m(2), ii)) ;
                    coeff  = @(omega1) abs(coeff1(omega1) + coeff2(omega1)).^2 ;

                    if map(1)==map(2)
                        omega_max = abs(E(i)-E(j))/ht;
                        omega_sign = -1 ;
                    else
                        omega_max = omegaDebye ;
                        omega_sign = +1 ;
                    end

                    func_Raman = @(omega1, omega2 ) 2*pi/ht^2 .* (map(1)+bose(omega1, T)) .* (map(2)+bose(omega2, T)) .* GammaValley_ph(omega1) .* GammaValley_ph(omega2) .* coeff(omega1) ;
                    func_Raman_singlevar =  @(omega1) func_Raman(omega1, (2*map(2)-1)*(E(i) - E(j))/ht + omega_sign.*omega1) ;

                    rate_second_order(i,j,ii) = Integrate_Raman(func_Raman_singlevar, omega_min, omega_max, abs(E(m(1))-E(i))/ht);

                else

                    coeff1 = @(omega1) valley_coupling(j,m(1)).*valley_coupling(m(1),i)./(E(m(1))-E(i)-sign(E(m(1))-E(i)).*ht.*omega1 + 0.5.*ht.*1i.*Gammac(m(1), ii)) ;
                    coeff2 = @(omega1) valley_coupling(j,m(2)).*valley_coupling(m(2),i)./(E(m(2))-E(i)-sign(E(m(2))-E(i)).*ht.*omega1 + 0.5.*ht.*1i.*Gammac(m(2), ii)) ;

                    func_Raman1 = @(omega1, omega2 ) 2*pi/ht^2 .* (map(1,1)+bose(omega1, T)) .* (map(1,2)+bose(omega2, T)) .* GammaValley_ph(omega1) .* GammaValley_ph(omega2) .* abs(coeff1(omega1)).^2 ;
                    func_Raman2 = @(omega1, omega2 ) 2*pi/ht^2 .* (map(2,1)+bose(omega1, T)) .* (map(2,2)+bose(omega2, T)) .* GammaValley_ph(omega1) .* GammaValley_ph(omega2) .* abs(coeff2(omega1)).^2 ;

                    if map(1,1)==map(1,2)
                        omega_max1 = abs(E(i)-E(j))/ht ;
                        omega1_sign1 = -1 ;
                    else
                        omega_max1 = omegaDebye ;
                        omega1_sign1 = +1 ;
                    end

                    if map(2,1)==map(2,2)
                        omega_max2 = abs(E(i)-E(j))/ht ;
                        omega1_sign2 = -1 ;
                    else
                        omega_max2 = omegaDebye ;
                        omega1_sign2 = +1 ;
                    end
                    
                    func_Raman1_singlevar =  @(omega1) func_Raman1(omega1,(2*map(1,2)-1)*(E(i) - E(j))/ht + omega1_sign1.*omega1) ;
                    func_Raman2_singlevar =  @(omega1) func_Raman2(omega1, (2*map(2,2)-1)*(E(i) - E(j))/ht + omega1_sign2.*omega1) ;

                    rate_second_order(i,j,ii) = Integrate_Raman(func_Raman1_singlevar, omega_min, omega_max1, abs(E(m(1))-E(i))/ht) + Integrate_Raman(func_Raman2_singlevar, omega_min, omega_max2, abs(E(m(2))-E(i))/ht) ;
                    if i==3 && j==1 && B==3
                        rate_second_order(i,j,ii) = integral(func_Raman1_singlevar, omega_min, omegaDebye) ;
                    end
                end
                if i==3 && j==1
                    coeff1 = @(omega1) 1./(E(5)-E(i)-ht.*omega1 + 0.5.*ht.*1i.*1e11) ;
                    func_Raman1 = @(omega1, omega2 ) 2*pi/ht^2 .* (bose(omega1, T)) .* (1+bose(omega2, T)) .* GammaValley_ph(omega1) .* GammaValley_ph(omega2) .* abs(coeff1(omega1)).^2 ;
                    func_Raman1_singlevar =  @(omega1) func_Raman1(omega1, abs(E(i) - E(j))/ht + omega1) ;
                end
                rate(i,j, ii) =  rate_second_order(i,j,ii) + rate_first_order(i,j,ii);

            end
        end
    end
    
end

%% Solving rate equations
resolution = temp_res ;

Rates_matrix = zeros(4,4, resolution) ;
eigenvalues_all = zeros(4,4) ;
eigenvectors_all = zeros(4,4) ;
eigenvalues = zeros(resolution,4) ;
eigenvectors = zeros(resolution,4) ;

for i=1:4
    for j=1:4
        if j~=i
            Rates_matrix(i,j,:) =  rate(j,i, :)  ;
        else
            Rates_matrix(i,i,:) = -sum(rate(i,:, :), 2) ;
        end
    end
end

for i=1:resolution
   [eigenvectors_all(:,:), eigenvalues_all(:,:)] = eig(reshape(Rates_matrix(:,:,i), [4,4])) ;
   
   tmp = sort(max(abs(eigenvalues_all(:,:)))) ;
   tmp2 = eigenvectors_all(:,find(max(abs(eigenvalues_all(:,:))) == tmp(2))) ;
   eigenvectors(i,:) = abs(tmp2(:)) ;
   
   eigenvalues(i,:) = sort(max(abs(eigenvalues_all(:,:))))  ;
end

%% Plotting
figure(1)
clf
semilogy(T_list ,abs(eigenvalues(:,2)))
xlabel('Temperature (mK)')
ylabel('Relaxation rate (s^-1)')
