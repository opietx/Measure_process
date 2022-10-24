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
r = 69e-9 ; %dipole matrix moment

Delta = 0.40e-9 ; %splitting at the anticrossing point
Rk = 25.813e3 ; %quantum resistance

GammaValley_ph = @(omega) (ht.^2./(2.*pi)).*(ht.*omega).^5 .* r.^2 ./ (4.*pi.*rho.*ht.^6) .* (It./vt.^7 + 2.*Il./vl.^7) .* e ;
GammaValley_jh = @(omega) 2.* R / Rk * omega * ht^2 * 3*r.^2 ./ l^2 ;

%% Simulation B dependence
B_res = 400 ;
B_list=linspace(0.5, 20, B_res) ;
T=0.1 ; %Temperature
rate = zeros(4,4, B_res) ;

for ii=1:length(B_list)
    fprintf('\nPercentage: %.2f', ii)
    B=B_list(ii);
    omega = g*ub*B/ht ;
    
    % Matrix of coefficients between states
    a1 = -(Evs-ht*omega) ./ sqrt((Evs-ht*omega)^2+Delta^2) ;
    a2 = -(Evs+ht*omega) ./ sqrt((Evs+ht*omega)^2+Delta^2) ;
    
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
    E = zeros(5, 1) ;
    E(1) = -(Evs + ht*omega)./2 ;
    E(2) = -sqrt((Evs - ht*omega).^2+Delta^2) / 2 ;
    E(3) = +sqrt((Evs - ht*omega)^2+Delta^2) / 2 ;
    E(4) = (Evs + ht*omega)./2 ;
       
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
                
                rate(i,j,ii) = rate_johnson +rate_phonon ;
                
            end
        end
    end
    
    
end

%% Solving rate equations
resolution = B_res ;% temp_res

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
figure
loglog(B_list, abs(eigenvalues(:,2)), '-')