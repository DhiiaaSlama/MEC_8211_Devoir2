
%% SCHEMA NUMERIQUE D ORDRE 2 AVEC TERME SOURCE DE PREMIER ORDRE 
%% NB: Il faut avoir installé la " Signal Processing ToolBox" pour rouler ce code
clear; clc; close('all');

% importer resultats COMSOL 
file ='C:\Users\houss\Desktop\DEVOIR 1 DHIA-HOUSSEM\Devoir 2\OrderCalc.xlsx';
COMSOL = table2array(readtable(file,'Sheet','COMSOL','PreserveVariableNames',true));
ne = height(COMSOL) ; 



% inputs
D = 1 ; % diametre m
ray = D/2 ;  %rayon m
Ce = 10 ; % Condition de Dirichlet mol/m3 
k = 4*10^(-9) ; % constante de la reaction de 1er ordre s^(-1) 
Deff = 10^(-10) ; % coefficient de diffusivité effective m2/s 
s = 10^(-8) ; % terme source analytique mol/m3/s 
Nr = [10,20,40]; % Nombre de noeuds; 
res{1,length(Nr)}=[]; % resultats ;
err{1,length(Nr)}=[]; % erreur L1,L2 et Linf

tsim = 10^10; % duree de la simulation secondes
Nt = 100; % nombre de pas temps ;  
annees = tsim/(86400*365); % Duree de simulation en annee 

for x =1:length(Nr)
% steps 
ri = 0 ; rf = D/2  ; dr = (rf - ri)/Nr(x); 
ti = 0;  tf = tsim ; dt = (tf-ti)/Nt ;

h = dr; 
R = ri:dr:rf; nr = numel(R) ; 
T = ti:dt:tf; nt = numel(T) ;

[r,t] = meshgrid(R,T); r = r'; t = t'; % x and t matrices




% initialisation de la concentration
c = zeros(nr,nt);
% Matrice de coefficients du systeme lineaire Ax=b
m = zeros(nr,nr);
m(1,1) = -3 ;
m(1,2) = 4 ; 
m(1,3) = -1 ;
m(nr,nr) = 1; 

% Boucle principale 
for i=2:Nr(x)  
   % coefficients de la matrice   
 alpha = ( -(Deff*dt/(2*R(i)*h))-(Deff*dt/h^2) ) / (1-k*dt);  % i+1
 beta = (1+2*Deff*dt/h^2)/ (1-k*dt);  % i 
 gamma = (Deff*dt/(2*R(i)*h)-(Deff*dt/h^2))/ (1-k*dt);  % i-1
 
    m(i,i) = beta ;
    m(i,i+1) = alpha ; 
    m(i,i-1) = gamma ;   
end



matcoeffs= [ alpha; beta ; gamma];


for t=1:nt-1
    K = [0 ; c(2:nr-1,t) ; Ce] ;  % Cold
    X = m\K ; 
    c(:,t+1) = X;
end     


% solution analytique avec terme source constant 
C_analytique = zeros(nr,1); 
for i=1:nr
C_analytique (i) =  0.25 * s/Deff * ray^2 * (R(i)^2/ray^2 -1) + Ce; 
end 

%% Resultats Animations

% % Animation de la solution numerique dans le temps
%     for j = 1:nt
%         plot(R,c(:,j)); xlabel('Distance [m]'); ylabel('Concentration [mol]');
%         title(['Concentration Versus la Distance au Temps = ' num2str(round(T(j),3)) ' s']) ; 
%         grid('on'); drawnow; %pause(1);
%     end

%     % Affichage Solution Analytique vs Solution Numerique
%     figure(2)
%     plot(R,C_analytique)
%     hold on 
%     plot(R,c(:,nt),'o')
%     grid on
%     xlabel('Distance (m)')
%     ylabel ('Concentration (mol/m3)')
%     title(['Concentration Versus la Distance au Temps = ' num2str(round(T(j),3)) ' s']) ;
%     legend('solution analytique', 'solution numérique')
% 
%     % Outputs 
%     % solution de reference (COMSOL)
    solref = downsample(COMSOL(:,2), fix(ne/Nr(x)));
    res{1,x}=[R(:) c(:,nt) solref];
    % Calcul d'erreur 
    err1 = 1/Nr(x)*sum(abs(c(:,nt)-solref));
    err2 = sqrt(1/Nr(x)*sum(abs(c(:,nt)-solref).^2));
    errinf = max(abs(c(:,nt)-solref));
    err{1,x} = [err1; err2;errinf]; 
    
end 

%  A = cell2mat(err); errmat = array2table( A, 'RowNames', ["err1"; "err2"; "errinf"],'VariableNames',["Nr10","Nr20","Nr40"]);


 

 
 
 
