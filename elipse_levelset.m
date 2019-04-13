function [mapa] = elipse_levelset(Vstup_obraz, pocatecni_maska, max_it)

% Aktiv_kontur je funkce/algoritmus pro segmetaèní metodu aktivních hranových kontur. 

% Vstup_obraz ... vstupní obraz, na kterém se provede segmentace 
% pocatecni_maska ... která je urèena manuálnì
% max_it ... maximální poèet iterací (èím víc tím pøesnìjší, ale pomalejší)
% mapa ... výstup funkce, tedy mapa, kde je ohranièena oblast zájmu


%% __________________________________________________________________________
% Flexibilní kontury - Aktivní hranové kontury

%************************************************************************** 
% Kód pøevzaný z: 
%     -> https://www.mathworks.com/matlabcentral/fileexchange/19567-active-contour-segmentation?s_tid=srchtitle
%     -> Shawn Lankton

% Tato technika deformuje poèáteèní køivku tak, že oddìluje popøedí od pozadí 
% na základì prostøedkù obou oblastí. Technika je velmi robustní k inicializaci 
% a dává velmi pìkné výsledky.

%Pøeddefinování promìnných
alpha = 0.2; 
display = true;
mapa = mask2phi(pocatecni_maska); % Vytvoøení zvýraznìné mapy z poèáteèní masky
  
%Hlavní cyklus
for its = 1:max_it   
    idx = find(mapa <= 1.2 & mapa >= -1.2);  %nalezneme køivku úzkého pásma
    
    %Nalezení vnítøního a vnìjšího prùmìru
    vnejsi_body = find(mapa<=0);           %vnìjší body
    vnitrni_body = find(mapa>0);           %vnítøní body    
    u = sum(Vstup_obraz(vnejsi_body))/(length(vnejsi_body)+ eps);   % vnítøní_prùmìr
    v = sum(Vstup_obraz(vnitrni_body))/(length(vnitrni_body)+ eps); % vnìjší prùmìr
                                                            % eps ... relativní pøesnost promìnlivého bodu    
    F = (Vstup_obraz(idx)-u).^2-(Vstup_obraz(idx)-v).^2;    % síla/váha obrazové informace
    curvature = get_curvature(mapa,idx);  % síla od zakøivení (penalizace)
    dphidt = F./max(abs(F)) + alpha*curvature;  %stoupání gradientu pro minimalizaci energie
    
    dt = .45/(max(abs(dphidt))+eps);  %Udržování Courant–Friedrichs–Lewy (CFL) stavu
    mapa(idx) = mapa(idx) + dt.*dphidt; %Rozvíjení køivky
    mapa = sussman(mapa, .5); %Udržuje SDF hladký 
                              % SDF =  the Shape Diameter Function
    %Prostøední výstup
%     if((display>0)&&(mod(its,20) == 0)) 
%       showCurveAndPhi(Vstup_obraz,mapa,its);  
%     end
end
  
  %Vytvoøení masky pro SDF
  seg = mapa<=0; % Získání masky z nastavené úrovnì
  

%__________________________POMOCNÉ FUNKCE _________________________________

%Zobrazení
% function showCurveAndPhi(Vystup, phi, i)
%   imshow(Vstup_obraz,'initialmagnification',200,'displayrange',[0 255]); hold on;
%   contour(phi, [0 0], 'g','LineWidth',4);
%   contour(phi, [0 0], 'k','LineWidth',2);
%   hold off; title([num2str(i) ' Iterations']); drawnow;
  
%Pøevede masku na SDF
function phi = mask2phi(init_a)
  phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
  
%Vypoèítá zakøivení podél SDF
function zakriveni = get_curvature(phi,idx)
    [dimy, dimx] = size(phi);        
    [y x] = ind2sub([dimy,dimx],idx);  % získání pøíslušných indexù

    %získání dolních indexù sousedù
    ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

    %Kontrola hranic
    ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;    

    %Zíkání indexù pro 8 sousedù
    idup = sub2ind(size(phi),yp1,x);    
    iddn = sub2ind(size(phi),ym1,x);
    idlt = sub2ind(size(phi),y,xm1);
    idrt = sub2ind(size(phi),y,xp1);
    idul = sub2ind(size(phi),yp1,xm1);
    idur = sub2ind(size(phi),yp1,xp1);
    iddl = sub2ind(size(phi),ym1,xm1);
    iddr = sub2ind(size(phi),ym1,xp1);
    
    %Získání centrálních derivátù SDF na x, y
    phi_x  = -phi(idlt)+phi(idrt);
    phi_y  = -phi(iddn)+phi(idup);
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
    phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
             +0.25*phi(iddr)+0.25*phi(idul);
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;
    
    %Výpoèet zakøivení (Kappa)
    zakriveni = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
              (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        

%Na úrovni re-inicializace pomocí sussman method
function D = sussman(D, dt)
  % dopøedu/dozadu rozdíly (diference)
  a = D - shiftR(D); % zpìt
  b = shiftL(D) - D; % vpøed
  c = D - shiftD(D); % zpìt
  d = shiftU(D) - D; % vpøed
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
  
%Celé maticové derivace
function shift = shiftD(M)
  shift = shiftR(M')';

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
  shift = shiftL(M')';
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    


