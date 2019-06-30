function [dq, F0, F, ncon] = dqdt_2D_net(t, q, params)
% Параметры сети 
n   =  params.n;
c   =  params.c;
k   =  params.k;
l0  =  params.l0;
masses = params.masses;
% Параметры тверого тела
m0  =  params.m0;
J0  =  params.J0; 
c0  =  params.c0;
k0  =  params.k0;
pp  =  params.pp;
g   =  params.g;

% Положение центра масс полигона
r0   = q(1:2); 
% Угол поворота полигона
phi0 = q(3); 
% Скорость центра масс полигона
v0   = q(2*n+4:2*n+5);
% Угловая скорость полигона
w0   = q(2*n+6);
% Фактическое положение точек с учетом положения и ориентации тела  
pp   = (repmat(r0,1,size(pp,1)) + [cos(phi0) -sin(phi0);sin(phi0) cos(phi0)]*pp')';

% Точки сети 
r  = reshape(q(4:3+n*2),2,n); 
% Скорости точек сети
v  = reshape(q(7+n*2:6+n*4),2,n);

%
% Внутренние силы сети
%
% dr_i = r_i - r_{i-1}
drp = r(:,2:n) - r(:,1:n-1);
% dv_i = v_i - v_{i-1}
dvp = v(:,2:n) - v(:,1:n-1);
% e_i = unit vector from i-1 to i
edrp = drp./repmat(sqrt(sum(drp.^2)),2,1);
% 
dvpm = sum(dvp.*edrp);

l = sqrt(sum(drp.^2));
% logic array 
tensioned = l > l0;

Fp = [tensioned.*((l-l0).*c + dvpm.*k) 0];
Fm = [0 Fp(1:end-1)];

eFp = [ edrp [0;0] ];
eFm = [ [0;0] -eFp(:,1:end-1)];

F   = repmat(Fp,2,1).*eFp + repmat(Fm,2,1).*eFm;
%
% Сила со стороны твердого тела (полигона)
%
F0 = [0;0]; 
M0 = 0;

[dri, i_inside] = points2polygon(r', pp);
ndri = sqrt(sum(dri.*dri,2));
ndri = ndri(i_inside~=0);
rcon = r(:,i_inside~=0);
vcon = v(:,i_inside~=0);
ncon = size(rcon,2);
dri  = dri(i_inside~=0,:);
if ncon > 0               
        rho= rcon - repmat(r0,1,ncon);
        vw = (rho'*[0 w0;-w0 0])';        
        vi = repmat(v0,1,ncon) + vw;
        % Проекция относительно скорости на нормаль  
        dv = repmat(dot(vcon - vi,(dri./repmat(ndri,1,2))'),2,1).*(dri./repmat(ndri,1,2))';
        dvt = (vcon - vi)-dv; 
        dvtn = sqrt(sum(dvt.*dvt,1));
        dvtn(dvtn<0.0005) = 1e10;
        % Вычисляем силу на узел
        Fb   =  (c0*dri' - dv*k0);  
        Ft   =  ((dvt./repmat(dvtn,2,1))).*repmat(sqrt(sum(Fb.*Fb,1)),2,1)*params.mu;        
        
        F(:,i_inside~=0) = F(:,i_inside~=0) + Fb - 0*Ft;         
        
        % На тело
        F0 = F0 - sum(Fb,2) + sum(Ft,2); 
        Mb = sum(cross([rho; zeros(1,ncon)],[-Fb + Ft; zeros(1,ncon)]),2);       
        
        M0 = M0 + Mb(3); 
end

a   = F./repmat(masses,2,1) + repmat(g,1,n);
a0  = F0/m0 + g;
e0  = M0/J0;

a(1:2) = 0;
a(end-1:end) = 0;

dq = [v0;w0;v(:);a0;e0;a(:)];

end

