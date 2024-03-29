clc
clear all
close all

% experimental sample parameter
b = 5*1e-3;     % width of the beam [m]
t = 0.5*1e-3;   % thickness of the beam [m]
l = 4*1e-3;     % hinge lengths [m]
H = 5*1e-3;     % tilted height [m]
L = 15*1e-3;    % span [m]
E = 0.5*1e6;    % modulus [Pa]


A = b*t;
II = b*t^3/12;

k = E*A/l;
kr = E*II/l;

syms W F d

lb = sqrt( (L-l)^2 + H^2  );
alpha = asin( H/lb  );
alpha_p = asin( (H-d)/lb  );

d_alpha = alpha-alpha_p;
dl = (L-l) - sqrt( (lb)^2 - (H-d)^2  );

W_left = 1/2*kr*d_alpha^2 + 1/2*k*dl^2;

dW_left = diff(W_left,d);
F_exp = 2*dW_left;

 

d_plot = 0:H/50:2.3*H;
F_plot = double( subs( F_exp,d,d_plot )   );
%%
fontsize = 20;
figure(1)
plot(d_plot*1000,F_plot,'k','LineWidth',2)
hold on
plot([0,12],[0,0],'k--')
xlim([0 12]);

xlabel('Displacement [mm]','FontSize',fontsize,'FontWeight','bold','FontName','Arial')
ylabel('Force [N]','FontSize',fontsize,'FontWeight','bold','FontName','Arial')
set(gca,'FontName','Calibri','FontSize',fontsize,'FontWeight',...
    'bold');
 
% find the stationary point
dy = diff(F_plot) ./ diff(d_plot);
signs = sign(dy);
dy_sign = diff(signs);
SP = find(abs(dy_sign)==2); % position of stationary point
vSP = F_plot(SP);
if isempty(vSP)==true
    type = 0;   % no snapping
    disp('Monotonic');
else
    if vSP(1)*vSP(2)<0
        type = 2;  % bistable where 2 stationary have different signs
        disp('Bistable')
    else
        type = 1; % monostable where 2 stationary have the same positive signs
        disp('Monostable')
    end
end


% generate energy plot
W_plot(1) = 0;
d_plot(1) = 0;
for i=2:length(d_plot)
   W_plot(i) = W_plot(i-1)+ (d_plot(i)-d_plot(i-1))*F_plot(i); 
end

figure(2)
plot(d_plot*1e3,W_plot,'k','LineWidth',2)
hold on
xlabel('Displacement [mm]','FontSize',fontsize,'FontWeight','bold','FontName','Arial')
ylabel('Energy [N.m]','FontSize',fontsize,'FontWeight','bold','FontName','Arial')
set(gca,'FontName','Calibri','FontSize',fontsize,'FontWeight',...
    'bold');
xlim([0 12]);
