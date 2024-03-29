clc 
clear all
close all

% fixed geometrical parameters
b = 5*1e-3;     % mm
E = 0.5*1e6;  % MPa
l = 4*1e-3;
L = 15*1e-3;

% changed parameters
H_range = (0.1:0.1:6)*1e-3;
t_range = (0.4:0.05:2)*1e-3;

 
syms d

% calculate the instability characteristics of different sets of parameters
for i = 1:length(H_range)
    for j = 1:length(t_range)
       H = H_range(i);
       t = t_range(j);
       d_plot = 0:H/50:2*H;
       F_exp =(E*b*t*(2*H - 2*d)*(l - L + (H^2 - (H - d)^2 + (L - l)^2)^(1/2)))/(l*(H^2 - (H - d)^2 + (L - l)^2)^(1/2)) + (E*b*t^3*(asin(H/(H^2 + (L - l)^2)^(1/2)) - asin((H - d)/(H^2 + (L - l)^2)^(1/2))))/(6*l*(1 - (H - d)^2/(H^2 + (L - l)^2))^(1/2)*(H^2 + (L - l)^2)^(1/2));
       F_plot = double( subs( F_exp,d,d_plot )   );
       
       % find stationary point
       dy = diff(F_plot) ./ diff(d_plot);
       signs = sign(dy);
       dy_sign = diff(signs);
       SP = find(abs(dy_sign)==2); % position of stationary point
       vSP = F_plot(SP);
       if isempty(vSP)==true
           type(i,j) = 0;   % no snapping
           map(i,j) = NaN;
       else
           if vSP(1)*vSP(2)<0
               type(i,j) = 2;  % bistable where 2 stationary have different signs
           else 
               type(i,j) = 1; % monostable where 2 stationary have the same positive signs
           end
           map(i,j) = vSP(2)/vSP(1);
       end
   
    end
end

% plot the phase map
contourf(H_range/L,t_range/L,map','--','LineColor','none');
grid off
shading flat
colorbar
hold on

% plot the critical line
contour(H_range/L,t_range/L,map',[0 0],'--','color','k','LineWidth',1.8)
grid off
shading flat


original_colormap = pink;  % Replace with your desired colormap
MapLen = size(original_colormap,1);
half_colormap = original_colormap(2*MapLen/10:7*MapLen/10, :);  % Select the first half of the colormap
colormap(half_colormap);
colorbar

fs = 20;
ylabel('Hinge thickness ratio, t/L','FontSize',fs,'FontWeight','bold','FontName','Calibri')
xlabel('Amplitude ratio, H/L','FontSize',fs,'FontWeight','bold','FontName','Calibri')
set(gca,'FontName','Calibri','FontSize',fs,'FontWeight',...
    'bold');
