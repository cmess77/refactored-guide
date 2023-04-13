delta_pr_fan = (2.2 - 1.5)/30;
delta_pr_comp = (28 - 20) / 30;
pr_fan = 1.5:delta_pr_fan:2.2;
pr_comp = 20:delta_pr_comp:28;
bypass_ratio = 5.9;
I = zeros(31,31);
TSFC = zeros(31,31);
eta_th = zeros(31,31);
eta_p = zeros(31,31);
eta_o = zeros(31,31);

[XX, YY] = meshgrid(pr_fan, pr_comp);
for ii = 1:numel(pr_fan)
  for j = 1:numel(pr_comp)
    trackXX = XX(ii,1);
    trackYY = YY(j,1);
    [I(ii,j),TSFC(ii,j),eta_th(ii,j),eta_o(ii,j),eta_p(ii,j)] = compute_values(XX(ii,j),YY(ii,j),bypass_ratio,0);
  endfor
endfor

surface(XX,YY,eta_o);
xlabel('Fan Pressure Ratio'); ylabel('Compressor Pressure Ratio'); zlabel('Overall Efficiency');
title('Overall Efficiency as a function of fan and compressor pressure ratios');

bypass_ratio_range = 1:0.1:30;
[~,~,~,~,~] = compute_values(2,24,bypass_ratio_range,1);