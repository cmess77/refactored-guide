function [I,TSFC,eta_th,eta_o,eta_p] = compute_values(pr_fan,pr_comp,bypass_ratios, bpr_range)
  Ta = 223.252; %K
  Pa = 2.65e4; %Pa
  Ma = 0.8; %mach
  g_a = 1.4;
  R = 287; %J/kg*K
  Ua = mach2vel(Ma,g_a,R,Ta); 
  g_fan = 1.4;
  g_to_burner = 1.4;
  g_post_burner = 1.35;

  cp_fan = cp_from_R(R,g_fan); % J/kg*K
  cp_to_burner = cp_from_R(R,g_to_burner); % J/kg*K
  cp_post_burner = cp_from_R(R,g_post_burner); % J/kg*K

  T_o4_max = 1500; %K
  T_o4 = 0; %K
  pr_burner = 0.97;
  eta_diffuser = 0.94;
  eta_compressor = 0.87;
  eta_fan = 0.92;
  eta_burner = 0.98;
  eta_turbine = 0.85;
  eta_core_nozzle = 0.97;
  eta_fan_nozzle = 0.98;
  h_c = 43000000; %J/kg
  f_st = 0.06;

  % Diffuser outlet conditions
  T_o2 = Ta*(1 + (((g_fan - 1)/2)*power(Ma,2)));
  P_o2 = Pa*power((1 + (eta_diffuser*((T_o2/Ta) - 1))),g_fan/(g_fan-1));

  %compressor outlet conditions
  T_o3 = T_o2*(1 + ((1/eta_compressor)*(power(pr_comp,(g_to_burner - 1)/g_to_burner) - 1)));
  P_o3 = P_o2 * pr_comp;
  W_comp_in = cp_to_burner*(T_o3 - T_o2);

  % burner outlet conditions
  f_max = ((T_o4_max/T_o3) - 1)/ ((h_c/(cp_post_burner*T_o3)) - (T_o4_max/T_o3));
  if f_max > f_st
    f = f_st;
    T_o4 = (1/(1+f))*(((f*eta_burner*h_c)/cp_post_burner) + T_o3);
  else
    f = f_max;
    T_o4 = T_o4_max;
  endif

  P_o4 = pr_burner * P_o3;

  % fan outlet conditions
  P_o8 = P_o2 * pr_fan;
  T_o8 = T_o2 * (1 + ((1/eta_fan)*(power(pr_fan,(g_fan-1)/g_fan) - 1)));
  u_ef = sqrt(2*eta_fan_nozzle*(g_fan/(g_fan - 1))*R*T_o8*(1 - power(Pa/P_o8,(g_fan - 1)/g_fan)));

  % Turbine outlet conditions
  W_fan_in = [];
  W_turb_out = [];
  T_o5 = [];
  P_o5 = [];
  u_e_core = [];
  for ii = 1:numel(bypass_ratios)
    fan_work_value = bypass_ratios(ii)*cp_fan*(T_o8 - T_o2);
    W_fan_in = [W_fan_in fan_work_value];
    
    turb_work_value = fan_work_value + W_comp_in;
    W_turb_out = [W_turb_out turb_work_value];
    
    temp_value = T_o4 - (turb_work_value/((1+f)*cp_post_burner));
    T_o5 = [T_o5 temp_value];
    
    press_value = P_o4 * (1 - ((1/eta_turbine)*(1 - (temp_value/T_o4))))^(g_post_burner/(g_post_burner-1));
    P_o5 = [P_o5 press_value];
    
    vel_value = sqrt(2*eta_core_nozzle*(g_post_burner/(g_post_burner-1))*R*temp_value*(1 - power(Pa/press_value,(g_post_burner-1)/g_post_burner)));
    u_e_core = [u_e_core vel_value];
  endfor

  for ii = 1:numel(u_e_core)
    temp_num = u_e_core(ii);
    if ~iscomplex(temp_num)
      u_e_core(ii) = temp_num;
    else
      u_e_core(ii) = 0;
    endif
  endfor

  u_e_core = u_e_core(u_e_core ~= 0);

  eta_th = [];
  eta_p = [];
  eta_o = [];
  I = [];
  TSFC = [];
  for ii = 1:numel(u_e_core)
    temp_eta_th = (((1+f)*(power(u_e_core(ii),2)/2)) + ...
             (bypass_ratios(ii)*(power(u_ef,2)/2)) - ...
             ((bypass_ratios(ii) + 1)*(power(Ua,2)/2))) / (f*h_c);
    eta_th = [eta_th temp_eta_th];
    
    temp_eta_p = (((bypass_ratios(ii)*(u_ef - Ua)) + (((1+f)*u_e_core(ii)) - Ua)) * Ua) ...
            / (((1+f)*(power(u_e_core(ii),2)/2)) + (bypass_ratios(ii)*(power(u_ef,2)/2)) - ...
             ((bypass_ratios(ii) + 1)*(power(Ua,2)/2)));
    eta_p = [eta_p temp_eta_p];
    
    temp_eta_o = temp_eta_th * temp_eta_p; 
    eta_o = [eta_o temp_eta_o];
    
    temp_I = (((1+f)*u_e_core(ii)) - Ua) + ((bypass_ratios(ii)+1)*(power(Ua,2)/2));
    I = [I temp_I];
    
    temp_TSFC = f/temp_I;
    TSFC = [TSFC temp_TSFC];
  endfor

  if(bpr_range)
    figure;
    plot(bypass_ratios(1:numel(eta_th)),eta_th,';eta\_th;');
    hold on;
    plot(bypass_ratios(1:numel(eta_p)),eta_p,';eta\_p;');
    plot(bypass_ratios(1:numel(eta_o)),eta_o,';eta\_o;');
    xlabel('Bypass Ratio'); ylabel('Efficiency');
    legend show;
    title('Efficiencies vs. Bypass Ratio');
    hold off;

    figure;
    plot(bypass_ratios(1:numel(I)),I);
    xlabel('Bypass Ratio'); ylabel('I');
    title('Specific Thrust vs. Bypass Ratio');
    figure;
    plot(bypass_ratios(1:numel(TSFC)),TSFC);
    xlabel('Bypass Ratio'); ylabel('TSFC');
    title('TSFC vs. Bypass Ratio');

    max_eta = 0;
    for ii = 1:numel(eta_o)
      if max_eta < eta_o(ii)
        max_eta = eta_o(ii);
        max_eta_index = ii;
      else
        continue;
      endif
    end  
    best_bypass_ratio = bypass_ratios(max_eta_index);  
    printf("Best bypass ratio is %f\n",best_bypass_ratio);
  endif  
endfunction  