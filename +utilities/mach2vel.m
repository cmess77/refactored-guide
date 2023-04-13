function vel = mach2vel(mach, gamma, R, T)
  
  vel = mach*sqrt(gamma*R*T);
endfunction