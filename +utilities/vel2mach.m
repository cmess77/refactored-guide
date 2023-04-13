function mach = vel2mach(vel,gamma,R,T)
  
  mach = vel / sqrt(gamma*R*T);
  
endfunction  