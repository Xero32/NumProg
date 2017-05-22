

  #ifndef LEAPFROG1D_H
  #define LEPAFROG1D_H
  
  #include "basic.h"
  #include "gridfunc1d.h"
  
  /*one step with Leapfrog*/
  void
  step_leapfrog1d_wave(pcgridfunc1d u_old, pcgridfunc1d v_old,
  pgridfunc1d u_new, pgridfunc1d v_new, double t, double delta, void* data);
  

  #endif