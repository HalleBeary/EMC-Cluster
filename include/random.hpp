#ifndef random_HPP
#define random_HPP

//Random number generator


/********************************************************************/
/* Long period (>2×1018) random number generator of L'Ecuyer with   */
/* Bays-Durham shuffle and added safeguards. Returns uniform random */ 
/* deviate between 0.0 and 1.0 (exclusive of the endpoint values).  */
/* Call with iso a negative int to initialize; thereafter, do not   */
/* alter iso between successive deviates in a sequence. RNMX should */ 
/* approximate the largest floating value that is less than 1.      */
/********************************************************************/

double Randomizer() // random generator
{
/* Initialized data */
  static int iso  = -212312;//-1345;
  static int iso2 = 54524;//123456789;
  static int ggiv[32] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static int iy = 0;
  static int j, k;
  
  /* System generated locals */
  int tmp;
  double ret_val, r2;

  tmp = iso;
  if (tmp <= 0)
  {
    tmp = (-tmp > 1 ? -tmp : 1); /* Computing MAX: iso = max(-iso,1) */
    iso2 = tmp;     
    for (j = 40; j >= 1; --j)
    {
      k    = tmp / 53668;
      tmp = (tmp - k * 53668) * 40014 - k * 12211;

      if (tmp < 0)
	tmp += 2147483563;
      if (j <= 32)
	ggiv[j - 1] = tmp;
    }
    iy = ggiv[0];
  }
  k = tmp / 53668;
  tmp = (tmp - k * 53668) * 40014 - k * 12211;
  if (tmp < 0)
    tmp += 2147483563;
  
  k = iso2 / 52774;
  iso2 = (iso2 - k * 52774) * 40692 - k * 3791;
  if (iso2 < 0)
    iso2 += 2147483399;
  
  j = iy / 67108862 + 1;
  iy = ggiv[j - 1] - iso2;
  ggiv[j - 1] = tmp;

  iso = tmp;
  if (iy < 1)
    iy += 2147483562;
  
  r2 = iy * 4.6566130573917691e-10; /* Computing MIN: rand = min(AM*iy,RNMX) */ 
  ret_val = (r2 < .99999987999999995 ? r2 : .99999987999999995);

  return ret_val;
}


#endif