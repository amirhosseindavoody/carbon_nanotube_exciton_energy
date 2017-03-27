#include <stdio.h>
#include <valarray>
#include <math.h>

#include <Eigen/Dense>

#include "my_math.h"
using namespace my_math;

int gcd ( int a, int b )
{
  int c;
  while ( a != 0 )
  {
     c = a;
     a = b%a;
     b = c;
  }
  return b;
}