#include "ibex.h"

int
main (int argc, char *argv[])
{
  ibex::Interval x(6,6), y(9-2,9+2);
  std::cout << "Product of 6 times 9 is " << x*y << std::endl;
  return 0;
}

