#include "ibex_Plugin.h"

int
main (int argc, char *argv[])
{
  ibex::Interval ans = answer_to_the_ultimate_question();
  std::cout << "Product of 6 times 9 is " << ans << std::endl;
  return 0;
}
