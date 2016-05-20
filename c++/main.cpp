#include "Histogram.h"
#include <iostream>
using std::cout;
using std::endl;
int main()
{
  std::string filepath     = "../../data/collected.out";
  std::string metafilepath = "../../data/meta.out";
  Histogram test(metafilepath, filepath);
  //test.Rapidity(100);
  //test.EtaMax();
  test.NPOM_count();
  //test.nfMax();
  //test.NFNB(80);
  //test.NPOMS_NPOMH(5,5,20);
  //test.NFNB(100);

  return 0;
}
