#pragma once
#include "globals.h"

class fps {
  public:
    bnum num;
    bnum exponent;

    fps();
    fps(bnum n, bnum ps);
    fps operator+(fps b);
    fps operator-(fps b);
    fps operator*(fps b);
    void operator=(bnum b);
    
    friend ostream & operator<<(ostream & os, const fps & a) {
      os << ll(a.num);
      return os;
    }
};
