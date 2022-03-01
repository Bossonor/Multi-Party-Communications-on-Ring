#pragma once
#include "globals.h"
#include "fps.h"
using namespace std;

class poly {
  public:
    bnum len; // highest degree
    vector<fps> data;

    poly();
    poly(bnum l, bnum ps);
    poly operator+(poly b);
    poly operator-(poly b);
    //poly operator*(poly b);
    //poly operator/(poly b);

    friend ostream & operator<<(ostream & os, const poly & a) {
      os << a.data;
      return os;
    }
}