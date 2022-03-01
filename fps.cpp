#include "fps.h"

fps::fps() {
  num = 0;
  exponent = 0;
}

fps::fps(bnum n, bnum ps) {
  exponent = ps;
  num = n % exponent;
  if(num < 0) num += exponent;
}

fps fps::operator+(fps b) {
  fps tmp;
  if(exponent != b.exponent) return tmp;
  tmp.exponent = exponent;
  tmp.num = (num + b.num) % exponent;
  if(tmp.num < 0) tmp.num += exponent;
  return tmp;
}

fps fps::operator-(fps b) {
  fps tmp;
  if(exponent != b.exponent) return tmp;
  tmp.exponent = exponent;
  tmp.num = (num - b.num) % exponent;
  if(tmp.num < 0) tmp.num += exponent;
  return tmp;
}

fps fps::operator*(fps b) {
  fps tmp;
  if(exponent != b.exponent) return tmp;
  tmp.exponent = exponent;
  tmp.num = num * b.num % exponent;
  if(tmp.num < 0) tmp.num += exponent;
  return tmp;
}

void fps::operator=(bnum b) {
  num = b % exponent;
  if(num < 0) num += exponent;
}

int main() {
  return 0;
}