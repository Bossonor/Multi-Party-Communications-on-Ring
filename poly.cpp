#include "poly.h"
using namespace std;

poly::poly() {
    len = 0;
}

poly::poly(bnum l, bnum ps) {
    len = l;
    data.resize(l + 1);
    fps num(0, ps);
    for(bnum i = 0; i <= l; i++) {
        if(i == l) num = bnum(rand() % (ps - 1) + 1);// The coefficient of the highest order term of the polynomial is not 0
        else num = bnum(rand() % ps);
        data[i] = num;
    }
}

poly poly::operator+(poly b) {
    poly tmp;
    tmp.len = max(len, b.len);
    tmp.data.resize(tmp.len + 1);
    for(bnum i = 0; i <= tmp.len; i++) {
        if(i <= len && i <= b.len) tmp.data[i] = data[i] + b.data[i];
        if(i <= len && i > b.len) tmp.data[i] = data[i];
        if(i > len && i <= b.len) tmp.data[i] = b.data[i];
    }
    for(bnum i = tmp.len; i >= 0; i--) {
        if(tmp.data[i].num == 0) {
            tmp.len--;
            tmp.data.pop_back();
        }
        else break;
    }
    return tmp;
}

poly poly::operator-(poly b) {
    poly tmp;
    tmp.len = max(len, b.len);
    tmp.data.resize(tmp.len + 1);
    for(bnum i = 0; i <= tmp.len; i++) {
        if(i <= len && i <= b.len) tmp.data[i] = data[i] - b.data[i];
        if(i <= len && i > b.len) tmp.data[i] = data[i];
        if(i > len && i <= b.len) tmp.data[i] = b.data[i];
    }
    for(bnum i = tmp.len; i >= 0; i--) {
        if(tmp.data[i].num == 0) {
            tmp.len--;
            tmp.data.pop_back();
        }
        else break;
    }
    return tmp;
}

