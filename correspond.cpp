#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
//#include <vector>
#include <NTL/vec_long.h>
#include <NTL/pair_ZZ_pX_long.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
//#include <algorithm>

using namespace std;
using namespace NTL;

const ZZ prime = conv<ZZ>(2); // p 
const long positive_integer = 3;//RandomBits_long(3) + 1; // 1 <= s <= 32bits
const ZZ exponent_ps = power(prime, positive_integer); // p^s

/* find all prime factor of p^r-1 */
vec_long find_all_prime_factor(ZZ number) {
    vec_long res;
    PrimeSeq all_prime; // <=30 bits
    long p = all_prime.next();
    while(p <= number) {
        if(number % p == 0) {
            res.append(p);
            number /= p;
        }
        p = all_prime.next();
    }
    return res;
}

/* test if poly is primitive */
bool test_primitive(ZZ_pX poly, vec_long p_array, ZZ degree) {
    cout << endl;
    cout << "test_poly: " << poly << endl;
    cout << "primitive test: " << endl;
    ZZ_pX tmp_poly;
    for(auto num : p_array) {
        cout << "  factor: " << degree / num << endl;
        tmp_poly.SetLength(conv<long>(degree / num + 1));
        SetCoeff(tmp_poly, tmp_poly.rep.length() - 1, 1);
        SetCoeff(tmp_poly, 0, -1);
        cout << "  tmp_poly % test_poly: " << tmp_poly % poly << endl;
        if((tmp_poly % poly).rep.length() == 0) {
            cout << "  test false" << endl;
            return false;
        } 
    }
    cout << "  test true" << endl;
    return true; 
}

/* initiate the data related with p^r-1 */
void initiate_pr_relate(long &degree, ZZ &pr, vec_long &pr_factor) {
    long pr_bit;
    do{
        degree = 5;//RandomBits_long(3) + 2; // 1 <= r <= 32bits; r = degree - 1
        pr = power(prime, degree - 1) - 1; // p^r - 1
        pr_bit = NumBits(pr);
    }while(pr_bit >= sizeof(long) * 8);
    cout << endl;
    cout << "degree_r: "<< degree - 1 << endl;
    cout << "exponent_pr-1: "<< pr << endl;
    cout << "all_prime_factor of p^r-1 : ";
    pr_factor = find_all_prime_factor(pr);
    cout << pr_factor << endl;
}

/* get inverse element from ring on Galois Ring by random*/
ZZ_pE inverse_elem(ZZ_pE ring, long degree) {
    //try {
        //return inv(ring);
    //}
    //catch(...) {
        //ZZ_p::init(prime);
        ZZ_pX tmp_poly;
        //ZZ_pX tmp_poly = conv<ZZ_pX>(inv(ring));
        //if(tmp_poly.rep.length() < degree + 1) 
        tmp_poly.SetLength(degree + 1);       
        //ZZ_p::init(exponent_ps);
        while(conv<ZZ_pE>(tmp_poly) * ring != conv<ZZ_pE>(1))
            for(int i = 0; i <= degree; i++)
                //SetCoeff(tmp_poly, i, tmp_poly[i] + conv<ZZ_p>(RandomBits_long(positive_integer - 1) * conv<long>(prime)));
                SetCoeff(tmp_poly, i, tmp_poly[i] + conv<ZZ_p>(RandomBits_long(positive_integer)));
        return conv<ZZ_pE>(tmp_poly);
    //}
}

/* get inverse element from ring on Galois Ring by euclidean if positive_integer = 1*/
ZZ_pX euclidean(ZZ_pX poly, ZZ_pX mod, ZZ_pX &s, ZZ_pX &t) {
    if(mod == conv<ZZ_pX>(0)) {
        s = conv<ZZ_pX>(1);
        t = conv<ZZ_pX>(0);
        return poly;
    }
    ZZ_pX gcd = euclidean(mod, poly % mod, s, t);
    ZZ_pX tmp_poly = s;
    s = t;
    t = tmp_poly - poly / mod * t;
    return gcd;
}
ZZ_pE extended_euclidean(ZZ_pX poly, ZZ_pX mod) {
    ZZ_pX s, t;
    ZZ_pX gcd = euclidean(poly, mod, s, t);
    if(gcd != conv<ZZ_pX>(1)) return conv<ZZ_pE>(0);
    return conv<ZZ_pE>(s % mod);
}

int main() {
/* initiate all data */
    ZZ_p::init(prime); 
    cout << "prime_p: " << prime << endl;
    cout << "positive_integer_s: " << positive_integer << endl;
    cout << "exponent_ps: "<< exponent_ps << endl;

    long degree_r;
    ZZ exponent_pr;
    vec_long pr_prime_factor;
    initiate_pr_relate(degree_r, exponent_pr, pr_prime_factor);

/* get fps_poly which is a degree of r, monic, primitive and irreducible poly */
    ZZ_pX Fps_poly;
    while(1) {
        ZZ_pX pr_poly;
        pr_poly.SetLength(conv<long>(exponent_pr + 1));
        SetCoeff(pr_poly, pr_poly.rep.length() - 1, 1);
        SetCoeff(pr_poly, 0, -1);
        vec_pair_ZZ_pX_long poly_factor;
        CanZass(poly_factor, pr_poly);
        cout << "all factors of poly which deg is p^r-1:" << endl;
        cout << poly_factor << endl;
        for(auto pr : poly_factor) {
            ZZ_pX tmp_poly = pr.a;
            if(LeadCoeff(tmp_poly) == 1 && tmp_poly.rep.length() == degree_r && DetIrredTest(tmp_poly) == 1) {
                if(test_primitive(tmp_poly, pr_prime_factor, exponent_pr)) {
                    Fps_poly = tmp_poly;
                    break;
                }
            }
        }
        if(Fps_poly.rep.length() != 0) break;
        initiate_pr_relate(degree_r, exponent_pr, pr_prime_factor);           
    }   
    cout << endl;
    cout <<"Fps_poly: " << Fps_poly << endl;

/* set Fps_poly and exponent_ps as environment*/
    ZZ_pE::init(Fps_poly);
    ZZ_p::init(exponent_ps);

/* get the set T of GR */
    cout << endl;
    vec_ZZ_pE set_T;
    ZZ_pX tmp_poly;
    tmp_poly.SetLength(1);
    set_T.append(conv<ZZ_pE>(tmp_poly));
    tmp_poly.SetLength(2);
    SetCoeff(tmp_poly, 1, 1);
    ZZ_pE tmp_ring(1);
    for(int i = 0; i < exponent_pr; i++) {
        set_T.append(tmp_ring);
        tmp_ring *= conv<ZZ_pE>(tmp_poly);
    }
    cout << "set T: " << endl;
    cout << set_T << endl;
    //ZZ_p::init(exponent_ps);
    //for(int i = 2; i < set_T.length(); i++) cout << inverse_elem(set_T[i], degree_r) << endl;

/* get n + 1 distinct elements from T */
    cout << endl;
    long count_n =  set_T.length() - 1;//RandomBnd(set_T.length());
    cout << "n: " << count_n << endl;
    vec_ZZ_pE elem_in_T = VectorCopy(set_T, count_n + 1);
    cout << "n + 1 elements in T: " << endl;
    cout << elem_in_T << endl;

/* get n + 1 units from GR, here are all "1" */
    //vector<ZZ_pE> tmp_std_set;
    //for(auto num : set_T) tmp_std_set.push_back(num);
    //random_shuffle(tmp_std_set.begin(), tmp_std_set.end());
    vec_ZZ_pE unit_in_GR;
    unit_in_GR.SetLength(count_n + 1, conv<ZZ_pE>(1));
    //for(int i = 0; i < count_n; i++) unit_in_GR.append(tmp_std_set[i + 1] - tmp_std_set[i]);
    //unit_in_GR.append(tmp_std_set[tmp_std_set.size() - 1] - tmp_std_set[0]);
    cout << "n + 1 units in GR: " << endl; 
    cout << unit_in_GR << endl;

/* set exponent_ps as environment */
    //ZZ_p::init(exponent_ps);

/* generate degree = t */
    cout << endl;
    long degree_t = count_n / 2;//RandomBnd(count_n);
    cout << "degree t: " << degree_t << endl;

/* generate degree = 2t */
    long degree_2t = count_n - 1;//RandomBnd(count_n);
    cout << "degree 2t: " << degree_2t << endl;

/* generate random GR[x] and degree = (t or 2t) && compute RS_code and RS_code[0] = c && compute [c]t and [c]2t */
    cout << endl;
    ZZ_pEX GR_X_t_c = random_ZZ_pEX(degree_t);
    cout << "c: " << conv<ZZ_pE>(unit_in_GR[0] * eval(GR_X_t_c, elem_in_T[0])) << endl;
    cout << endl;
    cout << "random GR[x]_t_c: " << endl; 
    cout << GR_X_t_c << endl;
    vec_ZZ_pE RS_code_t_c;
    for(int i = 0; i <= count_n; i++) RS_code_t_c.append(unit_in_GR[i] * eval(GR_X_t_c, elem_in_T[i]));
    cout << "RS_code_t_c: " << endl; 
    cout << RS_code_t_c << endl;
    vec_ZZ_pE c_t;
    for(int i = 1; i <= count_n; i++) c_t.append(unit_in_GR[i] * eval(GR_X_t_c, elem_in_T[i]));
    cout << "c_t: " << endl;
    cout << c_t << endl;

    cout << endl;
    ZZ_pEX GR_X_2t_c = random_ZZ_pEX(degree_2t);
    SetCoeff(GR_X_2t_c, 0, coeff(GR_X_t_c, 0));
    cout << "random GR[x]_2t_c: " << endl;
    cout << GR_X_2t_c << endl;
    vec_ZZ_pE RS_code_2t_c;
    for(int i = 0; i <= count_n; i++) RS_code_2t_c.append(unit_in_GR[i] * eval(GR_X_2t_c, elem_in_T[i]));
    cout << "RS_code_2t_c: " << endl; 
    cout << RS_code_2t_c << endl;
    vec_ZZ_pE c_2t;
    for(int i = 1; i <= count_n; i++) c_2t.append(unit_in_GR[i] * eval(GR_X_2t_c, elem_in_T[i]));
    cout << "c_2t: " << endl; 
    cout << c_2t << endl;

/* generate random GR[x] and degree = t && compute RS_code and RS_code[0] = a && compute [a]t */
    cout << endl;
    ZZ_pEX GR_X_t_a = random_ZZ_pEX(degree_t);
    cout << "a: " << conv<ZZ_pE>(unit_in_GR[0] * eval(GR_X_t_a, elem_in_T[0])) << endl;
    cout << "random GR[x]_t_a: " << endl;
    cout << GR_X_t_a << endl;
    vec_ZZ_pE RS_code_t_a;
    for(int i = 0; i <= count_n; i++) RS_code_t_a.append(unit_in_GR[i] * eval(GR_X_t_a, elem_in_T[i]));
    cout << "RS_code_t_a: " << endl;
    cout << RS_code_t_a << endl;
    vec_ZZ_pE a_t;
    for(int i = 1; i <= count_n; i++) a_t.append(unit_in_GR[i] * eval(GR_X_t_a, elem_in_T[i]));
    cout << "a_t: " << endl;
    cout << a_t << endl;

/* generate random GR[x] and degree = t && compute RS_code and RS_code[0] = b && compute [b]t*/
    cout << endl;
    ZZ_pEX GR_X_t_b = random_ZZ_pEX(degree_t);
    cout << "b: " << conv<ZZ_pE>(unit_in_GR[0] * eval(GR_X_t_b, elem_in_T[0])) << endl;
    cout << "random GR[x]_t_b: " << endl;
    cout << GR_X_t_b << endl;
    vec_ZZ_pE RS_code_t_b;
    for(int i = 0; i <= count_n; i++) RS_code_t_b.append(unit_in_GR[i] * eval(GR_X_t_b, elem_in_T[i]));
    cout << "RS_code_t_b: " << endl;
    cout << RS_code_t_b << endl;
    vec_ZZ_pE b_t;
    for(int i = 1; i <= count_n; i++) b_t.append(unit_in_GR[i] * eval(GR_X_t_b, elem_in_T[i]));
    cout << "b_t: " << endl;
    cout << b_t << endl;

/* compute (ab, [ab]2t) = (a, [a]t) * (b, [b]t) */
    cout << endl;
    vec_ZZ_pE RS_code_2t_ab;
    for(int i = 0; i <= count_n; i++) RS_code_2t_ab.append(RS_code_t_a[i] * RS_code_t_b[i]);
    cout << "RS_code_2t_ab: " << endl;
    cout << RS_code_2t_ab << endl;

/* compute [ab]2t = [a]t * [b]t */
    cout << endl;
    vec_ZZ_pE ab_2t;
    for(int i = 0; i < count_n; i++) ab_2t.append(a_t[i] * b_t[i]);
    cout << "ab_2t: " << endl;
    cout << ab_2t << endl;

/* compute (e, [e]2t) = (ab, [ab]2t) - (c, [c]2t) */
    cout << endl;
    vec_ZZ_pE RS_code_2t_e;
    for(int i = 0; i <= count_n; i++) RS_code_2t_e.append(RS_code_2t_ab[i] - RS_code_2t_c[i]);
    cout << "RS_code_2t_e: " << endl;
    cout << RS_code_2t_e << endl;

/* compute [e]2t = [ab]2t - [c]2t */
    cout << endl;
    vec_ZZ_pE e_2t;
    for(int i = 0; i < count_n; i++) e_2t.append(ab_2t[i] - c_2t[i]);
    cout << "e_2t: " << endl;
    cout << e_2t << endl;

/* compute e from [e]2t by coeff * A = [e]2t */
/*  
    cout << endl;
    mat_ZZ_pE A;
    A.SetDims(count_n - 1, count_n - 1);
    for(int j = 1; j < count_n; j++) {
        ZZ_pE tmp_ring(unit_in_GR[j + 1] * elem_in_T[1]);
        for(int i = 1; i < count_n; i++) {
            A[i - 1].put(j - 1, tmp_ring);
            tmp_ring *= elem_in_T[j + 1];
        }
    }
    cout << "matrix A: " << endl;
    cout << A << endl;

    cout << endl;
    ZZ_pE det_of_A = determinant(A);
    vec_ZZ_pE e_GR_tmp_coeff, tmp_set;
    for(int i = 1; i < count_n; i++) tmp_set.append(e_2t[i]);
    solve(det_of_A, e_GR_tmp_coeff, A, tmp_set); //coeff * A = [e]2t

    ZZ_pEX GR_X_2t_e;
    GR_X_2t_e.SetLength(e_GR_tmp_coeff.length());
    for(int i = 0; i < e_GR_tmp_coeff.length(); i++) SetCoeff(GR_X_2t_e, i, e_GR_tmp_coeff[i]);
    cout << "GR_X_2t_e: " << endl;
    cout << GR_X_2t_e << endl;

    cout << endl;
    cout << "e: " << conv<ZZ_pE>(unit_in_GR[0] * eval(GR_X_2t_e, elem_in_T[0])) << endl;
*/
/* compute e from [e]2t by langrange interpolation */
    cout << endl;
    ZZ_pEX GR_X_2t_e;
    for(int i = 1; i < count_n; i++) {
        ZZ_pEX tmp_poly1;
        tmp_poly1.SetLength(1);
        SetCoeff(tmp_poly1, 0, conv<ZZ_pE>(1));
        for(int j = 2; j <= count_n; j++) {
            if(i + 1 == j) continue;
            ZZ_pEX tmp_poly2;
            tmp_poly2.SetLength(2);
            SetCoeff(tmp_poly2, 1, conv<ZZ_pE>(1));
            SetCoeff(tmp_poly2, 0, conv<ZZ_pE>(-1) * elem_in_T[j]);
            //ZZ_p::init(prime);
            //tmp_ring = elem_in_T[i + 1] - elem_in_T[j];
            //ZZ_p::init(exponent_ps);
            tmp_poly2 *= inverse_elem(elem_in_T[i + 1] - elem_in_T[j], degree_r);
            //tmp_poly2 *= extended_euclidean(conv<ZZ_pX>(elem_in_T[i + 1] - elem_in_T[j]), Fps_poly);
            tmp_poly1 *= tmp_poly2;
        }
        GR_X_2t_e += tmp_poly1 * e_2t[i];
    }
    cout << "GR_X_2t_e: " << endl;
    cout << GR_X_2t_e << endl;

    cout << endl;
    cout << "e: " << conv<ZZ_pE>(unit_in_GR[0] * eval(GR_X_2t_e, elem_in_T[0])) << endl;

/* generate random GR[x] and degree = t && compute RS_code and RS_code[0] = e && compute [e]t */
    cout << endl;
    ZZ_pEX GR_X_t_e = random_ZZ_pEX(degree_t);
    SetCoeff(GR_X_t_e, 0, GR_X_2t_e[0]);
    cout << "random GR[x]_t_e: " << endl;
    cout << GR_X_t_e << endl;
    vec_ZZ_pE RS_code_t_e;
    for(int i = 0; i <= count_n; i++) RS_code_t_e.append(unit_in_GR[i] * eval(GR_X_t_e, elem_in_T[i]));
    cout << "RS_code_t_e: " << endl;
    cout << RS_code_t_e << endl;
    vec_ZZ_pE e_t;
    for(int i = 1; i <= count_n; i++) e_t.append(unit_in_GR[i] * eval(GR_X_t_e, elem_in_T[i]));
    cout << "e_t: " << endl;
    cout << e_t << endl;

/* compute (ab, [ab]t) = (e, [e]t) + (c, [c]t) */
    cout << endl;
    vec_ZZ_pE RS_code_t_ab;
    for(int i = 0; i <= count_n; i++) RS_code_t_ab.append(RS_code_t_e[i] + RS_code_t_c[i]);
    cout << "RS_code_t_ab: " << endl;
    cout << RS_code_t_ab << endl;

/* compute [ab]t = [e]t + [c]t */
    cout << endl;
    vec_ZZ_pE ab_t;
    for(int i = 0; i < count_n; i++) ab_t.append(e_t[i] + c_t[i]);
    cout << "ab_t: " << endl;
    cout << ab_t << endl;
  
}