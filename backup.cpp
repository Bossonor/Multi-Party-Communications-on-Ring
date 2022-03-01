#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <vector>
#include <algorithm>

using namespace std;
using namespace NTL;

const ZZ prime = conv<ZZ>(3); // p 
const long positive_integer = RandomBits_long(3) + 1; // 1 <= s <= 32bits
const ZZ exponent_ps = power(prime, positive_integer); // p^s

/* find all prime factor of p^r-1 */
vector<long> find_all_prime_factor(ZZ number) {
    vector<long> res;
    PrimeSeq all_prime; // <=30 bits
    long p = all_prime.next();
    while(p <= number) {
        if(number % p == 0) {
            res.push_back(p);
            number /= p;
        }
        p = all_prime.next();
    }
    return res;
}

/* test if poly is primitive */
bool test_primitive(ZZ_pX poly, vector<long> p_array, ZZ degree) {
    cout << endl << "primitive test: " << endl;
    cout << "poly: " << poly << endl;
    ZZ_pX tmp_poly;
    for(auto num : p_array) {
        cout << "factor: " << degree / num << endl;
        tmp_poly.SetLength(conv<long>(degree / num + 1));
        SetCoeff(tmp_poly, tmp_poly.rep.length() - 1, 1);
        SetCoeff(tmp_poly, 0, -1);
        cout << "tmp_poly % poly: " << tmp_poly % poly << endl;
        if((tmp_poly % poly).rep.length() == 0) {
            cout << "test false" << endl;
            return false;
        } 
    }
    cout << "test true" << endl;
    return true; 
}

/* initiate the data related with p^r-1 */
void initiate_pr_relate(long &degree, ZZ &pr, vector<long> &pr_factor) {
    long pr_bit;
    do{
        degree = 5;//RandomBits_long(3) + 2; // 1 <= r <= 32bits; r = degree - 1
        pr = power(prime, degree - 1) - 1; // p^r - 1
        pr_bit = NumBits(pr);
    }while(pr_bit >= sizeof(long) * 8);
    cout << "degree_r: "<< degree - 1 << endl;
    cout << "exponent_pr-1: "<< pr << endl;
    cout << "all_prime_factor of p^r-1 : ";
    pr_factor = find_all_prime_factor(pr);
    for(auto num : pr_factor) cout << num << " ";
    cout << endl;
}

int main() {
/* initiate all data */
    ZZ_p::init(prime); 
    cout << "prime_p: " << prime << endl;
    cout << "positive_integer_s: " << positive_integer << endl;
    cout << "exponent_ps: "<< exponent_ps << endl;

    long degree_r;
    ZZ exponent_pr;
    vector<long> pr_prime_factor;
    initiate_pr_relate(degree_r, exponent_pr, pr_prime_factor);

/* get fps_poly which is a degree of r, monic, primitive and irreducible poly */
    ZZ_pX Fps_poly;   
    while(1) {
        ZZ_pX pr_poly;
        pr_poly.SetLength(conv<long>(exponent_pr + 1));
        SetCoeff(pr_poly, pr_poly.rep.length() - 1, 1);
        SetCoeff(pr_poly, 0, -1);
        Vec<Pair<ZZ_pX, long>> poly_factor;
        CanZass(poly_factor, pr_poly);
        cout << "all factors of poly:" << endl << poly_factor << endl;
        for(int i = 0; i < poly_factor.length(); i++) {
            ZZ_pX tmp_poly = poly_factor[i].a;
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
    cout << endl << "Fps_poly: " << Fps_poly << endl;

/* reset */
    ZZ_p::init(exponent_ps);
    ZZ_pE::init(Fps_poly);

/* get the set T of GR */
    cout << endl << "set T: ";
    vector<ZZ_pE> set_T;
    ZZ_pX tmp_poly;
    tmp_poly.SetLength(1);
    set_T.push_back(conv<ZZ_pE>(tmp_poly));
    cout << set_T[0] << " ";
    tmp_poly.SetLength(2);
    SetCoeff(tmp_poly, 1, 1);
    ZZ_pE tmp_ring(1);
    for(int i = 0; i < exponent_pr; i++) {
        cout << tmp_ring << " ";
        set_T.push_back(tmp_ring);
        set_T.erase(unique(set_T.begin(), set_T.end()), set_T.end());
        tmp_ring *= conv<ZZ_pE>(tmp_poly);
    }
    cout << endl;

/* get n + 1 distinct elements from T */
    long count_n =  set_T.size() - 1;//RandomBnd(set_T.size());
    cout << endl << "n: " << count_n << endl;
    vector<ZZ_pE> tmp_set = set_T;
    random_shuffle(tmp_set.begin(), tmp_set.end());
    vector<ZZ_pE> elem_in_T;
    elem_in_T.assign(tmp_set.begin(), tmp_set.begin() + count_n + 1);
    cout << "n + 1 elements in T: ";
    for(auto num : elem_in_T) cout << num << " ";
    cout << endl;

/* get n + 1 units from GR, here from any two element's difference in T*/
    random_shuffle(tmp_set.begin(), tmp_set.end());
    vector<ZZ_pE> unit_in_GR;
    cout << endl << "n + 1 units in GR: ";
    for(int i = 0; i < count_n; i++) {
        unit_in_GR.push_back(tmp_set[i + 1] - tmp_set[i]);
        cout << unit_in_GR[i] << " ";
    }
    unit_in_GR.push_back(tmp_set[tmp_set.size() - 1] - tmp_set[0]);
    cout << unit_in_GR[count_n] << " " << endl;

/* generate random GR[x] && degree = k - 1 */
    long degree_k = count_n - 1;//RandomBnd(count_n);
    cout << endl << "degree k: " << degree_k << endl;
    ZZ_pEX GR_X = random_ZZ_pEX(degree_k);
    cout << "random GR[x]: " << GR_X << endl;

/* compute the RS_code */
    vector<ZZ_pE> RS_code;
    cout << endl << "RS_code: ";
    for(int i = 0; i <= count_n; i++) {
        RS_code.push_back(unit_in_GR[i] * eval(GR_X, elem_in_T[i]));
        cout << RS_code[i] << " ";
    }
    cout << endl;
}