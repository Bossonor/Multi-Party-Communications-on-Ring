#include <NTL/ZZ_pXFactoring.h>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;
using namespace NTL;

const ZZ prime = RandomPrime_ZZ(4); //p
const long positive_integer = 1;//RandomBits_long(3) + 1; // 1 <= s <= 32bits
const ZZ exponent_ps = power(prime, positive_integer); // p^s

/* test if poly is irreducible */
bool test_irred(ZZ_pX poly) {
    ZZ_pContext context;
    context.save();
    ZZ_p::init(prime);
    bool flag = DetIrredTest(poly) == 1 ? true : false;
    context.restore();
    cout << "flag:" << flag << endl;
    return flag;
}

/* find all prime factor of p^r-1 */
vector<long> find_all_prime_factor(ZZ number) {
    vector<long> res;
    PrimeSeq all_prime; // <=30 bits
    long p = all_prime.next();
    while(p <= SqrRoot(number)) {
        if(number % p == 0) res.push_back(p);
        p = all_prime.next();
    }
    return res;
}

/* test if poly is primitive */
bool test_primitive(ZZ_pX poly, vector<long> p_array, ZZ degree) {
    cout << "poly: " << poly << endl;
    ZZ_pX tmp_poly;
    for(auto num : p_array) {
        cout << "factor:" << degree / num << endl;
        tmp_poly.SetLength(conv<long>(degree / num + 1));
        SetCoeff(tmp_poly, tmp_poly.rep.length() - 1, 1);
        SetCoeff(tmp_poly, 0, -1);
        cout << "tmp_poly % poly: " << tmp_poly % poly << endl;
        if((tmp_poly % poly).rep.length() == 0) return false; 
    }
    return true; 
}

/* initiate the data related with p^r-1 */
void initiate_pr_relate(long &degree, ZZ &pr, vector<long> &pr_factor) {
    long pr_bit;
    do{
        degree = 2;//RandomBits_long(3) + 2; // 1 <= r <= 32bits; r = degree - 1
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

/* get an example of GR */ 
ZZ_p generate_random_GR(ZZ_pX poly, ZZ_p root, long degree) {
    ZZ_pX tmp_poly = random_ZZ_pX(RandomBits_long(degree * 2) + 1);
    tmp_poly %= poly;
    return eval(tmp_poly, root);
}

int main() {
/* initiate all data */
    ZZ_p::init(exponent_ps); 
    cout << "prime_p: " << prime << endl;
    cout << "positive_integer_s: " << positive_integer << endl;
    cout << "exponent_ps: "<< exponent_ps << endl;
    
    long degree_r;
    ZZ exponent_pr;
    vector<long> pr_prime_factor;
    initiate_pr_relate(degree_r, exponent_pr, pr_prime_factor);   

/* get fps_poly which is a monic primitive poly, here the poly is irreducible*/
    ZZ_pX Fps_poly;   
    while(1) {
        ZZ_pX pr_poly;
        pr_poly.SetLength(conv<long>(exponent_pr + 1));
        SetCoeff(pr_poly, pr_poly.rep.length() - 1, 1);
        SetCoeff(pr_poly, 0, -1);
        Vec<Pair<ZZ_pX, long>> poly_factor;
        CanZass(poly_factor, pr_poly);
        cout << "all factors of poly" << poly_factor << endl;
        for(int i = 0; i < poly_factor.length(); i++) {
            ZZ_pX tmp_poly = poly_factor[i].a;
            if(LeadCoeff(tmp_poly) == 1 && tmp_poly.rep.length() == degree_r) {
                if(test_irred(tmp_poly) && test_primitive(tmp_poly, pr_prime_factor, exponent_pr)) {
                    Fps_poly = tmp_poly;
                    break;
                }
            }
        }
        if(Fps_poly.rep.length() != 0) break;
        initiate_pr_relate(degree_r, exponent_pr, pr_prime_factor);           
    }   
    cout << "Fps_poly:" << Fps_poly << endl;

/* get root from fps_poly */
    ZZ_p Fps_root = FindRoot(Fps_poly);
    cout << "Fps_poly's root:" << Fps_root << endl;

/* get GR's nonzero element of order p^r-1 in GR*/
    ZZ_p ord_pr_elem = generate_random_GR(Fps_poly, Fps_root, degree_r);
    cout << "GR's nonzero element of order p^r-1 :" << ord_pr_elem << endl;

/* get the set T of GR */
    std::set<ZZ> tmp_set = {ZZ(0), ZZ(1)};
    ZZ_p tmp_elem(1);
    for(int i = 1; i < exponent_pr; i++) {
        tmp_elem *= ord_pr_elem;
        tmp_set.insert(conv<ZZ>(tmp_elem));
    }
    vector<ZZ_p> set_T;
    cout << "set T : ";
    for(auto num : tmp_set) {
        cout << num << " ";
        set_T.push_back(conv<ZZ_p>(num));
    }
    cout << endl;

/* get n + 1 distinct elements from T */
    long count_n =  set_T.size() - 1;//RandomBnd(set_T.size());
    cout << "n : " << count_n << endl;
    /*tmp_set.clear();
    while(tmp_set.size() < count_n + 1) {
        srand((unsigned)time(0));
        long x = rand() % set_T.size();
        tmp_set.insert(conv<ZZ>(set_T[x]));
    }*/
    vector<ZZ_p> elem_in_T = set_T;//;
    cout << "n + 1 elements in T: ";
    //for(auto num : tmp_set) elem_in_T.push_back(conv<ZZ_p>(num));
    random_shuffle(elem_in_T.begin(), elem_in_T.end());
    for(auto num : elem_in_T) cout << num << " ";
    cout << endl;

/* get n + 1 units from GR, here from any two element's difference in T*/
    tmp_set.clear();
    while(tmp_set.size() < count_n + 1) {
        srand((unsigned)time(0));
        long x = rand() % set_T.size();
        srand((unsigned)time(0));
        long y = rand() % set_T.size();
        tmp_set.insert(conv<ZZ>(set_T[x] - set_T[y]));
    }
    vector<ZZ_p> unit_in_GR;
    cout << "n + 1 units in GR: ";
    for(auto num : tmp_set) unit_in_GR.push_back(conv<ZZ_p>(num));
    random_shuffle(unit_in_GR.begin(), unit_in_GR.end());
    for(auto num : unit_in_GR) cout << num << " ";
    cout << endl;

/* generate random GR[x] && degree = k - 1 */
    long degree_k = count_n - 1;//RandomBnd(count_n);
    cout << "degree k: " << degree_k << endl;
    ZZ_pX GR_X;
    GR_X.SetLength(degree_k);
    for(int i = 0; i < degree_k; i++) SetCoeff(GR_X, i, generate_random_GR(Fps_poly, Fps_root, degree_r));
    cout << "random GR[x]: " << GR_X << endl;

/* compute the RS_code */
    vector<ZZ_p> RS_code;
    cout << "RS_code: ";
    for(int i = 0; i <= count_n; i++) {
        RS_code.push_back(unit_in_GR[i] * eval(GR_X, elem_in_T[i]));
        cout << RS_code[i] << " ";
    }
    cout << endl;
}