/*
* Creation date: 2014.04.16
* Creators: Youcef Lemsafer
* Authors: Youcef Lemsafer
*/
#include <string.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <chrono>



typedef long long int64;
typedef unsigned long long uint64;

#ifdef __cplusplus
extern "C" {
#endif

bool is_prime(unsigned long long p);
/*
* We assume that P = 1 and Q = (1 - D)/4
* where D is the first element in {5, -7, 9, -11, ...}
* such that the Jacobi symbol (D/n) = -1.
*/
bool is_slprp(uint64 n, int64 D);
int64 jacobi_symbol(int64 a, uint64 m);

#ifdef __cplusplus
}
#endif

void
print_primes_below(uint64 limit, std::ostream & out)
{
    out << std::endl;
    uintptr_t j = 1;
    for(uintptr_t n(0); n < limit; ++n)
    {
        if(is_prime(n)) {
            out.fill(' ');
            out.width(11);
            out << n;
            if( (j & 7) == 0 ) {
                out << std::endl << std::endl;
                j = 0;
            }
            ++j;
        }
    }
    out << std::endl;
}


template <typename V>
void
add_slprp_test(V & tests, uint64 n, int64 D, bool expected)
{
    tests.push_back(std::make_pair(std::make_pair(n, D), expected));
}

std::vector<std::pair<std::pair<uint64, int64>, bool> >
build_slprp_tests()
{
    decltype(build_slprp_tests()) testsVec;
    // Values below were taken from http://oeis.org/A217255
    // these are strong Lucas pseudo primes with parameters (P, Q)
    // selected using Selfridge's method A.
    // I used Wolfram alpha to find the second parameter:
    // e.g. JacobiSymbol[D, 324899] where D in {5, -7, 9, -11, 13, -15}
    // yields {1, 1, 1, 1, 1, -1} hence D = -15
    add_slprp_test(testsVec, 5459ULL,   -7LL,  true);
    add_slprp_test(testsVec, 5777ULL,   5LL,   true);
    add_slprp_test(testsVec, 10877ULL,  5LL,   true);
    add_slprp_test(testsVec, 16109ULL,  13LL,  true);
    add_slprp_test(testsVec, 18971ULL,  -11LL, true);
    add_slprp_test(testsVec, 22499ULL,  -15LL, true);
    add_slprp_test(testsVec, 24569ULL,  -7LL,  true);
    add_slprp_test(testsVec, 25199ULL,  -7LL,  true);
    add_slprp_test(testsVec, 40309ULL,  -7LL,  true);
    add_slprp_test(testsVec, 58519ULL,  -7LL,  true);
    add_slprp_test(testsVec, 75077ULL,  5LL,   true);
    add_slprp_test(testsVec, 97439ULL,  -7LL,  true);
    add_slprp_test(testsVec, 100127ULL, 5LL,   true);
    add_slprp_test(testsVec, 113573ULL, 5LL,   true);
    add_slprp_test(testsVec, 115639ULL, -7LL,  true);
    add_slprp_test(testsVec, 130139ULL, -15LL, true);
    add_slprp_test(testsVec, 155819ULL, -7LL,  true);
    add_slprp_test(testsVec, 158399ULL, -7LL,  true);
    add_slprp_test(testsVec, 161027ULL, 5LL,   true);
    add_slprp_test(testsVec, 162133ULL, 5LL,   true);
    add_slprp_test(testsVec, 176399ULL, -7LL,  true);
    add_slprp_test(testsVec, 176471ULL, -15LL, true);
    add_slprp_test(testsVec, 189419ULL, -7LL,  true);
    add_slprp_test(testsVec, 192509ULL, 13LL,  true);
    add_slprp_test(testsVec, 197801ULL, -11LL, true);
    add_slprp_test(testsVec, 224369ULL, -7LL,  true);
    add_slprp_test(testsVec, 230691ULL, -7LL,  true);
    add_slprp_test(testsVec, 231703ULL, 5LL,   true);
    add_slprp_test(testsVec, 243629ULL, -15LL, true);
    add_slprp_test(testsVec, 253259ULL, -7LL,  true);
    add_slprp_test(testsVec, 268349ULL, -15LL, true);
    add_slprp_test(testsVec, 288919ULL, 13LL,  true);
    add_slprp_test(testsVec, 313499ULL, -11LL, true);
    add_slprp_test(testsVec, 324899ULL, -15LL, true);
    // End of values taken from A217255

    // Begin Prime values of n
    add_slprp_test(testsVec, 61051ULL, -31LL, true);
    add_slprp_test(testsVec, 399499ULL, -43LL, true);
    add_slprp_test(testsVec, 4355311ULL, 61LL, true);
    add_slprp_test(testsVec, 5715319ULL, -67LL, true);
    add_slprp_test(testsVec, 4294967311ULL, -7LL, true);
    add_slprp_test(testsVec, 4294967357ULL, 5LL, true);
    add_slprp_test(testsVec, 4294967371ULL, -11LL, true);
    add_slprp_test(testsVec, 18536536771ULL, -103LL, true);
    add_slprp_test(testsVec, 263456789093ULL, 5LL, true);
    add_slprp_test(testsVec, 1099511627791ULL, -7LL, true);
    add_slprp_test(testsVec, 1099511627803ULL, 5LL, true);
    add_slprp_test(testsVec, 1099511627831ULL, 13LL, true);
    add_slprp_test(testsVec, 9007199254740997ULL, 5LL, true);
    add_slprp_test(testsVec, 9007199256950959ULL, -67LL, true);
    add_slprp_test(testsVec, 576460752303423619ULL, 13LL, true);
    add_slprp_test(testsVec, 576460752303423649ULL, -11LL, true);
    add_slprp_test(testsVec, 576460752303423733ULL, 5LL, true);
    add_slprp_test(testsVec, 576460752303423737ULL, 5LL, true);
    add_slprp_test(testsVec, 576460752303423749ULL, -7LL, true);
    add_slprp_test(testsVec, 576460752303423761ULL, 13LL, true);
    add_slprp_test(testsVec, 700000000000000289ULL, -11LL, true);
    add_slprp_test(testsVec, 2305843009213693921ULL, -7LL, true);
    add_slprp_test(testsVec, 2305843009213693951ULL, 17LL, true);
    add_slprp_test(testsVec, 2305843009213693967ULL, 5LL, true);
    add_slprp_test(testsVec, 2305843009220204479ULL, -71LL, true);
    add_slprp_test(testsVec, 9000000000000000041ULL, -11LL, true);
    add_slprp_test(testsVec, 9223372036854775507ULL, 5LL, true);
    add_slprp_test(testsVec, 9223372036866316891ULL, -83LL, true);
    add_slprp_test(testsVec, 9223372036918424191ULL, 97LL, true);
    add_slprp_test(testsVec, 9223372037009299789ULL, -107LL, true);
    add_slprp_test(testsVec, 9223372037410063339ULL, 97LL, true);
    add_slprp_test(testsVec, 10000000000000000051ULL, -7LL, true);
    add_slprp_test(testsVec, 12345678901234567891ULL, -11L, true);
    add_slprp_test(testsVec, 17000000000000000003ULL, 5LL, true);
    add_slprp_test(testsVec, 18446744064055665229ULL, 97LL, true);
    add_slprp_test(testsVec, 18446744064067846081ULL, -83LL, true);
    add_slprp_test(testsVec, uint64(-1) - 58, int64(5), true);
    // End Prime values of n

    // Begin Composite values of n
    // 2047 = 23*89 is a strong pseudo prime base 2
    add_slprp_test(testsVec, 2047ULL,       5LL,  false);
    add_slprp_test(testsVec, 4294967333ULL, 5LL,  false);
    add_slprp_test(testsVec, 4294967359ULL, 13LL, false);
    add_slprp_test(testsVec, 27300250277ULL, 5LL, true);   /* SLPSP 116833x233669 */
    add_slprp_test(testsVec, 27300250279ULL, -23LL, false);
    add_slprp_test(testsVec, 263456789039ULL, -7LL, false);
    add_slprp_test(testsVec, 263456789041ULL, -7LL, false);
    add_slprp_test(testsVec, 263456789047ULL, 5LL, false);
    add_slprp_test(testsVec, 263456789051ULL, -15LL, false);
    add_slprp_test(testsVec, 263456789081ULL, -7LL, false);
    add_slprp_test(testsVec, 263456789119ULL, -7LL, false);
    add_slprp_test(testsVec, 263456789129ULL, -11LL, false);
    add_slprp_test(testsVec, 1099511627801ULL, -7LL, false);
    add_slprp_test(testsVec, 1099511627813ULL, 5LL, false);
    add_slprp_test(testsVec, 39972590422099ULL, -7LL, true);    /* SLPSP 203419x196503721 */
    add_slprp_test(testsVec, 71015542332359ULL, -15LL, true);   /* SLPSP 5958839x11917681 */
    add_slprp_test(testsVec, 71026840741877ULL, 5LL, true);  /* SLPSP 5959313x11918629 */
    add_slprp_test(testsVec, 71027509003163ULL, 5LL, true);  /* SLPSP 3440627x20643769 */
    add_slprp_test(testsVec, 71028048949859ULL, -15LL, true);  /* SLPSP 8039x8835433381 */
    add_slprp_test(testsVec, 71028425928179ULL, -7LL, true);  /* SLPSP 5959381x11918759 */
    add_slprp_test(testsVec, 71028496858199ULL, -7LL, true);  /* SLPSP 2541089x27951991 */
    add_slprp_test(testsVec, 71029125005299ULL, -7LL, true);  /* SLPSP 13721x34301x150919 */
    add_slprp_test(testsVec, 71029563550243ULL, 5LL, true);  /* SLPSP 426763x166437961 */
    add_slprp_test(testsVec, 71030483737309ULL, -7LL, true);  /* SLPSP 7298827x9731767 */
    add_slprp_test(testsVec, 71035764857873ULL, 5LL, true);  /* SLPSP 2665253x26652541 */
    add_slprp_test(testsVec, 71036072324099ULL, -15LL, true);  /* SLPSP 8428289x8428291 */
    add_slprp_test(testsVec, 71040345630799ULL, -7LL, true);  /* SLPSP 3349163x21211373 */
    add_slprp_test(testsVec, 71040523518089ULL, -7LL, true);  /* SLPSP 878737x80843897 */
    add_slprp_test(testsVec, 71040611030759ULL, -7LL, true);  /* SLPSP 4214279x16857121 */
    add_slprp_test(testsVec, 71041841765927ULL, 5LL, true);  /* SLPSP 2665367x26653681 */
    add_slprp_test(testsVec, 71042052321701ULL, -7LL, true);  /* SLPSP 4214323x16857287 */
    add_slprp_test(testsVec, 71044690405037ULL, 5LL, true);  /* SLPSP 397337x178802101 */
    add_slprp_test(testsVec, 71045053572089ULL, -15LL, true);  /* SLPSP 5960077x11920157 */
    add_slprp_test(testsVec, 71045424333659ULL, -7LL, true);  /* SLPSP 1147021x61939079 */
    add_slprp_test(testsVec, 71047943319499ULL, -7LL, true);  /* SLPSP 1854131x38318729 */
    add_slprp_test(testsVec, 71051156822777ULL, 5LL, true);  /* SLPSP 5960333x11920669 */
    add_slprp_test(testsVec, 71051349391541ULL, -11LL, true);  /* SLPSP 1720603x41294447 */
    add_slprp_test(testsVec, 71054701394137ULL, 5LL, true);  /* SLPSP 4552393x15608209 */
    add_slprp_test(testsVec, 71055518472509ULL, -11LL, true);  /* SLPSP 6882611x10323919 */
    add_slprp_test(testsVec, 71055733818029ULL, -15LL, true);  /* SLPSP 580997x122299657 */
    add_slprp_test(testsVec, 71058959643359ULL, -11LL, true);  /* SLPSP 960647x73969897 */
    add_slprp_test(testsVec, 71061275485199ULL, -7LL, true);  /* SLPSP 1901x37380997099 */
    add_slprp_test(testsVec, 71061689478779ULL, -7LL, true);  /* SLPSP 5331479x13328701 */
    add_slprp_test(testsVec, 71062113271313ULL, 5LL, true);  /* SLPSP 126653x561077221 */
    add_slprp_test(testsVec, 71062298016079ULL, -7LL, true);  /* SLPSP 4866973x14600923 */
    add_slprp_test(testsVec, 71063077990277ULL, 5LL, true);  /* SLPSP 5960833x11921669 */
    add_slprp_test(testsVec, 71066892975077ULL, 5LL, true);  /* SLPSP 5960993x11921989 */
    add_slprp_test(testsVec, 83528108424479ULL, -7LL, true);    /* SLPSP 7290697x11456807 */
    add_slprp_test(testsVec, 83558429460899ULL, -11LL, true);    /* SLPSP 9141029x9141031 */
    add_slprp_test(testsVec, 9007199254741003ULL, 5LL, false);
    add_slprp_test(testsVec, 576460752303423623ULL, 5LL, false);
    add_slprp_test(testsVec, 576460752303423629ULL, -7LL, false);
    add_slprp_test(testsVec, 576460752303423641ULL, -7LL, false);
    add_slprp_test(testsVec, 576460752303423643ULL, 5LL, false);
    add_slprp_test(testsVec, 9223372036854775511ULL, -7LL, false);
    add_slprp_test(testsVec, 9223372036854775517ULL, 5LL, false);
    add_slprp_test(testsVec, 9223372036854775649ULL, -7LL, false);
    add_slprp_test(testsVec, 9223372036854775753ULL, 5LL, false);
    add_slprp_test(testsVec, 9223372036854775853ULL, 5LL, false);
    add_slprp_test(testsVec, 9223372036854776011ULL, -11LL, false);
    add_slprp_test(testsVec, 17000000000000000023ULL, 5LL, false);
    add_slprp_test(testsVec, 18446743979220271189ULL, -7LL, false);
    add_slprp_test(testsVec, uint64(-1) - 32, int64(5), false);
    // End Composite values of n

    return testsVec;
}

void
run_slprp_tests()
{
    auto const & tests = build_slprp_tests();
    
    for(auto const & tst : tests) {
        auto const n = tst.first.first;
        auto const D = tst.first.second;
        auto const expected = tst.second;
        auto const actual = is_slprp(n, D);
        std::ostringstream testRes;
        testRes << "is_slprp(" << n << ", " << D << ") = " << actual;
        for(auto i(testRes.str().size() + 1); i < 59; ++i) {
            testRes << " ";
        }
        if( expected != actual ) {
            testRes << " FAIL (Expected: " << expected << ", Actual:" << actual << ")" << std::endl;
        } else {
            testRes << " OK" << std::endl;
        }
        std::cout << testRes.str();
    }
}


void
run_is_prime_tests(std::vector<std::pair<uint64, bool>> const & is_prime_tests, std::ostream & out)
{
    uint64 count = 0, failures_count = 0;
    for(auto const & test : is_prime_tests)
    {
        auto result = is_prime(test.first);
        auto ok = (result == test.second);
        std::ostringstream testRes;
        testRes << "is_prime(" << test.first <<")";
        out << testRes.str();
        for(auto i(testRes.str().size() + 1); i < 60; ++i) {
            out << " ";
        }
        if(ok) {
            out << "OK";
        } else {
            out << "FAIL (Expected: " << test.second << " Actual: " << result << ")";
            ++failures_count;
        }
        out << std::endl;
        ++count;
    }
    if( failures_count != 0 )
    {
        out << failures_count << " tests failed (" << count << " run)." << std::endl;
    }
    else
    {
        out << "ALL " << count << " tests passed." << std::endl;
    }
}



std::vector<std::pair<uint64, bool>>
build_is_prime_tests()
{
    std::vector<std::pair<uint64, bool>> testVec;
    testVec.push_back(std::make_pair(0, false));
    testVec.push_back(std::make_pair(1, false));
    testVec.push_back(std::make_pair(2, true ));
    testVec.push_back(std::make_pair(3, true ));
    testVec.push_back(std::make_pair(4, false));
    testVec.push_back(std::make_pair(5, true ));
    testVec.push_back(std::make_pair(6, false));
    testVec.push_back(std::make_pair(7, true ));
    testVec.push_back(std::make_pair(8, false));
    testVec.push_back(std::make_pair(9, false));
    testVec.push_back(std::make_pair(499, true));
    testVec.push_back(std::make_pair(509, true));
    testVec.push_back(std::make_pair(511, false));
    testVec.push_back(std::make_pair(512, false));
    testVec.push_back(std::make_pair(751, true));
    testVec.push_back(std::make_pair(997, true));
    testVec.push_back(std::make_pair(2047, false));
    testVec.push_back(std::make_pair(9973, true));
    testVec.push_back(std::make_pair(10007, true));
    testVec.push_back(std::make_pair(15841, false));
    testVec.push_back(std::make_pair(16381, true));
    testVec.push_back(std::make_pair(29341, false));
    testVec.push_back(std::make_pair(42799, false));
    testVec.push_back(std::make_pair(49141, false));
    testVec.push_back(std::make_pair(65537, true));     // F4
    testVec.push_back(std::make_pair(1194649, false));   // 1093^2 is 2-sprp
    testVec.push_back(std::make_pair(12327121, false));  // 3511^2 is 2-sprp
    testVec.push_back(std::make_pair(25326001, false));  // 25326001 = 2251 * 11251, 2-sprp
    testVec.push_back(std::make_pair(((1ULL)<<31) - 1, true));  // M31
    testVec.push_back(std::make_pair(((1ULL)<<32) + 1, false)); // F5

    testVec.push_back(std::make_pair(39972590422099ULL, false));  /* SLPSP 203419x196503721 */
    testVec.push_back(std::make_pair(71015542332359ULL, false));  /* SLPSP 5958839x11917681 */
    testVec.push_back(std::make_pair(71026840741877ULL, false));  /* SLPSP 5959313x11918629 */
    testVec.push_back(std::make_pair(71027509003163ULL, false));  /* SLPSP 3440627x20643769 */
    testVec.push_back(std::make_pair(71028048949859ULL, false));  /* SLPSP 8039x8835433381 */
    testVec.push_back(std::make_pair(71028425928179ULL, false));  /* SLPSP 5959381x11918759 */
    testVec.push_back(std::make_pair(71028496858199ULL, false));  /* SLPSP 2541089x27951991 */
    testVec.push_back(std::make_pair(71029125005299ULL, false));  /* SLPSP 13721x34301x150919 */
    testVec.push_back(std::make_pair(71029563550243ULL, false));  /* SLPSP 426763x166437961 */
    testVec.push_back(std::make_pair(71030483737309ULL, false));  /* SLPSP 7298827x9731767 */
    testVec.push_back(std::make_pair(71035764857873ULL, false));  /* SLPSP 2665253x26652541 */
    testVec.push_back(std::make_pair(71036072324099ULL, false));  /* SLPSP 8428289x8428291 */
    testVec.push_back(std::make_pair(71040345630799ULL, false));  /* SLPSP 3349163x21211373 */
    testVec.push_back(std::make_pair(71040523518089ULL, false));  /* SLPSP 878737x80843897 */
    testVec.push_back(std::make_pair(71040611030759ULL, false));  /* SLPSP 4214279x16857121 */
    testVec.push_back(std::make_pair(71041841765927ULL, false));  /* SLPSP 2665367x26653681 */
    testVec.push_back(std::make_pair(71042052321701ULL, false));  /* SLPSP 4214323x16857287 */
    testVec.push_back(std::make_pair(71044690405037ULL, false));  /* SLPSP 397337x178802101 */
    testVec.push_back(std::make_pair(71045053572089ULL, false));  /* SLPSP 5960077x11920157 */
    testVec.push_back(std::make_pair(71045424333659ULL, false));  /* SLPSP 1147021x61939079 */
    testVec.push_back(std::make_pair(71047943319499ULL, false));  /* SLPSP 1854131x38318729 */
    testVec.push_back(std::make_pair(71051156822777ULL, false));  /* SLPSP 5960333x11920669 */
    testVec.push_back(std::make_pair(71051349391541ULL, false));  /* SLPSP 1720603x41294447 */
    testVec.push_back(std::make_pair(71054701394137ULL, false));  /* SLPSP 4552393x15608209 */
    testVec.push_back(std::make_pair(71055518472509ULL, false));  /* SLPSP 6882611x10323919 */
    testVec.push_back(std::make_pair(71055733818029ULL, false));  /* SLPSP 580997x122299657 */
    testVec.push_back(std::make_pair(71058959643359ULL, false));  /* SLPSP 960647x73969897 */
    testVec.push_back(std::make_pair(71061275485199ULL, false));  /* SLPSP 1901x37380997099 */
    testVec.push_back(std::make_pair(71061689478779ULL, false));  /* SLPSP 5331479x13328701 */
    testVec.push_back(std::make_pair(71062113271313ULL, false));  /* SLPSP 126653x561077221 */
    testVec.push_back(std::make_pair(71062298016079ULL, false));  /* SLPSP 4866973x14600923 */
    testVec.push_back(std::make_pair(71063077990277ULL, false));  /* SLPSP 5960833x11921669 */
    testVec.push_back(std::make_pair(71066892975077ULL, false));  /* SLPSP 5960993x11921989 */
    testVec.push_back(std::make_pair(83528108424479ULL, false));  /* SLPSP 7290697x11456807 */
    testVec.push_back(std::make_pair(83558429460899ULL, false));  /* SLPSP 9141029x9141031 */

    testVec.push_back(std::make_pair(((1ULL)<<61) - 1, true));  // M61
    testVec.push_back(std::make_pair(4611686014132420609ULL, false)); // M31^2
    testVec.push_back(std::make_pair((1ULL << 63) +  0, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  1, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  2, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  3, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  4, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  5, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  6, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  7, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  8, false));
    testVec.push_back(std::make_pair((1ULL << 63) +  9, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 10, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 11, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 12, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 13, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 14, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 15, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 16, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 17, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 18, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 19, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 20, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 21, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 22, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 23, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 24, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 25, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 26, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 27, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 28, false));
    testVec.push_back(std::make_pair((1ULL << 63) + 29, true));

    testVec.push_back(std::make_pair(uint64(-1)-58, true));
    testVec.push_back(std::make_pair(uint64(-1)-32, false));
    testVec.push_back(std::make_pair(uint64(-1)-2, false));
    testVec.push_back(std::make_pair(uint64(-1), false));

    return testVec;
}


template <typename T>
void
add_jacobi_test(T & container, int64 a, uint64 m, int64 r)
{
    container.push_back(std::make_pair(std::make_pair(a, m), r));
}

std::vector<std::pair<std::pair<int64, uint64>, int64>>
build_jacobi_symbol_tests()
{
    std::vector<std::pair<std::pair<int64, uint64>, int64>> tests;
    add_jacobi_test(tests,               0LL,                     1ULL,     1LL);
    add_jacobi_test(tests,               0LL,                     2ULL,     0LL);
    add_jacobi_test(tests,               0LL,                    17ULL,     0LL);
    add_jacobi_test(tests,               1LL,                     1ULL,     1LL);
    add_jacobi_test(tests,               1LL,                     2ULL,     1LL);
    add_jacobi_test(tests,               1LL,                     3ULL,     1LL);
    add_jacobi_test(tests,               2LL,                   127ULL,     1LL);
    add_jacobi_test(tests,               2LL,                   125ULL,    -1LL);
    add_jacobi_test(tests,              -5LL,                   127ULL,     1LL);
    add_jacobi_test(tests,               5LL,                   127ULL,    -1LL);
    add_jacobi_test(tests,             -11LL,                   127ULL,    -1LL);
    add_jacobi_test(tests,              14LL,                     7ULL,     0LL);
    add_jacobi_test(tests,            1212LL,            1200000007ULL,    -1LL);
    add_jacobi_test(tests,            1236LL,                 20003ULL,     1LL);
    add_jacobi_test(tests,            1411LL,                   317ULL,    -1LL);
    add_jacobi_test(tests,            1735LL,                   507ULL,     1LL);
    add_jacobi_test(tests,          222222LL,                304679ULL,    -1LL);
    add_jacobi_test(tests,          222222LL,                324899ULL,     1LL);
    add_jacobi_test(tests,        92177777LL,           12000008007ULL,    -1LL);
    add_jacobi_test(tests,      9000000007LL,         (1ULL << 63) + 1,    -1LL);
    add_jacobi_test(tests,     -9000000007LL,         (1ULL << 63) + 1,    -1LL);
    add_jacobi_test(tests,     12678676777LL,       222222222222223ULL,     1LL);
    add_jacobi_test(tests,     70001000387LL,               uint64(-1),     1LL);
    add_jacobi_test(tests,    -70001000387LL,               uint64(-1),    -1LL);
    add_jacobi_test(tests,  12299356789231LL,             678675867ULL,     1LL);
    add_jacobi_test(tests,  54665867345453LL,  17199909785980786011ULL,     0LL);
    add_jacobi_test(tests, int64(1ULL << 63),  17199909785980786011ULL,     1LL);
    add_jacobi_test(tests, int64(1ULL << 63),               uint64(-1),    -1LL);

    return tests;
}

void
run_jacobi_symbol_tests( 
      std::vector<std::pair<std::pair<int64, uint64>, int64>> const & tests
    , std::ostream & out
    )
{
    for(auto const & test : tests)
    {
        auto const & input = test.first;
        auto expected = test.second;
        auto actual = jacobi_symbol(input.first, input.second);
        std::ostringstream ostr;
        ostr << "jacobi_symbol(" << input.first 
                    << ", " << input.second << ") ";
        out << ostr.str();
        for(auto i(ostr.str().size()), i_max(decltype(ostr.str().size())(59)); i < i_max; ++i)
        {
            out << " ";
        }

        if( actual != expected ) {
            out << "FAIL (Expected: " << expected << "  Actual: " << actual << ")";
        } else {
            out << "OK";
        }
        out << std::endl;
    }
}


int main(int argc, char** argv)
{
    if( argc == 2 ) {
        if( !strcmp(argv[1], "-t") ) {
            run_jacobi_symbol_tests(build_jacobi_symbol_tests(), std::cout);
            run_slprp_tests();
            run_is_prime_tests(build_is_prime_tests(), std::cout);
        }
    } else {
        auto hasStart = false, hasLength = false, hasLimit = false, isCountMode = false;
        uint64 start, length, limit;
        for(auto pargv = argv + 1; *pargv; ++pargv) {
            if( ! strcmp(*pargv, "-c") ) {
                isCountMode = true;
            } else if( ! strcmp(*pargv, "-s") ) {
                if( *(pargv + 1) ) {
                    hasStart = true;
                    start =  _strtoui64(*(pargv + 1), nullptr, 10);
                    ++pargv;
                } else {
                    break;
                }
            } else if( ! strcmp(*pargv, "-l") ) {
                if( *(pargv + 1) ) {
                    hasLength = true;
                    length =  _strtoui64(*(pargv + 1), nullptr, 10);
                    ++pargv;
                } else {
                    break;
                }
            } else if( ! strcmp(*pargv, "-p") ) {
                if( *(pargv + 1) ) {
                    hasLimit = true;
                    limit =  _strtoui64(*(pargv + 1), nullptr, 10);
                    ++pargv;
                } else {
                    break;
                }
            }
        }
        if( !(hasStart && hasLength) && !hasLimit ) {
            std::cerr << "Invalid command line" << std::endl;
            return 1;
        }
        if( isCountMode ) {
            auto n0 = (hasStart && hasLength) ? start : 0ull;
            auto n1 = (hasStart && hasLength) ? start + length : limit;
            auto count = 0ull;
            auto startTime = std::chrono::high_resolution_clock::now();
            for(auto n = n0; n < n1; ++n) {
                if( is_prime(n) ) {
                    ++count;
                }
            }
            std::cout << "Seconds: " << 
                std::chrono::duration<double>(std::chrono::high_resolution_clock::now() 
                        - startTime).count() << '\n'
                      << "Primes: " << count << std::endl;
            return 0;
        }
        if( hasStart && hasLength ) {
            for(auto n = start; n < start + length; ++n) {
                if( is_prime(n) ) {
                    std::cout << n << "\n";
                }
            }
            return 0;
        }
        if( hasLimit ) {
            print_primes_below(limit, std::cout);
            return 0;
        }
    }
    return 1;
}

