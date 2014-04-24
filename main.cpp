/*
* Creation date: 2014.04.21
* Creator: Youcef Lemsafer
* Authors: Youcef Lemsafer
*/
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>

typedef long long int64;
typedef unsigned long long uint64;

/*
* We assume that P = 1 and Q = (1 - D)/4
* where D is the first element in {5, -7, 9, -11, ...}
* such that the Jacobi symbol (D/n) = -1.
*/
extern "C" bool is_slprp(uint64 n, int64 D);


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


int main(int argc, char** argv)
{
    if( argc == 2 ) {
        if( !strcmp(argv[1], "-t") ) {
            run_slprp_tests();
        }
    }
}

