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
    add_slprp_test(testsVec, 4294967311ULL, -7LL, true);
    add_slprp_test(testsVec, 4294967357ULL, 5LL, true);
    add_slprp_test(testsVec, 4294967371ULL, -11LL, true);
    add_slprp_test(testsVec, 263456789093ULL, 5LL, true);
    add_slprp_test(testsVec, 1099511627791ULL, -7LL, true);
    add_slprp_test(testsVec, 1099511627803ULL, 5LL, true);
    add_slprp_test(testsVec, 1099511627831ULL, 13LL, true);
    add_slprp_test(testsVec, 9007199254740997ULL, 5LL, true);
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
    add_slprp_test(testsVec, 9000000000000000041ULL, -11LL, true);
    add_slprp_test(testsVec, 10000000000000000051ULL, -7LL, true);
    add_slprp_test(testsVec, 12345678901234567891ULL, -11L, true);
    add_slprp_test(testsVec, 17000000000000000003ULL, 5LL, true);
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
    add_slprp_test(testsVec, 83528108424479ULL, -7LL, true);    /* SLPSP 7290697x11456807 */
    add_slprp_test(testsVec, 83558429460899ULL, -11LL, true);    /* SLPSP 9141029x9141031 */
    add_slprp_test(testsVec, 9007199254741003ULL, 5LL, false);
    add_slprp_test(testsVec, 576460752303423623ULL, 5LL, false);
    add_slprp_test(testsVec, 576460752303423629ULL, -7LL, false);
    add_slprp_test(testsVec, 576460752303423641ULL, -7LL, false);
    add_slprp_test(testsVec, 576460752303423643ULL, 5LL, false);
    add_slprp_test(testsVec, 17000000000000000023ULL, 5LL, false);
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

