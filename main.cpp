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
    add_slprp_test(testsVec, uint64(5459),   int64(-7),  true);
    add_slprp_test(testsVec, uint64(5777),   int64(5),   true);
    add_slprp_test(testsVec, uint64(10877),  int64(5),   true);
    add_slprp_test(testsVec, uint64(16109),  int64(13),  true);
    add_slprp_test(testsVec, uint64(18971),  int64(-11), true);
    add_slprp_test(testsVec, uint64(22499),  int64(-15), true);
    add_slprp_test(testsVec, uint64(24569),  int64(-7),  true);
    add_slprp_test(testsVec, uint64(25199),  int64(-7),  true);
    add_slprp_test(testsVec, uint64(40309),  int64(-7),  true);
    add_slprp_test(testsVec, uint64(58519),  int64(-7),  true);
    add_slprp_test(testsVec, uint64(75077),  int64(5),   true);
    add_slprp_test(testsVec, uint64(97439),  int64(-7),  true);
    add_slprp_test(testsVec, uint64(100127), int64(5),   true);
    add_slprp_test(testsVec, uint64(113573), int64(5),   true);
    add_slprp_test(testsVec, uint64(115639), int64(-7),  true);
    add_slprp_test(testsVec, uint64(130139), int64(-15), true);
    add_slprp_test(testsVec, uint64(155819), int64(-7),  true);
    add_slprp_test(testsVec, uint64(158399), int64(-7),  true);
    add_slprp_test(testsVec, uint64(161027), int64(5),   true);
    add_slprp_test(testsVec, uint64(162133), int64(5),   true);
    add_slprp_test(testsVec, uint64(176399), int64(-7),  true);
    add_slprp_test(testsVec, uint64(176471), int64(-15), true);
    add_slprp_test(testsVec, uint64(189419), int64(-7),  true);
    add_slprp_test(testsVec, uint64(192509), int64(13),  true);
    add_slprp_test(testsVec, uint64(197801), int64(-11), true);
    add_slprp_test(testsVec, uint64(224369), int64(-7),  true);
    add_slprp_test(testsVec, uint64(230691), int64(-7),  true);
    add_slprp_test(testsVec, uint64(231703), int64(5),   true);
    add_slprp_test(testsVec, uint64(243629), int64(-15), true);
    add_slprp_test(testsVec, uint64(253259), int64(-7),  true);
    add_slprp_test(testsVec, uint64(268349), int64(-15), true);
    add_slprp_test(testsVec, uint64(288919), int64(13),  true);
    add_slprp_test(testsVec, uint64(313499), int64(-11), true);
    add_slprp_test(testsVec, uint64(324899), int64(-15), true);
    // End of values taken from A217255

    // Begin Prime values of n
    add_slprp_test(testsVec, uint64(4294967311), int64(-7), true);
    add_slprp_test(testsVec, uint64(4294967357), int64(5), true);
    add_slprp_test(testsVec, uint64(4294967371), int64(-11), true);
    add_slprp_test(testsVec, uint64(1099511627791), int64(-7), true);
    add_slprp_test(testsVec, uint64(1099511627803), int64(5), true);
    add_slprp_test(testsVec, uint64(1099511627831), int64(13), true);
    add_slprp_test(testsVec, uint64(9007199254740997), int64(5), true);
    add_slprp_test(testsVec, uint64(576460752303423619), int64(13), true);
    add_slprp_test(testsVec, uint64(576460752303423649), int64(-11), true);
    add_slprp_test(testsVec, uint64(576460752303423733), int64(5), true);
    add_slprp_test(testsVec, uint64(576460752303423737), int64(5), true);
    add_slprp_test(testsVec, uint64(576460752303423749), int64(-7), true);
    add_slprp_test(testsVec, uint64(576460752303423761), int64(13), true);
    add_slprp_test(testsVec, uint64(700000000000000289), int64(-11), true);
    add_slprp_test(testsVec, uint64(2305843009213693921), int64(-7), true);
    add_slprp_test(testsVec, uint64(2305843009213693951), int64(17), true);
    add_slprp_test(testsVec, uint64(2305843009213693967), int64(5), true);
    add_slprp_test(testsVec, uint64(9000000000000000041), int64(-11), true);
    // End Prime values of n

    // Begin Composite values of n
    // 2047 = 23*89 is a strong pseudo prime base 2
    add_slprp_test(testsVec, uint64(2047),       int64(5),  false);
    add_slprp_test(testsVec, uint64(4294967333), int64(5),  false);
    add_slprp_test(testsVec, uint64(4294967359), int64(13), false);
    add_slprp_test(testsVec, uint64(1099511627801), int64(-7), false);
    add_slprp_test(testsVec, uint64(1099511627813), int64(5), false);
    add_slprp_test(testsVec, uint64(9007199254741003), int64(5), false);
    add_slprp_test(testsVec, uint64(576460752303423623), int64(5), false);
    add_slprp_test(testsVec, uint64(576460752303423629), int64(-7), false);
    add_slprp_test(testsVec, uint64(576460752303423641), int64(-7), false);
    add_slprp_test(testsVec, uint64(576460752303423643), int64(5), false);
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

