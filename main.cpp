/*
* Creation date: 2014.04.21
* Creator: Youcef Lemsafer
* Authors: Youcef Lemsafer
*/
#include <iostream>
#include <vector>

typedef long long int64;
typedef unsigned long long uint64;

/*
* We assume that P = 1 and Q = (1 - D)/4
* where D is the first element in {5, -7, 9, -11, ...}
* such that the Jacobi symbol (D/n) = -1.
*/
extern "C" bool is_slprp(uint64 n, int64 D);

std::vector<std::pair<uint64, int64> >
build_slprp_tests()
{
    std::vector<std::pair<uint64, int64> > testsVec;
    testsVec.push_back(std::make_pair(uint64(2047), int64(5)));
    testsVec.push_back(std::make_pair(uint64(5459), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(5777), int64(5)));
    testsVec.push_back(std::make_pair(uint64(10877), int64(5)));
    testsVec.push_back(std::make_pair(uint64(16109), int64(13)));
    testsVec.push_back(std::make_pair(uint64(18971), int64(-11)));
    testsVec.push_back(std::make_pair(uint64(22499), int64(-15)));
    testsVec.push_back(std::make_pair(uint64(24569), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(25199), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(40309), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(58519), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(75077), int64(5)));
    testsVec.push_back(std::make_pair(uint64(97439), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(100127), int64(5)));
    testsVec.push_back(std::make_pair(uint64(113573), int64(5)));
    testsVec.push_back(std::make_pair(uint64(115639), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(130139), int64(-15)));
    testsVec.push_back(std::make_pair(uint64(155819), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(158399), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(161027), int64(5)));
    testsVec.push_back(std::make_pair(uint64(162133), int64(5)));
    testsVec.push_back(std::make_pair(uint64(176399), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(176471), int64(-15)));
    testsVec.push_back(std::make_pair(uint64(189419), int64(-7)));
    testsVec.push_back(std::make_pair(uint64(192509), int64(13)));
    testsVec.push_back(std::make_pair(uint64(197801), int64(-11)));
    return testsVec;
}

void
run_slprp_tests()
{
    std::vector<std::pair<uint64, int64> > const & tests = build_slprp_tests();
    
    for(auto const & tst : tests) {
        std::cout << "is_slprp(" << tst.first 
            << ", " << tst.second << ") = " 
            << is_slprp(tst.first, tst.second) << std::endl;
    }
}

int main(int argc, char** argv)
{
    run_slprp_tests();
}

