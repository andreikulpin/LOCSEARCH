/* 
 * File:   testfunc.hpp
 * Author: kulpin
 */

#ifndef TESTFUNC_HPP
#define TESTFUNC_HPP

#include <mpproblem.hpp>
#include <box/box.hpp>
#include <common/utilmacro.hpp>

namespace OPTITEST {

    class TestObjective : public COMPI::Functor <double> {
    public:

        TestObjective(int n) : mN(n) {
        }

        double func(const double* x) const {
//            return 4 * SGSQR(x[0] - 5) + SGSQR(x[1] - 6);
            return 100 * SGSQR(x[1] - x[0] * x[0]) + SGSQR(1 - x[1]);
        }

    private:
        int mN;
    };

    /**
     * Factory to produce instances of test problem
     */
    class TestProblemFactory {
    public:

        /**
         * Constructor
         * @param n problem dimension
         * @param a left border for a box
         * @param b right border for a box
         */
        TestProblemFactory(int n, double a, double b) : mN(n), mA(a), mB(b) {
        }

        COMPI::MPProblem<double>* getProblem() const {
            COMPI::MPProblem<double>* prob = new COMPI::MPProblem<double>();
            prob->mVarTypes.assign(mN, COMPI::MPProblem<double>::GENERIC);
            prob->mObjectives.push_back(std::make_shared<TestObjective>(mN));
            prob->mBox = new snowgoose::Box<double>(mN);
            for (int i = 0; i < mN; i++) {
                prob->mBox->mA[i] = mA;
                prob->mBox->mB[i] = mB;
            }
            return prob;
        }

    private:
        int mN;
        double mA;
        double mB;
    };

}

#endif /* TESTFUNC_HPP */

