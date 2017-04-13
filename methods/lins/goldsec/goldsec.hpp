/* 
 * File:   dichotls.hpp
 * Author: medved
 *
 * Created on February 26, 2016, 6:45 PM
 */

#ifndef GOLDENSEC_HPP
#define GOLDENSEC_HPP

#include <vector>

#include <common/lineseach.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>


namespace LOCSEARCH {

    /**
     * Golden section line search
     */
    template <class FT> class GoldenSecLS : public LineSearch<FT> {
    public:

        struct Options {
            /**
             * Initial step
             */
            FT mSInit = 0.01;

            /**
             * Max steps forward to estimate minimum
             */
            FT mMaxForwardSteps = 8;

            /**
             * Max steps to backtrack when estimating the interval
             */
            FT mMaxBackSteps = 2;
            /**
             * Delta when to stop search
             */
            FT mDelta = 0.1;

        };

        /**
         * The constructor
         *
         * @param prob optimization problem
         */
        GoldenSecLS(const COMPI::MPProblem<FT>& prob) :
        mProblem(prob) {
            mXnew.resize(prob.mVarTypes.size());
            mDir.resize(prob.mVarTypes.size());
        }

        bool search(const FT* d, FT* x, FT& v) override {
            FT lp = 0;
            FT lpp = 0;
            FT l = mOptions.mSInit;
            FT fcur = v;
            constexpr FT rho = 0.382;
            bool forward = true;
            bool goodDir = false;
            bool rv = false;
            const int n = mProblem.mVarTypes.size();
            FT* xnew = mXnew.data();
            const snowgoose::Box<FT>& box = *(mProblem.mBox);
            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            for (int i = 0; i < mOptions.mMaxForwardSteps; i++) {
                std::cout << " i = " << i << ", " << (forward ? "forward" : "back") << "\n";
                snowgoose::VecUtils::vecSaxpy(n, x, d, l, xnew);
                if (!snowgoose::BoxUtils::isIn(xnew, box)) {
                    std::cout << "OUT OF BOX\n";
                    break;
                }
                FT fn = obj->func(xnew);
                if (fn < fcur) {
                    goodDir = true;
                    fcur = fn;
                    std::cout << "S: " << fn << ", l = " << l << "\n";
                    if (forward) {
                        lpp = lp;
                        lp = l;
                        l = lp + (lp - lpp) * ((1 - rho) / rho);
                        std::cout << "lpp = " << lpp << ", lp = " << lp << ", l = " << l << "\n";
                    } else {
                        break;
                    }
                } else {
                    std::cout << "F: " << fn << ", l = " << l << "\n";
                    if (i == 0) {
                        forward = false;
                        // TMP
                        //break;
                    }
                    if (forward) {
                        break;
                    } else {
                        if (i > mOptions.mMaxBackSteps) {
                            break;
                        }
                        lp = l;
                        l *= rho;
                        continue;
                    }
                }
            }

            if (goodDir && (mOptions.mDelta > 0)) {
                std::cout << "Start golden search\n";
                if (forward)
                    snowgoose::VecUtils::vecSaxpy(n, x, d, lpp, x);
                else
                    snowgoose::VecUtils::vecSaxpy(n, x, d, lp, xnew);
                FT *ndir = mDir.data();
                snowgoose::VecUtils::vecSaxpy(n, xnew, x, -1., ndir);
                FT a = 0;
                FT b = 1;
                FT L = rho;
                FT R = 1 - rho;

                auto getv = [&](FT coe) {
                    snowgoose::VecUtils::vecSaxpy(n, x, ndir, coe, xnew);
                    return obj->func(xnew);
                };
                FT fL = fcur;
                FT fR = getv(R);
                std::cout << "Before loop: fL = " << fL << ", fR = " << fR << "\n";
                while (b - a > mOptions.mDelta) {
                    if (fL <= fR) {
                        b = R;
                        R = L;
                        fR = fL;
                        L = a + rho * (b - a);
                        fL = getv(L);
                    } else {
                        a = L;
                        L = R;
                        fL = fR;
                        R = a + (1 - rho) * (b - a);
                        fR = getv(R);
                    }
                    std::cout << "In loop: D = " << b - a << ", fL = " << fL << ", fR = " << fR << "\n";
                }
                if (fL <= fR) {
                    if (fL < fcur) {
                        snowgoose::VecUtils::vecSaxpy(n, x, ndir, L, x);
                        fcur = fL;
                    }
                } else if (fR < fcur) {
                    snowgoose::VecUtils::vecSaxpy(n, x, ndir, R, x);
                    fcur = fR;
                }
            }

            if (fcur < v) {
                v = fcur;
                return true;
            } else
                return false;
        }

        std::string about() const {
            return "Golden section line search.";
        }

        /**
         * Retrieve options reference
         * @return options reference
         */
        Options& getOptions() {
            return mOptions;
        }

    private:
        std::vector<FT> mXnew;
        std::vector<FT> mDir;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
    };
}

#endif /* GOLDENSEC_HPP */

