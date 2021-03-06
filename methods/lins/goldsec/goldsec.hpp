/* 
 * File:   dichotls.hpp
 * Author: medved
 *
 * Created on February 26, 2016, 6:45 PM
 */

#ifndef GOLDENSEC_HPP
#define GOLDENSEC_HPP

#include <vector>
#include <sstream>
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
            FT mSInit = 0.1;
            /**
             * Max steps forward to estimate minimum
             */
            FT mMaxForwardSteps = 16;
            /**
             * Max steps to backtrack when estimating the interval
             */
            FT mMaxBackSteps = 2;
            /**
             * Delta when to stop search
             */
            FT mDelta = 0.1;
            /**
             * Trace on/off
             */
            bool mDoTracing = false;

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
            FT lbest = 0;
            FT fbest = v;
            constexpr FT rho = 0.382;
            bool forward = true;
            bool goodDir = false;
            bool doGoldenSearch = false;
            const int n = mProblem.mVarTypes.size();
            FT* xnew = mXnew.data();
            const snowgoose::Box<FT>& box = *(mProblem.mBox);
            auto obj = mProblem.mObjectives.at(0);
            int i = 0;
            snowgoose::VecUtils::vecCopy(n, d, mDir.data());
            FT* dir = mDir.data();

            for (;;) {
                if (mOptions.mDoTracing) {
                    std::cout << " i = " << i << ", " << (forward ? "forward" : "back")
                            << ", l = " << l << ", fbest = " << fbest << "\n";
                }
                // TMP FIX!!!
                //snowgoose::BoxUtils::projectDirection(dir, x, box);

                snowgoose::VecUtils::vecSaxpy(n, x, dir, l, xnew);
                FT fn;
                if (!snowgoose::BoxUtils::isIn(xnew, box)) {
                    if (mOptions.mDoTracing)
                        std::cout << "OUT OF BOX\n";
                    FT alpha = l;
                    for (int j = 0; j < n; j++) {
                        FT nalpha = alpha;
                        if (dir[j] > 0) {
                            nalpha = (box.mB[j] - x[j]) / dir[j];
                        } else if (dir[j] < 0) {
                            nalpha = (box.mA[j] - x[j]) / dir[j];
                        }
                        alpha = SGMIN(alpha, nalpha);
                    }
                    l = alpha;
                    snowgoose::VecUtils::vecSaxpy(n, x, dir, l, xnew);
                    fn = obj->func(xnew);
                    if (mOptions.mDoTracing)
                        std::cout << "alpha = " << alpha << ", fn = " << fn << "\n";
                    if (fn <= fbest) {
                        fbest = fn;
                        lbest = l;
                        break;
                    }
                } else {
                    fn = obj->func(xnew);
                }
                if (fn < fbest) {
                    goodDir = true;
                    fbest = fn;
                    lbest = l;
                    if (mOptions.mDoTracing)
                        std::cout << "S: " << fn << ", l = " << l << "\n";
                    if (forward) {
                        if (i > mOptions.mMaxForwardSteps)
                            break;
                        lpp = lp;
                        lp = l;
                        l = lp + (lp - lpp) * ((1 - rho) / rho);
                        if (mOptions.mDoTracing)
                            std::cout << "lpp = " << lpp << ", lp = " << lp << ", l = " << l << "\n";
                    } else {
                        doGoldenSearch = true;
                        break;
                    }
                } else {
                    if (mOptions.mDoTracing)
                        std::cout << "F: " << fn << ", l = " << l << "\n";
                    if (i == 0)
                        forward = false;

                    if (forward) {
                        doGoldenSearch = true;
                        break;
                    } else {
                        if (i >= mOptions.mMaxBackSteps)
                            break;
                        lp = l;
                        l *= rho;
                    }
                }
                i++;
            }

            if (doGoldenSearch && (mOptions.mDelta > 0)) {
                if (mOptions.mDoTracing)
                    std::cout << "Start golden search\n";
                const FT beg = forward ? lpp : 0;
                const FT stretch = forward ? (l - lpp) : lp;

                auto remap = [&] (FT coe) {
                    return beg + stretch * coe;
                };

                FT a = 0;
                FT b = 1;
                FT L = rho;
                FT R = 1 - rho;


                auto getv = [&](FT coe) {
                    snowgoose::VecUtils::vecSaxpy(n, x, dir, remap(coe), xnew);
                    return obj->func(xnew);
                };
                FT fL = fbest;
                FT fR = getv(R);
                if (mOptions.mDoTracing) {
                    std::cout << "Before loop: fL = " << fL << ", fR = " << fR << "\n";
                    //std::cout << "Check  fL = " << getv(L) << ", fR = " << getv(R) << "\n";

                }
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
                    if (mOptions.mDoTracing) {
                        std::cout << "In loop: D = " << b - a << ", fL = " << fL << ", fR = " << fR << "\n";
                        //std::cout << "Check  fL = " << getv(L) << ", fR = " << getv(R) << "\n";
                    }
                }
                if (fL <= fR) {
                    if (fL < fbest) {
                        fbest = fL;
                        lbest = remap(L);
                    }
                } else {
                    if (fR < fbest) {
                        fbest = fR;
                        lbest = remap(R);
                    }
                }
            }

            if (fbest < v) {
                snowgoose::VecUtils::vecSaxpy(n, x, dir, lbest, x);
                v = fbest;
                return true;
            } else
                return false;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Golden Section Line Search with \n";
            os << "Initial step = " << mOptions.mSInit << "\n";
            os << "Maximal reduction steps to locate minimum = " << mOptions.mMaxBackSteps << "\n";
            os << "Maximal increase steps to locate minimum = " << mOptions.mMaxForwardSteps << "\n";
            os << "Minimal interval for minimum size = " << mOptions.mDelta << "\n";
            return os.str();
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

