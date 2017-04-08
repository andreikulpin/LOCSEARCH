/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef VARCOORGRAD_HPP
#define  VARCOORGRAD_HPP

#include <sstream>
#include <vector>
#include <functional>
#include <solver.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>


namespace LOCSEARCH {

    /**
     * Coordinate descent for box constrained problems with unequal dynamically adjacent
     * steps along directions combined with gradient descent
     */
    template <typename FT> class VarCoorGrad : public COMPI::Solver<FT> {
    public:

        /**
         * Determines stopping conditions
         * @param fval current best value found
         * @param x current best point
         * @param stpn step number
         * @return true if the search should stop
         */
        using Stopper = std::function<bool(FT fval, const FT* x, int stpn) >;

        /**
         * Watches the current step
         * @param fval current best value found
         * @param current best point
         * @param stpn step number
         * @param gran - current granularity vector
         */
        using Watcher = std::function<void(FT fval, const FT* x, const std::vector<FT>& gran, int stpn) >;

        /**
         * Options for Gradient Box Descent method
         */
        struct Options {
            /**
             * Initial value of granularity
             */
            FT mHInit = 0.01;

            /**
             * Increase in the case of success
             */
            FT mInc = 1.75;
            /**
             * Decrease in the case of failure
             */
            FT mDec = 0.5;
            /**
             * Lower bound for granularity
             */
            FT mHLB = 1e-08;
            /**
             * Upper bound on granularity
             */
            FT mHUB = 1e+02;
            /**
             * Gradient descent multiplier (if <= 0 then don't do gradient step)
             */
            FT mGradStep = 1;
            /**
             * Maximal number of consecutive gradient steps
             */
            FT mGradMaxSteps = 8;
            /**
             * Gradient search speedup parameter
             */
            FT mGradSpeedup = 2;
            /**
             * Golden search delta (set <= 0 to cancel golden search)
             */
            FT mGoldenSearchDelta = 0.1;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        VarCoorGrad(const COMPI::MPProblem<FT>& prob) :
        mProblem(prob) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
        }

        /**
         * Perform search
         * @param x start point and result
         * @param v  the resulting value
         * @return true if search converged and false otherwise
         */
        bool search(FT* x, FT& v) override {
            bool rv = false;
            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT fcur = obj->func(x);
            int n = mProblem.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            std::vector<FT> sft;
            sft.assign(n, mOptions.mHInit);
            FT* xold = new FT[n];
            FT* xnew = new FT[n];
            FT* grad = new FT[n];
            FT* ndir = new FT[n];

            auto inc = [this] (FT h) {
                FT t = h * mOptions.mInc;
                t = SGMIN(t, mOptions.mHUB);
                return t;
            };

            auto dec = [this](FT h) {
                FT t = h * mOptions.mDec;
                t = SGMAX(t, mOptions.mHLB);
                return t;
            };

            int dir = 1;


            auto step = [&] () {
                for (int i = 0; i < n;) {
                    const FT h = sft[i];
                    const FT dh = h;
                    FT y = x[i] + dir * dh;
                    y = SGMAX(y, box.mA[i]);
                    y = SGMIN(y, box.mB[i]);
                    const FT tmp = x[i];
                    x[i] = y;
                    const FT fn = obj->func(x);
                    const FT dx = y - tmp;
                    grad[i] = (dx == 0) ? 0 : (fn - fcur) / dx;
                    if (fn >= fcur) {
                        x[i] = tmp;
                        if (dir == 1)
                            dir = -1;
                        else {
                            sft[i] = dec(h);
                            i++;
                            dir = 1;
                        }
                    } else {
                        sft[i] = inc(h);
                        fcur = fn;
                        dir = 1;
                        i++;
                    }
                }
            };

            int sn = 0;
            bool br = false;
            while (!br) {
                sn++;
                FT fold = fcur;
                step();
                if (fcur == fold) {
                    FT H = snowgoose::VecUtils::maxAbs(n, sft.data(), nullptr);
                    if (H <= mOptions.mHLB)
                        break;
                } else
                    rv = true;
                if (mOptions.mGradStep > 0) {
                    snowgoose::VecUtils::vecCopy(n, x, xold);
                    FT l = mOptions.mGradStep * snowgoose::VecUtils::maxAbs(n, sft.data(), nullptr);
                    //FT l = mOptions.mGradStep;
                    for (int i = 0; i < mOptions.mGradMaxSteps; i++) {
                        snowgoose::VecUtils::vecSaxpy(n, x, grad, -l, xnew);
                        if (!snowgoose::BoxUtils::isIn(xnew, box)) {
                            std::cout << "OUT OF BOX\n";
                            break;
                        }
                        snowgoose::BoxUtils::project(xnew, box);
                        FT fn = obj->func(xnew);
                        if (fn < fcur) {
                            snowgoose::VecUtils::vecCopy(n, xnew, x);
                            fcur = fn;
                            std::cout << "S: " << fn << ", l = " << l << "\n";
                            l *= mOptions.mGradSpeedup;
                        } else {
                            std::cout << "F: " << fn << ", l = " << l << "\n";
                            if (mOptions.mGoldenSearchDelta > 0) {
                                std::cout << "Start golden search\n";
                                constexpr FT rho = 0.382;
                                FT a = 0;
                                FT b = 1;
                                FT L = a + rho * (b - a);
                                FT R = a + (1 - rho) * (b - a);
                                snowgoose::VecUtils::vecSaxpy(n, xnew, xold, -1., ndir);
                                auto getv = [&](FT coe) {
                                    snowgoose::VecUtils::vecSaxpy(n, xold, ndir, coe, xnew);
                                    return obj->func(xnew);
                                };
                                FT fL = getv(L);
                                FT fR = getv(R);
                                std::cout << "Befor loop: fcur = " << fcur <<"fL = " << fL << ", fR = " << fR << "\n";
                                while (b - a > mOptions.mGoldenSearchDelta) {
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
                                    std::cout << "In loop: fL = " << fL << ", fR = " << fR << "\n";
                                }
                                if (fL <= fR) {
                                    if (fL < fcur) {
                                        snowgoose::VecUtils::vecSaxpy(n, xold, ndir, L, x);
                                        fcur = fL;
                                    }
                                } else if(fR < fcur) {
                                        snowgoose::VecUtils::vecSaxpy(n, xold, ndir, R, x);
                                        fcur = fR;
                                }

                            }
                            break;
                        }
                    }
                }

                for (auto w : mWatchers) {
                    w(fcur, x, sft, sn);
                }
                for (auto s : mStoppers) {
                    if (s(fcur, x, sn)) {
                        br = true;
                        break;
                    }
                }
            }
            v = fcur;
            delete [] xnew;
            delete [] grad;
            delete [] xold;
            delete [] ndir;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Coordinate descent method with variable adaptation combined with gradient descent\n";
            os << "Initial step = " << mOptions.mHInit << "\n";
            os << "Increment multiplier = " << mOptions.mInc << "\n";
            os << "Decrement multiplier = " << mOptions.mDec << "\n";
            os << "Upper bound on the step = " << mOptions.mHUB << "\n";
            os << "Lower bound on the step = " << mOptions.mHLB << "\n";
            os << "Gradient step = " << mOptions.mGradStep << "\n";
            os << "Gradient speedup = " << mOptions.mGradSpeedup << "\n";
            os << "Gradient max steps = " << mOptions.mGradMaxSteps << "\n";
            return os.str();
        }

        /**
         * Retrieve options
         * @return options
         */
        Options & getOptions() {
            return mOptions;
        }

        /**
         * Retrieve stoppers vector reference
         * @return stoppers vector reference
         */
        std::vector<Stopper>& getStoppers() {
            return mStoppers;
        }

        /**
         * Get watchers' vector
         * @return watchers vector
         */
        std::vector<Watcher>& getWatchers() {
            return mWatchers;
        }

    private:

        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        std::vector<Stopper> mStoppers;
        std::vector<Watcher> mWatchers;

    };
}

#endif /* GFSDESC_HPP */

