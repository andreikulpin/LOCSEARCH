/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef GFSDESC_HPP
#define  GFSDESC_HPP

#include <sstream>
#include <common/locsearch.hpp>
#include <common/lineseach.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>


namespace LOCSEARCH {

    /**
     * Modified gradient descent for gradient free box constrained optimization
     */
    template <typename FT> class GFSDesc : public LocalSearch <FT> {
    public:

        /**
         * Determines stopping conditions
         */
        class Stopper {
        public:
            /**
             * Returns true when the search should stop
             * @param xdiff difference between old and new x
             * @param fdiff difference between old and new f value
             * @param gmin gradient minimal component
             * @param fval function value
             * @param n current step number 
             */
            virtual bool stopnow(FT xdiff, FT fdiff, FT gmin, FT fval, int n) = 0;
        };

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
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        GFSDesc(const COMPI::MPProblem<FT>& prob, Stopper& stopper, LineSearch<FT>* ls = nullptr) :
        mProblem(prob),
        mStopper(stopper),
        mLS(ls) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
        }

        /**
         * Perform search
         * @param x start point and result
         * @param v  the resulting value
         * @return true if search converged and false otherwise
         */
        bool search(FT* x, FT& v) {

            bool rv;

            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT uold = obj->func(x);
            int n = mProblem.mVarTypes.size();
            FT sft[n];
            FT g[n];
            FT h = mOptions.mHInit;
            FT gmin;
            int imin;
            const snowgoose::Box<double>& box = *(mProblem.mBox);

            // Computes pseudo gradient
            auto getg = [&] {
                FT fv = obj->func(x);
                FT sum = 0;
                gmin = 0;
                imin = -1;
                for (;;) {
                    for (int i = 0; i < n; i++) {
                        FT sf, sb;
                        FT y = x[i] - h;
                        if (y < box.mA[i]) {
                            y = box.mA[i];
                        }
                        sb = y - x[i];
                        FT tmp = x[i];
                        x[i] = y;
                        FT fb = obj->func(x);
                        x[i] = tmp;

                        y = x[i] + h;
                        if (y > box.mB[i]) {
                            y = box.mB[i];
                        }
                        sf = y - x[i];
                        tmp = x[i];
                        x[i] = y;
                        FT ff = obj->func(x);
                        x[i] = tmp;

                        FT db = fb - fv;
                        FT df = ff - fv;
                        FT gi = 0;
                        if ((df >= 0) && (db >= 0)) {
                            g[i] = 0;
                            sft[i] = 0;
                        } else if (df < db) {
                            sft[i] = sf;
                            g[i] = df;
                            gi = df;
                        } else {
                            sft[i] = sb;
                            g[i] = -db;
                            gi = db;
                        }
                        if (gi < gmin) {
                            imin = i;
                            gmin = gi;
                        }
                        sum += SGABS(gi);
                    }
                    if (imin == -1) {
                        if (h > mOptions.mHLB)
                            h *= mOptions.mDec;
                        else break;
                    } else
                        break;
                }
                return sum;
            };



            int k = 0;
            FT xk[n];



            for (;;) {
                k++;
                double s = getg();
                if (imin == -1) {
                    rv = true;
                    break;
                }
                FT u;



                if (mLS == nullptr) {
                    snowgoose::VecUtils::vecSaxpy(n, x, g, -h / s, xk);
                    //VecUtils::vecSaxpy(n, x, g, h * 1. / gmin, xk);
                    //VecUtils::vecSaxpy(n, x, g, 1. / gmin, xk);
                    snowgoose::BoxUtils::project(xk, *(mProblem.mBox));
                } else {
                    FT rg[n];
                    FT vv;
                    snowgoose::VecUtils::vecCopy(n, x, xk);
                    snowgoose::VecUtils::vecMult(n, g, -h / s, rg);
                    mLS->search(rg, xk, vv);
                }
                u = obj->func(xk);

                if (u > uold + gmin) {
                    snowgoose::VecUtils::vecCopy(n, x, xk);
                    xk[imin] = x[imin] + sft[imin];
                    u = uold + gmin;
                }
                FT xdiff = snowgoose::VecUtils::vecDist(n, x, xk);
                FT fdiff = u - uold;
                if (mStopper.stopnow(xdiff, fdiff, gmin, uold, k)) {
                    rv = true;
                    if (u < uold) {
                        snowgoose::VecUtils::vecCopy(n, xk, x);
                        uold = u;
                    }
                    break;
                }
                if (u >= uold) {
                    if (h > mOptions.mHLB)
                        h *= mOptions.mDec;
                } else {
                    if (h < mOptions.mHUB)
                        h *= mOptions.mInc;
                    snowgoose::VecUtils::vecCopy(n, xk, x);
                    uold = u;
                }
            }
            v = uold;

            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Gradient free descent method\n";
            os << "Initial step = " << mOptions.mHInit << "\n";
            os << "Increment multiplier = " << mOptions.mInc << "\n";
            os << "Decrement multiplier = " << mOptions.mDec << "\n";
            os << "Upper bound on the step = " << mOptions.mHUB << "\n";
            os << "Lower bound on the step = " << mOptions.mHLB << "\n";
            if(mLS != nullptr)
                os << "Line search: " << mLS->about() << "\n";
            return os.str();
        }

        /**
         * Retrieve options
         * @return options
         */
        Options& getOptions() {
            return mOptions;
        }

    private:

        /*
        void project(FT* x) {
            int n = LocalOptimizer<FT>::getObjective()->getDim();
            for (int i = 0; i < n; i++) {
                if (x[i] < LocalBoxOptimizer<FT>::mBox.mA[i])
                    x[i] = LocalBoxOptimizer<FT>::mBox.mA[i];
                else if (x[i] > LocalBoxOptimizer<FT>::mBox.mB[i])
                    x[i] = LocalBoxOptimizer<FT>::mBox.mB[i];
            }
        }
         */

        Stopper &mStopper;
        const COMPI::MPProblem<FT>& mProblem;
        LineSearch<FT>* mLS;
        Options mOptions;
    };
}

#endif /* GFSDESC_HPP */

