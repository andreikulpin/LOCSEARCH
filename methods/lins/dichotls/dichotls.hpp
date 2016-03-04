/* 
 * File:   dichotls.hpp
 * Author: medved
 *
 * Created on February 26, 2016, 6:45 PM
 */

#ifndef DICHOTLS_HPP
#define DICHOTLS_HPP

#include <common/lineseach.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
#include <common/vec.hpp>

namespace LOCSEARCH {

    /**
     * Dichotomic line search
     */
    template <class FT> class DichotLS : public LineSearch<FT> {
    public:

        /**
         * Determines stopping conditions
         */
        class Stopper {
        public:
            /**
             * Stop iterations when return true
             * @param s current step length
             * @param k number of moves done
             * @param vold previous value of objective
             * @param vnew new value of objective
             * @return true to stop iterations, false otherwise
             */
            virtual bool stopnow(FT s, int k, FT vold, FT vnew) = 0;
        };

        struct Options {
            /**
             * Initial step
             */
            FT mSInit = 0.1;

            /**
             * Accelerate search when objective is improved
             */
            FT mAccelerate = 1.1;

            /**
             * Decrease step when objective get worse
             */
            FT mSlowDown = 0.5;
        };

        /**
         * The constructor
         *
         * @param box bounding box
         */
        DichotLS(const COMPI::MPProblem<FT>& prob, Stopper& stopper) :
        mProblem(prob),
        mStopper(stopper) {
        }

        bool search(const FT* d, FT* x, FT& v) {
            FT s = mOptions.mSInit;
            int n = mProblem.mVarTypes.size();
            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            FT xk[n];
            v = obj->func(x);
            FT vn = v;
            FT vo = v;
            int k = 0;
            bool rv = false;
            FT S = 0;
            auto trystep = [&] () {
                snowgoose::VecUtils::vecSaxpy(n, x, d, S + s, xk);
                if (COMPI::MPUtils::isFeasible(mProblem, xk)) {
                    vn = obj->func(xk);
                    if (vn < v)
                        return true;
                    else
                        return false;
                } else
                    return false;
            };

            for (;; k++) {
                if (mStopper.stopnow(s, k, vo, vn))
                    break;
                vo = v;
                if (!trystep()) {
                    s *= mOptions.mSlowDown;
                    continue;
                } else {
                    rv = true;
                    S += s;
                    s *= mOptions.mAccelerate;
                    v = vn;
                }
            }
            if (rv)
                snowgoose::VecUtils::vecCopy(n, xk, x);
            return rv;
        }

        std::string about() const {
            return "Dichotomic line search";
        }

        /**
         * Retrieve options reference
         * @return options reference
         */
        Options& getOptions() {
            return mOptions;
        }

    private:

        Stopper &mStopper;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
    };
}

#endif /* DICHOTLS_HPP */

