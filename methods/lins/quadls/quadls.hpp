/* 
 * File:   dichotls.hpp
 * Author: medved
 *
 * Created on February 26, 2016, 6:45 PM
 */

#ifndef QUADLS_HPP
#define	QUADLS_HPP

#include <common/lineseach.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
#include <common/vec.hpp>

namespace LOCSEARCH {

    /**
     * Line search based on quadratic interpolation
     */
    template <class FT> class QuadLS : public LineSearch<FT> {
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
        QuadLS(const COMPI::MPProblem<FT>& prob, Stopper& stopper) :
        mProblem(prob),
        mStopper(stopper) {
        }

        bool search(const FT* d, FT* x, FT& v) {
            FT s = mOptions.mSInit;
            const int n = mProblem.mVarTypes.size();
            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            FT xk[n];
            FT vo = obj->func(x);
            v = vo;
            FT vn = vo;
            FT S = 0;
            int k = 0;
            bool rv = false;
            const int hsize = 3;
            FT sh[hsize];
            FT fh[hsize];

            int fill = 0;
            auto pushhist = [&] (FT sv, FT fv) {
                if (fill == hsize) {
                    for (int i = 0; i < hsize - 1; i++) {
                        sh[i] = sh[i + 1];
                        fh[i] = fh[i + 1];
                    }
                    sh[hsize - 1] = sv;
                    fh[hsize - 1] = fv;
                } else {
                    sh[fill] = sv;
                    fh[fill] = fv;
                    fill++;
                }
            };
            auto compsq = [&] () {
                FT a, b, c, ss;
                if (fill == hsize) {
                    a = (1. / (sh[0] - sh[2])) * ((fh[0] - fh[1]) / (sh[0] - sh[1]) - (fh[1] - fh[2]) / (sh[1] - sh[2]));
                    b = (fh[0] - fh[1]) / (sh[0] - sh[1]) - a * (sh[0] + sh[1]);
                    c = fh[0] - a * sh[0] * sh [0] - b * sh[0];
                    if (a > 0) {
                        ss = -b / (2. * a);
                    } else {
                        ss = -1.;
                    }
                } else {
                    ss = -1.;
                }
                return ss;
            };
            auto trystep = [&] (FT stp) {
                if (stp <= 0)
                    return false;
                //snowgoose::VecUtils::vecSaxpy(n, x, d, stp, xk);
                if (COMPI::MPUtils::isFeasible(mProblem, xk)) {
                    vn = obj->func(xk);
                    if (vn < v)
                        return true;
                    else
                        return false;
                } else
                    return false;
            };
            pushhist(S, vo);
            for (;; k++) {
                if (mStopper.stopnow(s, k, vo, vn))
                    break;
                vo = v;
                FT sq = compsq();
                bool succ = trystep(sq);
                if (succ) {
                    S = sq;
                } else {
                    succ = trystep(S + s);
                    if (succ) {
                        S += s;
                        s *= mOptions.mAccelerate;
                    } else {
                        s *= mOptions.mSlowDown;
                    }
                }
                if (succ) {
                    pushhist(S, vn);
                    rv = true;
                    v = vn;
                }
            }
            if (rv) {
                snowgoose::VecUtils::vecSaxpy(n, x, d, S, x);
            }
            return rv;
        }

        
        std::string about() const {
            return "Quadratic approximation grad free line search";
        }
        
        /**
         * Retrieve options reference
         * @return options reference
         */
        Options & getOptions() {
            return mOptions;
        }

    private:

        Stopper &mStopper;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
    };
}

#endif	/* QUADLS_HPP */

