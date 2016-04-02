/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef VARCOORDESC_HPP
#define  VARCOORDESC_HPP

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
     * Coordinate descent for box constrained problems with variing steps along directions
     */
    template <typename FT> class VarCoorDesc : public LocalSearch <FT> {
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
             * @param gran current granularity
             * @param fval function value
             * @param n current step number 
             */
            virtual bool stopnow(FT xdiff, FT fdiff, const FT* gran, FT fval, int n) = 0;
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
        VarCoorDesc(const COMPI::MPProblem<FT>& prob, Stopper& stopper) :
        mProblem(prob),
        mStopper(stopper) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
            int n = prob.mBox->mDim;
            mSft = new FT [n];
        }

        ~VarCoorDesc() {
            delete [] mSft;
        }

        /**
         * Perform search
         * @param x start point and result
         * @param v  the resulting value
         * @return true if search converged and false otherwise
         */
        bool search(FT* x, FT& v) {
            bool rv = false;

            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT fcur = obj->func(x);
            int n = mProblem.mVarTypes.size();
            for (int i = 0; i < n; i++) {
                mSft[i] = mOptions.mHInit;
            }
            const snowgoose::Box<double>& box = *(mProblem.mBox);

            // One step
            auto step = [&] () {
                FT xd = 0;
                for (int i = 0; i < n; i++) {
                    FT h = mSft[i];
                    FT y = x[i] - h;
                    if (y < box.mA[i]) {
                        y = box.mA[i];
                    }
                    FT tmp = x[i];
                    x[i] = y;
                    FT fn = obj->func(x);
                    if (fn >= fcur) {
                        x[i] = tmp;
                    } else {
                        FT t = h * mOptions.mInc;                        
                        t = SGMIN(t, mOptions.mHUB);
                        mSft[i] = t;
                        fcur = fn;
                        xd += h * h;
                        continue;
                    }

                    y = x[i] + h;
                    if (y > box.mB[i]) {
                        y = box.mB[i];
                    }
                    tmp = x[i];
                    x[i] = y;
                    fn = obj->func(x);
                    if (fn >= fcur) {
                        FT t = h * mOptions.mDec;                        
                        t = SGMAX(t, mOptions.mHLB);
                        mSft[i] = t;
                        x[i] = tmp;
                    } else {
                        FT t = h * mOptions.mInc;                        
                        t = SGMIN(t, mOptions.mHUB);
                        mSft[i] = t;
                        fcur = fn;
                        xd += h * h;
                        continue;
                    }
                }
                return xd;
            };

            int sn = 0;
            for (;;) {
                sn++;
                FT fold = fcur;
                FT xdiff = step();
                FT fdiff = fold - fcur;
                if (fcur < fold) {
                    rv = true;
                } else {
                    FT H = snowgoose::VecUtils::maxAbs(n, mSft, nullptr);
                    if (H <= mOptions.mHLB)
                        break;                    
                }
                if (mStopper.stopnow(xdiff, fdiff, mSft, fcur, sn)) {
                    break;
                }

            }
            v = fcur;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Coordinate descent method with variable adaptation\n";
            os << "Initial step = " << mOptions.mHInit << "\n";
            os << "Increment multiplier = " << mOptions.mInc << "\n";
            os << "Decrement multiplier = " << mOptions.mDec << "\n";
            os << "Upper bound on the step = " << mOptions.mHUB << "\n";
            os << "Lower bound on the step = " << mOptions.mHLB << "\n";
            return os.str();
        }

        /**
         * Retrieve options
         * @return options
         */
        Options & getOptions() {
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
        FT* mSft;
    };
}

#endif /* GFSDESC_HPP */

