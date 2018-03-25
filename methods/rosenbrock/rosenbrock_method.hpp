/* 
 * File:   rosenbrock_method.hpp
 * Author: kulpin
 */

#ifndef ROSENBROCK_METHOD_HPP
#define  ROSENBROCK_METHOD_HPP

#include <sstream>
#include <vector>
#include <functional>
#include <memory>
#include <solver.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
#include <common/lineseach.hpp>


namespace LOCSEARCH {
    /**
     * Rosenbrock method implementation
     */
    template <typename FT> class RosenbrockMethod : public COMPI::Solver<FT> {
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
         * Options for Rosenbrock method
         */
        struct Options {
            /**
             * Initial value of granularity
             */
            FT* mHInit = nullptr;
            /**
             * Increase in the case of success
             */
            FT mInc = 2;
            /**
             * Decrease in the case of failure
             */
            FT mDec =  - 0.5;
            /**
             * Lower bound for granularity
             */
            FT mHLB = - 1e+02;
            /**
             * Upper bound on granularity
             */
            FT mHUB = 1e+02;
            /**
             * Trace on/off
             */
            bool mDoTracing = false;
            /**
             * Max unsuccessful steps number
             */
            FT mEps = 0.6;
            /**
             * Max unsuccessful steps number
             */
            int maxUnsuccessStepsNumber = 10;
            /**
             * Max steps number
             */
            int maxStepsNumber = 100;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        RosenbrockMethod(const COMPI::MPProblem<FT>& prob) :
        mProblem(prob) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
        }

        /**
         * Retrieve the pointer to the line search
         * @return 
         */
        std::unique_ptr<LineSearch<FT>>&getLineSearch() {
            return mLS;
        }

        /**
         * Perform search
         * @param x start point and result
         * @param v  the resulting value
         * @return true if search converged and false otherwise
         */
        bool search(FT* x, FT& v) override {
            bool rv = false;
            auto obj = mProblem.mObjectives.at(0);
            int n = mProblem.mVarTypes.size();
            int n_sqr = n * n;
            
            FT fcur = obj->func(x);
            
            FT xOld[n]; 
            snowgoose::VecUtils::vecCopy(n, x, xOld);
            FT fOld = fcur; 
            
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            
            std::vector<FT> sft(mOptions.mHInit, mOptions.mHInit + n);
            FT stepLenghts[n];
            snowgoose::VecUtils::vecSet(n, 0., stepLenghts);

            FT * dirs = new FT[n_sqr];
            snowgoose::VecUtils::vecSet(n * n, 0., dirs);
            for (int i = 0; i < n; i++) {
                dirs[i * n + i] = 1;
            }
           
            FT * a = new FT[n_sqr];
            FT * b = new FT[n_sqr];
            FT * d = new FT[n_sqr];
            
            int unsuccessSteps = 0;
            int sterNum = 1;
            bool br = false;
            
            if (mOptions.mDoTracing) {
                printArray("x", n, x);
                printVector("sft", n, sft);
                printArray("stepLenghts", n, stepLenghts);
                printMatrix("dirs", n, n, dirs);
            }

            auto inc = [this] (FT h) {
                FT t = h;
                t = h * mOptions.mInc;
                t = SGMIN(t, mOptions.mHUB);
                return t;
            };

            auto dec = [this](FT h) {
                FT t = h;
                t = h * mOptions.mDec;
                t = SGMAX(t, mOptions.mHLB);
                return t;
            };

            /*
             * Attepmt yielding new minimum along each base direction.
             * @return true if step along at least one direction was successful
             */
            auto step = [&] () {
                if (mOptions.mDoTracing) {
                    std::cout << "\n*** Step " << sterNum << " ***\n";
                }
                bool isStepSuccessful = false;
                FT xn[n]; 
                snowgoose::VecUtils::vecCopy(n, x, xn);
                FT fn = fcur;

                for (int i = 0; i < n; i++) {
                    const FT h = sft[i];
                    FT xtmp[n]; 
                    snowgoose::VecUtils::vecSaxpy(n, xn, &(dirs[i * n]), h, xtmp);
                    FT ftmp = obj->func(xtmp);
                    
                    if (mOptions.mDoTracing) {
                        printArray("xtmp", n, xtmp);
                        // std::cout << "ftmp = " << ftmp << "\n";
                    }

                    if (ftmp < fn) {
                        isStepSuccessful = true;
                        stepLenghts[i] += h;
                        sft[i] = inc(h);
                        snowgoose::VecUtils::vecCopy(n, xtmp, xn);
                        fcur = ftmp;

                    } else {
                        sft[i] = dec(h);
                    }
                }

                snowgoose::VecUtils::vecCopy(n, xn, x);
                return isStepSuccessful;
            };

            auto ortogonalize = [&] () {
                if (mOptions.mDoTracing) {
                    std::cout << "\n** Ortogonalize **" << "\n";
                }
                
                snowgoose::VecUtils::vecSet(n * n, 0., a);

                for (int i = 0; i < n; i++) {
                    if (stepLenghts[i] == 0) {
                        snowgoose::VecUtils::vecCopy(n, &(dirs[i * n]), &(a[i * n]));
                    
                    } else {
                        for (int j = i; j < n; j++) {
                            // printArray("dirs ", n, dirs[j]);
                            snowgoose::VecUtils::vecSaxpy(n, &(a[i * n]), &(dirs[j * n]), stepLenghts[j], &(a[i * n]));
                            // printArray("a ", n, a[i]);
                        }
                    }

                    snowgoose::VecUtils::vecCopy(n, &(a[i * n]), &(b[i * n]));
                    
                    for (int j = 0; j < i; j++) {
                        FT scalarMlp = snowgoose::VecUtils::vecScalarMult(n, &(a[i * n]), &(d[j * n]));
                        snowgoose::VecUtils::vecSaxpy(n, &(b[i * n]), &(d[j * n]), - scalarMlp, &(b[i * n]));
                    }
                    
                    FT norm = snowgoose::VecUtils::vecNormTwo(n, &(b[i * n]));
                    snowgoose::VecUtils::vecMult(n, &(b[i * n]), 1 / norm, &(d[i * n]));
                }

                for (int i = 0; i < n; i++) {
                   snowgoose::VecUtils::vecCopy(n, &(d[i * n]), &(dirs[i * n]));
                }
                
                if (mOptions.mDoTracing) {
                    printMatrix("a", n, n, a);
                    printMatrix("b", n, n, b);
                    printMatrix("d", n, n, d);
                    std::cout << "*************" << "\n";
                }
            };
            
            while (!br) {
                bool success = step();
                unsuccessSteps = success ? 0 : unsuccessSteps + 1;
                sterNum++;
                
                if (mOptions.mDoTracing) {
                    std::cout << (success ? "Success" : "Not success") << "\n";
                    printArray("x", n, x);
                    printVector("sft", n, sft);
                }

                if (!success) {
                    if (fcur < fOld || unsuccessSteps >= mOptions.maxUnsuccessStepsNumber) {
                        FT dist = snowgoose::VecUtils::vecDist(n, xOld, x); 
                        
                        if (dist > mOptions.mEps) {
                            fOld = fcur;
                            snowgoose::VecUtils::vecCopy(n, x, xOld);
                            ortogonalize();
                            sft.assign(mOptions.mHInit, mOptions.mHInit + n);
                            snowgoose::VecUtils::vecSet(n, 0., stepLenghts);
                            
                        } else {
                            br = true;
                            
                            if (unsuccessSteps >= mOptions.maxUnsuccessStepsNumber) {
                                std::cout << "Stopped as number of unsuccessful steps reached its limit of " 
                                        << mOptions.maxUnsuccessStepsNumber << "\n";
                                
                            } else {
                                std::cout << "Stopped as last jump length was less than eps" << "\n";
                                std::cout << "Eps = " << mOptions.mEps << " Dist = " << dist << "\n";
                            }
                        }
                    
                    } else {
                        br = true;
                        for (int i = 0; i < n; i++) {
                            br = SGABS(sft[i]) <= mOptions.mEps;
                        }
                        
                        if (br) {
                            std::cout << "Stopped as all step lengths was less than eps\n";
                        }
                    }
                }
                
                if (sterNum >= mOptions.maxStepsNumber) {
                    br = true;
                    std::cout << "Stopped as number of steps was too big\n";
                }

                for (auto w : mWatchers) {
                    w(fcur, x, sft, sterNum);
                }
                for (auto s : mStoppers) {
                    if (s(fcur, x, sterNum)) {
                        br = true;
                        break;
                    }
                }
            }
            v = fcur;

            delete [] dirs;
            delete [] a;
            delete [] b;
            delete [] d;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
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
        std::unique_ptr<LineSearch<FT>> mLS;
        
        void printMatrix(const char * name, int n, int m, FT * matrix) {
            std::cout << name << " =\n";
            for (int i = 0; i < n; i++) {
                std::cout << snowgoose::VecUtils::vecPrint(m, &(matrix[i + n])) << "\n";
            }
        }
        
        void printArray(const char * name, int n, FT * array) {
            std::cout << name << " = ";
            std::cout << snowgoose::VecUtils::vecPrint(n, array) << "\n";
        }
        
        void printVector(const char * name, int n, std::vector<FT> vector) {
            std::cout << name << " = ";
            std::cout << "[ ";
            for (int i = 0; i < n; i++) {
                std::cout << vector[i] << ", ";
            }
            std::cout << " ]" << "\n";
        }
    };
}

#endif 

