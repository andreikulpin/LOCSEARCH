/* 
 * File:   lineseach.hpp
 * Author: medved
 *
 * Created on February 26, 2016, 5:54 PM
 */

#ifndef LINESEACH_HPP
#define LINESEACH_HPP

namespace LOCSEARCH {

    /**
     * Performs line search in a give segment
     */
    template <class FT> class LineSearch {
    public:

        /**
         * Search along a direction
         * @param d direction
         * @param x starting point and the result
         * @param v resulting value
         * @return true if success
         */
        virtual bool search(const FT* d, FT* x, FT& v) = 0;

        /**
         * Information about the solver
         * @return information
         */
        virtual std::string about() const {
            return "Unknown line search";
        }
    };
}

#endif /* LINESEACH_HPP */

