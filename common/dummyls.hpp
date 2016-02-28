/* 
 * File:   dummyls.hpp
 * Author: medved
 *
 * Created on February 27, 2016, 10:15 AM
 */

#ifndef DUMMYLS_HPP
#define	DUMMYLS_HPP

#include <common/sgerrcheck.hpp>
#include "lineseach.hpp"

namespace LOCSEARCH {

    /**
     * Performs line search in a give segment
     */
    template <class FT> class DummyLineSearch : public LineSearch <FT> {
    public:

        /**
         * Constructor
         * @param reperr true if we want search call reports an error
         */
        DummyLineSearch(bool reperr = false) : mReportError(reperr) {
        }

        bool search(const FT* d, FT* x, FT& v) {
            if(mReportError){
                SG_ERROR_REPORT("Search method in LineSearch does nothing.")
            }
            return true;
        }

    private:
        bool mReportError;
    };
}

#endif	/* DUMMYLS_HPP */

