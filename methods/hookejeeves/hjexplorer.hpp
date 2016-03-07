/* 
 * File:   hjexplorer.hpp
 * Author: medved
 *
 * Created on March 5, 2016, 3:23 PM
 */

#ifndef HJEXPLORER_HPP
#define	HJEXPLORER_HPP

#include <string>

namespace LOCSEARCH {

    /**
     * Local exploration of the vicinity to be used in Hooke-Jeeves Method
     */
    template <class FT> class HJExplorer {
    public:

        /**
         * Explore the vicinity of the 'x' point to find better value
         * @param x start vector on entry, resulting vector on exit 
         * @return new obtained value 
         */
        virtual FT explore(FT* x) = 0;

        /**
         * The text description of the method
         * @return string description
         */
        virtual std::string about() const {
            return "Unnamed exporer";
        }
    };

};


#endif	/* HJEXPLORER_HPP */

