/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   locsearch.hpp
 * Author: mikhail
 *
 * Created on 19 февраля 2016 г., 13:30
 */

#ifndef LOCSEARCH_HPP
#define LOCSEARCH_HPP

#include <string>

namespace LOCSEARCH {

    /**
     * Abstract class for local search methods
     */
    template <class FT> class LocalSearch {
    public:

        /**
         * Searches for the optimial value
         * @param x starting point at the beginning and the result at the end
         * @param v on the exit - found value
         * @return true if the search was successul, false if 
         *         no improvement was obtained
         */
        virtual bool search(FT* x, FT& v) = 0;

        /**
         * Information about the  solver
         * @return information
         */
        virtual std::string about() const {
            return "No information about solver";
        }
    };


}


#endif /* LOCSEARCH_HPP */

