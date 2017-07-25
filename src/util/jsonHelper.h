/*!
 * \file jsonHelper.h
 *
 * \brief Utilities for reading json file information in the environment
 *
 * \date Started: 2017
 *
 * \authors Mahdi Soltan Mohammadi, Michelle Strout
 *
 * Copyright (c) 2017, University of Arizona <br>
 * All rights reserved. <br>
 * See ../../COPYING for details. <br>
 */

#ifndef JSONHELPER_H_
#define JSONHELPER_H_


#include <iostream>
#include <iegenlib.h>
#include <parser/jsoncons/json.hpp>

using jsoncons::json;
using namespace iegenlib;

// Reads a list of UFCs from a json structure and stores them in the environment  
void addUFCs(json &ufcs);


// Reads a list of universially quantified constraints from a json structure
// and stores them in the environment
void adduniQuantRules(json &uqCons);

// Reads iterators that we should not project from a json sructure
void notProjectIters(Relation* rel, std::set<int> &parallelTvs, json &np);

#endif
