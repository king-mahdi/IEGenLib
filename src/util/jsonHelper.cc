/*!
 * \file jsonHelper.cc
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

#include "jsonHelper.h"

//using jsoncons::json;
//using namespace iegenlib;

// Reads a list of UFCs from a json structure and stores them in the environment  
void addUFCs(json &ufcs){

  for (size_t j = 0; j < ufcs["UFS"].size(); ++j){

    bool bijective = false;
    if( ufcs["UFS"][j]["Bijective"].as<string>() == string("true") ){
      bijective = true;
    }
    iegenlib::MonotonicType monotonicity = iegenlib::Monotonic_NONE;
    if(ufcs["UFS"][j]["Monotonicity"].as<string>() == 
                            string("Monotonic_Nondecreasing")){
      monotonicity = iegenlib::Monotonic_Nondecreasing;
    } else if(ufcs["UFS"][j]["Monotonicity"].as<string>() == 
                                  string("Monotonic_Increasing")){
      monotonicity = iegenlib::Monotonic_Increasing;
    }

    iegenlib::appendCurrEnv(ufcs["UFS"][j]["Name"].as<string>(),// Name
         new Set(ufcs["UFS"][j]["Domain"].as<string>()),   // Domain 
         new Set(ufcs["UFS"][j]["Range"].as<string>()),    // Range
         bijective,                                              // Bijective?
         monotonicity                                            // Monotonicity?
                                );
  }
}


// Reads a list of universially quantified constraints from a json structure
// and stores them in the environment
void adduniQuantRules(json &uqCons){

  UniQuantRule *uqRule;
  for (size_t j = 0; j < uqCons.size(); ++j){

    // forall e1, e2, ... : p => q
    delete uqRule;
    uqRule = new UniQuantRule(uqCons[j]["Type"].as<string>(), 
                uqCons[j]["UniQuantVar"].as<string>(), 
                uqCons[j]["p"].as<string>(), uqCons[j]["q"].as<string>());
    currentEnv.addUniQuantRule( uqRule );
  }
}

// Reads iterators that we should not project from a json sructure
void notProjectIters(Relation* rel, std::set<int> &parallelTvs, json &np){

  iegenlib::TupleDecl td = rel->getTupleDecl();
  for (size_t j = 0; j < np.size(); ++j){
    string tvS = np[j].as<string>();
    int tvN = -1;
    for (unsigned int c = 0 ; c < td.getSize() ; c++){
      if( tvS == td.elemToString(c) ){
        tvN = c;
        break;
      }
    }
    parallelTvs.insert( tvN );
  }
}

