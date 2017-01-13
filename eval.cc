/*!
 * \file simplifyForPartialParallelDriver.cc
 *
 * This file is a driver for using simplifyForPartialParallel 
 * and pertaining functions in IEgenLib.
 * These functions can be used to carry out simplification of 
 * data dependence relation in non-affine relations containing
 * uninterpreted function calls.
 *
 * To compile, after building IEgenLib with:

       ./configure    
       make

    Run the following command with YOUR OWN ADDRESSES:

        g++ -o EXECUTABLE simplifyDriver.cc 
            -I IEGENLIB_HOME/src IEGENLIB_HOME/build/src/libiegenlib.a -lisl -std=c++11

    You can also run following command after running "make install" for IEgenLIB:

       g++ -o EXECUTABLE simplifyDriver.cc 
           -I INSTALLATION_FOLDER/include/iegenlib cpp_api_example.cc
              INSTALLATION_FOLDER/lib/libiegenlib.a -lisl -std=c++11
   
    IEGENLIB_HOME indicates where you have your copy of IEgenLIB.
    For example if you are compiling this file in its original location that is 
    IEGENLIB_HOME, then you can run following to compile:

    g++ -o eval eval.cc -I src build/src/libiegenlib.a -lisl -std=c++11

 * Now to run the driver, you should put your dependence relations inside
 * JSON files and give them as inputs to the driver, one or more files at a time.
 * JSON file format is demonstrated by examples gs_csr.json and ilu_csr.json that
 * can be found in same directory as this driver. The json input files include commnet
 * fields that have the application code (Gauss-Seidel and ILU) in them, and says where
 * in the code each dependence relation is getting extract. So you can run following
 * afetr compiling the driver to get the simplified relations:
 
   ./eval eval.list 

 *     
 *
 * \date Date Started: 3/21/2016
 *
 * \authors Michelle Strout, Mahdi Soltan Mohammadi
 *
 * Copyright (c) 2016, University of Arizona <br>
 * All rights reserved. <br>
 * See COPYING for details. <br>
 * 
 *
 *
 */


#include <iostream>
#include <fstream>
#include "iegenlib.h"
#include "parser/jsoncons/json.hpp"

using jsoncons::json;
using iegenlib::Set;
using iegenlib::Relation;
using namespace std;

std::ofstream rout("result.txt", std::ofstream::out);
std::ofstream rnout("result_notes.txt", std::ofstream::out);

#define n_combo 10
enum methods {
  Super_Affine_Set,
  Approximation,
  Domain_Specific,
  Monotonicity,
  Projection,
  Simplification
};
enum result {
  Unsat,
  Projected
};

string method_name(int m)
{
  if( m == Approximation)  return string("Approximation");
  else if( m == Domain_Specific)  return string("Domain_Specific");
  else if( m == Monotonicity)  return string("Monotonicity");
  else if( m == Super_Affine_Set)  return string("Super_Affine_Set");
  else if( m == Projection)  return string("Projection");
  else if( m == Simplification)  return string("Simplification");

  return string("Not Specified");
}

int res_tot[n_combo][2];
int n_t_conjs;
int n_t_iters;

void simplify(string inputFile);

// Utility function
bool printRelation(string msg, Relation *rel);
void EXPECT_EQ(string a, string b);
void EXPECT_EQ(Relation *a, Relation *b);
int str2int(string str);


//----------------------- MAIN ---------------
int main(int argc, char **argv)
{
  std::string jname; 
  if (argc != 2){
    cout<<"\n\nYou need to specify the input file containing name of input JSON files:"
          "\n./simplifyDriver eval.list\n\n";
  } else {
     
    ifstream in(argv[1]);
    
    while(in>>jname){
      simplify(string( jname ));
    }
  }

  rout<<"\n---- Final results:"<<"\n\n";
  rout<<"Number of cunjunctions: "<<n_t_conjs<<"\n";
  rout<<"Number of iterators: "<<n_t_iters<<"\n\n";
  rout<<"Name of method     , UnSat Conj., % of Total, Projected TVs, % of Total"<<"\n";
  rout<<"Super_Affine_Set   , "<<res_tot[Super_Affine_Set][Unsat]<<", "
      <<(int)(((double)res_tot[Super_Affine_Set][Unsat]/(double)n_t_conjs)*100)<<", "
      <<"NA, " //<<res_tot[Super_Affine_Set][Projected]<<", "
      <<"NA, " //<<(int)(((double)res_tot[Super_Affine_Set][Projected]/(double)n_t_iters)*100)
      <<"\n";
  rout<<"Domain_Specific    , "<<res_tot[Domain_Specific][Unsat]<<", "
      <<(int)(((double)res_tot[Domain_Specific][Unsat]/(double)n_t_conjs)*100)<<", "
      <<"NA, " //<<res_tot[Domain_Specific][Projected]<<", "
      <<"NA, " //<<(int)(((double)res_tot[Domain_Specific][Projected]/(double)n_t_iters)*100)
      <<"\n";
  rout<<"Monotonicity       , "<<res_tot[Monotonicity][Unsat]<<", "
      <<(int)(((double)res_tot[Monotonicity][Unsat]/(double)n_t_conjs)*100)<<", "
      <<"NA, " //<<res_tot[Monotonicity][Projected]<<", "
      <<"NA, " //<<(int)(((double)res_tot[Monotonicity][Projected]/(double)n_t_iters)*100)
      <<"\n";
  rout<<"Projection         , "<<res_tot[Projection][Unsat]<<", "
      <<(int)(((double)res_tot[Projection][Unsat]/(double)n_t_conjs)*100)<<", "
      <<res_tot[Projection][Projected]<<", "
      <<(int)(((double)res_tot[Projection][Projected]/(double)n_t_iters)*100)
      <<"\n";
  rout<<"Simplification     , "<<res_tot[Simplification][Unsat]<<", "
      <<(int)(((double)res_tot[Simplification][Unsat]/(double)n_t_conjs)*100)<<", "
      <<res_tot[Simplification][Projected]<<", "
      <<(int)(((double)res_tot[Simplification][Projected]/(double)n_t_iters)*100)
      <<"\n\n\n";

    return 0;
}

// Reads information from a JSON file (inputFile), and applies
// the simplification algorithm to the sets found in the file. 
void simplify(string inputFile)
{
  iegenlib::setCurrEnv();
  std::set<int> parallelTvs;
  // (0)
  // Read the data from inputFile
  ifstream in(inputFile);
  json data;
  in >> data;

  for(size_t p = 0; p < data.size(); ++p){    // Dependence relations (DR) found in the file

    int res_c[n_combo][2] = {0};
    int n_c_conjs = 0;
    int n_c_iters = 0;

    n_c_conjs += data[p].size();
    n_t_conjs += data[p].size();


    rnout<<"\n\n---- "<<data[p][0]["Name"].as<string>()<<":\n\n";

    for (size_t i = 0; i < data[p].size(); ++i){// Conjunctions found for one DR in the file

      bool UNS = false;

      // (1) Introduce the uninterpreted function symbols to environment, and 
      //     indicate their domain, range, whether they are bijective, or monotonic.
      if( i == 0 ){  // Read these data only once. 
                     // They are stored in the first conjunction.
        for (size_t j = 0; j < data[p][i]["UFS"].size(); ++j){
 
          bool bijective = false;
          if( data[p][i]["UFS"][j]["Bijective"].as<string>() == string("true") ){
            bijective = true;
          }
          iegenlib::MonotonicType monotonicity = iegenlib::Monotonic_NONE;
          if(data[p][i]["UFS"][j]["Monotonicity"].as<string>() == 
                                       string("Monotonic_Nondecreasing")){
            monotonicity = iegenlib::Monotonic_Nondecreasing;
          } else if(data[p][i]["UFS"][j]["Monotonicity"].as<string>() == 
                                          string("Monotonic_Increasing")){
            monotonicity = iegenlib::Monotonic_Increasing;
          }

          iegenlib::appendCurrEnv(data[p][i]["UFS"][j]["Name"].as<string>(),// Name
              new Set(data[p][i]["UFS"][j]["Domain"].as<string>()),   // Domain 
              new Set(data[p][i]["UFS"][j]["Range"].as<string>()),    // Range
              bijective,                                              // Bijective?
              monotonicity                                            // Monotonicity?
                                  );
        }  
      }

      // (2) Putting constraints in an iegenlib::Relation
      Relation* rel = new Relation(data[p][i]["Relation"].as<string>());

      Relation* rel_copy = new Relation( rel->inArity() , rel->outArity() );
      *rel_copy = *rel;

      int orig_arity = rel->arity();
      int new_arity = rel->arity();

      n_c_iters += rel->arity();
      n_t_iters += rel->arity();
      
      // (3) Specify loops that are going to be parallelized, 
      //     so we are not going to project them out.
      if( i == 0 ){  // Read these data only once. 
                     // They are stored in the first conjunction.
  
        for (size_t j = 0; j < data[p][0]["Do Not Project Out"].size(); ++j){
          string tvS = data[p][0]["Do Not Project Out"][j].as<string>();
          int tvN = 0;
          iegenlib::TupleDecl td = rel->getTupleDecl();
          for (unsigned int c = 0 ; c < td.getSize() ; c++){
            if( tvS == td.elemToString(c) ){
              tvN = c;
              break;
            }
          }
          parallelTvs.insert( tvN );
        }
      }


      //################# Testing Super_Affine_Set #################
      rel_copy->normalize();
     
      if( rel_copy->isDefault()) {
        res_c[Super_Affine_Set][Unsat]++;
        res_tot[Super_Affine_Set][Unsat]++;
        UNS = true;

        rnout<<"  The relation "<<i<<" is unsat: Super_Affine_Set"<<"\n";
      }
      else {
        //find # of new founded equalities;
        //new_arity = rel_copy->arity();
        //res_c[Super_Affine_Set][Projected] += orig_arity - new_arity;
        //res_c[Super_Affine_Set][Projected] += orig_arity - new_arity;
      } 


      //################# Testing domain specific information #################
      if( ! UNS ){  //Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        *rel_copy = *rel;
        Relation *rel_extend;
        for (size_t j = 0; j < data[p][0]["User Defined"].size(); ++j){
   
          rel_extend = rel_copy->addUFConstraints(
                        data[p][0]["User Defined"][j]["Func1"].as<string>(),
                        data[p][0]["User Defined"][j]["operator"].as<string>(),
                        data[p][0]["User Defined"][j]["Func2"].as<string>()
                                              );
          *rel_copy = *rel_extend;
          delete rel_extend;
        }

        rel_copy->normalize();
        if( rel_copy->isDefault()) {
          res_c[Domain_Specific][Unsat]++;
          res_tot[Domain_Specific][Unsat]++;

          rnout<<"  The relation "<<i<<" is unsat: Domain_Specific"<<"\n";
        }
        else {
          //find # of new founded equalities;
        } 
      }


      //################# Testing Monotonicity #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        *rel_copy = *rel;
        Relation *rel_extend;
        rel_copy->normalize();
        rel_extend = rel_copy->addConstraintsDueToMonotonicity();
        *rel_copy = *rel_extend;
        delete rel_extend;

        rel_copy->normalize();
        if( rel_copy->isDefault()) {

          res_c[Monotonicity][Unsat]++;
          res_tot[Monotonicity][Unsat]++;

          rnout<<"  The relation "<<i<<" is unsat: Monotonicity"<<"\n";
        }
        else {
          //find # of new founded equalities;
        } 
      }


      //################# Testing Projection #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy

        Relation *temp;
        int lastTV = rel->arity()-1;
        *rel_copy = *rel;
        for (int tv = lastTV ; tv >= 0 ; tv-- ) {

          if ( parallelTvs.find(tv) != parallelTvs.end() ){
            continue;
          }

          // Project out if it is not an UFCall argument
          temp = rel_copy->projectOut(tv); 

          if ( temp ){
            delete rel_copy;        
            rel_copy = temp;
          }
          if( rel_copy->isDefault() ){
            UNS = true;
            break;
          }
        }

        if( UNS ) {
          res_c[Projection][Unsat]++;
          res_tot[Projection][Unsat]++;

          rnout<<"  The relation "<<i<<" is unsat: Projection"<<"\n";
        }
        else {
          int pr = rel->arity() - rel_copy->arity();
          res_c[Projection][Projected] += pr;
          res_tot[Projection][Projected] += pr;

          rnout<<"  Projected "<<pr<<" TVs from relation "<<i<<": Projection"<<"\n";
        } 
      }

      //################# Testing simplifyForPartialParallel #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        // Applying heuristic for removing expensive iterators
        int numConstToRemove = str2int(data[p][0]["Remove Constraints"].as<string>());
        rel->RemoveExpensiveConsts(parallelTvs, numConstToRemove );
        // Add user defined constraints
        Relation *rel_extend;
        for (size_t j = 0; j < data[p][0]["User Defined"].size(); ++j){
          rel_extend = rel->addUFConstraints(
                       data[p][0]["User Defined"][j]["Func1"].as<string>(),
                       data[p][0]["User Defined"][j]["operator"].as<string>(),
                       data[p][0]["User Defined"][j]["Func2"].as<string>()
                                            );
          *rel = *rel_extend;
          delete rel_extend;
        }
        // Simplifyng the constraints relation
        if(rel_copy) delete rel_copy;

        //rel->normalize();

        rel_copy = rel->simplifyForPartialParallel(parallelTvs);

        //if( rel_copy ) 
          //rel_copy->normalize();

//if( inputFile == string("ilu_csr.json") && (i == 6 || i == 7) )
  //printRelation(string("Norm 6/7: ") , rel_copy);

        if( !rel_copy ){ //|| rel_copy->isDefault() ) {
          res_c[Simplification][Unsat]++;
          res_tot[Simplification][Unsat]++;

          rnout<<"  The relation "<<i<<" is unsat: Simplification"<<"\n";
        }
        else {
          int pr = rel->arity() - rel_copy->arity();
          res_c[Simplification][Projected] += pr;
          res_tot[Simplification][Projected] += pr;

          rnout<<"  Projected "<<pr<<" TVs from relation "<<i<<": Simplification"<<"\n";
        } 
      }

      delete rel;
    }

    rout<<"---- "<<data[p][0]["Name"].as<string>()<<":\n\n";
    rout<<"Number of cunjunctions: "<<n_c_conjs<<"\n";
    rout<<"Number of iterators: "<<n_c_iters<<"\n\n";
    rout<<"Name of method     , UnSat Conj., % of Total, Projected TVs, % of Total"<<"\n";
    rout<<"Super_Affine_Set   , "<<res_c[Super_Affine_Set][Unsat]<<", "
        <<(int)(((double)res_c[Super_Affine_Set][Unsat]/(double)n_c_conjs)*100)<<", "
        <<"NA, " //<<res_c[Super_Affine_Set][Projected]<<", "
        <<"NA, " //<<(int)(((double)res_c[Super_Affine_Set][Projected]/(double)n_c_iters)*100)
        <<"\n";
    rout<<"Domain_Specific    , "<<res_c[Domain_Specific][Unsat]<<", "
        <<(int)(((double)res_c[Domain_Specific][Unsat]/(double)n_c_conjs)*100)<<", "
        <<"NA, " //<<res_c[Domain_Specific][Projected]<<", "
        <<"NA, " //<<(int)(((double)res_c[Domain_Specific][Projected]/(double)n_c_iters)*100)
        <<"\n";
    rout<<"Monotonicity       , "<<res_c[Monotonicity][Unsat]<<", "
        <<(int)(((double)res_c[Monotonicity][Unsat]/(double)n_c_conjs)*100)<<", "
        <<"NA, " //<<res_c[Monotonicity][Projected]<<", "
        <<"NA, " //<<(int)(((double)res_c[Monotonicity][Projected]/(double)n_c_iters)*100)
        <<"\n";
    rout<<"Projection         , "<<res_c[Projection][Unsat]<<", "
        <<(int)(((double)res_c[Projection][Unsat]/(double)n_c_conjs)*100)<<", "
        <<res_c[Projection][Projected]<<", "
        <<(int)(((double)res_c[Projection][Projected]/(double)n_c_iters)*100)
        <<"\n";
    rout<<"Simplification     , "<<res_c[Simplification][Unsat]<<", "
        <<(int)(((double)res_c[Simplification][Unsat]/(double)n_c_conjs)*100)<<", "
        <<res_c[Simplification][Projected]<<", "
        <<(int)(((double)res_c[Simplification][Projected]/(double)n_c_iters)*100)
        <<"\n";
    rout<<"\n\n";


  } // End of p loop

}







bool printRelation(string msg, Relation *rel){

    if ( rel ) {

        cout<<"\n\n"<<msg<<rel->toISLString()<<"\n\n";
    } else {

        cout<<"\n\n"<<msg<<"Not Satisfiable"<<"\n\n";
    }

    return true;
}

void EXPECT_EQ(string a, string b){

    if( a != b ){

        cout<<"\n\nExpected: "<<a;
        cout<<"\n\nActual:"<< b <<"\n\n";
    }
}

void EXPECT_EQ(Relation *a, Relation *b){

    if( a != b ){
        cout<<"\n\nIncorrect results: Expected or Actual is NULL.\n\n";
    }
}

int str2int(string str){
  int i;
  sscanf (str.c_str(),"%d",&i);
  return i;
}

