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

    g++ -o build/bin/eval src/drivers/eval.cc -I src build/src/libiegenlib.a -lisl -std=c++11

 * Now to run the driver, you should put your dependence relations inside
 * JSON files and give them as inputs to the driver, one or more files at a time.
 * JSON file format is demonstrated by examples gs_csr.json and ilu_csr.json that
 * can be found in same directory as this driver. The json input files include commnet
 * fields that have the application code (Gauss-Seidel and ILU) in them, and says where
 * in the code each dependence relation is getting extract. So you can run following
 * afetr compiling the driver to get the simplified relations:
 
   build/bin/eval data/eval.list 

 *     
 *
 * \date Date Started: 1/10/2017
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
using iegenlib::Exp;
using iegenlib::Conjunction;
using iegenlib::Set;
using iegenlib::Relation;
using namespace std;

#define n_combo 10
enum methods {
  Super_Affine_Set,
  Approximation,
  Domain_Specific,
  Monotonicity,
  DS_And_Mon,
  Projection,
  App_Proj,
  Simplification
};
enum result {
  conj,
  Unsat,
  nit,
  Projected
};

int res_tot[n_combo][4];
std::ofstream rout("data/result.txt", std::ofstream::out);
std::ofstream rnout("data/result_notes.txt", std::ofstream::out);
void eval(string inputFile);

// Utility function
void printSetExp(string msg, std::set<Exp> se, iegenlib::TupleDecl td);
void printSetExpTF(string msg, std::set<Exp> se, iegenlib::TupleDecl td);
string method_name(int m);
bool printRelation(string msg, Relation *rel);
bool printRelationTF(string msg, Relation *rel);
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

    rout<<"Glossary:"<<"\n\n";
    rout<<"CA = Conjunctions available to determing satisfiability by this method\n"
        <<"CU = Number of  Conjunctions determined as unsatisfiable\n"
        <<"P  = Percentage of unsatisfiable conjunctions of available ones\n"
        <<"IA = Iterators available to determing satisfiability by this method\n"
        <<"IP = Number of  projected iterators\n"
        <<"P  = Percentage of projected iterators of available ones\n\n";

    ifstream in(argv[1]);
    while(in>>jname){
      eval(string( jname ));
    }
  }

  double ac,uc,ai,ui;
  rout<<"\n---- Final results:"<<"\n\n";
  rout<<"Name of method     | CA | CU | P | IA | IP | P\n";

  ac = res_tot[Super_Affine_Set][conj]; uc = res_tot[Super_Affine_Set][Unsat];
  rout<<"Super_Affine_Set   | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
      <<"NA | "<<"NA | "<<"NA | "<<"\n";

  ac = res_tot[Domain_Specific][conj]; uc = res_tot[Domain_Specific][Unsat];
  rout<<"Domain_Specific    | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
      <<"NA | "<<"NA | "<<"NA | "<<"\n";

  ac = res_tot[Monotonicity][conj]; uc = res_tot[Monotonicity][Unsat];
  rout<<"Monotonicity       | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
      <<"NA | "<<"NA | "<<"NA | "<<"\n";

  ac = res_tot[DS_And_Mon][conj]; uc = res_tot[DS_And_Mon][Unsat];
  rout<<"DS_And_Mon         | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
      <<"NA | "<<"NA | "<<"NA | "<<"\n";

  ac = res_tot[Projection][conj]; uc = res_tot[Projection][Unsat];
  ai = res_tot[Projection][nit]; ui = res_tot[Projection][Projected];
  rout<<"Projection         | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
      <<ai<<" | "<<ui<<" | "<<(int)((ui/ai)*100)<<" |\n";

  ac = res_tot[App_Proj][conj]; uc = res_tot[App_Proj][Unsat];
  ai = res_tot[App_Proj][nit]; ui = res_tot[App_Proj][Projected];
  rout<<"App_Proj           | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
      <<ai<<" | "<<ui<<" | "<<(int)((ui/ai)*100)<<" |\n";

  ac = res_tot[Simplification][conj]; uc = res_tot[Simplification][Unsat];
  ai = res_tot[Simplification][nit]; ui = res_tot[Simplification][Projected];
  rout<<"Simplification     | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
      <<ai<<" | "<<ui<<" | "<<(int)((ui/ai)*100)<<" |\n";

  return 0;
}

// Reads information from a JSON file (inputFile), and applies
// the simplification algorithm to the sets found in the file. 
void eval(string inputFile)
{
  iegenlib::setCurrEnv();
  std::set<int> parallelTvs;
  // (0)
  // Read the data from inputFile
  ifstream in(inputFile);
  json data;
  in >> data;

  for(size_t p = 0; p < data.size(); ++p){    // Dependence relations (DR) found in the file

    std::set<Exp> addConst;
    int res_c[n_combo][4] = {0};
    int n_c_conjs = 0;
    n_c_conjs += data[p].size();

    rnout<<"\n\n---- "<<data[p][0]["Name"].as<string>()<<":\n\n";

    for (size_t i = 0; i < data[p].size(); ++i){// Conjunctions found for one DR in the file

      int n_c_iters = 0;
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

      char msg[100];
      sprintf(msg, "@@@ Conjunction No. %d: ", int(i) );
      printRelationTF( string(msg) , rel);

      Relation* rel_copy = new Relation( rel->inArity() , rel->outArity() );
      *rel_copy = *rel;

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

      n_c_iters += rel->arity();
      n_c_iters -= parallelTvs.size();


      Relation* rel_norm = new Relation( rel->inArity() , rel->outArity() );
      *rel_norm = *rel;
      rel_norm->normalize();

      //################# Testing Super_Affine_Set #################
      res_c[Super_Affine_Set][conj]++;
      res_tot[Super_Affine_Set][conj]++;

      rel_copy->normalize();

      Relation* rel_t = new Relation( rel->inArity() , rel->outArity() );
      *rel_t = *rel;
      addConst = rel_t->boundDomainRange();
      printRelationTF(string("@@@ Conjunction with Super Affine"
                             " Set (Domain/Range) constraints: ") , rel_t);
      printSetExpTF(string("@@@ Added constraints by Super Affine"
                           " Set: ") , addConst, rel->getTupleDecl());
      delete rel_t;

      if( rel_copy->isDefault()) {
        res_c[Super_Affine_Set][Unsat]++;
        res_tot[Super_Affine_Set][Unsat]++;
        UNS = true;

        rnout<<"\n### The relation "<<i<<" is unsat: Super_Affine_Set"<<"\n";
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
        res_c[Domain_Specific][conj]++;
        res_tot[Domain_Specific][conj]++;

        *rel_copy = *rel;
        Relation *rel_extend;

        // Doing Super Affine Set before adding domain specific constraints
        rel_copy->normalize();

        for (size_t j = 0; j < data[p][0]["User Defined"].size(); ++j){
   
          rel_extend = rel_copy->addUFConstraints(
                        data[p][0]["User Defined"][j]["Func1"].as<string>(),
                        data[p][0]["User Defined"][j]["operator"].as<string>(),
                        data[p][0]["User Defined"][j]["Func2"].as<string>()
                                              );
          *rel_copy = *rel_extend;
          delete rel_extend;
        }

        addConst = constraintsDifference( rel_copy->mConjunctions.front() , rel_norm->mConjunctions.front() );
        printRelationTF(string("@@@ Conjunction with Domain Specific constraints: ") , rel_copy);
        printSetExpTF(string("@@@ Added constraints by Domain Specific: ") , addConst, rel->getTupleDecl());

        rel_copy->normalize();
        if( rel_copy->isDefault()) {
          res_c[Domain_Specific][Unsat]++;
          res_tot[Domain_Specific][Unsat]++;

          rnout<<"\n### The relation "<<i<<" is unsat: Domain_Specific"<<"\n";
        }
        else {
          //find # of new founded equalities;
        } 
      }


      //################# Testing Monotonicity #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        res_c[Monotonicity][conj]++;
        res_tot[Monotonicity][conj]++;

        *rel_copy = *rel;
        Relation *rel_extend;

        // Doing Super Affine Set before adding domain specific constraints
        rel_copy->normalize();

        rel_extend = rel_copy->addConstraintsDueToMonotonicity();

        *rel_copy = *rel_extend;
        delete rel_extend;

        addConst = constraintsDifference( rel_copy->mConjunctions.front() , rel_norm->mConjunctions.front() );
        printRelationTF(string("@@@ Conjunction with Monotonicity constraints: ") , rel_copy);
        printSetExpTF(string("@@@ Added constraints by Monotonicity: ") , addConst, rel->getTupleDecl());

        rel_copy->normalize();
        if( rel_copy->isDefault()) {

          res_c[Monotonicity][Unsat]++;
          res_tot[Monotonicity][Unsat]++;

          rnout<<"\n### The relation "<<i<<" is unsat: Monotonicity"<<"\n";
        }
        else {
          //find # of new founded equalities;
        } 
      }


      //################# Testing DS_And_Mon    #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        res_c[DS_And_Mon][conj]++;
        res_tot[DS_And_Mon][conj]++;

        *rel_copy = *rel;
        Relation *rel_extend;

        // Doing Super Affine Set before adding domain specific constraints
        rel_copy->normalize();

        for (size_t j = 0; j < data[p][0]["User Defined"].size(); ++j){
   
          rel_extend = rel_copy->addUFConstraints(
                        data[p][0]["User Defined"][j]["Func1"].as<string>(),
                        data[p][0]["User Defined"][j]["operator"].as<string>(),
                        data[p][0]["User Defined"][j]["Func2"].as<string>()
                                              );
          *rel_copy = *rel_extend;
          delete rel_extend;
        }

        rel_extend = rel_copy->addConstraintsDueToMonotonicity();
        *rel_copy = *rel_extend;
        delete rel_extend;

        addConst = constraintsDifference( rel_copy->mConjunctions.front() , rel->mConjunctions.front() );
        printRelationTF(string("@@@ Conjunction with DS_And_Mon constraints: ") , rel_copy);
        printSetExpTF(string("@@@ Added constraints by Monotonicity: ") , addConst, rel->getTupleDecl());

        rel_copy->normalize();
        if( rel_copy->isDefault()) {

          res_c[DS_And_Mon][Unsat]++;
          res_tot[DS_And_Mon][Unsat]++;

          rnout<<"\n### The relation "<<i<<" is unsat: DS_And_Mon"<<"\n";
        }
        else {
          //find # of new founded equalities;
        } 
      }


      //################# Testing Projection #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        // Keeping track of available iterators for projection
        res_c[Projection][conj]++;
        res_tot[Projection][conj]++;

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

          rnout<<"\n### The relation "<<i<<" is unsat: Projection"<<"\n";
        }
        else {

          res_c[Projection][nit] += n_c_iters;
          res_tot[Projection][nit] += n_c_iters;

          int pr = rel->arity() - rel_copy->arity();
          res_c[Projection][Projected] += pr;
          res_tot[Projection][Projected] += pr;

          rnout<<"\n### Projected "<<pr<<" TVs from relation "<<i<<": Projection"<<"\n";
        } 
      }

      //################# Testing Projection+Approaximation #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        // Keeping track of available iterators for projection
        res_c[App_Proj][conj]++;
        res_tot[App_Proj][conj]++;

        Relation *temp;
        int lastTV = rel->arity()-1;
        *rel_copy = *rel;

        // Applying heuristic for removing expensive iterators
        int numConstToRemove = str2int(data[p][0]["Remove Constraints"].as<string>());
        rel_copy->RemoveExpensiveConsts(parallelTvs, numConstToRemove );

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
          res_c[App_Proj][Unsat]++;
          res_tot[App_Proj][Unsat]++;

          rnout<<"\n### The relation "<<i<<" is unsat: App_Proj"<<"\n";
        }
        else {

          res_c[App_Proj][nit] += n_c_iters;
          res_tot[App_Proj][nit] += n_c_iters;

          int pr = rel->arity() - rel_copy->arity();
          res_c[App_Proj][Projected] += pr;
          res_tot[App_Proj][Projected] += pr;

          rnout<<"\n### Projected "<<pr<<" TVs from relation "<<i<<": App_Proj"<<"\n";
        } 
      }

      //################# Testing simplifyForPartialParallel #################
      if( ! UNS ){  // Right now only testing if Super_Affine_Set
                    // has not already determined unsatisfiabiltiy
        res_c[Simplification][conj]++;
        res_tot[Simplification][conj]++;

        // Applying heuristic for removing expensive iterators
        int numConstToRemove = str2int(data[p][0]["Remove Constraints"].as<string>());
        rel->RemoveExpensiveConsts(parallelTvs, numConstToRemove );

        rel->normalize();

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

        rel_copy = rel->simplifyForPartialParallel(parallelTvs);

        if( rel_copy ){
          addConst = constraintsDifference( rel_copy->mConjunctions.front() , 
                                            rel_norm->mConjunctions.front() );
          printRelationTF(string("@@@ Conjunction after simplification algorithm: ")
                          , rel_copy);
          printSetExpTF(string("@@@ Added constraints by Monotonicity: ") , 
                        addConst, rel->getTupleDecl());
           
          rel_copy->normalize();
        } else {
          printRelationTF(string("@@@ Conjunction after simplification algorithm: ")
                          , rel_copy);
          addConst.clear();
          printSetExpTF(string("@@@ Added constraints by Monotonicity: ") , 
                        addConst, rel->getTupleDecl());
        }

        if( !rel_copy || rel_copy->isDefault() ) {
          res_c[Simplification][Unsat]++;
          res_tot[Simplification][Unsat]++;

          rnout<<"\n### The relation "<<i<<" is unsat: Simplification"<<"\n";
        }
        else {

          res_c[Simplification][nit] += n_c_iters;
          res_tot[Simplification][nit] += n_c_iters;

          int pr = rel->arity() - rel_copy->arity();
          res_c[Simplification][Projected] += pr;
          res_tot[Simplification][Projected] += pr;

          rnout<<"\n### Projected "<<pr<<" TVs from relation "<<i<<": Simplification"<<"\n";
        } 
      }

      delete rel;
    }

    double ac,uc,ai,ui;
    rout<<"---- "<<data[p][0]["Name"].as<string>()<<":\n\n";
    rout<<"Name of method     | CA | CU | P | IA | IP | P\n";

    ac = res_c[Super_Affine_Set][conj]; uc = res_c[Super_Affine_Set][Unsat];
    rout<<"Super_Affine_Set   | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
        <<"NA | "<<"NA | "<<"NA | "<<"\n";

    ac = res_c[Domain_Specific][conj]; uc = res_c[Domain_Specific][Unsat];
    rout<<"Domain_Specific    | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
        <<"NA | "<<"NA | "<<"NA | "<<"\n";

    ac = res_c[Monotonicity][conj]; uc = res_c[Monotonicity][Unsat];
    rout<<"Monotonicity       | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
        <<"NA | "<<"NA | "<<"NA | "<<"\n";

    ac = res_c[DS_And_Mon][conj]; uc = res_c[DS_And_Mon][Unsat];
    rout<<"DS_And_Mon         | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
        <<"NA | "<<"NA | "<<"NA | "<<"\n";

    ac = res_c[Projection][conj]; uc = res_c[Projection][Unsat];
    ai = res_c[Projection][nit]; ui = res_c[Projection][Projected];
    rout<<"Projection         | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
        <<ai<<" | "<<ui<<" | "<<(int)((ui/ai)*100)<<" |\n";

    ac = res_c[App_Proj][conj]; uc = res_c[App_Proj][Unsat];
    ai = res_c[App_Proj][nit]; ui = res_c[App_Proj][Projected];
    rout<<"App_Proj           | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
        <<ai<<" | "<<ui<<" | "<<(int)((ui/ai)*100)<<" |\n";

    ac = res_c[Simplification][conj]; uc = res_c[Simplification][Unsat];
    ai = res_c[Simplification][nit]; ui = res_c[Simplification][Projected];
    rout<<"Simplification     | "<<ac<<" | "<<uc<<" | "<<(int)((uc/ac)*100)<<" | "
        <<ai<<" | "<<ui<<" | "<<(int)((ui/ai)*100)<<" |\n";
    rout<<"\n\n";

  } // End of p loop

}


void printSetExp(string msg, std::set<Exp> se, iegenlib::TupleDecl td){

  cout<<"\n"<<msg<<" { ";
  for (std::set<Exp>::iterator it=se.begin(); it!=se.end(); it++){
    cout<<(*it).prettyPrintString( td );
    Exp te = (Exp)(*it);
    if( te.isEquality() ) cout<<" = 0 ";
    else cout<<" >= 0 ";
    if( std::next(it) != se.end()) cout<<" && ";
  }
  cout<<"}\n";
}


void printSetExpTF(string msg, std::set<Exp> se, iegenlib::TupleDecl td){

  rnout<<"\n"<<msg<<" { ";
  for (std::set<Exp>::iterator it=se.begin(); it!=se.end(); it++){
    rnout<<(*it).prettyPrintString( td );
    Exp te = (Exp)(*it);
    if( te.isEquality() ) rnout<<" = 0 ";
    else rnout<<" >= 0 ";
    if( std::next(it) != se.end()) rnout<<" && ";
  }
  rnout<<"}\n";
}

bool printRelationTF(string msg, Relation *rel){

    if ( rel ){
        rnout<<"\n"<<msg<<rel->toISLString()<<"\n";
    } else {
        rnout<<"\n"<<msg<<"Not Satisfiable"<<"\n";
    }

    return true;
}

string method_name(int m)
{
  if( m == Approximation)  return string("Approximation");
  else if( m == Domain_Specific)  return string("Domain_Specific");
  else if( m == Monotonicity)  return string("Monotonicity");
  else if( m == Super_Affine_Set)  return string("Super_Affine_Set");
  else if( m == Projection)  return string("Projection");
  else if( m == App_Proj)  return string("App_Proj");
  else if( m == Simplification)  return string("Simplification");

  return string("Not Specified");
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

