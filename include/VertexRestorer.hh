#ifndef VertexRestorer_h
#define VertexRestorer_h 1

#include "marlin/Processor.h"
#include <iostream>
#include <stdlib.h>
#include "lcio.h"
#include "TFile.h"
#include "TTree.h"
#include <string>

using namespace lcio ;
using namespace marlin ;


/** 
 *  <h4> ==== VertexRestorer Processor ==== </h4>
 *  
 * Used to bla bla
 * 
 * @param ROOTFileName: Name of the output ROOT file
 *   
 * @author Bilokin Sviatoslav, LAL
 * @version $Id: VertexRestorer.h,v 1.0 2010/02/18 Exp $ 
 */


class VertexRestorer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new VertexRestorer ; }
  
  
  VertexRestorer() ;
  
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  
  
 protected:

  std::string _colName ;
  std::string _clusterName ;
  

  int _nRun ;
  int _nEvt ;
  
  TFile* hfile ;
  std::string _hfilename ;
  TTree* _Tree;
} ;

#endif



