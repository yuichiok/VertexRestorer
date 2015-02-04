#include "VertexRestorer.hh"

using namespace lcio ;
using namespace marlin ;


VertexRestorer aVertexRestorer ;


VertexRestorer::VertexRestorer() : Processor("VertexRestorer") {
  
  // modify processor description
  _description = "VertexRestorer";
  

  // register steering parameters: name, description, class-variable, default value
  
  registerProcessorParameter( "ROOTFileName",
                              "Output ROOT File Name",
                              _hfilename,
                              std::string("VertexRestorer.root") );
  
/*  registerInputCollection( LCIO::CALORIMETERHIT , 
			   "CalorimeterHitCollection",
			   "Name of the Calorimeter hit collection"  ,
			   _colName ,
			   std::string("ConvCalorimeterHits") ) ;
*/
  
}


void VertexRestorer::init() 
{ 
   // usually a good idea to
   printParameters() ;
 
   _nRun = 0 ;
   _nEvt = 0 ;   
   
   hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;
 
   _Tree = new TTree( "tree", "tree" );
   
   _Tree->Branch("evtn", &_nEvt, "evtn/I");
   _Tree->Branch("runN", &_nRun, "runN/I");

}

void VertexRestorer::processRunHeader( LCRunHeader* run) 
{ 
   _nRun++ ;
} 

void VertexRestorer::processEvent( LCEvent * evt ) 
{ 
	_nEvt ++ ;
  
	try 
	{
		
	}
	catch( DataNotAvailableException &e)
	{
		std::cout << "Whoops!....\n";
	}
    
}


void VertexRestorer::check( LCEvent * evt ) {}



void VertexRestorer::end()
{   
  hfile->Write();
  hfile->Close();
}

