#ifndef VertexRestorer_h
#define VertexRestorer_h 1

#include "marlin/Processor.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include "lcio.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include "MathOperator.hh"
#include "CheatOperator.hh"
#include "MyReconstructedParticle.hh"
#include "TrackOperator.hh"
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include <EVENT/LCCollection.h>
#include <EVENT/Vertex.h>
#include <EVENT/Track.h>
#include <EVENT/ParticleID.h>
#include <IMPL/ParticleIDImpl.h>
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

namespace TTbarAnalysis
{
	class VertexRestorer : public Processor {
	  
	 public:
	  
	  virtual Processor*  newProcessor() { return new VertexRestorer ; }
	  
	  
	  VertexRestorer() ;
	  
	  virtual void init() ;
	  virtual void processRunHeader( LCRunHeader* run ) ;
	  virtual void processEvent( LCEvent * evt ) ; 
	  virtual void check( LCEvent * evt ) ; 
	  virtual void end() ;
	  
	  void PrintTrack(EVENT::Track * track); 
	  void PrintParticle(EVENT::ReconstructedParticle * particle);
	  float GetDeviation(const EVENT::ReconstructedParticle * particle, double * pos);
	  float GetError(const EVENT::ReconstructedParticle * particle);

	  std::vector< MyReconstructedParticle * > * RefineVertices(EVENT::LCCollection * main, EVENT::LCCollection * sec, EVENT::LCCollection * v0col = NULL);
	  std::vector< MyReconstructedParticle * > AddParticles(const std::vector< EVENT::ReconstructedParticle * > & pri, EVENT::Vertex * sec, const std::vector< EVENT::ReconstructedParticle * > * toCompare = NULL);

	  void CompareCollections(std::vector< MyReconstructedParticle * > * detected, EVENT::LCCollection * missed, EVENT::LCCollection * bstar, EVENT::LCCollection * missedvtx = NULL);
	  void CompareCollectionsRel(std::vector< MyReconstructedParticle * > * detected, EVENT::LCCollection * prongs, EVENT::LCCollection * bstar, EVENT::LCCollection *rel);

	  bool CompareParticles(const EVENT::ReconstructedParticle * particle1, const EVENT::ReconstructedParticle * particle2, bool verbose = false); 
	  bool CompareParticles(const EVENT::MCParticle * particle1, const EVENT::ReconstructedParticle * particle2); 
	  std::vector< EVENT::ReconstructedParticle * > * RestoreVerticesPFO(EVENT::LCCollection * main, EVENT::LCCollection * sec); 
	  std::vector< EVENT::ReconstructedParticle * > * RestoreVertices(EVENT::Vertex * primary, EVENT::LCCollection * sec); 
	  std::vector< EVENT::ReconstructedParticle * > * RestoreVertices(EVENT::LCCollection * sec, EVENT::LCCollection * rel); 
	  std::vector< EVENT::ReconstructedParticle * > * RestoreOneTrackVertices(EVENT::LCCollection * jet, EVENT::LCCollection * rel, EVENT::LCCollection * sec); 
	  std::vector< EVENT::ReconstructedParticle * > GetAdditionalParticles(const std::vector< EVENT::ReconstructedParticle * > & pri, EVENT::Vertex * sec, const std::vector< EVENT::ReconstructedParticle * > * toCompare = NULL);
	  std::vector< EVENT::ReconstructedParticle * > GetAdditionalParticles(const std::vector< EVENT::ReconstructedParticle * > & pri, const std::vector< EVENT::ReconstructedParticle * > * toCompare = NULL);
	  
	  EVENT::ReconstructedParticle * CopyParticle(const EVENT::ReconstructedParticle * particle, const EVENT::Vertex * vertex, float theta = 0.0);
	  bool IsDublicate(const EVENT::ReconstructedParticle * particle, std::vector< EVENT::ReconstructedParticle * > & data);
	  bool IsDublicate(MyReconstructedParticle * particle, std::vector< MyReconstructedParticle * > & data);
	  void WritePrimaries(TrackOperator & opera, EVENT::LCCollection * pricol, EVENT::LCCollection * sec = NULL, EVENT::LCCollection * rel = NULL);
	  void AddParticleID(EVENT::ReconstructedParticle * particle, float theta);
	  void ClearVariables();
	  void AnalyseSecondaries(EVENT::LCCollection * main, EVENT::LCCollection * missed);
	  bool IsParticleFromIP(const EVENT::ReconstructedParticle * particle);
	  bool TakeParticle(EVENT::ReconstructedParticle * primary, const EVENT::Vertex * vertex, std::vector<double> * parameters = NULL);

	  EVENT::MCParticle * CompareToCollection(EVENT::ReconstructedParticle * particle, EVENT::LCCollection * missed, EVENT::LCCollection * rel);
	  void TestMethod(EVENT::LCCollection * prongs, EVENT::LCCollection * pfo, EVENT::LCCollection * rel);
	  void TestCollections(EVENT::LCCollection * tagged, EVENT::LCCollection * rel);
	 protected:

	  std::string _colName ;
	  std::string _colRelName ;
	  std::string _colPriName ;
	  std::string _colSecName ;
	  std::string _colJetName ;
	  std::string _colJetRelName ;
	  std::string _clusterName ;
	  std::string _colMCMissedName ;
	  std::string _colMissedName ;
	  std::string _colMissedVtxName ;
	  std::string _colBStarName ;
	  std::string _colV0Name ;
	  std::string _colEGPName;
	 
	  static const int MAXN = 30;
	  float _fakeMomentum[MAXN];
	  float _bstarMomentum[MAXN];
	  float _missedMomentum[MAXN];
	  //float _v0Momentum[MAXN];
	  
	  float _fakeDeviation[MAXN];
	  float _fakeSecDeviation[MAXN];
	  //float _v0Deviation[MAXN];
	  float _bstarDeviation[MAXN];
	  float _missedDeviation[MAXN];
	  float _missedSecDeviation[MAXN];

	  float _fakeError[MAXN];
	  float _bstarError[MAXN];
	  float _missedError[MAXN];
	  //float _v0Error[MAXN];
	  float _fakeObs[MAXN];
	  float _missedObs[MAXN];

	  float _fakeCostheta[MAXN];
	  float _missedCostheta[MAXN];

	  float _allmissedMomentum[MAXN];
	  
	  float _missedAngle[MAXN];
	  float _v0Angle[MAXN];
	  float _bstarAngle[MAXN];
	  float _fakeAngle[MAXN];
	  static const int MAXP = 200;
	  int _primariesTotal;
	  float _primeOffset[MAXP];
	  float _primeDeviation[MAXP];
	  float _primeMomentum[MAXP];
	  float _primeError[MAXP];
		
	  double ip[3];
	  float _aParameter;
	  float _bParameter;
	  
	  int _missedDetected;
	  int _bstarDetected;
	  int _fakeDetected;

	  int _missedTotal;
	  int _bstarTotal;
	  int _detectedTotal;
	  //int _v0Total;
	  int _missedRecoTotal;
	  float _missedRecoMomentum[MAXN];
	  float _missedRecoTheta[MAXN];
	  float _missedRecoDeviation[MAXN];
	  int _nRun ;
	  int _nEvt ;

	  int _particlesTotal;
	  float _costhetaTestedTotal[MAXP];
	  float _costhetaTestedLost[MAXP];
	  float _costhetaTestedGood[MAXP];
	  float _costhetaTestedMultiple[MAXP];

	  int _ntotaltest;
	  int _nlosttest;
	  int _ngoodtest;
	  int _nmultitest;

	  float _offsetCut;
	  float _angleCut;

	  TFile* hfile ;
	  std::string _hfilename ;
	  TTree* _Tree;
	  TTree* _TestTree;
	  TTree* _primaryTree;
	  TrackOperator myTrackOperator;
	  EVENT::Vertex * myPrimary;

	} ;
}
	#endif



