#include "VertexRestorer.hh"

using namespace lcio ;
using namespace marlin ;
using std::vector;
using std::string;
using std::abs;
using EVENT::Track;
using IMPL::ReconstructedParticleImpl;
using IMPL::VertexImpl;
using EVENT::ParticleID;
using IMPL::ParticleIDImpl;
namespace TTbarAnalysis
{
  VertexRestorer aVertexRestorer ;


  VertexRestorer::VertexRestorer() : Processor("VertexRestorer") {
	  
    // modify processor description
    _description = "VertexRestorer";
	  

    // register steering parameters: name, description, class-variable, default value
	  
    registerProcessorParameter( "ROOTFileName",
				"Output ROOT File Name",
				_hfilename,
				string("VertexRestorer.root") );
	  
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE , 
			     "PFOCollection",
			     "Name of the Calorimeter hit collection"  ,
			     _colName ,
			     string("PandoraPFOs") ) ;
    registerInputCollection( LCIO::MCPARTICLE ,                                                                                                                                             
                             "MCCollectionName",                                                                                    
			     "Name of the MC Collection"  ,   
			     _MCcolName,
			     string("MCParticle") 
			     ) ;
    registerInputCollection( LCIO::VERTEX,
			     "PrimaryCollectionName" , 
			     "Name of the PrimaryVertex collection"  ,
			     _colPriName ,
			     std::string("PrimaryVertex")
			     );

    registerInputCollection( LCIO::VERTEX,
			     "FinalVertexCollectionName" ,
			     "what is this?? Only used for purgatory writting?"  ,
			     _colFinalVertex ,
			     std::string("RefinedVertex")
			     );

    registerOutputCollection( LCIO::VERTEX,
			      "OutputCollectionName" , 
			      "Name of the Vertex collection"  ,
			      _outputcolName ,
			      std::string("RecoveredJets_vtx")
			      );
    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			      "OutputJetCollectionName" , 
			      "Name of the Jet collection"  ,
			      _outputjetcolName ,
			      std::string("RecoveredJets")
			      );
    registerOutputCollection( LCIO::LCRELATION,
			      "OutputRelCollectionName" , 
			      "Name of the jet vertex relation collection"  ,
			      _outputjetrelcolName ,
			      std::string("RecoveredJets_rel")
			      );
    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			      "OutputRelRPCollectionName" , 
			      "Name of the jet vertex RP collection"  ,
			      _outputjetrelRPcolName ,
			      std::string("RecoveredJets_vtx_RP")
			      );
    registerInputCollection( LCIO::VERTEX,
			     "SecondaryCollectionName" , 
			     "Name of the BuildUpVertex collection"  ,
			     _colSecName ,
			     std::string("FinalJets_vtx")
			     );
    registerInputCollection( LCIO::VERTEX,
			     "SecondaryRPCollectionName" , 
			     "Name of the BuildUpVertex collection"  ,
			     _colSecRPName ,
			     std::string("FinalJets_vtx_RP")
			     );
    registerInputCollection( LCIO::VERTEX,
			     "V0CollectionName" , 
			     "Name of the BuildUpVertex_V0 collection"  ,
			     _colV0Name ,
			     std::string("BuildUpVertex_V0")
			     );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "JetCollectionName" , 
			     "Name of the Jet collection"  ,
			     _colJetName ,
			     std::string("FinalJets")
			     );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "MissedCollectionName" , 
			     "Name of the Jet collection"  ,
			     _colMissedName ,
			     std::string("MissedParticles")
			     );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "MCMissedCollectionName" , 
			     "Name of the Jet collection"  ,
			     _colMCMissedName ,
			     std::string("MCMissedParticles")
			     );
    registerInputCollection( LCIO::TRACK,
			     "NotUsedTracksCollectionName" , 
			     "Name of tracks not passed to PFA collection"  ,
			     _colTrashCanName ,
			     std::string("TracksFailBothCanFormPfoFlags")
			     );
    registerInputCollection( LCIO::LCRELATION,
			     "RecoMcTruthCollectionName" , 
			     "Name of the RecoMcTruth collection"  ,
			     _colRelName ,
			     std::string("RecoMCTruthLink")
			     );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "MissedVtxCollectionName" , 
			     "Name of the Jet collection"  ,
			     _colMissedVtxName ,
			     std::string("MissedParticlesVtx")
			     );
    registerInputCollection( LCIO::MCPARTICLE,
			     "BStarCollectionName" , 
			     "Name of the BStar collection"  ,
			     _colBStarName ,
			     std::string("BStar")
			     );
    registerInputCollection( LCIO::VERTEX,
			     "MCVertexCollectionName",
			     "Name of the MCVertex collection",
			     _colMCvertexName,
			     std::string("MCVertex")
			     );
    registerInputCollection( LCIO::MCPARTICLE,
			     "EGProngsCollectionName" , 
			     "Name of the EGProngs collection"  ,
			     _colEGPName ,
		    
			     //std::string("MissedParticles")
			     std::string("EGProngs")
			     );
    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "JetRelCollectionName",
			     "Name of the Jet relation collection",
			     _colJetRelName,
			     std::string("FinalJets_rel")
			     );
    _useTracks = 0;
    registerProcessorParameter("UseTracks" , 
			       "Use tracks for association"  ,
			       _useTracks,
			       _useTracks
			       );
    _Test = 0;
    registerProcessorParameter("TestEfficiency" , 
			       "Use tracks for association"  ,
			       _Test,
			       _Test
			       );
    _Bfield=3.5;
    registerProcessorParameter( "Bfield",
                             "ILD Magnetic field", 
                             _Bfield,
                             _Bfield
                             );
    _angleCut = 0.05;
    registerProcessorParameter("angleCut" , 
			       "angle cut for association"  ,
			       _angleCut,
			       _angleCut
			       );
    _offsetCut = 0.05;
    registerProcessorParameter("offsetCut" , 
			       "Offset cut for association"  ,
			       _offsetCut,
			       _offsetCut
			       );
    _aParameter = 0.005;
    registerProcessorParameter("a" , 
			       "a parameter of accuracy in mm"  ,
			       _aParameter,
			       _aParameter
			       );
    _bParameter = 0.01;
    registerProcessorParameter("b" , 
			       "b parameter of accuracy in mm"  ,
			       _bParameter,
			       _bParameter
			       );
	    
    ip[0] = 0.0;
    ip[1] = 0.0;
    ip[2] = 0.0;
	  
  }


  void VertexRestorer::init() 
  { 
    // usually a good idea to
    printParameters() ;
	 
    _nRun = 0 ;
    _nEvt = 0 ;   
	   
    hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;
	 
    _Tree = new TTree( "Stats", "tree" );
	   
    _Tree->Branch("detectedTotal", &_detectedTotal, "detectedTotal/I");
    _Tree->Branch("missedTotal", &_missedTotal, "missedTotal/I");
    _Tree->Branch("allmissedMomentum", _allmissedMomentum, "allmissedMomentum[missedTotal]/F");
    _Tree->Branch("bstarTotal", &_bstarTotal, "bstarTotal/I");
	   
    _Tree->Branch("missedDetected", &_missedDetected, "missedDetected/I");
    _Tree->Branch("missedMomentum", _missedMomentum, "missedMomentum[missedDetected]/F");
    _Tree->Branch("missedError", _missedError, "missedError[missedDetected]/F");
    _Tree->Branch("missedOffset", _missedOffset, "missedOffset[missedDetected]/F");
    _Tree->Branch("missedAngle", _missedAngle, "missedAngle[missedDetected]/F");
    _Tree->Branch("missedObs", _missedObs, "missedObs[missedDetected]/F");
    _Tree->Branch("missedZ", _missedZ, "missedZ[missedDetected]/F");
    _Tree->Branch("missedVtxHits", _missedVtxHits, "missedVtxHits[missedDetected]/I");
    _Tree->Branch("missedFtdHits", _missedFtdHits, "missedFtdHits[missedDetected]/I");
    _Tree->Branch("missedDeviation", _missedDeviation, "missedDeviation[missedDetected]/F");
    _Tree->Branch("missedD0Deviation", _missedD0Deviation, "missedD0Deviation[missedDetected]/F");
    _Tree->Branch("missedZ0Deviation", _missedZ0Deviation, "missedZ0Deviation[missedDetected]/F");
    _Tree->Branch("missedZ0", _missedZ0, "missedZ0[missedDetected]/F");
    _Tree->Branch("missedD0", _missedD0, "missedD0[missedDetected]/F");
    _Tree->Branch("missedZ0Error", _missedZ0Error, "missedZ0Error[missedDetected]/F");
    _Tree->Branch("missedD0Error", _missedD0Error, "missedD0Error[missedDetected]/F");
    _Tree->Branch("missedPvtx", _missedSecDeviation, "missedPvtx[missedDetected]/F");
    _Tree->Branch("missedCostheta", _missedCostheta, "missedCostheta[missedDetected]/F");
    _Tree->Branch("missedCosthetaVtx", _missedCosthetaVtx, "missedCosthetaVtx[missedDetected]/F");
	   
    /*	   _Tree->Branch("bstarDetected", &_bstarDetected, "bstarDetected/I");
	   _Tree->Branch("bstarMomentum", _bstarMomentum, "bstarMomentum[bstarDetected]/F");
	   _Tree->Branch("bstarError", _bstarError, "bstarError[bstarDetected]/F");
	   _Tree->Branch("bstarAngle", _bstarAngle, "bstarAngle[bstarDetected]/F");
	   _Tree->Branch("bstarDeviation", _bstarDeviation, "bstarDeviation[bstarDetected]/F");
	   
	   _Tree->Branch("v0Detected", &_v0Detected, "v0Detected/I");
	   _Tree->Branch("v0Momentum", _v0Momentum, "v0Momentum[v0Detected]/F");
	   _Tree->Branch("v0Error", _v0Error, "v0Error[v0Detected]/F");
	   _Tree->Branch("v0Angle", _v0Angle, "v0Angle[v0Detected]/F");
	   _Tree->Branch("v0Deviation", _v0Deviation, "v0Deviation[v0Detected]/F");*/
	   
    _Tree->Branch("fakeDetected", &_fakeDetected, "fakeDetected/I");
    _Tree->Branch("fakeMomentum", _fakeMomentum, "fakeMomentum[fakeDetected]/F");
    _Tree->Branch("fakeError", _fakeError, "fakeError[fakeDetected]/F");
    _Tree->Branch("fakeAngle", _fakeAngle, "fakeAngle[fakeDetected]/F");
    _Tree->Branch("fakeObs", _fakeObs, "fakeObs[fakeDetected]/F");
    _Tree->Branch("fakeZ", _fakeZ, "fakeZ[fakeDetected]/F");
    _Tree->Branch("fakeOffset", _fakeOffset, "fakeOffset[fakeDetected]/F");
    _Tree->Branch("fakeVtxHits", _fakeVtxHits, "fakeVtxHits[fakeDetected]/I");
    _Tree->Branch("fakeFtdHits", _fakeFtdHits, "fakeFtdHits[fakeDetected]/I");
    _Tree->Branch("fakeDeviation", _fakeDeviation, "fakeDeviation[fakeDetected]/F");
    _Tree->Branch("fakeD0Deviation", _fakeD0Deviation, "fakeD0Deviation[fakeDetected]/F");
    _Tree->Branch("fakeZ0Deviation", _fakeZ0Deviation, "fakeZ0Deviation[fakeDetected]/F");
    _Tree->Branch("fakeZ0", _fakeZ0, "fakeZ0[fakeDetected]/F");
    _Tree->Branch("fakeD0", _fakeD0, "fakeD0[fakeDetected]/F");               
    _Tree->Branch("fakeZ0Error", _fakeZ0Error, "fakeZ0Error[fakeDetected]/F");
    _Tree->Branch("fakeD0Error", _fakeD0Error, "fakeD0Error[fakeDetected]/F");
    _Tree->Branch("fakePvtx", _fakeSecDeviation, "fakePvtx[fakeDetected]/F");
    _Tree->Branch("fakeCostheta", _fakeCostheta, "fakeCostheta[fakeDetected]/F");
    _Tree->Branch("fakeCosthetaVtx", _fakeCosthetaVtx, "fakeCosthetaVtx[fakeDetected]/F");
	   
    _Tree->Branch("purDetected", &_purDetected, "purDetected/I");
    _Tree->Branch("purOffset", _purOffset, "purOffset[purDetected]/F");
    _Tree->Branch("purAngle", _purAngle, "purAngle[purDetected]/F");
    _Tree->Branch("purCostheta", _purCostheta, "purCostheta[purDetected]/F");
    _Tree->Branch("purMomentum", _purMomentum, "purMomentum[purDetected]/F");
    _Tree->Branch("missedRecoTotal", &_missedRecoTotal, "missedRecoTotal/I");
    _Tree->Branch("missedRecoMomentum", _missedRecoMomentum, "missedRecoMomentum[missedRecoTotal]/F");
    _Tree->Branch("missedRecoTheta", _missedRecoTheta, "missedRecoTheta[missedRecoTotal]/F");


    _primaryTree =  new TTree( "Primaries", "tree" );
    _primaryTree->Branch("primeZ", &_primeZ, "primeZ/F");
    _primaryTree->Branch("primeT", &_primeT, "primeT/F");
    _primaryTree->Branch("primariesTotal", &_primariesTotal, "primariesTotal/I");
    _primaryTree->Branch("primeType", &_primeType, "primeType[primariesTotal]/I");
    _primaryTree->Branch("primeOffset", &_primeOffset, "primeOffset[primariesTotal]/F");
    _primaryTree->Branch("primeMomentum", &_primeMomentum, "primeMomentum[primariesTotal]/F");
    _primaryTree->Branch("primeDeviation", &_primeDeviation, "primeDeviation[primariesTotal]/F");
    _primaryTree->Branch("primeError", &_primeError, "primeError[primariesTotal]/F");
    _primaryTree->Branch("primeTeorError", &_primeTeorError, "primeTeorError[primariesTotal]/F");
    _primaryTree->Branch("primeChi2", &_primeChi2, "primeChi2[primariesTotal]/F");
    _primaryTree->Branch("primeCostheta", &_primeCostheta, "primeCostheta[primariesTotal]/F");
    _primaryTree->Branch("primeVtxHits", &_primeVtxHits, "primeVtxHits[primariesTotal]/I");
    _primaryTree->Branch("primeFtdHits", &_primeFtdHits, "primeFtdHits[primariesTotal]/I");
    _primaryTree->Branch("primeEtdHits", &_primeEtdHits, "primeEtdHits[primariesTotal]/I");
	   
    _purgatoryTree =  new TTree( "Purgatory", "tree" );
    _purgatoryTree->Branch("purgatoryTotal", &_purgatoryTotal, "purgatoryTotal/I");
    _purgatoryTree->Branch("purgatoryOffset", _purgatoryOffset, "purgatoryOffset[purgatoryTotal]/F");
    _purgatoryTree->Branch("purgatoryMomentum", _purgatoryMomentum, "purgatoryMomentum[purgatoryTotal]/F");
    _purgatoryTree->Branch("purgatoryDeviation", _purgatoryDeviation, "purgatoryDeviation[purgatoryTotal]/F");
    _purgatoryTree->Branch("purgatoryError", _purgatoryError, "purgatoryError[purgatoryTotal]/F");
    _purgatoryTree->Branch("purgatoryTeorError", _purgatoryTeorError, "purgatoryTeorError[purgatoryTotal]/F");
    _purgatoryTree->Branch("purgatoryCostheta", _purgatoryCostheta, "purgatoryCostheta[purgatoryTotal]/F");
    _purgatoryTree->Branch("purgatoryVtxHits", _purgatoryVtxHits, "purgatoryVtxHits[purgatoryTotal]/I");
    _purgatoryTree->Branch("purgatoryFtdHits", _purgatoryFtdHits, "purgatoryFtdHits[purgatoryTotal]/I");
    _purgatoryTree->Branch("purgatoryParent", _purgatoryParent, "purgatoryParent[purgatoryTotal]/I");
	   
    _trashTree =  new TTree( "Zombies", "tree" );
    _trashTree->Branch("trashTotal", &_trashTotal, "trashTotal/I");
    _trashTree->Branch("trashVtxHits", _trashVtxHits, "trashVtxHits[trashTotal]/I");
    _trashTree->Branch("trashTpcHits", _trashTpcHits, "trashTpcHits[trashTotal]/I");
    _trashTree->Branch("trashHasMC", _trashHasMCParticle, "trashHasMC[trashTotal]/I");
    _trashTree->Branch("trashIsReco", _trashIsReco, "trashIsReco[trashTotal]/I");
    _trashTree->Branch("trashIsDublicate", _trashIsDublicate, "trashIsDublicate[trashTotal]/I");
    _trashTree->Branch("trashCostheta", _trashCostheta, "trashCostheta[trashTotal]/F");

    _TestTree =  new TTree( "Test", "tree" );
    /*_TestTree->Branch("ntotaltest", &_ntotaltest, "ntotaltest/I");
      _TestTree->Branch("ngoodtest", &_ngoodtest, "ngoodtest/I");
      _TestTree->Branch("nlosttest", &_nlosttest, "nlosttest/I");
      _TestTree->Branch("nmultitest", &_nmultitest, "nmultitest/I");
      _TestTree->Branch("costhetaTestedTotal", _costhetaTestedTotal, "costhetaTestedTotal[ntotaltest]/F");
      _TestTree->Branch("costhetaTestedGood", _costhetaTestedGood, "costhetaTestedGood[ngoodtest]/F");
      _TestTree->Branch("costhetaTestedLost", _costhetaTestedLost, "costhetaTestedLost[nlosttest]/F");
      _TestTree->Branch("costhetaTestedMultiple", _costhetaTestedMultiple, "costhetaTestedMultiple[nmultitest]/F");*/
    _TestTree->Branch("ntest1", &_ntest, "ntest/I");
    _TestTree->Branch("parameter1", _testparameter, "parameter1[ntest]/F");
    _TestTree->Branch("test1", _testDistance, "test1[ntest]/F");

  }

  void VertexRestorer::processRunHeader( LCRunHeader* run) 
  { 
    _nRun++ ;
  } 
  float VertexRestorer::GetDeviation(const EVENT::ReconstructedParticle * particle, double * pos)
  {
    vector<float> direction = MathOperator::getDirection(particle->getMomentum());
    float accuracy = GetError(particle); // sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p * pow(sin(angles[1]), 4.0/3.0)) );
    float offset = MathOperator::getDistanceTo(ip, direction, pos);
    return offset / accuracy;
  }
  void  VertexRestorer::PrintTrack(Track * track)
  {
    if (!track) 
      {
	return;
      }
    streamlog_out(DEBUG) << "|\t" << track->getD0() 
	      << "|\t" << track->getZ0() 
	      << "|\t" << track->getPhi() 
	      << "|\t" << track->getOmega() 
	      << "|\t" << track->getTanLambda()
	      << '\n';
  }

  bool VertexRestorer::IsParticleFromIP(const EVENT::ReconstructedParticle * particle)
  {
    if (!particle) 
      {
	return false;
      }
    const Vertex * vertex = particle->getStartVertex();
    if (vertex) 
      {
	return vertex->isPrimary();
      }
    streamlog_out(DEBUG) << "Particle has no vertex!\n";
    return false;
  }
  void  VertexRestorer::CopyParameters(EVENT::LCParameters & from, EVENT::LCParameters & to)
  {
    vector<string> intkeys;
    from.getIntKeys(intkeys);
    vector<string> floatkeys;
    from.getFloatKeys(floatkeys);
    vector<string> stringkeys;
    from.getStringKeys(stringkeys);
    //streamlog_out(DEBUG) << "Int keys size: " << intkeys.size() << "\n";
    for (unsigned int i = 0; i < intkeys.size(); i++) 
      {
	vector<int> intvals;
	from.getIntVals(intkeys[i], intvals);
	//streamlog_out(DEBUG) << "Int key "<< intkeys[i] << " has " << intvals.size() << " values\n";
	to.setValues(intkeys[i], intvals);
      }
    //streamlog_out(DEBUG) << "Float keys size: " << floatkeys.size() << "\n";
    for (unsigned int i = 0; i < floatkeys.size(); i++) 
      {
	vector<float> floatvals;
	from.getFloatVals(floatkeys[i], floatvals);
	/*streamlog_out(DEBUG) << "Float key "<< floatkeys[i] << " has " << floatvals.size() << " values\n";
	  for (int j = 0; j < floatvals.size(); j++) 
	  {
	  streamlog_out(DEBUG) << "\t" << floatvals[j] <<"\n";
	  }*/
	to.setValues(floatkeys[i], floatvals);
      }
    //streamlog_out(DEBUG) << "String keys size: " << stringkeys.size() << "\n";
    for (unsigned int i = 0; i < stringkeys.size(); i++) 
      {
	vector<string> stringvals;
	from.getStringVals(stringkeys[i], stringvals);
	//streamlog_out(DEBUG) << "String key "<< stringkeys[i] << " has " << stringvals.size() << " values\n";
	//	streamlog_out(DEBUG) << "\t" << from.getStringVal(stringkeys[0]) <<"\n";
	to.setValues(stringkeys[i], stringvals);
      }
  }
  std::vector< Vertex * > * VertexRestorer::convert(const std::vector< LCObject * > & objs)
  {
    std::vector< Vertex * > * result = new std::vector< Vertex * >();
    for (int i = 0; i < objs.size(); i++) 
      {
	result->push_back(dynamic_cast< Vertex * >(objs[i]));
      }
    return result;
  }
		
  void VertexRestorer::PrintParticle(ReconstructedParticle * particle)
  {
    if (!particle) 
      {
	return;
      }
    int vtxhits = 0;
    int tpchits = 0;
    int ftdhits = 0;
    streamlog_out(DEBUG) << std::fixed << std::setw( 6 ) << std::setprecision( 3 ) << std::setfill( ' ' );
    //int id =  particle->getTracks()[0]->getSubdetectorHitNumbers()[4];
    if (abs(particle->getCharge()) > 0.5) 
      {
	vtxhits =  particle->getTracks()[0]->getSubdetectorHitNumbers()[0];
	tpchits =  particle->getTracks()[0]->getSubdetectorHitNumbers()[6];
	ftdhits =  particle->getTracks()[0]->getSubdetectorHitNumbers()[5];
      }
    vector<float> direction = MathOperator::getDirection(particle->getMomentum());
    float costheta = std::cos(MathOperator::getAngles(direction)[1]);
    //float chi2 = (particle->getTracks().size() > 0)? particle->getTracks()[0]->getChi2() / (float)particle->getTracks()[0]->getNdf() : -1.0;
    streamlog_out(DEBUG)<<"|"<<vtxhits << ":" <<ftdhits<<":" <<tpchits <<"\t\t|"<<particle->getMass()<<"\t\t|"<<particle->getCharge() <<"\t\t|"<< costheta <<"\t\t|"<<particle->getEnergy() <<"\t\t|\n";
  }
  void VertexRestorer::processEvent( LCEvent * evt ) 
  { 
    _nEvt ++ ;
		  
    try 
      {
	streamlog_out(DEBUG) << "********VertexRecovery*"<<_nEvt << "**********\n";
	TrackOperator opera;
	//opera.test();
	LCCollection* trackrelcol = evt->getCollection( "MarlinTrkTracksMCTruthLink" );
	LCCollection* pricol = evt->getCollection( _colPriName );
	Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
	myPrimary = primary;
	LCCollection* v0col = evt->getCollection( _colV0Name );
	LCCollection* maincol = evt->getCollection( _colName );
	LCCollection* seccol = evt->getCollection( _colSecName );
	LCCollection* relation = evt->getCollection(_colRelName);
	LCCollection* jetcol = evt->getCollection( _colJetName );
	LCCollection* relcol = evt->getCollection( _colJetRelName );
	//LCCollection* genvertexcol = evt->getCollection( _colMCvertexName );


				
	RecoveryOperator recovery(primary, maincol,_Bfield);
	//streamlog_out(DEBUG) << "HERE!";
	LCCollection* trashcol = (_useTracks)? evt->getCollection( _colTrashCanName ) : NULL;
	/*vector< Vertex * > newvertices = recovery.RecoverBuildVertices(seccol,trashcol);
	  IMPL::LCCollectionVec * newseccol = new IMPL::LCCollectionVec ( EVENT::LCIO::VERTEX ) ;
	  IMPL::LCCollectionVec * newvtxRPcol = new IMPL::LCCollectionVec ( EVENT::LCIO::RECONSTRUCTEDPARTICLE ) ;
	  for (int i = 0; i < newvertices.size(); i++) 
	  {
	  newseccol->addElement(newvertices[i]);
	  newvtxRPcol->addElement(newvertices[i]->getAssociatedParticle());
	  }
	  CopyParameters(seccol->parameters (), newseccol->parameters ());
	  evt->addCollection( newseccol, _outputcolName ) ;
	  CopyParameters(secRPcol->parameters (), newvtxRPcol->parameters ());
	  evt->addCollection( newvtxRPcol, _outputjetrelRPcolName ) ;
	*/
	//vector < ReconstructedParticle * > * result = RestoreVerticesPFO(maincol, seccol);
				
	//LCCollection* trashcol = evt->getCollection( _colTrashCanName );
	//IMPL::LCCollectionVec * newjetcol = new IMPL::LCCollectionVec ( (const IMPL::LCCollectionVec&)(*jetcol) ) ;
	std::cout<<_Test<<std::endl;
	if (!_Test) 
	  {
	    LCCollection* secRPcol = evt->getCollection( _colSecRPName );
	    IMPL::LCCollectionVec * newjetcol = new IMPL::LCCollectionVec ( EVENT::LCIO::RECONSTRUCTEDPARTICLE ) ;
	    newjetcol->setSubset();
	    IMPL::LCCollectionVec * newjetrelcol = new IMPL::LCCollectionVec ( EVENT::LCIO::LCRELATION ) ;
	    IMPL::LCCollectionVec * newseccol = new IMPL::LCCollectionVec ( EVENT::LCIO::VERTEX ) ;
	    IMPL::LCCollectionVec * newvtxRPcol = new IMPL::LCCollectionVec ( EVENT::LCIO::RECONSTRUCTEDPARTICLE ) ;
	    vector< Vertex * > newvertices = recovery.RecoverJetVertices(jetcol,relcol,seccol,trashcol, newjetrelcol);
	    for (unsigned int i = 0; i < newvertices.size(); i++) 
	      {
		newseccol->addElement(newvertices[i]);
		newvtxRPcol->addElement(newvertices[i]->getAssociatedParticle());
		for (int j = 0; j < newvertices[i]->getAssociatedParticle()->getParticles().size(); j++) 
		  {
		    if (newvertices[i]->getAssociatedParticle()->getParticles()[j]->getClusters().size() > 0) 
		      {
			continue;
		      }
		    //streamlog_out(DEBUG) << "HERE!!!\n" ;
		    //newvtxRPcol->addElement(newvertices[i]->getAssociatedParticle()->getParticles()[j]);
		  }
	      }
	    LCRelationNavigator navigator(newjetrelcol);
	    CopyParameters(jetcol->parameters (), newjetcol->parameters ());
	    PIDHandler pidh(newjetcol);
	    int alid = -1;
	    try
	      {
		alid = pidh.getAlgorithmID("lcfiplus");
	      }
	    catch(UTIL::UnknownAlgorithm &e)
	      {
		streamlog_out(DEBUG) << "No algorithm lcfiplus!\n";
		streamlog_out(DEBUG) << "Jet number: " << jetcol->getNumberOfElements() << "\n";
		alid = -1;
	      }
	    int newalid = (alid < 0)? 0 : pidh.addAlgorithm("vtxrec", pidh.getParameterNames (alid));
	    for (int i = 0; i < jetcol->getNumberOfElements(); i++) 
	      {
		newjetcol->addElement(jetcol->getElementAt(i));
		ReconstructedParticle * jetpart = dynamic_cast< ReconstructedParticle * >(newjetcol->getElementAt(i));
		vector< Vertex * > * vertices = convert(navigator.getRelatedToObjects(jetpart));
		vector<float> params;
		float btag = 0.;
		if (alid > -1) 
		  {
		    const ParticleID& pid =  pidh.getParticleID(jetpart,alid);
		    params = pid.getParameters();
		    btag = params[pidh.getParameterIndex(alid,"BTag")];
		  }
		bool changed = false;
		float bmomentum = 0.0;
		for (int j = 0; j < vertices->size(); j++) 
		  {
		    bmomentum += MathOperator::getModule(vertices->at(j)->getAssociatedParticle()->getMomentum());
		    if (vertices->at(j)->getParameters()[0] > 0.0) 
		      {
			changed = true;
		      }
							
		  }
		float btag2 =  MathOperator::getRandom(btag , 1.0);//(btag + (1.0 - btag)/2.0
		if (bmomentum > 15 && btag > 0.6 && changed) 
		  {
		    params[pidh.getParameterIndex(alid,"BTag")] = btag2;
		  }
		streamlog_out(DEBUG) << "BTAG: " << btag << " BTAG NEW: " << btag2 << " p: " << bmomentum  << "\n";
		pidh.setParticleID(jetpart, 42, 9999, 0.0, newalid, params);
	      }
	    newjetcol->setSubset( true);
	    evt->addCollection( newjetcol, _outputjetcolName ) ;
	    CopyParameters(seccol->parameters (), newseccol->parameters ());
	    evt->addCollection( newseccol, _outputcolName ) ;
	    CopyParameters(relcol->parameters (), newjetrelcol->parameters ());
	    evt->addCollection( newjetrelcol, _outputjetrelcolName ) ;
	    CopyParameters(secRPcol->parameters (), newvtxRPcol->parameters ());
	    evt->addCollection( newvtxRPcol, _outputjetrelRPcolName ) ;
	    _totaltracks += recovery.GetStatistics();
	    streamlog_out(DEBUG) << "Total added tracks: " << _totaltracks << "\n";
	  }
	else 
	  {

	    std::cout<<" save histo "<<std::endl;                                                                                                
	    ///*/
	    //vector < MyReconstructedParticle * > * result = RefineVertices(maincol, seccol, v0col);
	    vector < MyReconstructedParticle * > * result = RestoreVertices(jetcol, relcol, seccol, trashcol);
	    //vector < MyReconstructedParticle * > * result =__CheatRestoreVertices(jetcol, genvertexcol);
	    LCCollection* bstarcol = evt->getCollection( _colBStarName );
	    CompareCollectionsRel(result, evt->getCollection(_colEGPName), bstarcol , trackrelcol); // relation
	    _Tree->Fill();
	    //TestMethod2(evt->getCollection("MissedParticles"), seccol, evt->getCollection(_colRelName));
	    //_TestTree->Fill();
				
	    WritePrimaries(opera, pricol, evt->getCollection(_colEGPName), relation);
	    //WriteZombies(trashcol, maincol, evt->getCollection("MCParticlesSkimmed"), relation, trackrelcol);
	    WritePurgatory(pricol, evt->getCollection( _colFinalVertex ), maincol, evt->getCollection(_colEGPName), relation, evt->getCollection(_MCcolName));
	    _primaryTree->Fill();
	    _purgatoryTree->Fill();
	    _trashTree->Fill();//*/
	  }
	ClearVariables();
      }
    catch( DataNotAvailableException &e)
      {
	streamlog_out(DEBUG) << "Whoops!....\n";
	streamlog_out(DEBUG) << e.what();
      }
		    
  }
  void VertexRestorer::ClearVariables()
  {
    _primariesTotal = 0;
    _purgatoryTotal = 0;
    for (int i = 0; i < MAXP; i++) 
      {
	_primeOffset[i] = -1.0;
	_primeType[i] = 0;
	_primeDeviation[i] = -1.0;
	_primeMomentum[i] = -1.0;
	_primeError[i] = -1.0;
	_primeTeorError[i] = -1.0;
	_primeChi2[i] = -1.0;
	_primeCostheta[i] = -2.0;
	_primeVtxHits[i] = -1.0;
	_primeFtdHits[i] = -1.0;
      }
    /*for (int i = 0; i < MAXN; i++) 
      {
      _missedAngle[i] = -1.0;
      _fakeAngle[i] = -1.0;
      }*/
  }
  void VertexRestorer::WriteZombies(LCCollection * trashcol, LCCollection * pfos, LCCollection * mccol, LCCollection * rel, LCCollection * trackrel)
  {
    LCRelationNavigator navigator(rel);
    _trashTotal = trashcol->getNumberOfElements();
    vector<MCParticle *> taken;
    for (int i = 0; i < _trashTotal; i++) 
      {
	Track * trashtrack = dynamic_cast< Track * >(trashcol->getElementAt(i)); 
	_trashVtxHits[i] = trashtrack->getSubdetectorHitNumbers()[0];
	_trashTpcHits[i] = trashtrack->getSubdetectorHitNumbers()[6];
	ReconstructedParticle * trashreco = myTrackOperator.ReconstructParticle(trashtrack,_Bfield);
	PrintParticle(trashreco);
	MCParticle * trashmc = CompareToCollection(trashreco, mccol, trackrel, true);
	if (!trashmc) 
	  {
	    streamlog_out(DEBUG) << "No mc particle!!\n";
	    _trashHasMCParticle[i] = 0;
	    continue;
	  }
	_trashIsDublicate[i] = 0;
	for (int j = 0; j < taken.size(); j++) 
	  {
	    if (taken[j] == trashmc) 
	      {
		streamlog_out(DEBUG) << "Particle is dublicate\n";
		_trashIsDublicate[i] = 1;
	      }
	  }
	taken.push_back(trashmc);
	_trashHasMCParticle[i] = 1;
	streamlog_out(DEBUG) << "MC p: " << MathOperator::getModule(trashmc->getMomentum())
		  << " q: " << trashmc->getCharge()
		  << "\n";
	vector<float> direction = MathOperator::getDirection(trashmc->getMomentum());
	_trashCostheta[i] = std::cos( MathOperator::getAngles(direction)[1]);
	vector< LCObject * > obj = navigator.getRelatedFromObjects(trashmc);
	vector< float > weights = navigator.getRelatedFromWeights(trashmc);
	float maxweight = 0.5;
	ReconstructedParticle * trashreco2 = NULL;
	for (int j = 0; j < obj.size(); j++) 
	  {
	    if (weights[j] > maxweight) 
	      {
		maxweight = weights[j];
		trashreco2 = dynamic_cast< ReconstructedParticle * >(obj[j]);
	      }
	  }
	if (trashreco2) 
	  {
	    streamlog_out(DEBUG) << "Reco particle found:\n";
	    _trashIsReco[i] = 1;
	    PrintParticle(trashreco2);
	  }
	else 
	  {
	    _trashIsReco[i] = 0;
	    streamlog_out(DEBUG) << "Reco particle NOT found!\n";
	  }
      }
    streamlog_out(DEBUG) << "Total: " << _trashTotal << "\n";
  }

  void VertexRestorer::WritePrimaries(TrackOperator & opera, EVENT::LCCollection * pricol, EVENT::LCCollection * sec, EVENT::LCCollection * rel)
  {
    Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
    double * pos = MathOperator::toDoubleArray(primary->getPosition(), 3);

    vector < ReconstructedParticle * > primaries = primary->getAssociatedParticle()->getParticles();
    _primariesTotal = primaries.size();
    _primeT = sqrt(primary->getPosition()[0] * primary->getPosition()[0] + primary->getPosition()[1] * primary->getPosition()[1]);
    _primeZ= primary->getPosition()[2];
    for (unsigned int i = 0; i < primaries.size(); i++)
      {
	//opera.test(primaries[i]);
	if (sec && rel && CompareToCollection(primaries[i], sec,rel)) 
	  {
	    streamlog_out(DEBUG) << "EXCLUDED!\n";
	    continue;
	  }
	vector<float> direction = MathOperator::getDirection(primaries[i]->getMomentum());
	double * start = opera.GetStartPoint(primaries[i]);
	//PrintParticle(primaries[i]);
	_primeType[i] = primaries[i]->getType();
	_primeOffset[i] = MathOperator::getDistanceTo(pos, direction, start);//opera.GetOffset(primaries[i]);
	//_primeOffset[i] = opera.GetOffset(primaries[i]);
	_primeMomentum[i] = MathOperator::getModule(primaries[i]->getMomentum());
				
	_primeError[i] = opera.GetOffsetErrorSimple(primaries[i]);
	_primeTeorError[i] = GetError(primaries[i]);
	_primeChi2[i] = primaries[i]->getTracks()[0]->getChi2()/(float)primaries[i]->getTracks()[0]->getNdf();
	//_primeChi2[i] = primaries[i]->getTracks()[0]->getChi2();
	_primeCostheta[i] = std::cos(MathOperator::getAngles(direction)[1]);
	_primeVtxHits[i] = primaries[i]->getTracks()[0]->getSubdetectorHitNumbers()[0];
	_primeFtdHits[i] = primaries[i]->getTracks()[0]->getSubdetectorHitNumbers()[5];
	_primeEtdHits[i] = primaries[i]->getTracks()[0]->getSubdetectorHitNumbers()[10];
	//streamlog_out(DEBUG) << "Offset: " << _primeOffset[i] << "; error^2: " << _primeError[i] << '\n'; 
	_primeDeviation[i] = myTrackOperator.GetOffsetSignificance(primaries[i],pos); //_primeOffset[i] / _primeError[i]; 
	//streamlog_out(DEBUG) << "Offset: " << _primeOffset[i] << " deviation: " << _primeDeviation[i] << '\n';
				
      }

  }
  MCParticle * VertexRestorer::GetInterestingParent(MCParticle * daughter)
  {
    if (!daughter || daughter->getParents().size() < 1) 
      {
	//streamlog_out(DEBUG) << "\tno parents\n";
	return NULL;
      }
    MCParticle * parent = daughter->getParents()[0];
    int parentpdg = abs(parent->getPDG());
    if ((parentpdg > 200 && parentpdg < 300) ||
	(parentpdg > 10100 && parentpdg < 10225) ||
	(parentpdg > 110 && parentpdg < 120 )|| 
	parentpdg == 2212 || 
	parentpdg == 2112 || 
	parentpdg == 11 ||
	parentpdg == 13) 
      {
	//streamlog_out(DEBUG) << "\tcontinue to search with " << parentpdg <<'\n';
	return GetInterestingParent(parent);
      }
    //streamlog_out(DEBUG) << "Satisfied with " << parentpdg <<'\n';
    return parent;
  }
  void VertexRestorer::WritePurgatory(LCCollection * pricol, LCCollection * sec, LCCollection * pfos, LCCollection *egprongs, LCCollection *rel, LCCollection *mccol)
  {
    Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
    double * pos = MathOperator::toDoubleArray(primary->getPosition(), 3);
    vector < ReconstructedParticle * > primaries = primary->getAssociatedParticle()->getParticles();
    vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
    particles->reserve(particles->size() + primaries.size());
    particles->insert(particles->end(),primaries.begin(),primaries.end());
    int snumber = sec->getNumberOfElements();
    //int v0number = v0col->getNumberOfElements();
    for (int i = 0; i < snumber; i++) 
      {
	Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
	vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
	particles->reserve(particles->size() + secondaries.size());
	particles->insert(particles->end(),secondaries.begin(),secondaries.end());
      }
    int pfocharge = 0;
    for (int i = 0; i < pfos->getNumberOfElements(); i++) 
      {
	ReconstructedParticle * purparticle = dynamic_cast< ReconstructedParticle * > (pfos->getElementAt(i));
	if (std::abs(purparticle->getCharge()) < 0.9) 
	  {
	    continue;
	  }
	pfocharge++;
	if (IsDublicate(purparticle, *particles)) 
	  {
	    continue;
	  }
	if (sec && rel && CompareToCollection(purparticle, egprongs,rel)) 
	  {
	    streamlog_out(DEBUG) << "EXCLUDED2!\n";
	    continue;
	  }
	MCParticle * mcp = CompareToCollection(purparticle, mccol,rel);
	if (mcp) 
	  {
	    //streamlog_out(DEBUG) << "Start looking for parent of " << abs(mcp->getPDG()) << "\n";
	    MCParticle * parent = GetInterestingParent(mcp);
	    _purgatoryParent[_purgatoryTotal] = (parent)? abs(parent->getPDG()): 0; 
	    if (abs(mcp->getPDG()) == 13 && mcp->getParents().size() > 0 && _purgatoryParent[_purgatoryTotal] == 92 ) 
	      {
		MCParticle * parent = mcp->getParents()[0];
		if (abs(parent->getPDG()) == 211) 
		  {
		    _purgatoryParent[_purgatoryTotal] = 50000;
		  }
	      }
	  }
	else 
	  {
	    _purgatoryParent[_purgatoryTotal] = 0;
	  }

	vector<float> direction = MathOperator::getDirection(purparticle->getMomentum());
	double * start = myTrackOperator.GetStartPoint(purparticle);
	_purgatoryOffset[_purgatoryTotal] = MathOperator::getDistanceTo(pos, direction, start);//opera.GetOffset(primaries[i]);
	_purgatoryMomentum[_purgatoryTotal] = MathOperator::getModule(purparticle->getMomentum());
	_purgatoryTeorError[_purgatoryTotal] = GetError(purparticle);
	_purgatoryError[_purgatoryTotal] = std::sqrt(myTrackOperator.GetOffsetError(purparticle, myTrackOperator.GetStartPoint(purparticle), primary, _purgatoryOffset[_purgatoryTotal]));
	_purgatoryCostheta[_purgatoryTotal] = std::cos(MathOperator::getAngles(direction)[1]);
	_purgatoryVtxHits[_purgatoryTotal] = purparticle->getTracks()[0]->getSubdetectorHitNumbers()[0];
	_purgatoryFtdHits[_purgatoryTotal] = purparticle->getTracks()[0]->getSubdetectorHitNumbers()[5];
	_purgatoryDeviation[_purgatoryTotal] = _purgatoryOffset[_purgatoryTotal] / _purgatoryError[_purgatoryTotal]; 
	if ( purparticle->getStartVertex()) 
	  {
	    //streamlog_out(DEBUG) << "Particle has start vertex  "<< purparticle->getStartVertex()->isPrimary() << "\n";
	  }
	else 
	  {
	    //streamlog_out(DEBUG) << "PRIMARY PARTICLE HAS NO VERTEX!!!\n";
	  }
	_purgatoryTotal++;

      }
    streamlog_out(DEBUG) << "Primaries: " << primaries.size() 
	      << " primaries+secondaries: " << particles->size() 
	      << " purgatories: " << _purgatoryTotal
	      << " PFOs: " << pfocharge 
	      << "\n";

  }
  MCParticle * VertexRestorer::CompareToCollection(EVENT::ReconstructedParticle * particle, EVENT::LCCollection * prongs, EVENT::LCCollection * rel, bool byTrack)
  {
    LCRelationNavigator navigator(rel);
    vector< LCObject * > obj;// = navigator.getRelatedToObjects((byTrack)? particle->getTracks()[0]: particle);
    vector< float > weights;// = navigator.getRelatedToWeights((byTrack)? particle->getTracks()[0]: particle);
    if (byTrack) 
      {
	obj = navigator.getRelatedToObjects(particle->getTracks()[0]);
	weights = navigator.getRelatedToWeights(particle->getTracks()[0]);
      }
    else 
      {
	obj = navigator.getRelatedToObjects(particle);
	weights = navigator.getRelatedToWeights(particle);
      }
    if (obj.size() < 1)
      {
	streamlog_out(DEBUG) << "Particle was not generated!\n";
	return false;
      }
    MCParticle * winner = NULL;
    float maxweight = 0.90;
    for (unsigned int i = 0; i < obj.size(); i++) 
      {
	if (weights[i] > maxweight) 
	  {
	    winner = dynamic_cast< MCParticle * >(obj[i]);
	    maxweight = weights[i];
	  }
      }
    int numberOfProngs = prongs->getNumberOfElements();
    MCParticle * mcprong = NULL;
    for (int i = 0; i < numberOfProngs; i++) 
      {
	mcprong = dynamic_cast< MCParticle * >(prongs->getElementAt(i));
	if (mcprong == winner) 
	  {
	    return mcprong;
	  }
      }
    return NULL;
  }
  void VertexRestorer::CompareCollectionsRel(std::vector< MyReconstructedParticle * > * detected, EVENT::LCCollection * prongs, EVENT::LCCollection * bstar, EVENT::LCCollection *  rel)
  {
    _detectedTotal = detected->size();
		double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
    _missedTotal = prongs->getNumberOfElements();
    for (int i = 0; i < _missedTotal; i++) 
      {
	MCParticle * particle = dynamic_cast< MCParticle * >(prongs->getElementAt(i));
	_allmissedMomentum[i] = MathOperator::getModule(particle->getMomentum());
      }
    _bstarTotal = bstar->getNumberOfElements();
    _missedDetected = 0;
    _bstarDetected = 0;
    _fakeDetected = 0;
    _purDetected = 0;
    bool byTracks = true;
    for (int i = 0; i < _detectedTotal; i++) 
      {
	float taken = false;
	ReconstructedParticle * particle = detected->at(i)->Get();
	if (CompareToCollection(particle, prongs, rel, byTracks)) 
	  {
	    taken = true;
	    streamlog_out(DEBUG) << "Found a prong!\n";
	    //_missedError[_missedDetected] = GetError(particle);
	    _missedError[_missedDetected] = detected->at(i)->GetAccuracy();
	    _missedMomentum[_missedDetected] = MathOperator::getModule(particle->getMomentum());
	    _missedDeviation[_missedDetected] =  myTrackOperator.GetOffsetSignificance(particle,primaryPosition);// detected->at(i)->GetOffset() / _missedError[_missedDetected];
	    _missedD0Deviation[_missedDetected] =  std::abs(particle->getTracks()[0]->getD0() / std::sqrt(particle->getTracks()[0]->getCovMatrix()[0])); //
	    _missedZ0Deviation[_missedDetected] =  std::abs(particle->getTracks()[0]->getZ0() / std::sqrt(particle->getTracks()[0]->getCovMatrix()[9])); //
	    _missedZ0[_missedDetected] =  particle->getTracks()[0]->getZ0(); //
	    _missedD0[_missedDetected] =  particle->getTracks()[0]->getD0(); 
	    _missedZ0Error[_missedDetected] =  particle->getTracks()[0]->getCovMatrix()[9];
	    _missedD0Error[_missedDetected] =  particle->getTracks()[0]->getCovMatrix()[0];
	    _missedOffset[_missedDetected] =  detected->at(i)->GetOffset() ;
	    _missedAngle[_missedDetected] = detected->at(i)->GetAngle();
	    _missedObs[_missedDetected] = detected->at(i)->GetObservable();
	    _missedZ[_missedDetected] = detected->at(i)->GetZ();
	    _missedVtxHits[_missedDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[0];
	    _missedFtdHits[_missedDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[4];
	    _missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset();
	    _missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
	    _missedCosthetaVtx[_missedDetected] = detected->at(i)->GetCosthetaVtx();
	    _missedDetected++;
	  }
	if (CompareToCollection(particle, bstar, rel, byTracks)) 
	  {
	    taken = true;
	    _bstarError[_bstarDetected] = GetError(particle);
	    _bstarMomentum[_bstarDetected] = MathOperator::getModule(particle->getMomentum());
	    _bstarDetected++;
	  }
	if (!taken) 
	  {
	    _fakeError[_fakeDetected] = detected->at(i)->GetAccuracy();
	    //_fakeError[_fakeDetected] = GetError(particle);
	    _fakeMomentum[_fakeDetected] = MathOperator::getModule(particle->getMomentum());
	    _fakeDeviation[_fakeDetected] =   myTrackOperator.GetOffsetSignificance(particle,primaryPosition); //detected->at(i)->GetOffset() / _fakeError[_fakeDetected];
	    _fakeD0Deviation[_fakeDetected] =  std::abs(particle->getTracks()[0]->getD0() / std::sqrt(particle->getTracks()[0]->getCovMatrix()[0])); //
	    _fakeZ0Deviation[_fakeDetected] =  std::abs(particle->getTracks()[0]->getZ0() / std::sqrt(particle->getTracks()[0]->getCovMatrix()[9])); //
	    _fakeZ0[_fakeDetected] =  particle->getTracks()[0]->getZ0(); //
	    _fakeD0[_fakeDetected] =  particle->getTracks()[0]->getD0();               
	    _fakeZ0Error[_fakeDetected] =  particle->getTracks()[0]->getCovMatrix()[9];                                                                      
	    _fakeD0Error[_fakeDetected] =  particle->getTracks()[0]->getCovMatrix()[0];
	    _fakeSecDeviation[_fakeDetected] =  detected->at(i)->GetSecOffset();
	    _fakeAngle[_fakeDetected] = detected->at(i)->GetAngle(); 
	    _fakeObs[_fakeDetected] = detected->at(i)->GetObservable(); 
	    _fakeZ[_fakeDetected] = detected->at(i)->GetZ(); 
	    _fakeOffset[_fakeDetected] = detected->at(i)->GetOffset(); 
	    _fakeCostheta[_fakeDetected] = detected->at(i)->GetCostheta(); 
	    _fakeCosthetaVtx[_fakeDetected] = detected->at(i)->GetCosthetaVtx(); 
	    _fakeVtxHits[_fakeDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[0];
	    _fakeFtdHits[_fakeDetected] = particle->getTracks()[0]->getSubdetectorHitNumbers()[4];
	    _fakeDetected++;
	  }
	if (!taken && !particle->getStartVertex())//&& !IsDublicate(particle, myPrimary->getAssociatedParticle()->getParticles())) 
	  {
	    _purMomentum[_purDetected] = MathOperator::getModule(particle->getMomentum());
	    _purOffset[_purDetected] =  detected->at(i)->GetOffset();
	    _purAngle[_purDetected] =  detected->at(i)->GetAngle();
	    _purCostheta[_purDetected] =  detected->at(i)->GetCostheta();
	    streamlog_out(DEBUG) << "Purgatory detected!\n";
	    _purDetected++;
	  }
      }
    streamlog_out(DEBUG) << "Total detected: " << _detectedTotal
	      << "; missed total: " << _missedTotal
	      << "; missed detected: " << _missedDetected
	      << "; B* detected: " << _bstarDetected
	      << '\n';
  }

  void VertexRestorer::CompareCollections(vector< MyReconstructedParticle * > * detected, LCCollection * missed, LCCollection * bstar, LCCollection * missedvtx)
  {
    _missedTotal = missed->getNumberOfElements();
    for (int i = 0; i < _missedTotal; i++) 
      {
	ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * >(missed->getElementAt(i));
	_allmissedMomentum[i] = MathOperator::getModule(particle->getMomentum());
      }
    _bstarTotal = bstar->getNumberOfElements();
    _detectedTotal = detected->size();

    _missedDetected = 0;
    _bstarDetected = 0;
    _fakeDetected = 0;

    for (int i = 0; i < _detectedTotal; i++) 
      {
	float taken = false;
	ReconstructedParticle * detectedParticle = detected->at(i)->Get();
	for (int j = 0; j < _missedTotal; j++) 
	  {
	    ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * > (missed->getElementAt(j));
	    if (CompareParticles(particle, detectedParticle)) 
	      {
		taken = true;
		_missedError[_missedDetected] = detected->at(i)->GetAccuracy();
		_missedMomentum[_missedDetected] = MathOperator::getModule(particle->getMomentum());
		_missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() / _missedError[_missedDetected];
		_missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset();
		_missedAngle[_missedDetected] = detected->at(i)->GetAngle();
		_missedObs[_missedDetected] = detected->at(i)->GetObservable();
		_missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
		_missedDetected++;
	      }
	  }
	if (missedvtx && !taken)
	  {
	    int missedvtxTotal = missedvtx->getNumberOfElements();
	    for (int j = 0; j < missedvtxTotal; j++) 
	      {
		ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * > (missedvtx->getElementAt(j));
		if (CompareParticles(particle, detectedParticle)) 
		  {
		    taken = true;
		    //_missedDetected++;
		    _missedError[_missedDetected] = detected->at(i)->GetAccuracy();
		    //_missedError[_missedDetected] = GetError(detectedParticle);
		    _missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() / _missedError[_missedDetected];
		    _missedAngle[_missedDetected] = detected->at(i)->GetAngle();
		    _missedObs[_missedDetected] = detected->at(i)->GetObservable();
		    _missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
		    _missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset();
		    _missedMomentum[_missedDetected++] = MathOperator::getModule(particle->getMomentum());
		  }
	      }
	  }
	for (int j = 0; j < _bstarTotal; j++) 
	  {
	    MCParticle * particle = dynamic_cast< MCParticle * > (bstar->getElementAt(j));
	    if (CompareParticles(particle, detectedParticle)) 
	      {
		_bstarError[_bstarDetected] = GetError(detectedParticle);
		_bstarMomentum[_bstarDetected] = MathOperator::getModule(particle->getMomentum());
		//_bstarAngle[_bstarDetected] =  detectedParticle->getParticleIDs()[0]->getLikelihood();
		_bstarDetected++;
		//_bstarDeviation[_bstarDetected] = GetDeviation(detectedParticle, MathOperator::toDoubleArray(detectedParticle->getStartVertex()->getPosition(), 3));
		taken = true;
	      }
	  }
	if (!taken) 
	  {
	    _fakeError[_fakeDetected] = detected->at(i)->GetAccuracy();
	    _fakeDeviation[_fakeDetected] =  detected->at(i)->GetOffset() / _fakeError[_fakeDetected];
	    _fakeSecDeviation[_fakeDetected] =  detected->at(i)->GetSecOffset();
	    _fakeAngle[_fakeDetected] = detected->at(i)->GetAngle(); 
	    _fakeObs[_fakeDetected] = detected->at(i)->GetObservable(); 
	    _fakeOffset[_fakeDetected] = detected->at(i)->GetOffset(); 
	    _fakeCostheta[_fakeDetected] = detected->at(i)->GetCostheta(); 
	    _fakeMomentum[_fakeDetected++] = MathOperator::getModule(detectedParticle->getMomentum());
	  }
      }
    //_fakeDetected = _detectedTotal - _missedDetected - _bstarDetected;
    streamlog_out(DEBUG) << "Total detected: " << _detectedTotal
	      << "; missed total: " << _missedTotal
	      << "; missed detected: " << _missedDetected
	      << "; B* detected: " << _bstarDetected
	      << '\n';

  }

  vector< MyReconstructedParticle * > * VertexRestorer::RefineVertices(LCCollection * main, LCCollection * sec, LCCollection * v0col)
  {
    vector < MyReconstructedParticle * > * result = new vector< MyReconstructedParticle * >();
    vector< ReconstructedParticle * > primaries;
    int mnumber = main->getNumberOfElements();
    for (int i = 0; i < mnumber; i++) 
      {
	primaries.push_back(dynamic_cast< ReconstructedParticle * > (main->getElementAt(i)));
      }
    vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
    int snumber = sec->getNumberOfElements();
    int v0number = v0col->getNumberOfElements();
    for (int i = 0; i < snumber; i++) 
      {
	Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
	vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
	particles->reserve(particles->size() + secondaries.size());
	particles->insert(particles->end(),secondaries.begin(),secondaries.end());
      }
    for (int i = 0; i < v0number; i++) 
      {
	Vertex * secondary = dynamic_cast< Vertex * >(v0col->getElementAt(i));
	vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
	particles->reserve(particles->size() + secondaries.size());
	//particles->insert(particles->end(),secondaries.begin(),secondaries.end());
      }
    for (int i = 0; i < snumber; i++) 
      {
	Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
	//ParametrizeVertex(secondary);
	vector< MyReconstructedParticle * > additional = AddParticles(primaries, secondary, particles);
	for (unsigned int j = 0; j < additional.size(); j++) 
	  {
	    if (!IsDublicate(additional[j], *result)) 
	      {
		result->push_back(additional[j]);
	      }
	  }
      }
    return result;
  }
  vector< float > VertexRestorer::ParametrizeVertex(const Vertex * sec)
  {
    const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
    vector< float > result;
    if (secondaries.size() < 2) 
      {
	return result;
      }
    /*float mindprime = 1000;
      float maxdprime = 0;
      double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
      //streamlog_out(DEBUG) << "Vertex position: " << MathOperator::getModule(sec->getPosition()) << " Chi2: " << sec->getChi2() << "\n";
      for (int i = 0; i < secondaries.size(); i++) 
      {
      double * secPos = myTrackOperator.GetStartPoint(secondaries[i]);
      vector< float > secDir = MathOperator::getDirection(secondaries[i]->getMomentum());
      double * vtxPos = MathOperator::toDoubleArray(sec->getPosition(),3);
      vector<float> diff = MathOperator::getDirection(vtxPos, secPos);
      double diif[3];
      for (int m = 0; m < 3; m++) 
      {
      diif[m] = diff[m];
      }
      float distance = MathOperator::getAngle(diif, secondaries[i]->getMomentum());
      if (distance < mindprime) 
      {
      mindprime = distance;
      }
      if (distance > maxdprime && distance < 2.0) 
      {
      maxdprime = distance;
      }

      }
      //streamlog_out(DEBUG) << "\tdmin'= " << mindprime << "\tdmax'= " << maxdprime << "\n";
      result.push_back(mindprime);
      result.push_back(maxdprime);
      return result;*/
    float meanOffset = 0.0;
    float sigmaOffset = 0.0;
    float meanAngle = 0.0;
    float sigmaAngle = 0.0;
    float meanVtx = 0.0;
    float sigmaVtx = 0.;
    float meanObs = 0.0;
    float sigmaObs = 0.0;
    vector<float> offsets;
    vector<float> angles;
    vector<float> observables;
    for (int i = 0; i < secondaries.size(); i++) 
      {
	ReconstructedParticle * primary = secondaries[i];
	vector< float > direction = MathOperator::getDirection(primary->getMomentum());
	double * trackPosition = myTrackOperator.GetStartPoint(primary);
	double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
	double * secondaryPosition = MathOperator::toDoubleArray(sec->getPosition(),3);
	float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
	vector<float> diff = MathOperator::getDirection(secondaryPosition, trackPosition);
	double diif[3];
	for (int m = 0; m < 3; m++) 
	  {
	    diif[m] = diff[m];
	  }
	float angle = MathOperator::getAngle(diif, primary->getMomentum());
	int vtxhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[0];
	//streamlog_out(DEBUG) << "\t\toffset: " << primaryOffset / GetError(primary) << "\n";
	streamlog_out(DEBUG) << "\t\tangle: " << angle << "\n";
	float momentum = MathOperator::getModule(primary->getMomentum());
	meanVtx += vtxhits;
	meanAngle += angle;
	angles.push_back(angle);
	offsets.push_back(primaryOffset / GetError(primary));
	meanOffset += primaryOffset / GetError(primary);
	//observables.push_back(primaryOffset / GetError(primary) / angle);
	//meanObs += primaryOffset / GetError(primary) / angle;
	observables.push_back(momentum);
	meanObs += momentum;

      }
    meanObs/= secondaries.size();
    meanVtx /= secondaries.size();
    meanAngle /= secondaries.size();
    meanOffset /= secondaries.size();
    for (int i = 0; i < secondaries.size(); i++) 
      {
	ReconstructedParticle * primary = secondaries[i];
	int vtxhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[0];
	sigmaOffset += (meanOffset - offsets[i])*(meanOffset - offsets[i]);
	sigmaAngle  += std::pow(meanAngle - angles[i], 2);
	sigmaVtx += std::pow(vtxhits - meanVtx, 2);
	sigmaObs += std::pow( observables[i] - meanObs, 2);
      }
    sigmaOffset = std::sqrt(sigmaOffset / secondaries.size());
    sigmaAngle = std::sqrt(sigmaAngle / secondaries.size());
    sigmaVtx = std::sqrt(sigmaVtx / secondaries.size());
    sigmaObs = std::sqrt(sigmaObs / secondaries.size());

    result.push_back(meanOffset);
    result.push_back(sigmaOffset);
    result.push_back(meanAngle);
    result.push_back(sigmaAngle);
    result.push_back(meanVtx);
    result.push_back(sigmaVtx);
    result.push_back(meanObs);
    result.push_back(sigmaObs);
    return result;
  }
  float VertexRestorer::GetMinDiffDistance(const ReconstructedParticle * candidate, const Vertex * sec, float & obs)
  {
    const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
    double * candidatePos = myTrackOperator.GetStartPoint(candidate);
    vector< float > candidateDir = MathOperator::getDirection(candidate->getMomentum());
    float alpha = 0.0;
    float maxdistance = 0.0;
    float mindistance = 100.0;
    double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
    for (unsigned int i = 0; i < secondaries.size(); i++) 
      {
	double * secPos = myTrackOperator.GetStartPoint(secondaries[i]);
	vector< float > secDir = MathOperator::getDirection(secondaries[i]->getMomentum());
	float distance = MathOperator::getDistanceBtw(candidatePos, candidateDir, secPos, secDir);

	//float angle = MathOperator::getAngle(candidate->getMomentum(), secondaries[i]->getMomentum());
	if (distance < mindistance) 
	  {
	    mindistance = distance;
	    obs = myTrackOperator.GetDprime(candidate, secondaries[i], primaryPosition);
	  }
      }

    return mindistance;
  }
  //bool VertexRestorer::TakeParticle(EVENT::ReconstructedParticle * primary, const EVENT::Vertex * sec, vector<double> * parameters)
  MyReconstructedParticle * VertexRestorer::TakeParticle(ReconstructedParticle * primary, const Vertex * sec, ReconstructedParticle * jet)
  {
    if (std::abs(primary->getCharge() ) < 0.9) 
      {
	return NULL;
      }
    const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
    double * trackPosition = myTrackOperator.GetStartPoint(primary);
    float trackDistance = MathOperator::getModule(trackPosition);
    if (trackDistance > 20.0) 
      {
	return NULL;
      }
    //return true;
    //double observable = MathOperator::getModule(primary->getMomentum()) / MathOperator::getModule(sec->getAssociatedParticle()->getMomentum());
    vector< float > direction = MathOperator::getDirection(primary->getMomentum());
    vector< float > directionTrack = MathOperator::getDirection(trackPosition);
    double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
    double * secondaryPosition = MathOperator::toDoubleArray(sec->getPosition(),3);
    float sin = std::sin(MathOperator::getAngles(direction)[1]);
    float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
    //float accuracy = GetError(primary);
    //float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primary, trackPosition, myPrimary, primaryOffset));
    float accuracy = (myTrackOperator.GetOffsetErrorSimple(primary));
    float secondaryOffset = MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
    vector<float> diff = MathOperator::getDirection(secondaryPosition, trackPosition);
    vector< float > directionVtx = MathOperator::getDirection(secondaryPosition);
    double diif[3];
    double aright[3];
    //trackPosition[2] *= sin * sin;
    float primaryOffsetTest = MathOperator::getDistanceTo(ip, direction, trackPosition);

    float observable =0; //MathOperator::getModule(*MathOperator::vectorProduct(directionVtx, direction));//MathOperator::getAngle(primary->getMomentum(), secondaryPosition);
    float cosbeta = std::cos(MathOperator::getAngle(trackPosition, primary->getMomentum()));
    for (int m = 0; m < 3; m++) 
      {
	diif[m] = diff[m];
	aright[m] = cosbeta * MathOperator::getModule(trackPosition) * direction[m] - trackPosition[m];
	observable += aright[m] * directionVtx[m];
      }
    float angle = MathOperator::getAngle(diif, primary->getMomentum());
    float angleError = 1.0;
    if (angle > 0.0) 
      {
	angleError = std::sqrt(myTrackOperator.GetAngleError(angle, sec, primary));
      }
    float secOffset =  MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
    float costheta = std::cos(MathOperator::getAngles(direction)[1]); //std::abs( std::cos(MathOperator::getAngles(secDiraction)[1] ) );
    float costhetaVtx = std::cos(MathOperator::getAngles(directionVtx)[1]); //std::abs( std::cos(MathOperator::getAngles(secDiraction)[1] ) );
    //float dprime = myTrackOperator.GetDprime(primary, secondaries[0], primaryPosition);
    float pvtx = MathOperator::getModule(sec->getAssociatedParticle()->getMomentum());
    //float l =  GetMinDiffDistance(primary, sec, dprime);// std::abs( distance / secondaryOffset -0.5 ) * cos;
    //vector< float > limits = ParametrizeVertex(sec);
    int vtxhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[0];
    int sithits = primary->getTracks()[0]->getSubdetectorHitNumbers()[2];
    int ftdhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[5];
    if ( abs(costheta) > 0.9 ) 
      {
	accuracy = 2*GetError(primary);
      }
    float radius = std::abs(1. / primary->getTracks()[0]->getOmega());
    float phi = 1.571 - primary->getTracks()[0]->getPhi();
    float distance = MathOperator::getModule(secondaryPosition);
    //float observable = MathOperator::getAngle(primary->getMomentum(),sec->getAssociatedParticle()->getMomentum());
    //float observable = std::atan((radius - std::sqrt(radius*radius - distance*distance)) / distance);
    //float observable = (std::pow(secondaryPosition[0] - radius * std::cos(phi) - trackPosition[0],2) + std::pow(secondaryPosition[1] - radius * std::sin(phi) - trackPosition[1],2) - radius * radius) / (radius * radius);
    //float chi2 = getChi2(primary, sec);
    //vector< float > parametersVertex = ParametrizeVertex(sec);
    if (abs(costheta) < 0.9) 
      {
	//return NULL;
      }
    int ntracks = sec->getAssociatedParticle()->getParticles().size();
    float z = MathOperator::getModule(primary->getMomentum()) / MathOperator::getModule(jet->getMomentum());
    float anglecut = 0.9;//(abs(costheta) > 0.9)? 0.05: 0.1;
    float deviation = myTrackOperator.GetOffsetSignificance(primary,primaryPosition);
    observable = - observable / abs(observable) * deviation;  //MathOperator::getModule(aright) / primaryOffsetTest;//trackDistance;//MathOperator::getAngle(aright, primary->getMomentum());
    //float deviation = primaryOffset /accuracy;
    //bool result = //(primaryOffset /accuracy  > 2.5 *std::atan(angle * 100) || angle < 0.005)   &&
    //!(angle > 0.04 && z < 0.06) &&
    //bool result = primaryOffset /accuracy  > 20.0 * sqrt( angle )   &&
    bool result = //(deviation  > 2.0 + 50.*angle || angle < 0.005) &&
      //bool result = (primaryOffset /accuracy  > 0.5 + 90.*angle || angle < 0.005) && // best for dst
      //bool result = //(primaryOffset/accuracy >0.002+ 0.015 *angle / angleError ) &&
      //(((ftdhits > 0) && angle < anglecut/2) || (vtxhits > 1 && angle < anglecut));
      //(ftdhits + vtxhits > 1) &&
      //bool result = //(primaryOffset/accuracy >0.2+ 15 * sqrt(angle - 0.001 * (ntracks - 2)) ||  angle < 0.005) &&
      //(vtxhits > 2 || angle < 0.001) &&
      //(vtxhits> 3 && sithits > -1  || angle < 0.001) && // ||( abs(costheta) < 0.9 && vtxhits+ftdhits > 5 ) || angle < 0.001)
      angle < anglecut;// + 0.03 * sine;
			
    if (result) 
      {
	streamlog_out(DEBUG) << "Found a track with offset " << primaryOffset
		  << ", error " << accuracy//GetError(primary) 
		  << ", Angle " << angle //GetError(primary) 
		  << ", ntracks " << ntracks //distance//GetError(primary) 
		  << ", costheta " << costheta
		  << ", vtxhits " << vtxhits
		  << ", radius " << radius
	  //<< ", ftdhits " << ftdhits
		  <<" :\n";
	MyReconstructedParticle * newparticle = new MyReconstructedParticle(primary);
	newparticle->SetAngle(angle);
	newparticle->SetObservable(observable);
	newparticle->SetAccuracy(accuracy);
	newparticle->SetOffset(primaryOffset);
	newparticle->SetSecOffset(trackPosition[2]);
	newparticle->SetCostheta(costheta);
	newparticle->SetCosthetaVtx(costhetaVtx);
	if (jet) 
	  {
	    newparticle->SetZ(z);
	  }
	return newparticle;
      }
    else 
      {
	return NULL;
      }
  }

  vector< MyReconstructedParticle * > VertexRestorer::AddParticles(const vector< ReconstructedParticle * > & pri, Vertex * sec, const vector< ReconstructedParticle * > * toCompare, ReconstructedParticle * jet)
  {
    vector< MyReconstructedParticle * > result;
    if (sec->getAssociatedParticle()->getParticles().size() < 2) 
      {
	return result;
      }
    const vector< ReconstructedParticle * > * particles = (toCompare)? toCompare : &(sec->getAssociatedParticle()->getParticles());
    for (unsigned int i = 0; i < pri.size(); i++) 
      {
	vector<double> * parameters = new vector<double> ();
	ReconstructedParticle * primary = pri[i];
	MyReconstructedParticle * newparticle = TakeParticle(primary, sec, jet);
	if (newparticle) 
	  {
	    PrintParticle(primary);
	    /*bool found = false;
	      for (unsigned int j = 0; j < particles->size(); j++) 
	      {
	      ReconstructedParticle * secondary = particles->at(j);
	      if (CompareParticles(secondary, primary)) 
	      {
	      found = true;
	      break;
	      }
	      }*/
	    //if (!found) 
	    if (!IsDublicate(primary, *particles))
	      {
		streamlog_out(DEBUG) << "The particle is unique!\n";
		//float errorvtx = (sec->getChi2() < 0.00001)? 0.00001: sec->getChi2();
		result.push_back(newparticle);
	      }
	  }
      }
    return result;
			
  }
  void VertexRestorer::AddParticleID(EVENT::ReconstructedParticle * particle, float theta)
  {
    ParticleIDImpl * pid = new ParticleIDImpl();
    pid->setAlgorithmType(42);
    pid->setLikelihood(theta);
    particle->addParticleID(pid);
  }
  vector< ReconstructedParticle * > * VertexRestorer::RestoreVerticesPFO(LCCollection * main, EVENT::LCCollection * sec)
  {
    vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
    vector< ReconstructedParticle * > primaries;
    int mnumber = main->getNumberOfElements();
    for (int i = 0; i < mnumber; i++) 
      {
	primaries.push_back(dynamic_cast< ReconstructedParticle * > (main->getElementAt(i)));
      }
    int snumber = sec->getNumberOfElements();
    vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
    for (int i = 0; i < snumber; i++) 
      {
	Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
	vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
	particles->reserve(particles->size() + secondaries.size());
	particles->insert(particles->end(),secondaries.begin(),secondaries.end());
      }
    for (int i = 0; i < snumber; i++) 
      {
	Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
	vector< ReconstructedParticle * > additional = GetAdditionalParticles(primaries, secondary, particles);
	for (unsigned int j = 0; j < additional.size(); j++) 
	  {
	    if (!IsDublicate(additional[j], *result)) 
	      {
		result->push_back(additional[j]);
	      }
	  }
      }
    return result;

  }

  void VertexRestorer::AnalyseSecondaries(EVENT::LCCollection * main, EVENT::LCCollection * missed)
  {
    int mnumber = main->getNumberOfElements();
    _missedTotal = missed->getNumberOfElements();
    vector< ReconstructedParticle * > lost;
    for ( int i = 0; i < _missedTotal; i++) 
      {
	ReconstructedParticle * secondary = dynamic_cast< ReconstructedParticle * > (missed->getElementAt(i));
	bool found = false;
	for (int j = 0; j < mnumber; j++) 
	  {
	    ReconstructedParticle * primary = dynamic_cast< ReconstructedParticle * > (main->getElementAt(j));
	    if (CompareParticles(secondary, primary)) 
	      {
		found = true;
		break;
	      }
	  }
	if (!found) 
	  {
	    lost.push_back(secondary);
	  }
      }
    _missedRecoTotal = lost.size();
    for (int i = 0; i < _missedRecoTotal; i++) 
      {
	_missedRecoMomentum[i] = MathOperator::getModule(lost[i]->getMomentum());
	vector<float> direction = MathOperator::getDirection(lost[i]->getMomentum());
	_missedRecoTheta[i] = MathOperator::getAngles(direction)[1];
      }
  }
  bool VertexRestorer::CompareParticles(const ReconstructedParticle * particle1, const ReconstructedParticle * particle2, bool module)
  {
    if (particle1 == particle2) 
      {
	//streamlog_out(DEBUG) << "Equal pointers!\n";
	return true;
      }
    //return false;
    if (particle1->getCharge() * particle2->getCharge() < 0.0) 
      {
	return false;
      }
    float angle = MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum());
    if (angle > 0.02) 
      {
	return false;
      }
    float recomodule = MathOperator::getModule(particle2->getMomentum());
    float mcmodule = MathOperator::getModule(particle1->getMomentum());
    float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
    if (abs(particle2->getCharge()) > 0.9 && abs(particle1->getCharge())  > 0.9) 
      {
	int vtxhits2 =  particle2->getTracks()[0]->getSubdetectorHitNumbers()[0];
	int vtxhits1 =  particle1->getTracks()[0]->getSubdetectorHitNumbers()[0];
	if ((vtxhits1 && !vtxhits2) || (vtxhits2 && !vtxhits1)) 
	  {
	    return angle < 0.01;
	  }
      }
    if (ratio > 0.1) 
      {
	return false;
      }
    return true;
  }
  bool VertexRestorer::CompareParticles(const MCParticle * particle1, const ReconstructedParticle * particle2)
  {
    if (particle1->getCharge() * particle2->getCharge() < 0.0) 
      {
	return false;
      }
    float angle = MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum());
    if (angle > 0.02) 
      {
	return false;
      }
    float recomodule = MathOperator::getModule(particle2->getMomentum());
    float mcmodule = MathOperator::getModule(particle1->getMomentum());
    float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
    if (ratio > 0.02) 
      {
	return false;
      }
    return true;
  }
  /*vector < EVENT::ReconstructedParticle * > * VertexRestorer::RestoreVertices(Vertex * primary, LCCollection * sec)
    {
    //vector < EVENT::Vertex * > * result = new vector< EVENT::Vertex* >();
    vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
    vector< ReconstructedParticle * > primaries = primary->getAssociatedParticle()->getParticles();
    int number = sec->getNumberOfElements();
    vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
    for (int i = 0; i < number; i++) 
    {
    Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
    vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
    particles->reserve(particles->size() + secondaries.size());
    particles->insert(particles->end(),secondaries.begin(),secondaries.end());
    }

    for (int i = 0; i < number; i++) 
    {
    Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
    vector< ReconstructedParticle * > additional = GetAdditionalParticles(primaries, secondary, particles);
    for (unsigned int j = 0; j < additional.size(); j++) 
    {
    if (!IsDublicate(additional[j], *result)) 
    {
    result->push_back(additional[j]);
    }
    }
    }
    return result;
    }*/
  vector< ReconstructedParticle * > * VertexRestorer::RestoreOneTrackVertices(LCCollection * jetcol, LCCollection * rel, LCCollection * sec)
  {
    vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
    int number = jetcol->getNumberOfElements();
    LCRelationNavigator navigator(rel);
    PIDHandler pidh(jetcol);
			
    /*			int alid = -1;
			try
			{
			alid = pidh.getAlgorithmID("lcfiplus");
			}
			catch(UTIL::UnknownAlgorithm &e)
			{
			streamlog_out(DEBUG) << "No algorithm lcfiplus!\n";
			alid = 0;
			}
			vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
			for (int i = 0; i < sec->getNumberOfElements(); i++) 
			{
			Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
			vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
			particles->reserve(particles->size() + secondaries.size());
			particles->insert(particles->end(),secondaries.begin(),secondaries.end());
			}

			for (int i = 0; i < number; i++)
			{
			ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
			int nvtx = navigator.getRelatedToObjects(jet).size();
			const ParticleID& pid = pidh.getParticleID(jet,alid);
			vector<float> params = pid.getParameters();
			float btag = params[pidh.getParameterIndex(alid,"BTag")];
			float ctag = params[pidh.getParameterIndex(alid,"CTag")];
			if (btag < 0.3 || nvtx == 2 || ctag > 0.3) 
			{
			continue;
			}
			streamlog_out(DEBUG) << "Jet energy: " << jet->getEnergy()
			<< " b-tag: " << btag 
			<< " c-tag: " << ctag 
			<< " # of vtx: " << nvtx
			<<  '\n';
			vector< ReconstructedParticle * > additional = GetAdditionalParticles(jet->getParticles(), particles);
			result->reserve(result->size() + additional.size());
			result->insert(result->end(), additional.begin(),additional.end());

			}*/
    return result;
  }
  vector < MyReconstructedParticle * > * VertexRestorer::RestoreVertices(LCCollection * jetcol, LCCollection * rel, LCCollection * secvtx,  LCCollection * trash)
  {
    vector < MyReconstructedParticle * > * result = new vector< MyReconstructedParticle * >();
    int number = jetcol->getNumberOfElements();
    LCRelationNavigator navigator(rel);
    int snumber = secvtx->getNumberOfElements();
    vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
    vector< Vertex * > tagged;

    for (int i = 0; i < snumber; i++) 
      {
	Vertex * secondary = dynamic_cast< Vertex * >(secvtx->getElementAt(i));
	tagged.push_back(secondary);
	vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
	particles->reserve(particles->size() + secondaries.size());
	particles->insert(particles->end(),secondaries.begin(),secondaries.end());
      }
    vector< ReconstructedParticle * >  resurrected;
    if (trash) 
      {
	//streamlog_out(DEBUG) << "Resurrected particles: \n";
	int tracknumber = trash->getNumberOfElements();
	vector< ReconstructedParticle * > jetparticles;
	for (int i = 0; i < number; i++) 
	  {
	    ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
	    jetparticles.reserve(jet->getParticles().size() + jetparticles.size());
	    jetparticles.insert(jetparticles.end(), jet->getParticles().begin(), jet->getParticles().end());

	  }
	for (int i = 0; i < tracknumber; i++) 
	  {
	    Track * track = dynamic_cast<Track *>(trash->getElementAt(i));
	    ReconstructedParticle * particle = myTrackOperator.ReconstructParticle(track,_Bfield);
	    if (!IsDublicate(particle, *particles) && !IsDublicate(particle, jetparticles)) 
	      {
		//PrintParticle(particle);
		resurrected.push_back(particle);
	      }
	  }
	for (int i = 0; i < resurrected.size(); i++) 
	  {
	    for (int j = 0; j < i; j++) 
	      {
		if (j != i && CompareParticles(resurrected[i],resurrected[j])) 
		  {
		    streamlog_out(DEBUG) <<"Found same zombies i: "<< i << " j: " << j <<"!\n";
		    //PrintParticle(resurrected[i]);
		    //PrintParticle(resurrected[j]);
		  }
	      }
	  }
	streamlog_out(DEBUG) << "\n";
      }
    for (int i = 0; i < number; i++) 
      {
	ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
	int nvtx = navigator.getRelatedToObjects(jet).size();
	streamlog_out(DEBUG) << "Jet energy: " << jet->getEnergy() 
	  // << " b-tag: " << btag 
	  // << " c-tag: " << ctag 
		  << " # of vtx: " << nvtx 
		  <<  '\n';
	vector< LCObject * > objs = navigator.getRelatedToObjects(jet);
	for (int j = 0; j < nvtx; j++) 
	  {
	    Vertex * vertex = dynamic_cast< Vertex * >(objs[j]);
	    bool tag = false;
	    for (int i = 0; i < tagged.size(); i++) 
	      {
		if (tagged[i] == vertex) 
		  {
		    tag = true;
		  }
	      }
	    if (!tag) 
	      {
		continue;
	      }
	    vector< ReconstructedParticle * > toInject;
	    toInject.reserve(jet->getParticles().size() + resurrected.size());
	    toInject.insert(toInject.end(), jet->getParticles().begin(), jet->getParticles().end());
	    toInject.insert(toInject.end(), resurrected.begin(), resurrected.end());

	    vector< MyReconstructedParticle * > additional = AddParticles(toInject, vertex, particles, jet);
	    for (unsigned int j = 0; j < additional.size(); j++) 
	      {
		if (!IsDublicate(additional[j], *result)) 
		  {
		    result->push_back(additional[j]);
		  }
	      }
	  }
      }
    return result;
  }
  vector< MyReconstructedParticle * > * VertexRestorer::__CheatRestoreVertices(LCCollection * jetcol, LCCollection * genvertexcol, LCCollection * trashcol)
  {
    vector < MyReconstructedParticle * > * result = new vector< MyReconstructedParticle * >();
    int jetnumber = jetcol->getNumberOfElements();
    int vtxnumber = genvertexcol->getNumberOfElements();
    vector< ReconstructedParticle * > toInject;
    for (int i = 0; i < jetnumber; i++) 
      {
	ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
	toInject.reserve(jet->getParticles().size()+toInject.size());
	toInject.insert(toInject.end(), jet->getParticles().begin(), jet->getParticles().end());
      }
    for (int i = 0; i < vtxnumber; i++) 
      {
	Vertex * vertex = dynamic_cast<Vertex * > (genvertexcol->getElementAt(i));
	vector< MyReconstructedParticle * > additional = AddParticles(toInject, vertex);
	for (unsigned int j = 0; j < additional.size(); j++) 
	  {
	    if (!IsDublicate(additional[j], *result)) 
	      {
		result->push_back(additional[j]);
	      }
	  }
      }
    return result;
  }

  vector< ReconstructedParticle * > VertexRestorer::GetAdditionalParticles(const vector< ReconstructedParticle * > & pri, Vertex * sec, const vector< ReconstructedParticle * > * toCompare)
  {
    vector< ReconstructedParticle * > result;
    if (sec->getAssociatedParticle()->getParticles().size() < 2) 
      {
	return result;
      }
    const vector< ReconstructedParticle * > * particles = (toCompare)? toCompare : &(sec->getAssociatedParticle()->getParticles());
    for (unsigned int i = 0; i < pri.size(); i++) 
      {
	ReconstructedParticle * primaryTrack = pri[i];
	if (abs(primaryTrack->getCharge())  < 0.9) 
	  {
	    continue;
	  }
	vector< float > direction = MathOperator::getDirection(primaryTrack->getMomentum());
	double * pos = MathOperator::toDoubleArray(sec->getPosition(), 3);
	double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
	double * trackPosition = myTrackOperator.GetStartPoint(primaryTrack);
	float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
	float secondaryOffset = MathOperator::getDistanceTo(pos, direction, trackPosition);
	float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primaryTrack, trackPosition, myPrimary, primaryOffset));
	vector<float> diff = MathOperator::getDirection(pos, trackPosition);
	double diif[3];
	for (int m = 0; m < 3; m++) 
	  {
	    diif[m] = diff[m];
	  }
	float tan  = secondaryOffset / MathOperator::getModule(diif);
	float angle = MathOperator::getAngle(primaryTrack->getMomentum(), diif);
	if (primaryOffset / accuracy < 3.0
	    //&& primaryOffset / accuracy < 3.0
				    
	    //&& secondaryOffset < 5.0
	    && angle < 0.03)
	  {
	    streamlog_out(DEBUG) << "Found a track with offset " << primaryOffset
		      << ", deviation " << primaryOffset / accuracy 
		      << ", angle " << angle
		      << ", distance " << MathOperator::getModule(trackPosition)
		      << ", error*e-4 " << accuracy
		      <<" :\n";
	    PrintParticle(primaryTrack);
	    bool found = false;
	    for (unsigned int j = 0; j < particles->size(); j++) 
	      {
		ReconstructedParticle * secondary = particles->at(j);
		if (CompareParticles(secondary, primaryTrack)) 
		  {
		    found = true;
		    break;
		  }
	      }
	    if (!found) 
	      {
		streamlog_out(DEBUG) << "The particle is unique!\n";
		//	if (!IsDublicate(primaryTrack,result)) 
		{
		  result.push_back(CopyParticle(primaryTrack,sec, tan));
		}
	      }
	  }
      }
    return result;
  }
		
  ReconstructedParticle * VertexRestorer::CopyParticle(const ReconstructedParticle * particle, const Vertex * vertex, float theta)
  {
    ReconstructedParticleImpl * nouveau = new ReconstructedParticleImpl();
    VertexImpl * vimp = new VertexImpl();
    nouveau->setMomentum(particle->getMomentum());
    nouveau->setMass(particle->getMass());
    nouveau->setCharge(particle->getCharge());
    nouveau->setEnergy(particle->getEnergy());
    vimp->setPosition(vertex->getPosition());
    vimp->setPrimary(vertex->isPrimary());
    nouveau->setStartVertex(vimp);
    nouveau->addTrack(particle->getTracks()[0]); //CRUNCH
    AddParticleID(nouveau, theta);
    return nouveau;
  }
  bool VertexRestorer::IsDublicate(const ReconstructedParticle * particle, const vector< ReconstructedParticle * > & data)
  {
    bool dublicate = false;
    for (unsigned int j = 0; j < data.size(); j++) 
      {
	if (CompareParticles(particle, data[j])) 
	  {
	    //streamlog_out(DEBUG) << "Dublicate found!!!!\n";
	    dublicate = true;
	    break;
	  }
      }
    return dublicate;
  }
  bool VertexRestorer::IsDublicate( MyReconstructedParticle * particle, vector< MyReconstructedParticle * > & data)
  {
    bool dublicate = false;
    for (unsigned int j = 0; j < data.size(); j++) 
      {
	if (CompareParticles(particle->Get(), data[j]->Get())) 
	  {
	    //streamlog_out(DEBUG) << "Dublicate found!!!!\n";
	    dublicate = true;
	    break;
	  }
      }
    return dublicate;
  }
  float VertexRestorer::GetError(const ReconstructedParticle * particle)
  {
    if (!particle || particle->getTracks().size() < 1) 
      {
	streamlog_out(DEBUG) << "The particle is null or 0 tracks!\n";
	return 0.0;
      }
    float p = MathOperator::getModule(particle->getMomentum());
    vector<float> direction = MathOperator::getDirection(particle->getMomentum());
    vector<float> angles = MathOperator::getAngles(direction);
    float accuracy = sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p * pow(sin(angles[1]), 4.0/3.0)) );
    if (abs(cos(angles[1])) > 0.95) 
      {
	//accuracy = sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p * pow(cos(angles[1]), 4.0/3.0)) );
      }
    return accuracy;
  }
  void VertexRestorer::TestMethod(EVENT::LCCollection * prongs, EVENT::LCCollection * pfo, EVENT::LCCollection * rel)
  {
    _nlosttest = 0;
    _ngoodtest = 0;
    _nmultitest = 0;
    _ntotaltest = 0;
    vector< MCParticle * > verifiedProngs;
    for (int i = 0; i < pfo->getNumberOfElements(); i++) 
      {
	ReconstructedParticle * recoparticle = dynamic_cast<ReconstructedParticle *>(pfo->getElementAt(i));
	if (std::abs(recoparticle->getCharge()) < 0.01) 
	  {
	    continue;
	  }
	int times = 0;
	MCParticle * mcprong = CompareToCollection(recoparticle, prongs, rel);
	if (mcprong) 
	  {
	    verifiedProngs.push_back(mcprong);
	    for (unsigned int j = 0; j < verifiedProngs.size(); j++) 
	      {
		if (verifiedProngs[j] == mcprong) 
		  {
		    times++;
		  }
	      }
	  }
	vector<float> direction = MathOperator::getDirection(recoparticle->getMomentum());
	float cos = std::cos(MathOperator::getAngles(direction)[1]);
	if (times == 1) 
	  {
	    streamlog_out(DEBUG) << "A good one found\n";
	    _costhetaTestedGood[_ngoodtest] = cos;
	    _ngoodtest++;
	  }
	if (times == 0) 
	  {
	    streamlog_out(DEBUG) << "A lost one found\n";

	    _costhetaTestedLost[_nlosttest] = cos;
	    _nlosttest++;
	  }
	if (times>1) 
	  {
	    streamlog_out(DEBUG) << "Several tracks found\n";
	    _costhetaTestedMultiple[_nmultitest] = cos;
	    _nmultitest++;
	  }

      }
    for (int i = 0; i < prongs->getNumberOfElements(); i++) 
      {
	MCParticle * mcparticle = dynamic_cast<MCParticle *>(prongs->getElementAt(i));
	vector<float> direction = MathOperator::getDirection(mcparticle->getMomentum());
	float cos = std::cos(MathOperator::getAngles(direction)[1]);
	_costhetaTestedTotal[_ntotaltest] = cos;
	_ntotaltest++;
      }
    return;
  }
  void VertexRestorer::TestMethod2(EVENT::LCCollection * prongs, EVENT::LCCollection * pfo, EVENT::LCCollection * rel)
  {
    double momentum[3];
    momentum[0] = 30.0;
    momentum[1] = 30.0;
    momentum[2] = 0.0;
    double zero[3];
    zero[0] = 0.0;
    zero[1] = 0.0;
    zero[2] = 0.0;
    double point[3];
    point[0] = 6.0;
    point[1] = 1.0;
    point[2] = 0.0;
    vector< float > direction = MathOperator::getDirection(momentum);
    float offset = MathOperator::getDistanceTo(point, direction, zero);
    streamlog_out(DEBUG) << "Our offset is " << offset << "mm\n";
    int vtxnumber = pfo->getNumberOfElements();

    _ntest = 0;
    for (int i = 0; i < vtxnumber; i++) 
      {
	Vertex * recovtx = dynamic_cast< Vertex * > (pfo->getElementAt(i));
	const vector< ReconstructedParticle * > particles = recovtx->getAssociatedParticle()->getParticles();
	float minDistance = 10.0;
	for (int j = 0; j < particles.size(); j++) 
	  {
	    double * secPos = myTrackOperator.GetStartPoint(particles[j]);
	    vector< float > secDir = MathOperator::getDirection(particles[j]->getMomentum());
	    double * vtxPos = MathOperator::toDoubleArray(recovtx->getPosition(),3);
	    for (int k = 0; k < j; k++) 
	      {
		if (k == j)
		  {
		    continue;
		  }
		double * secPos1 = myTrackOperator.GetStartPoint(particles[k]);
		vector< float > secDir1 = MathOperator::getDirection(particles[k]->getMomentum());
		//float distance = MathOperator::getAngle(particles[j]->getMomentum(), particles[i]->getMomentum());
		float distance = MathOperator::getDistanceBtw(secPos, secDir, secPos1, secDir1);
		streamlog_out(DEBUG) << "i: " << i << " j: " << j << " distance: " << distance << "\n";
		if (distance < minDistance) 
		  {
		    minDistance = distance;
		  }
	      }
	  }
	_testparameter[_ntest] = recovtx->getChi2();
	_testDistance[_ntest++] = minDistance;
      }//*/
    /*int prongnumber = pfo->getNumberOfElements();
      vector< ReconstructedParticle * > particles;
      for (int i = 0; i < prongnumber; i++) 
      {
      MCParticle * mcparticle = dynamic_cast<MCParticle *>(prongs->getElementAt(i));
      LCRelationNavigator navigator(rel);
      vector< LCObject * > obj = navigator.getRelatedFromObjects(mcparticle);
      vector< float > weights = navigator.getRelatedFromWeights(mcparticle);
      if (obj.size() < 1)
      {
      streamlog_out(DEBUG) << "Particle was not generated!\n";
      return false;
      }
      MCParticle * winner = NULL;
      float maxweight = 0.0;
      for (unsigned int j = 0; j < obj.size(); j++) 
      {
      if (std::abs(dynamic_cast< ReconstructedParticle * >(obj[j])->getCharge()) < 0.9) 
      {
      continue;
      }
      if (weights[j] > maxweight) 
      {
      winner = dynamic_cast< ReconstructedParticle * >(obj[j]);
      maxweight = weights[j];
      }
      }
      if (winner) 
      {
					
      }
    *///}
  }//*/
  void VertexRestorer::TestCollections(EVENT::LCCollection * tagged, EVENT::LCCollection * pfo)
  {
    int pfonumber = pfo->getNumberOfElements();
    const vector< ReconstructedParticle * > primaries = myPrimary->getAssociatedParticle()->getParticles();
    for (int i = 0; i < pfonumber; i++) 
      {
				
      }
  }

  void VertexRestorer::check( LCEvent * evt ) {}


  float VertexRestorer::getChi2(EVENT::ReconstructedParticle * primary, const EVENT::Vertex * sec)
  {
    double * trackPosition = myTrackOperator.GetStartPoint(primary);
    vector< float > direction = MathOperator::getDirection(primary->getMomentum());
    double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
    double * secondaryPosition = MathOperator::toDoubleArray(sec->getPosition(),3);
    float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
    float accuracy = GetError(primary);
    //float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primary, trackPosition, myPrimary, primaryOffset));
    //float accuracy = std::sqrt(myTrackOperator.GetOffsetErrorSimple(primary));
    //float secondaryOffset = MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
    vector<float> diff = MathOperator::getDirection(secondaryPosition, trackPosition);
    double diif[3];
    for (int m = 0; m < 3; m++) 
      {
	diif[m] = diff[m];
      }
    float momentum = MathOperator::getModule(primary->getMomentum());
    float angle = MathOperator::getAngle(diif, primary->getMomentum());
    //float angleError = 1.0;
    //if (angle > 0.0) 
    //{
    //	angleError = std::sqrt(myTrackOperator.GetAngleError(angle, sec, primary));
    //}
    int vtxhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[0];
    int sithits = primary->getTracks()[0]->getSubdetectorHitNumbers()[2];
    int ftdhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[5];
    int ntracks = sec->getAssociatedParticle()->getParticles().size();
    vector< float > parametersVertex = ParametrizeVertex(sec);
    float angleSigma = 0.02;
    float meanOffset = parametersVertex[0];
    float sigmaOffset = parametersVertex[1] / (ntracks);

    float meanAngle = parametersVertex[2];
    float sigmaAngle = (parametersVertex[3] < 0.005)? 0.005 :parametersVertex[3];

    float meanVtx = parametersVertex[4];
    float sigmaVtx = (parametersVertex[5] < 0.00001)? 1. / ntracks:parametersVertex[5];
			
    float meanObs = parametersVertex[6];
    float sigmaObs = (parametersVertex[7] < 0.00001)? 1. / ntracks : parametersVertex[7];
    streamlog_out(DEBUG) << "\tmean offset: " << meanOffset << " sigma: "  << sigmaOffset << "\n";
    streamlog_out(DEBUG) << "\tmean angle: " << meanAngle << " sigma: "  << sigmaAngle << "\n";
    streamlog_out(DEBUG) << "\tmean nvtxhits: " << meanVtx << " sigma: "  << sigmaVtx << "\n";
    streamlog_out(DEBUG) << "\tmean obs: " << meanObs << " sigma: "  << sigmaObs << "\n";
			
    //return 0.0;
    float chi2 = std::pow(primaryOffset / accuracy - meanOffset, 2) /std::pow( sigmaOffset, 2) +
      std::pow(angle - meanAngle, 2) /std::pow( sigmaAngle, 2) +
      std::pow(vtxhits - meanVtx, 2) /std::pow( sigmaVtx, 2 ); 
    std::pow(momentum - meanObs, 2) /std::pow( sigmaObs, 2);

    return chi2;
  }

  void VertexRestorer::end()
  {   
    hfile->Write();
    hfile->Close();
  }
}
