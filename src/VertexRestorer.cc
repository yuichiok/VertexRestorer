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

	  registerInputCollection( LCIO::VERTEX,
	      "PrimaryCollectionName" , 
	      "Name of the PrimaryVertex collection"  ,
	      _colPriName ,
	      std::string("PrimaryVertex")
	   );
	    registerInputCollection( LCIO::VERTEX,
		    "SecondaryCollectionName" , 
		    "Name of the BuildUpVertex collection"  ,
		    _colSecName ,
		    std::string("FinalJets_vtx")
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
	    registerInputCollection( LCIO::MCPARTICLE,
		    "EGProngsCollectionName" , 
		    "Name of the EGProngs collection"  ,
		    _colEGPName ,
		    std::string("EGProngs")
	    );
	    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
	    	"JetRelCollectionName",
		"Name of the Jet relation collection",
		_colJetRelName,
	    	std::string("FinalJets_rel")
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
	   _Tree->Branch("missedAngle", _missedAngle, "missedAngle[missedDetected]/F");
	   _Tree->Branch("missedObs", _missedObs, "missedObs[missedDetected]/F");
	   _Tree->Branch("missedDeviation", _missedDeviation, "missedDeviation[missedDetected]/F");
	   _Tree->Branch("missedSecDeviation", _missedSecDeviation, "missedSecDeviation[missedDetected]/F");
	   _Tree->Branch("missedCostheta", _missedCostheta, "missedCostheta[missedDetected]/F");
	   
	   _Tree->Branch("bstarDetected", &_bstarDetected, "bstarDetected/I");
	   _Tree->Branch("bstarMomentum", _bstarMomentum, "bstarMomentum[bstarDetected]/F");
	   _Tree->Branch("bstarError", _bstarError, "bstarError[bstarDetected]/F");
	   _Tree->Branch("bstarAngle", _bstarAngle, "bstarAngle[bstarDetected]/F");
	   _Tree->Branch("bstarDeviation", _bstarDeviation, "bstarDeviation[bstarDetected]/F");
	   
	   /*_Tree->Branch("v0Detected", &_v0Detected, "v0Detected/I");
	   _Tree->Branch("v0Momentum", _v0Momentum, "v0Momentum[v0Detected]/F");
	   _Tree->Branch("v0Error", _v0Error, "v0Error[v0Detected]/F");
	   _Tree->Branch("v0Angle", _v0Angle, "v0Angle[v0Detected]/F");
	   _Tree->Branch("v0Deviation", _v0Deviation, "v0Deviation[v0Detected]/F");*/
	   
	   _Tree->Branch("fakeDetected", &_fakeDetected, "fakeDetected/I");
	   _Tree->Branch("fakeMomentum", _fakeMomentum, "fakeMomentum[fakeDetected]/F");
	   _Tree->Branch("fakeError", _fakeError, "fakeError[fakeDetected]/F");
	   _Tree->Branch("fakeAngle", _fakeAngle, "fakeAngle[fakeDetected]/F");
	   _Tree->Branch("fakeObs", _fakeObs, "fakeObs[fakeDetected]/F");
	   _Tree->Branch("fakeDeviation", _fakeDeviation, "fakeDeviation[fakeDetected]/F");
	   _Tree->Branch("fakeSecDeviation", _fakeSecDeviation, "fakeSecDeviation[fakeDetected]/F");
	   _Tree->Branch("fakeCostheta", _fakeCostheta, "fakeCostheta[fakeDetected]/F");
	   
	   _Tree->Branch("missedRecoTotal", &_missedRecoTotal, "missedRecoTotal/I");
	   _Tree->Branch("missedRecoMomentum", _missedRecoMomentum, "missedRecoMomentum[missedRecoTotal]/F");
	   _Tree->Branch("missedRecoTheta", _missedRecoTheta, "missedRecoTheta[missedRecoTotal]/F");


	   _primaryTree =  new TTree( "Primaries", "tree" );
	   _primaryTree->Branch("primariesTotal", &_primariesTotal, "primariesTotal/I");
	   _primaryTree->Branch("primeOffset", &_primeOffset, "primeOffset[primariesTotal]/F");
	   _primaryTree->Branch("primeMomentum", &_primeMomentum, "primeMomentum[primariesTotal]/F");
	   _primaryTree->Branch("primeDeviation", &_primeDeviation, "primeDeviation[primariesTotal]/F");
	   _primaryTree->Branch("primeError", &_primeError, "primeError[primariesTotal]/F");
	   
	   _TestTree =  new TTree( "Test", "tree" );
	   _TestTree->Branch("ntotaltest", &_ntotaltest, "ntotaltest/I");
	   _TestTree->Branch("ngoodtest", &_ngoodtest, "ngoodtest/I");
	   _TestTree->Branch("nlosttest", &_nlosttest, "nlosttest/I");
	   _TestTree->Branch("nmultitest", &_nmultitest, "nmultitest/I");
	   _TestTree->Branch("costhetaTestedTotal", _costhetaTestedTotal, "costhetaTestedTotal[ntotaltest]/F");
	   _TestTree->Branch("costhetaTestedGood", _costhetaTestedGood, "costhetaTestedGood[ngoodtest]/F");
	   _TestTree->Branch("costhetaTestedLost", _costhetaTestedLost, "costhetaTestedLost[nlosttest]/F");
	   _TestTree->Branch("costhetaTestedMultiple", _costhetaTestedMultiple, "costhetaTestedMultiple[nmultitest]/F");

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
		std::cout << "|\t" << track->getD0() 
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
			std::cout << "Particle has no vertex!\n";
			return false;
		}
		void VertexRestorer::processEvent( LCEvent * evt ) 
		{ 
			_nEvt ++ ;
		  
			try 
			{
				std::cout << "********VertexRecovery*"<<_nEvt << "**********\n";
				TrackOperator opera;
				opera.test();
				LCCollection* maincol = evt->getCollection( _colName );
				LCCollection* seccol = evt->getCollection( _colSecName );
				LCCollection* v0col = evt->getCollection( _colV0Name );
				LCCollection* relation = evt->getCollection(_colRelName);
				//vector < ReconstructedParticle * > * result = RestoreVerticesPFO(maincol, seccol);
				LCCollection* pricol = evt->getCollection( _colPriName );
				Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
				myPrimary = primary;
				vector < MyReconstructedParticle * > * result = RefineVertices(maincol, seccol, v0col);
				//LCCollection* jetcol = evt->getCollection( _colJetName );
				//LCCollection* relcol = evt->getCollection( _colJetRelName );
				//vector < ReconstructedParticle * > * result = RestoreOneTrackVertices(jetcol, relcol, seccol);
				
				//vector < ReconstructedParticle * > * result = RestoreVertices(primary, seccol);
				//vector < EVENT::Vertex * > * result = RestoreVertices(jetcol, relcol);
				LCCollection* misscol = evt->getCollection( _colMissedName );
				LCCollection* missvtxcol =evt->getCollection( _colMissedVtxName );
				LCCollection* bstarcol = evt->getCollection( _colBStarName );
				//AnalyseSecondaries(maincol, misscol);
				CompareCollectionsRel(result, evt->getCollection(_colEGPName), bstarcol , relation);
				//CompareCollections(result, misscol, bstarcol, missvtxcol);//*/
				//TestMethod(evt->getCollection(_colEGPName), maincol, evt->getCollection(_colRelName));
				_Tree->Fill();
				//_TestTree->Fill();
				//WritePrimaries(opera, pricol, evt->getCollection(_colEGPName), relation);
				//_primaryTree->Fill();
				ClearVariables();
			}
			catch( DataNotAvailableException &e)
			{
				std::cout << "Whoops!....\n";
			}
		    
		}
		void VertexRestorer::ClearVariables()
		{
			_primariesTotal = 0;
			for (int i = 0; i < MAXP; i++) 
			{
				_primeOffset[i] = -1.0;
				_primeDeviation[i] = -1.0;
				_primeMomentum[i] = -1.0;
				_primeError[i] = -1.0;
			}
			/*for (int i = 0; i < MAXN; i++) 
			{
				_missedAngle[i] = -1.0;
				_fakeAngle[i] = -1.0;
			}*/
		}

		void VertexRestorer::WritePrimaries(TrackOperator & opera, EVENT::LCCollection * pricol, EVENT::LCCollection * sec, EVENT::LCCollection * rel)
		{
			Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
			double * pos = MathOperator::toDoubleArray(primary->getPosition(), 3);
			vector < ReconstructedParticle * > primaries = primary->getAssociatedParticle()->getParticles();
			_primariesTotal = primaries.size();
			for (unsigned int i = 0; i < primaries.size(); i++)
			{
				//opera.test(primaries[i]);
				if (sec && rel && CompareToCollection(primaries[i], sec,rel)) 
				{
					std::cout << "EXCLUDED!\n";
					continue;
				}
				vector<float> direction = MathOperator::getDirection(primaries[i]->getMomentum());
				double * start = opera.GetStartPoint(primaries[i]);
			        //PrintParticle(primaries[i]);
				_primeOffset[i] = MathOperator::getDistanceTo(pos, direction, start);//opera.GetOffset(primaries[i]);
				_primeMomentum[i] = MathOperator::getModule(primaries[i]->getMomentum());

				_primeError[i] = std::sqrt(opera.GetOffsetError(primaries[i], opera.GetStartPoint(primaries[i]), primary, _primeOffset[i]));
				//std::cout << "Offset: " << _primeOffset[i] << "; error^2: " << _primeError[i] << '\n'; 
				_primeDeviation[i] = _primeOffset[i] / _primeError[i]; 
			        //std::cout << "Offset: " << _primeOffset[i] << " deviation: " << _primeDeviation[i] << '\n';
				
			}

		}
		MCParticle * VertexRestorer::CompareToCollection(EVENT::ReconstructedParticle * particle, EVENT::LCCollection * prongs, EVENT::LCCollection * rel)
		{
			LCRelationNavigator navigator(rel);
			vector< LCObject * > obj = navigator.getRelatedToObjects(particle);
			vector< float > weights = navigator.getRelatedToWeights(particle);
			if (obj.size() < 1)
			{
				std::cout << "Particle was not generated!\n";
				return false;
			}
			MCParticle * winner = NULL;
			float maxweight = 0.0;
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
			/*vector< MCParticle * > generated;
			for (unsigned int i = 0; i < obj.size(); i++) 
			{
				MCParticle * mcpart = dynamic_cast< MCParticle * >(obj[i]);
				generated.push_back(mcpart);
			}
			//MCParticle * mcpart = dynamic_cast< MCParticle * >(obj[0]); // CRUNCH!

			int numberOfProngs = prongs->getNumberOfElements();
			for (int i = 0; i < numberOfProngs; i++) 
			{
				MCParticle * mcprong = dynamic_cast< MCParticle * >(prongs->getElementAt(i));
				for (unsigned int j = 0; j < obj.size(); j++) 
				{
					if (mcprong == generated[j]) 
					{
						return true;
					}
				}
			}

			return false;// */
		}
		void VertexRestorer::CompareCollectionsRel(std::vector< MyReconstructedParticle * > * detected, EVENT::LCCollection * prongs, EVENT::LCCollection * bstar, EVENT::LCCollection *  rel)
		{
			_detectedTotal = detected->size();
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

			for (int i = 0; i < _detectedTotal; i++) 
			{
				float taken = false;
				ReconstructedParticle * particle = detected->at(i)->Get();
				if (CompareToCollection(particle, prongs, rel)) 
				{
					taken = true;
					std::cout << "Found a prong!\n";
					_missedError[_missedDetected] = GetError(particle);
					_missedMomentum[_missedDetected] = MathOperator::getModule(particle->getMomentum());
					_missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() / _missedError[_missedDetected];
					_missedAngle[_missedDetected] = detected->at(i)->GetAngle();
					_missedObs[_missedDetected] = detected->at(i)->GetObservable();
					_missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset() / _missedError[_missedDetected];
					_missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
					_missedDetected++;
				}
				if (CompareToCollection(particle, bstar, rel)) 
				{
					taken = true;
					_bstarError[_bstarDetected] = GetError(particle);
					_bstarMomentum[_bstarDetected] = MathOperator::getModule(particle->getMomentum());
					_bstarDetected++;
				}
				if (!taken) 
				{
					_fakeError[_fakeDetected] = GetError(particle);
					_fakeMomentum[_fakeDetected] = MathOperator::getModule(particle->getMomentum());
					_fakeDeviation[_fakeDetected] =  detected->at(i)->GetOffset() / _fakeError[_fakeDetected];
					_fakeSecDeviation[_fakeDetected] =  detected->at(i)->GetSecOffset() / _fakeError[_fakeDetected];
					_fakeAngle[_fakeDetected] = detected->at(i)->GetAngle(); 
					_fakeObs[_fakeDetected] = detected->at(i)->GetObservable(); 
					_fakeCostheta[_fakeDetected] = detected->at(i)->GetCostheta(); 
					_fakeDetected++;
				}
			}
			std::cout << "Total detected: " << _detectedTotal
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
						_missedError[_missedDetected] = GetError(detectedParticle);
						_missedMomentum[_missedDetected] = MathOperator::getModule(particle->getMomentum());
						_missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() / _missedError[_missedDetected];
						_missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset() / _missedError[_missedDetected];
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
							_missedError[_missedDetected] = GetError(detectedParticle);
							_missedDeviation[_missedDetected] =  detected->at(i)->GetOffset() / _missedError[_missedDetected];
							_missedAngle[_missedDetected] = detected->at(i)->GetAngle();
							_missedObs[_missedDetected] = detected->at(i)->GetObservable();
							_missedCostheta[_missedDetected] = detected->at(i)->GetCostheta();
							_missedSecDeviation[_missedDetected] =  detected->at(i)->GetSecOffset() / _missedError[_missedDetected];
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
					_fakeError[_fakeDetected] = GetError(detectedParticle);
					_fakeDeviation[_fakeDetected] =  detected->at(i)->GetOffset() / _fakeError[_fakeDetected];
					_fakeSecDeviation[_fakeDetected] =  detected->at(i)->GetSecOffset() / _fakeError[_fakeDetected];
					_fakeAngle[_fakeDetected] = detected->at(i)->GetAngle(); 
					_fakeObs[_fakeDetected] = detected->at(i)->GetObservable(); 
					_fakeCostheta[_fakeDetected] = detected->at(i)->GetCostheta(); 
					_fakeMomentum[_fakeDetected++] = MathOperator::getModule(detectedParticle->getMomentum());
				}
			}
			//_fakeDetected = _detectedTotal - _missedDetected - _bstarDetected;
			std::cout << "Total detected: " << _detectedTotal
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
				particles->insert(particles->end(),secondaries.begin(),secondaries.end());
			}
			for (int i = 0; i < snumber; i++) 
			{
				Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
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
		bool VertexRestorer::TakeParticle(EVENT::ReconstructedParticle * primary, const EVENT::Vertex * sec, vector<double> * parameters)
		{
				if (std::abs(primary->getCharge() ) < 0.9) 
				{
					return false;
				}
				double * trackPosition = myTrackOperator.GetStartPoint(primary);
				float trackDistance = MathOperator::getModule(trackPosition);
				if (trackDistance > 60.0) 
				{
					return false;
				}
				//return true;
				vector< float > direction = MathOperator::getDirection(primary->getMomentum());
				double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
				double * secondaryPosition = MathOperator::toDoubleArray(sec->getPosition(),3);
				float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
				float accuracy = GetError(primary);
				//float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primary, trackPosition, myPrimary, primaryOffset));
				float secondaryOffset = MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
				vector<float> diff = MathOperator::getDirection(secondaryPosition, trackPosition);
				double diif[3];
				for (int m = 0; m < 3; m++) 
				{
					diif[m] = diff[m];
				}
				//float theta = secondaryOffset / MathOperator::getModule(secondaryPosition); // MathOperator::getAngle(diif, primary->getMomentum());
				float angle = MathOperator::getAngle(diif, primary->getMomentum());
				float secOffset =  MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
				vector< float > secDiraction = MathOperator::getDirection(secondaryPosition);
				float costheta = std::cos(MathOperator::getAngles(direction)[1]); //std::abs( std::cos(MathOperator::getAngles(secDiraction)[1] ) );
				float sine = (sec->getChi2() > 10.)?std::sin(MathOperator::getAngle(secondaryPosition, primary->getMomentum())) : 0;
				float cos = (sec->getChi2() > 10.)?std::cos(MathOperator::getAngle(secondaryPosition, trackPosition)) : 0;
				bool result = primaryOffset / accuracy > 3.0
					&& angle < 0.5;// + 0.03 * sine;

				if (result) 
				{
					float distance = myTrackOperator.GetDistanceBtw(myPrimary, sec, primary);
					float observable = std::abs( distance / secondaryOffset -0.5 ) * cos;
					std::cout << "Found a track with offset " << primaryOffset
						<< ", error " << accuracy//GetError(primary) 
						<< ", angle " << angle//GetError(primary) 
						<< ", DISTANCE " << distance//GetError(primary) 
						<< ", position " << trackDistance
						<<" :\n";
					if (parameters) 
					{
						parameters->push_back(angle);
						parameters->push_back(observable);
						parameters->push_back(primaryOffset);
						parameters->push_back(secOffset);
						parameters->push_back(costheta);
					}
				}
				else 
				{
					if (parameters) 
					{
						delete parameters;
					}
				}
				return result;
		}

		vector< MyReconstructedParticle * > VertexRestorer::AddParticles(const std::vector< EVENT::ReconstructedParticle * > & pri, EVENT::Vertex * sec, const std::vector< EVENT::ReconstructedParticle * > * toCompare)
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
				if (TakeParticle(primary, sec, parameters)) 
				{
					PrintParticle(primary);
					bool found = false;
					for (unsigned int j = 0; j < particles->size(); j++) 
					{
						ReconstructedParticle * secondary = particles->at(j);
						if (CompareParticles(secondary, primary)) 
						{
							found = true;
							break;
						}
					}
					
					if (!found)
					{
						//std::cout << "The particle is unique!\n";
						//float errorvtx = (sec->getChi2() < 0.00001)? 0.00001: sec->getChi2();
						MyReconstructedParticle * newparticle = new MyReconstructedParticle(primary);
						newparticle->SetAngle(parameters->at(0));
						newparticle->SetObservable(parameters->at(1));
						newparticle->SetOffset(parameters->at(2));
						newparticle->SetSecOffset(parameters->at(3));
						newparticle->SetCostheta(parameters->at(4));
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

		void VertexRestorer::PrintParticle(ReconstructedParticle * particle)
		{
			if (!particle) 
			{
				return;
			}
			std::cout << std::fixed << std::setw( 6 ) << std::setprecision( 3 ) << std::setfill( ' ' );
			int id = 0;
			if (particle->getParticleIDUsed()) 
			{
				std::cout << "Type " << particle->getParticleIDUsed()->getType() << '\n';
				id = particle->getParticleIDs()[0]->getPDG(); 
			}
			std::cout<<"|"<< id <<"\t\t|"<<particle->getMass()<<"\t\t|"<<particle->getCharge()  <<"\t\t|"<<particle->getEnergy() <<"\t\t|\n";
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
		bool VertexRestorer::CompareParticles(const ReconstructedParticle * particle1, const ReconstructedParticle * particle2, bool verbose)
		{
			if (particle1 == particle2) 
			{
				//std::cout << "Equal pointers!\n";
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
			if (ratio > 0.04) 
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
		vector < EVENT::ReconstructedParticle * > * VertexRestorer::RestoreVertices(Vertex * primary, LCCollection * sec)
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
		}
		vector< ReconstructedParticle * > * VertexRestorer::RestoreOneTrackVertices(LCCollection * jetcol, LCCollection * rel, LCCollection * sec)
		{
			vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
			int number = jetcol->getNumberOfElements();
			LCRelationNavigator navigator(rel);
			PIDHandler pidh(jetcol);
			int alid = pidh.getAlgorithmID("lcfiplus");
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
				std::cout << "Jet energy: " << jet->getEnergy()
					  << " b-tag: " << btag 
					  << " c-tag: " << ctag 
					  << " # of vtx: " << nvtx
					  <<  '\n';
				vector< ReconstructedParticle * > additional = GetAdditionalParticles(jet->getParticles(), particles);
				result->reserve(result->size() + additional.size());
				result->insert(result->end(), additional.begin(),additional.end());

			}
			return result;
		}
		vector< ReconstructedParticle * > VertexRestorer::GetAdditionalParticles(const vector< ReconstructedParticle * > & pri, const vector< ReconstructedParticle * > * particles )
		{
			vector< ReconstructedParticle * > result;
			for (unsigned int i = 0; i < pri.size(); i++) 
			{
				ReconstructedParticle * primaryTrack = pri[i];
				if (abs(primaryTrack->getCharge())  < 0.9) 
				{
					continue;
				}
				vector< float > direction = MathOperator::getDirection(primaryTrack->getMomentum());
				double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
				double * trackPosition = myTrackOperator.GetStartPoint(primaryTrack);
				float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
				float accuracy = std::sqrt(myTrackOperator.GetOffsetError(primaryTrack, trackPosition, myPrimary, primaryOffset));
				if (primaryOffset / accuracy > 5.0
				    && primaryOffset < 3.0)
				{
					std::cout << "Found a track with offset " << primaryOffset
						<< ", deviation " << primaryOffset / accuracy 
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
						std::cout << "The particle is unique!\n";
						result.push_back(primaryTrack);
					}
				}
			}
			return result;	    //&& primaryOffset / accuracy < 3.0
			
		}
		vector < EVENT::ReconstructedParticle * > * VertexRestorer::RestoreVertices(LCCollection * jetcol, LCCollection * rel)
		{
			vector < ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
			int number = jetcol->getNumberOfElements();
			LCRelationNavigator navigator(rel);
			
			for (int i = 0; i < number; i++) 
			{
				ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
				int nvtx = navigator.getRelatedToObjects(jet).size();
				std::cout << "Jet energy: " << jet->getEnergy() 
					 // << " b-tag: " << btag 
					 // << " c-tag: " << ctag 
					  << " # of vtx: " << nvtx 
					  <<  '\n';
				vector< LCObject * > objs = navigator.getRelatedToObjects(jet);
				for (int j = 0; j < nvtx; j++) 
				{
					Vertex * vertex = dynamic_cast< Vertex * >(objs[j]);
					vector< ReconstructedParticle * > additional = GetAdditionalParticles(jet->getParticles(), vertex);
					result->reserve(result->size() + additional.size());
					result->insert(result->end(), additional.begin(),additional.end());
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
				if (primaryOffset / accuracy < 5.0
				    //&& primaryOffset / accuracy < 3.0
				    
				    //&& secondaryOffset < 5.0
				    && angle < 0.03)
				{
					std::cout << "Found a track with offset " << primaryOffset
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
						std::cout << "The particle is unique!\n";
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
		bool VertexRestorer::IsDublicate(const ReconstructedParticle * particle, vector< ReconstructedParticle * > & data)
		{
			bool dublicate = false;
			for (unsigned int j = 0; j < data.size(); j++) 
			{
				if (CompareParticles(particle, data[j])) 
				{
					//std::cout << "Dublicate found!!!!\n";
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
					//std::cout << "Dublicate found!!!!\n";
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
				std::cout << "The particle is null or 0 tracks!\n";
				return 0.0;
			}
			float p = MathOperator::getModule(particle->getMomentum());
			vector<float> direction = MathOperator::getDirection(particle->getMomentum());
			vector<float> angles = MathOperator::getAngles(direction);
			float accuracy = sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p * pow(sin(angles[1]), 4.0/3.0)) );
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
					std::cout << "A good one found\n";
					_costhetaTestedGood[_ngoodtest] = cos;
					_ngoodtest++;
				}
				if (times == 0) 
				{
					std::cout << "A lost one found\n";

					_costhetaTestedLost[_nlosttest] = cos;
					_nlosttest++;
				}
				if (times>1) 
				{
					std::cout << "Several tracks found\n";
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
			/*for (int i = 0; i < prongs->getNumberOfElements(); i++) 
			{
				MCParticle * mcparticle = dynamic_cast<MCParticle *>(prongs->getElementAt(i));
				int times = 0;
				for (int j = 0; j < pfo->getNumberOfElements(); j++) 
				{
					ReconstructedParticle * recoparticle = dynamic_cast<ReconstructedParticle *>(pfo->getElementAt(j));
					if (CompareParticles(mcparticle, recoparticle)) 
					{
						times++;
					}
				}
				vector<float> direction = MathOperator::getDirection(mcparticle->getMomentum());
				float cos = std::cos(MathOperator::getAngles(direction)[1]);
				if (times == 1) 
				{
					std::cout << "A good one found\n";
					_costhetaTestedGood[_ngoodtest] = cos;
					_ngoodtest++;
				}
				if (times == 0) 
				{
					std::cout << "A lost one found\n";

					_costhetaTestedLost[_nlosttest] = cos;
					_nlosttest++;
				}
				if (times>1) 
				{
					std::cout << "Several tracks found\n";
					_costhetaTestedMultiple[_nmultitest] = cos;
					_nmultitest++;
				}
			}*/
		}
		void VertexRestorer::TestCollections(EVENT::LCCollection * tagged, EVENT::LCCollection * pfo)
		{
			int pfonumber = pfo->getNumberOfElements();
			const vector< ReconstructedParticle * > primaries = myPrimary->getAssociatedParticle()->getParticles();
			for (int i = 0; i < pfonumber; i++) 
			{
				
			}
		}

		void VertexRestorer::check( LCEvent * evt ) {}



		void VertexRestorer::end()
		{   
			hfile->Write();
			hfile->Close();
		}
}
