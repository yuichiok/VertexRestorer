#include "VertexRestorer.hh"

using namespace lcio ;
using namespace marlin ;
using std::vector;
using std::string;
using std::abs;
using EVENT::Track;
using IMPL::ReconstructedParticleImpl;
using IMPL::VertexImpl;
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
		    std::string("BuildUpVertex")
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
	    registerInputCollection( LCIO::MCPARTICLE,
		    "BStarCollectionName" , 
		    "Name of the BStar collection"  ,
		    _colBStarName ,
		    std::string("BStar")
	    );
	    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
	    	"JetRelCollectionName",
		"Name of the Jet relation collection",
		_colJetRelName,
	    	std::string("FinalJets_rel")
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
	   _Tree->Branch("bstarTotal", &_bstarTotal, "bstarTotal/I");
	   
	   _Tree->Branch("missedDetected", &_missedDetected, "missedDetected/I");
	   _Tree->Branch("missedMomentum", _missedMomentum, "missedMomentum[missedDetected]/F");
	   _Tree->Branch("missedError", _missedError, "missedError[missedDetected]/F");
	   _Tree->Branch("missedDeviation", _missedDeviation, "missedDeviation[missedDetected]/F");
	   
	   _Tree->Branch("bstarDetected", &_bstarDetected, "bstarDetected/I");
	   _Tree->Branch("bstarMomentum", _bstarMomentum, "bstarMomentum[bstarDetected]/F");
	   _Tree->Branch("bstarError", _bstarError, "bstarError[bstarDetected]/F");
	   _Tree->Branch("bstarDeviation", _bstarDeviation, "bstarDeviation[bstarDetected]/F");
	   
	   _Tree->Branch("fakeDetected", &_fakeDetected, "fakeDetected/I");
	   _Tree->Branch("fakeMomentum", _fakeMomentum, "fakeMomentum[fakeDetected]/F");
	   _Tree->Branch("fakeError", _fakeError, "fakeError[fakeDetected]/F");
	   _Tree->Branch("fakeDeviation", _fakeDeviation, "fakeDeviation[fakeDetected]/F");

	/*	std::cout << "|\tD0" bstar
		       << "|\tZ0"  
		       << "|\tPhi" 
		       << "|\tOmega" 
		       << "|\tTanl|" 
		       << '\n';*/
	}

	void VertexRestorer::processRunHeader( LCRunHeader* run) 
	{ 
	   _nRun++ ;
	} 
	float VertexRestorer::GetDeviation(const EVENT::ReconstructedParticle * particle, double * pos)
	{
		float p = MathOperator::getModule(particle->getMomentum());
		vector<float> direction = MathOperator::getDirection(particle->getMomentum());
		vector<float> angles = MathOperator::getAngles(direction);
		float accuracy = sqrt(_aParameter*_aParameter + _bParameter*_bParameter /( p * p * pow(sin(angles[1]), 4.0/3.0)) );
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

		void VertexRestorer::processEvent( LCEvent * evt ) 
		{ 
			_nEvt ++ ;
		  
			try 
			{
				std::cout << "********VertexRecovery*"<<_nEvt << "**********\n";
				//LCCollection* maincol = evt->getCollection( _colName );
				LCCollection* seccol = evt->getCollection( _colSecName );
				//vector < ReconstructedParticle * > * result = RestoreVerticesPFO(maincol, seccol);
				LCCollection* pricol = evt->getCollection( _colPriName );
				//LCCollection* jetcol = evt->getCollection( _colJetName );
				//LCCollection* relcol = evt->getCollection( _colJetRelName );
				Vertex * primary = dynamic_cast< Vertex * >(pricol->getElementAt(0));
				vector < ReconstructedParticle * > * result = RestoreVertices(primary, seccol);
				//vector < EVENT::Vertex * > * result = RestoreVertices(jetcol, relcol);
				LCCollection* misscol = evt->getCollection( _colMissedName );
				LCCollection* bstarcol = evt->getCollection( _colBStarName );
				CompareCollections(result, misscol, bstarcol);
				_Tree->Fill();
			}
			catch( DataNotAvailableException &e)
			{
			std::cout << "Whoops!....\n";
			}
		    
		}

	  	void VertexRestorer::CompareCollections(vector< ReconstructedParticle * > * detected, LCCollection * missed, LCCollection * bstar)
		{
			_missedTotal = missed->getNumberOfElements();
			_bstarTotal = bstar->getNumberOfElements();
			_detectedTotal = detected->size();

			_missedDetected = 0;
			_bstarDetected = 0;
			_fakeDetected = 0;

			for (int i = 0; i < _detectedTotal; i++) 
			{
				float taken = false;
				ReconstructedParticle * detectedParticle = detected->at(i);
				for (int j = 0; j < _missedTotal; j++) 
				{
					ReconstructedParticle * particle = dynamic_cast< ReconstructedParticle * > (missed->getElementAt(j));
					if (CompareParticles(particle, detectedParticle)) 
					{
						taken = true;
						//_missedDetected++;
						_missedDeviation[_missedDetected] = GetDeviation(detectedParticle, MathOperator::toDoubleArray(detectedParticle->getStartVertex()->getPosition(), 3));
						_missedError[_missedDetected] = GetError(detectedParticle);
						_missedMomentum[_missedDetected++] = MathOperator::getModule(particle->getMomentum());
						
					}
				}
				for (int j = 0; j < _bstarTotal; j++) 
				{
					MCParticle * particle = dynamic_cast< MCParticle * > (bstar->getElementAt(j));
					if (CompareParticles(particle, detectedParticle)) 
					{
						_bstarError[_bstarDetected] = GetError(detectedParticle);
						//_bstarDetected++;
						_bstarDeviation[_bstarDetected] = GetDeviation(detectedParticle, MathOperator::toDoubleArray(detectedParticle->getStartVertex()->getPosition(), 3));
						taken = true;
						_bstarMomentum[_bstarDetected++] = MathOperator::getModule(particle->getMomentum());
					}
				}
				if (!taken) 
				{
					_fakeError[_fakeDetected] = GetError(detectedParticle);
					_fakeDeviation[_fakeDetected] = GetDeviation(detectedParticle, MathOperator::toDoubleArray(detectedParticle->getStartVertex()->getPosition(), 3));
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
				result->reserve(result->size() + additional.size());
				result->insert(result->end(), additional.begin(),additional.end());
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
		bool VertexRestorer::CompareParticles(const ReconstructedParticle * particle1, const ReconstructedParticle * particle2)
		{
			if (particle1->getCharge() != particle2->getCharge()) 
			{
				return false;
			}
			float angle = MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum());
			if (angle > 0.01) 
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
		bool VertexRestorer::CompareParticles(const MCParticle * particle1, const ReconstructedParticle * particle2)
		{
			if (particle1->getCharge() != particle2->getCharge()) 
			{
				return false;
			}
			float angle = MathOperator::getAngle(particle1->getMomentum(), particle2->getMomentum());
			if (angle > 0.01) 
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

			for (int i = 0; i < number; i++) 
			{
				Vertex * secondary = dynamic_cast< Vertex * >(sec->getElementAt(i));
				vector< ReconstructedParticle * > additional = GetAdditionalParticles(primaries, secondary);
				result->reserve(result->size() + additional.size());
				result->insert(result->end(), additional.begin(),additional.end());
			}
			return result;
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
				vector< float > direction = MathOperator::getDirection(primaryTrack->getMomentum());
				double * pos = MathOperator::toDoubleArray(sec->getPosition(), 3);
				float offset = MathOperator::getDistanceTo(ip, direction, pos);
				float deviation = GetDeviation(primaryTrack, pos);
				float angle = MathOperator::getAngle(primaryTrack->getMomentum(), pos);
				if (abs(primaryTrack->getCharge()) > 0.0 
				    && offset < _offsetCut  
				    && angle < 0.05)
				    //&& deviation < 10
				    //&& !(MathOperator::getModule(primaryTrack->getMomentum()) > 5 && deviation > 6)) 
				{
					std::cout << "Found a track with offset " << offset 
						<< ", deviation " << deviation 
						<< ", angle " << angle
						<< ", error*e-4 " << GetError(primaryTrack) * 10000
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
						if (!IsDublicate(primaryTrack,result)) 
						{
							result.push_back(CopyParticle(primaryTrack,sec));
						}
					}
				}
			}
			return result;
		}
		
		ReconstructedParticle * VertexRestorer::CopyParticle(const ReconstructedParticle * particle, const Vertex * vertex)
		{
			ReconstructedParticleImpl * nouveau = new ReconstructedParticleImpl();
			VertexImpl * vimp = new VertexImpl();
			nouveau->setMomentum(particle->getMomentum());
			nouveau->setMass(particle->getMass());
			nouveau->setCharge(particle->getCharge());
			nouveau->setEnergy(particle->getEnergy());
			vimp->setPosition(vertex->getPosition());
			nouveau->setStartVertex(vimp);
			nouveau->addTrack(particle->getTracks()[0]); //CRUNCH
			return nouveau;
		}
		float VertexRestorer::IsDublicate(const ReconstructedParticle * particle, vector< ReconstructedParticle * > & data)
		{
			bool dublicate = false;
			for (unsigned int j = 0; j < data.size(); j++) 
			{
				if (CompareParticles(particle, data[j])) 
				{
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
			Track * track = particle->getTracks()[0];
			return MathOperator::getTrace(track->getCovMatrix(), 5);
		}

void VertexRestorer::check( LCEvent * evt ) {}



void VertexRestorer::end()
{   
hfile->Write();
hfile->Close();
}
}
