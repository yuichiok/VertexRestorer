#include "RecoveryOperator.hh"
using IMPL::ReconstructedParticleImpl;
using IMPL::VertexImpl;
using IMPL::LCRelationImpl;
using UTIL::LCRelationNavigator;
using EVENT::ReconstructedParticle;
using EVENT::LCCollection;
using EVENT::LCRelation;
using EVENT::Vertex;
using EVENT::Track;
using EVENT::LCObject;
using std::vector;
using std::string;
namespace TTbarAnalysis
{
	RecoveryOperator:: RecoveryOperator(EVENT::Vertex * vertex, EVENT::LCCollection *pfos)
	{
		_aParameter = 0.005;
		_bParameter = 0.01;
		myTotalTracksCounter = 0;
		myPrimary = vertex;
		myPFOs = pfos;
	}
	void RecoveryOperator::PrintParticle(ReconstructedParticle * particle)
	{
		if (!particle) 
		{
			return;
		}
		//streamlog_out(DEBUG) << std::fixed << std::setw( 6 ) << std::setprecision( 3 ) << std::setfill( ' ' );
		int id =  particle->getTracks()[0]->getSubdetectorHitNumbers()[4];
		if (particle->getParticleIDUsed()) 
		{
			streamlog_out(DEBUG) << "Type " << particle->getParticleIDUsed()->getType() << '\n';
			id = particle->getParticleIDs()[0]->getPDG(); 
		}
		float chi2 = (particle->getTracks().size() > 0)? particle->getTracks()[0]->getChi2() / (float)particle->getTracks()[0]->getNdf() : -1.0;
		streamlog_out(DEBUG)<<"|"<< id <<"\t\t|"<<particle->getMass()<<"\t\t|"<<particle->getCharge()  <<"\t\t|"<<particle->getEnergy() <<"\t\t|"<< chi2 <<"\t\t|\n";
	}
	vector< ReconstructedParticle * > RecoveryOperator::getPFOParticles()
	{
		vector< ReconstructedParticle * > result;
		int pfonumber = myPFOs->getNumberOfElements();
		for (int i = 0; i < pfonumber; i++) 
		{
			ReconstructedParticle * pfo = dynamic_cast< ReconstructedParticle * > (myPFOs->getElementAt(i));
			if (std::abs(pfo->getCharge() ) > 0.9) 
			{
				result.push_back(pfo);
			}
		}
		return result;
	}
	
	vector< ReconstructedParticle * > * RecoveryOperator::getVertexParticles(EVENT::LCCollection * secvtx, vector< Vertex * > * tagged)
	{
		vector< ReconstructedParticle * > * particles = new vector< ReconstructedParticle * >();
		int snumber = secvtx->getNumberOfElements();
		for (int i = 0; i < snumber; i++)
		{
		        Vertex * secondary = dynamic_cast< Vertex * >(secvtx->getElementAt(i));
		        tagged->push_back(secondary);
		        vector< ReconstructedParticle * > secondaries = secondary->getAssociatedParticle()->getParticles();
		        particles->reserve(particles->size() + secondaries.size());
		        particles->insert(particles->end(),secondaries.begin(),secondaries.end());
		}
		return particles;
	}
	vector< ReconstructedParticle * > RecoveryOperator::getTrackParticles(LCCollection * trash, vector< ReconstructedParticle * > * secparticles)
	{
		vector< ReconstructedParticle * > result;
		if (!trash) 
		{
			return result;
		}
		int tracknumber = trash->getNumberOfElements();
		for (int i = 0; i < tracknumber; i++) 
		{
			Track * track = dynamic_cast<Track *>(trash->getElementAt(i));
			ReconstructedParticle * particle = myTrackOperator.ReconstructParticle(track);
			if (!secparticles) 
			{
				result.push_back(particle);
				continue;
			}
			if (!IsDublicate(particle, *secparticles)) 
			{
				result.push_back(particle);
			}
			else 
			{
				//streamlog_out(DEBUG) << "Disposing created particle!\n";
				delete particle;
			}
		}
		streamlog_out(DEBUG) << "\n";
		return result;
	}
	vector< Vertex * > RecoveryOperator::RecoverJetVertices(LCCollection * jetcol, LCCollection * jetrelcol, LCCollection * secvtx, LCCollection * damagedcol, LCCollection * newjetrelcol)
	{
		vector< Vertex * > result;
		int jnumber = jetcol->getNumberOfElements();
		LCRelationNavigator navigator(jetrelcol);
		vector< Vertex * > * tagged = new vector< Vertex * > ();
		vector< ReconstructedParticle * > * particles = getVertexParticles(secvtx, tagged);
		vector< ReconstructedParticle * > zombies = getTrackParticles(damagedcol, particles); // CHANGED
		vector< ReconstructedParticle * > pfo = getPFOParticles();
		vector< ReconstructedParticle * > taken;
		vector< ReconstructedParticle * > toDispose;
		for (int i = 0; i < jnumber; i++) 
		{
			ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle *>(jetcol->getElementAt(i));
			int nvtx = navigator.getRelatedToObjects(jet).size();
			streamlog_out(DEBUG) << "Jet energy: " << jet->getEnergy() << "\n";
			vector< LCObject * > objs = navigator.getRelatedToObjects(jet);
			vector< Vertex * > oldvtx;
			float pHadron = getHadronMomentum(objs);
			float pcut = 80;
			for (int j = 0; j < nvtx; j++) 
			{
				Vertex * vertex = dynamic_cast< Vertex * >(objs[j]);
				streamlog_out(DEBUG) << "Vertex distance: " << MathOperator::getModule(vertex->getPosition()) 
					<< " tracks " << vertex->getAssociatedParticle()->getParticles().size() 
					<< " charge " << vertex->getAssociatedParticle()->getCharge() 
					<< " chi2: " << vertex->getChi2()
					<< " prob: " << vertex->getProbability()
					<< " algo: " << vertex->getAlgorithmType()
					<< " |p|: " << pHadron
					<< "\n";
				vector< ReconstructedParticle * > toInject;
				//toInject.reserve(jet->getParticles().size()+ zombies.size());
				toInject.reserve(pfo.size()+ zombies.size());
				
				//toInject.insert(toInject.end(), jet->getParticles().begin(), jet->getParticles().end());
				toInject.insert(toInject.end(), pfo.begin(), pfo.end());
				toInject.insert(toInject.end(), zombies.begin(), zombies.end());

				vector< ReconstructedParticle * > additional;
				if (pHadron < pcut) 
				{
					additional= AddParticles(toInject, vertex, particles, tagged);
				}
				vector< ReconstructedParticle * > filtered;
				for (unsigned int k = 0; k < additional.size(); k++) 
				{
					if (!IsDublicate(additional[k], taken)) 
					{
						PrintParticle(additional[k]);
						myTotalTracksCounter++;
						filtered.push_back(additional[k]);
						taken.push_back(additional[k]);
					}
				}
				Vertex * newvertex = CreateRecoveredVertex(filtered, vertex, true);
				streamlog_out(DEBUG) << "Vertex distance: " << MathOperator::getModule(newvertex->getPosition()) 
					<< " tracks " << newvertex->getAssociatedParticle()->getParticles().size() 
					<< " charge " << newvertex->getAssociatedParticle()->getCharge()
					<< "\n"
					<< "\n";
				//result.push_back(newvertex);
				oldvtx.push_back(newvertex);
				if (newjetrelcol) 
				{
					//newjetrelcol->addElement(CreateNewRelation(newvertex, jet));
				}
			}
			vector< Vertex * > newvtx = oldvtx;//RefineJetVertices(oldvtx);
			for (int j = 0; j < newvtx.size(); j++) 
			{ 
				result.push_back(newvtx[j]);	
				if (newjetrelcol) 
				{
					newjetrelcol->addElement(CreateNewRelation(newvtx[j], jet));
				}
			}
		}
		for (unsigned int k = 0; k < zombies.size(); k++) 
		{
			if (!IsDublicate(zombies[k], taken)) 
			{
				//streamlog_out(DEBUG) << "Disposing created particle!\n";
				ReconstructedParticle * zombie = zombies[k];
				delete zombie;
			}
		}
		return result;
	}
	vector< Vertex * > RecoveryOperator::RecoverBuildVertices(LCCollection * secvtx,LCCollection * damagedcol)
	{
		vector< Vertex * > result;
		int secnumber = secvtx->getNumberOfElements();
		vector< Vertex * > * tagged = new vector< Vertex * > ();
		vector< ReconstructedParticle * > * particles = getVertexParticles(secvtx, tagged);
		vector< ReconstructedParticle * > zombies = getTrackParticles(damagedcol, particles);
		vector< ReconstructedParticle * > pfo = getPFOParticles();
		vector< ReconstructedParticle * > taken;
		for (int i = 0; i < secnumber; i++) 
		{
			Vertex * vertex = dynamic_cast< Vertex * >(secvtx->getElementAt(i));
			streamlog_out(DEBUG) << "Vertex distance: " << MathOperator::getModule(vertex->getPosition()) 
				<< " tracks " << vertex->getAssociatedParticle()->getParticles().size() 
				<< " charge " << vertex->getAssociatedParticle()->getCharge() 
				<< " chi2: " << vertex->getChi2()
				<< " prob: " << vertex->getProbability()
				<< "\n";
			vector< ReconstructedParticle * > toInject;
			toInject.reserve(pfo.size() + zombies.size());
			toInject.insert(toInject.end(), zombies.begin(), zombies.end());
			toInject.insert(toInject.end(), pfo.begin(), pfo.end());

			vector< ReconstructedParticle * > additional = AddParticles(toInject, vertex, particles, tagged);
			vector< ReconstructedParticle * > filtered;
			for (unsigned int k = 0; k < additional.size(); k++) 
			{
				if (!IsDublicate(additional[k], taken)) 
				{
					PrintParticle(additional[k]);
					myTotalTracksCounter++;
					filtered.push_back(additional[k]);
					taken.push_back(additional[k]);
				}
			}
			Vertex * newvertex = CreateRecoveredVertex(filtered, vertex);
			streamlog_out(DEBUG) << "Vertex distance: " << MathOperator::getModule(newvertex->getPosition()) 
				<< " tracks " << newvertex->getAssociatedParticle()->getParticles().size() 
				<< " charge " << newvertex->getAssociatedParticle()->getCharge()
				<< "\n"
				<< "\n";
			result.push_back(newvertex);
		}
		return result;
	}
	vector< ReconstructedParticle * > RecoveryOperator::AddParticles(const vector< ReconstructedParticle * > & pri, Vertex * sec, const vector< ReconstructedParticle * > * toCompare, vector< Vertex * > * allVtx)
	{
		vector< ReconstructedParticle * > result;
		if (sec->getAssociatedParticle()->getParticles().size() < 1) 
		{
			return result;
		}
		const vector< ReconstructedParticle * > * particles = (toCompare)? toCompare : &(sec->getAssociatedParticle()->getParticles());
		for (unsigned int i = 0; i < pri.size(); i++) 
		{
			ReconstructedParticle * candidate = pri[i];
			//if (TakeParticle(candidate, sec) && IsMinimalAngle(candidate, sec, allVtx) && !IsDublicate(candidate, *particles)) 
			if (TakeParticle(candidate, sec) && !IsDublicate(candidate, *particles)) 
			//if (TakeParticle(candidate, sec)) 
			{
				if (!IsDublicate(candidate, result)) 
				{
					result.push_back(candidate);
				}
			}
		}
		return result;
		
	}
	bool RecoveryOperator::TakeParticle(EVENT::ReconstructedParticle * primary, const EVENT::Vertex * sec)
	{
		if (std::abs(primary->getCharge() ) < 0.9) 
		{
			return false;
		}
		const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
		double * trackPosition = myTrackOperator.GetStartPoint(primary);
		float trackDistance = MathOperator::getModule(trackPosition);
		if (trackDistance > 20.0) 
		{
			return false;
		}
		//return true;
		vector< float > direction = MathOperator::getDirection(primary->getMomentum());
		double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
		double * secondaryPosition = MathOperator::toDoubleArray(sec->getPosition(),3);
		float primaryOffset = MathOperator::getDistanceTo(primaryPosition, direction, trackPosition);
		//float accuracy = GetError(primary);
		float accuracy = myTrackOperator.GetOffsetErrorSimple(primary);
		//float secondaryOffset = MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
		vector<float> diff = MathOperator::getDirection(secondaryPosition, trackPosition);
		double diif[3];
		double aright[3];
		float observable =0; //MathOperator::getModule(*MathOperator::vectorProduct(directionVtx, direction));//MathOperator::getAngle(primary->getMomentum(), secondaryPosition);
		float cosbeta = std::cos(MathOperator::getAngle(trackPosition, primary->getMomentum()));
		vector< float > directionVtx = MathOperator::getDirection(secondaryPosition);
		for (int m = 0; m < 3; m++) 
		{
			diif[m] = diff[m];
			aright[m] = cosbeta * MathOperator::getModule(trackPosition) * direction[m] - trackPosition[m];
			observable += aright[m] * directionVtx[m];
		}
		float angle = MathOperator::getAngle(diif, primary->getMomentum());
		//float secOffset =  MathOperator::getDistanceTo(secondaryPosition, direction, trackPosition);
		//float costheta = std::cos(MathOperator::getAngles(direction)[1]); //std::abs( std::cos(MathOperator::getAngles(secDiraction)[1] ) );
		//float dprime = myTrackOperator.GetDprime(primary, secondaries[0], primaryPosition);
		//float l = 0.0;// GetMinDiffDistance(primary, sec, dprime);// std::abs( distance / secondaryOffset -0.5 ) * cos;
		//vector< float > limits = ParametrizeVertex(sec);
		int vtxhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[0];
		int ftdhits = primary->getTracks()[0]->getSubdetectorHitNumbers()[5];
		float costheta = std::cos(MathOperator::getAngles(direction)[1]); //std::abs( std::cos(MathOperator::getAngles(secDiraction)[1] ) );
		float angleError = 1.0;
		if (angle > 0.0) 
		{
			angleError = std::sqrt(myTrackOperator.GetAngleError(angle, sec, primary));
		}
		if ( std::abs(costheta) > 0.9 && vtxhits == 0 ) 
		{
			//accuracy = 4*GetError(primary);
		}
		//float p = MathOperator::getModule(sec->getAssociatedParticle()->getMomentum());
		//float anglecut = 0.1 - 0.1/250 * p;
		//float anglecut = 0.075 - 0.1*std::atan(p/5.0 - 10.0)/2.0/3.14;
		//float anglecut = 0.1/1.78 - 0.1*std::atan(p/40.0 - 1.0)/1.7814;
		float p = MathOperator::getModule(primary->getMomentum());
		float anglecut = 0.08; // 0.08
		float deviation = myTrackOperator.GetOffsetSignificance(primary);
		//float deviation = primaryOffset /accuracy;
		//bool result = (primaryOffset /accuracy  > 110.0 * angle + 0.2  || angle < 0.005) &&
		//bool result = (primaryOffset /accuracy  > 80.0 * angle + 0.5 || angle < 0.005) && // BEST
		//bool result = (primaryOffset /accuracy  > 50.0 * angle + 2. || angle < 0.005) && // BETTER
		//bool result = ( deviation > 50.0 * angle + 2. || angle < 0.005) &&
		bool result = ( deviation > 25*sqrt(angle)+1.0 || angle < 0.001) &&
				observable < 0.01 &&
		//bool result = (primaryOffset /accuracy  > 15.0 * sqrt( angle )  || angle < 0.001) &&
		//bool result = (primaryOffset /accuracy  > 2.2 *std::atan(angle * 100) || angle < 0.005)   &&
		//bool result = (primaryOffset /accuracy  > 17.0 * sqrt( angle )+0.4  || angle < 0.001) &&
			 //(vtxhits > 2 || angle < 0.001) &&
			//(((ftdhits > 0) && angle < anglecut) || (vtxhits > 1 && angle < anglecut)) &&
			//p > 0.3 &&

			//(ftdhits > 0) && angle < anglecut/2;
		//bool result =  (primaryOffset /accuracy  > 0.7 + 1.5 * angle / angleError) &&
			 angle < anglecut;// + 0.03 * sine;
		if (result) 
		{
			/*streamlog_out(DEBUG) << "Found a track with offset " << primaryOffset
				<< ", error " << accuracy//GetError(primary) 
				<< ", Angle " << angle //GetError(primary) 
				//<< ", DISTANCE " << dprime - MathOperator::getModule(secondaryPosition) //distance//GetError(primary) 
				<< ", position " << trackDistance
				<<" :\n";//*/
		}
		delete primaryPosition;
		delete secondaryPosition;
		delete trackPosition;
		return result;
	}
	Vertex * RecoveryOperator::CreateRecoveredVertex(vector< ReconstructedParticle * > & newtracks, Vertex * oldvertex, bool add)
	{
		//streamlog_out(DEBUG) << "Start to create new vertex\n";
		VertexImpl * newvertex = new VertexImpl();
		newvertex->setPrimary(false);
		newvertex->setChi2(oldvertex->getChi2());
		newvertex->setPosition(oldvertex->getPosition());
		newvertex->setProbability(oldvertex->getProbability());
		newvertex->setCovMatrix(oldvertex->getCovMatrix());
		ReconstructedParticle * oldparticle = oldvertex->getAssociatedParticle();
		ReconstructedParticleImpl * newparticle = new ReconstructedParticleImpl();
		double newmomentum[3];
		float newenergy = (!add)?0:oldparticle->getEnergy();
		float newcharge = (!add)?0:oldparticle->getCharge();
		newmomentum[0] = (!add)?0:oldparticle->getMomentum()[0];
		newmomentum[1] = (!add)?0:oldparticle->getMomentum()[1];
		newmomentum[2] = (!add)?0:oldparticle->getMomentum()[2];
		//streamlog_out(DEBUG) << "\tAssembling all prongs together\n";
		if (add) 
		{
			for (unsigned int i = 0; i < oldparticle->getParticles().size(); i++) 
			{
				newparticle->addParticle(oldparticle->getParticles()[i]);
			}
		}
		for (unsigned int i = 0; i < newtracks.size(); i++) 
		{
			newparticle->addParticle(newtracks.at(i));
			newmomentum[0] += newtracks[i]->getMomentum()[0];
			newmomentum[1] += newtracks[i]->getMomentum()[1];
			newmomentum[2] += newtracks[i]->getMomentum()[2];
			newenergy += newtracks[i]->getEnergy();
			newcharge += newtracks[i]->getCharge();
			//Add tracks
		}
		float newmass = std::sqrt(newenergy * newenergy - newmomentum[0]*newmomentum[0] -newmomentum[1]*newmomentum[1]-newmomentum[2]*newmomentum[2] );
		newparticle->setMomentum(newmomentum);
		newparticle->setCharge(newcharge);
		newparticle->setMass(newmass);
		newparticle->setEnergy(newenergy);
		//streamlog_out(DEBUG) << "\tSetting associated particle\n";
		newvertex->setAssociatedParticle(newparticle);
		newvertex->setAlgorithmType("lcfiplus");
		newvertex->addParameter (newtracks.size());
		streamlog_out(DEBUG) << "Finished vertex with " << newtracks.size() << " new tracks, "
			<< newparticle->getMass() - oldparticle->getMass()  << " mass gain, "
			<< newparticle->getCharge() - oldparticle->getCharge() << " charge difference\n";
		return newvertex;	
	}
	LCRelation * RecoveryOperator::CreateNewRelation(Vertex * newvertex, ReconstructedParticle * oldjet)
	{
		LCRelationImpl * relation = new LCRelationImpl();
		relation->setFrom(oldjet);
		relation->setTo(newvertex);
		relation->setWeight(1.0);
		return relation;
	}
	float RecoveryOperator::GetError(const ReconstructedParticle * particle)
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
		return accuracy;
	}
	bool RecoveryOperator::IsDublicate(const ReconstructedParticle * particle, const vector< ReconstructedParticle * > & data)
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
	bool RecoveryOperator::CompareParticles(const ReconstructedParticle * particle1, const ReconstructedParticle * particle2)
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
				return angle < 0.005;
			}
		}
		if (ratio > 0.1) 
		{
			return false;
		}
		return true;
	}
	vector< float > RecoveryOperator::ParametrizeVertex(const Vertex * sec)
	{
		const vector< ReconstructedParticle * > secondaries = sec->getAssociatedParticle()->getParticles();
		vector< float > result;
		if (secondaries.size() < 2) 
		{
			return result;
		}
		float mindprime = 1000;
		float maxdprime = 0;
		//double * primaryPosition = MathOperator::toDoubleArray(myPrimary->getPosition(),3);
		//streamlog_out(DEBUG) << "Vertex position: " << MathOperator::getModule(sec->getPosition()) << " Chi2: " << sec->getChi2() << "\n";
		for (unsigned int i = 0; i < secondaries.size(); i++) 
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
		result.push_back(mindprime);
		result.push_back(maxdprime);
		return result;
		
	}
	bool RecoveryOperator::IsMinimalAngle(ReconstructedParticle * candidate, Vertex * chosen, vector< Vertex * > * vertices)
	{
		double * canPos = myTrackOperator.GetStartPoint(candidate);
		double * chosenPos = MathOperator::toDoubleArray(chosen->getPosition(),3);
		vector<float> cdiff = MathOperator::getDirection(chosenPos, canPos);
		double cdiif[3];
		for (int m = 0; m < 3; m++)
		{
		        cdiif[m] = cdiff[m];
		}
		float chosenangle =  MathOperator::getAngle(cdiif, candidate->getMomentum());
		float minangle = 1.51;
		for (unsigned int i = 0; i < vertices->size(); i++) 
		{
			Vertex * sec = vertices->at(i);
			if (sec == chosen || sec->getAssociatedParticle()->getParticles().size() == 1) 
			{
				//streamlog_out(DEBUG) << "Skipping\n";
				continue;
			}
			double * vtxPos = MathOperator::toDoubleArray(sec->getPosition(),3);
			vector<float> diff = MathOperator::getDirection(vtxPos, canPos);
			double diif[3];
			for (int m = 0; m < 3; m++)
			{
			        diif[m] = diff[m];
			}
			float angle =  MathOperator::getAngle(diif, candidate->getMomentum());
			if (angle < minangle) 
			{
				minangle = angle;
			}
			delete vtxPos;
		}
		delete canPos;
		delete chosenPos;
		//streamlog_out(DEBUG) << "Chosen angle: " << chosenangle << " minangle: " << minangle << "\n";
		return chosenangle < minangle;
	}
	int RecoveryOperator::GetStatistics()
	{
		return myTotalTracksCounter;
	}
	vector< Vertex * > RecoveryOperator::RefineJetVertices(vector< Vertex * > & oldvtx)//ReconstructedParticle * oldjet, EVENT::LCCollection * jetrelcol)
	{
		//LCRelationNavigator navigator(jetrelcol);
		//streamlog_out(DEBUG) << "Jet energy: " << jet->getEnergy() << "\n";
		//vector< LCObject * > objs = navigator.getRelatedToObjects(jet);
		//vector< Vertex * > oldvtx;
		int nvtx = oldvtx.size();
		vector<ReconstructedParticle *> trackpool;
		for (int i = 0; i < nvtx; i++) 
		{
			Vertex * vertex = oldvtx[i];
			trackpool.reserve(trackpool.size() + vertex->getAssociatedParticle()->getParticles().size());
			trackpool.insert(trackpool.end(), vertex->getAssociatedParticle()->getParticles().begin(), vertex->getAssociatedParticle()->getParticles().end());
		}	
		if (nvtx != 2) 
		{
			return oldvtx;
		}
		vector<ReconstructedParticle *> newvtx1particles;
		vector<ReconstructedParticle *> newvtx2particles;
		double * secondaryPosition1 = MathOperator::toDoubleArray(oldvtx[0]->getPosition(),3);
		double * secondaryPosition2 = MathOperator::toDoubleArray(oldvtx[1]->getPosition(),3);
		for (int i = 0; i < trackpool.size(); i++) 
		{
			ReconstructedParticle * particle = trackpool[i];
			double * trackPosition = myTrackOperator.GetStartPoint(particle);
			vector<float> diff1 = MathOperator::getDirection(secondaryPosition1, trackPosition);
			vector<float> diff2 = MathOperator::getDirection(secondaryPosition2, trackPosition);
			double diif1[3];
			double diif2[3];
			for (int m = 0; m < 3; m++) 
			{
				diif1[m] = diff1[m];
				diif2[m] = diff2[m];
			}
			float angle1 = MathOperator::getAngle(diif1, particle->getMomentum());
			float angle2 = MathOperator::getAngle(diif2, particle->getMomentum());
			if (angle1 < angle2) 
			{
				newvtx1particles.push_back(particle);
			}
			else 
			{
				newvtx2particles.push_back(particle);
			}
		}
		vector<Vertex * > result;
		result.push_back(CreateRecoveredVertex(newvtx1particles, oldvtx[0], false));
		result.push_back(CreateRecoveredVertex(newvtx2particles, oldvtx[1], false));
		streamlog_out(DEBUG) << "Momentum: " << MathOperator::getModule(result[0]->getAssociatedParticle()->getMomentum()) << "\n";
		return result;
	}
	float RecoveryOperator::getHadronMomentum(vector< LCObject * > & objs)
	{
		float momentum[3];
		momentum[0] = 0;
		momentum[1] = 0;
		momentum[2] = 0;
		for (unsigned int i = 0; i < objs.size(); i++) 
		{
			Vertex * vertex = dynamic_cast< Vertex * >(objs[i]);
			for (unsigned int j = 0; j < 3; j++) 
			{
				momentum[j] += vertex->getAssociatedParticle()->getMomentum()[j];
			}
		}
		float module = MathOperator::getModule(momentum);
		return module;
	}


}
