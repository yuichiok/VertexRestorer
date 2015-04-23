#include "CheatOperator.hh"
using std::vector;
using std::string;
using std::abs;
using IMPL::ReconstructedParticleImpl;
using EVENT::ReconstructedParticle;
using IMPL::VertexImpl;
using UTIL::LCRelationNavigator;
using EVENT::Vertex;
using EVENT::MCParticle;
using EVENT::LCObject;
using EVENT::LCCollection;
namespace TTbarAnalysis
{
	CheatOperator:: CheatOperator(LCCollection * prongs, LCCollection * rel)
	{
		myProngs = prongs;
		myNavigator = new LCRelationNavigator(rel);
	}
	
	ReconstructedParticle * CheatOperator::_GetParticle(MCParticle * particle)
	{
		vector< LCObject * > obj = myNavigator->getRelatedFromObjects(particle);
		vector< float > weights = myNavigator->getRelatedToWeights(particle);
		ReconstructedParticle * winner = NULL;
		float maxweight = 0.0;
		for (unsigned int i = 0; i < obj.size(); i++) 
		{
			ReconstructedParticle * candidate = dynamic_cast< ReconstructedParticle * >(obj[i]);
			if (weights[i] > maxweight && std::abs(candidate->getCharge()) > 0.1) 
			{
				winner = candidate;
				maxweight = weights[i];
			}
		}
		return winner;
	}

	vector< Vertex * > * CheatOperator::_Construct()
	{
		vector< Vertex * > * result = new vector< Vertex * > ();
		vector< MCParticle * > genprongs;
		vector< int > trackIDs;
		myProngs->getParameters().getIntVals("trackIDs",trackIDs);
		vector< int > uniqueIDs;
		getUniqueIDs(trackIDs, uniqueIDs);
		for (int i = 0; i < myProngs->getNumberOfElements(); i++) 
		{
			genprongs.push_back(dynamic_cast< MCParticle * >(myProngs->getElementAt(i)));
		}
		for (unsigned int i = 0; i < uniqueIDs.size(); i++) 
		{
			std::cout << "Unique id: " << uniqueIDs[i];

		}
		return result;
	}

	void CheatOperator::getUniqueIDs(std::vector< int > & total, std::vector< int > & unique)
	{
		for (unsigned int i = 0; i < total.size(); i++) 
		{
			if (unique.size() == 0) 
			{
				unique.push_back(total[i]);
			}
			bool found = false;
			for (unsigned int j = 0; j < unique.size(); j++) 
			{
				if (unique[i] == total[i]) 
				{
					found = true;
					break;
				}
			}
			if (!found) 
			{
				unique.push_back(total[i]);
			}
		}
	}
	Vertex * construct(vector< MCParticle * > & genprongs, vector< int > & ids, int type)
	{
		VertexImpl * result = new VertexImpl();
		result->setPrimary(false);
		result->setAlgorithmType("cheatgenvertex");

		for (int i = 0; i < ids.size(); i++) 
		{
			if (ids[i] == type) 
			{
				MCParticle * particle = genprongs[i];
				
			}
		}
	}
}
