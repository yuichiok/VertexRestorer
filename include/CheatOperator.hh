#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCCollection.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/VertexImpl.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#ifndef CheatOperator_hh
#define CheatOperator_hh 1
namespace TTbarAnalysis
{
	class CheatOperator 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			CheatOperator (EVENT::LCCollection * prongs, EVENT::LCCollection * rel);
			virtual ~CheatOperator () 
			{
				delete myNavigator;
			};
		//
		//	Methods
		//
			std::vector< EVENT::Vertex * > * _Construct();	
			EVENT::ReconstructedParticle * _GetParticle(EVENT::MCParticle * particle);
		private:
		//
		//	Data
		//
			EVENT::LCCollection * myProngs;
			UTIL::LCRelationNavigator * myNavigator;
		//
		//	Private methods
		//i
			void getUniqueIDs(std::vector< int > & total, std::vector< int > & unique);
			EVENT::Vertex * construct(std::vector< EVENT::MCParticle * > & genprongs, std::vector< int > & ids, int type);
	};
}
#endif
