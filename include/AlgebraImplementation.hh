#ifndef ALGEBRAIMPLEMENTATION_H_1
#define ALGEBRAIMPLEMENTATION_H_1


#include "lcio.h"
#include <EVENT/ReconstructedParticle.h>
#include "TMatrixD.h"

namespace TTbarAnalysis
{
	class AlgebraImplementation 
	{
		public:
		//
		//	Constants
		//
	
		//
		//	Constructors
		//
			AlgebraImplementation (){};
			virtual ~AlgebraImplementation () {};
		//
		//	Methods
		//
			static int GetCovMatrixMomenta(EVENT::ReconstructedParticle const *, TMatrixD &);
			static int GetCovMatrixMomenta(EVENT::ReconstructedParticle const *, EVENT::FloatVec &);
		private:
		//
		//	Data
		//
			
		//
		//	Private methods
		//
	};
}
#endif
