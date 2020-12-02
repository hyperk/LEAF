#ifndef RootPosition_H
#define RootPosition_H


#include "TObject.h"
#include "Environments.h"

/**
* \class RootPosition
 *
 * This class holds Position information (aims to replace Pos3D)
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/12/01 17:12:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */
 

class RootPosition : public TObject {
	public:
		RootPosition();
		~RootPosition();
		
		
		double &operator[](int n) 		{ return Vtx[n];			}
		double &x()				{ return Vtx[0];			}
		double &y()				{ return Vtx[1];			}
		double &z()				{ return Vtx[2];			}
		double &t()				{ return Vtx[3];			}
		
		double R()				{ return Astro_GetLength(Vtx); 	} 
		double Radius()			{ return Astro_GetR(Vtx); 		} 
		
		void clear()				{ Vtx.clear(); Vtx.resize(4,0.);	}
		
		
	private:
		
		
		std::vector<double> Vtx;
		
		ClassDef(RootPosition,1) //EventRootInfo structure
};

#if !defined(__CLING__)
ClassImp(RootPosition)
#endif

#endif
