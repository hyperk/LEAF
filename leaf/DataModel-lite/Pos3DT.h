#ifndef Pos3DT_H
#define Pos3DT_H


#include "Environments.h"

/**
* \class Pos3DT
 *
 * This class holds Position information (aims to replace Pos3DT)
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/12/01 17:12:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */
 

class Pos3DT {
	public:
		Pos3DT();
		~Pos3DT();
		
		
		double &operator[](int n) 		{ return Vtx[n];			}
		double &x()				{ return Vtx[0];			}
		double &y()				{ return Vtx[1];			}
		double &z()				{ return Vtx[2];			}
		double &t()				{ return Vtx[3];			}
		
		double R()				{ return Astro_GetLength(Vtx); 	} 
		double Radius()			{ return Astro_GetR(Vtx); 		} 
		
		void clear()				{ Vtx.clear(); Vtx.resize(4,0.);	}
		
		std::vector<double> Vtx;
	
};

#endif
