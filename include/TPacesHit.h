#ifndef PACESHIT_H
#define PACESHIT_H

#include <cstdio>
#include <cmath>

#include "TFragment.h"
#include "TChannel.h"

#include "TVector3.h"

#include "TGRSIDetectorHit.h"


class TPacesHit : public TGRSIDetectorHit {
	public:
		TPacesHit();
		TPacesHit(const TPacesHit&);
		virtual ~TPacesHit();

	private:
      Int_t filter;

	public:
		/////////////////////////  Setters	/////////////////////////////////////
      inline void SetFilterPattern(const int &x)   { filter = x;   }                  //! 
      //void SetHit();

      TVector3 GetPosition(Double_t dist = 0.0) const; //!

		/////////////////////////  Getters	/////////////////////////////////////
      inline Int_t    GetFilterPattern() const         {   return filter;   }          //!

		/////////////////////////  TChannel Helpers /////////////////////////////////////
      bool   InFilter(Int_t);  //!

	public:
		virtual void Clear(Option_t *opt = "");		 //!
		virtual void Print(Option_t *opt = "") const; //!
      virtual void Copy(TPacesHit&) const;        //!

	ClassDef(TPacesHit,3);
};




#endif
