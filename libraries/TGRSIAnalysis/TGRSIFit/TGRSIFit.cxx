#include "TGRSIFit.h"

ClassImp(TGRSIFit);

TString TGRSIFit::fDefaultFitType("");

TGRSIFit::TGRSIFit(){
   this->Clear();
}

TGRSIFit::TGRSIFit(const TGRSIFit &copy) : TF1(copy){
   ((TGRSIFit&)copy).Copy(*this);
}

TGRSIFit::~TGRSIFit(){
   this->AddToGlobalList(kFALSE);
}

void TGRSIFit::Copy(TObject &obj) const{
   ((TGRSIFit&)obj).init_flag   = init_flag;
   ((TGRSIFit&)obj).goodfit_flag= goodfit_flag;
   TF1::Copy(obj);
}

void TGRSIFit::Print(Option_t *opt) const {
   if(strchr(opt,'+') != NULL){
      printf("Params Init: %d\n", init_flag);
      printf("Good Fit:    %d\n", goodfit_flag);
      TNamed::Print(opt);
   }
}

void TGRSIFit::Clear(Option_t *opt) {
   init_flag = false;
   goodfit_flag = false;
   fDefaultFitType.Clear();
}

void TGRSIFit::ClearParameters(Option_t *opt){
   for(int i =0; i< GetNpar();++i){
      SetParameter(i,0);
   }
}

Bool_t TGRSIFit::AddToGlobalList(Bool_t on){
   // Add to global list of functions (gROOT->GetListOfFunctions() )
   // return previous status (true of functions was already in the list false if not)
   
   if (!gROOT) return false;

   if (on )  {
      if(gROOT->GetListOfFunctions()->FindObject(this) != nullptr){
         return true;
      }
      gROOT->GetListOfFunctions()->Add(this);
      // do I need to delete previous one with the same name ???
      //TF1 * old = dynamic_cast<TF1*>( gROOT->GetListOfFunctions()->FindObject(GetName()) );
      //if (old) gROOT->GetListOfFunctions()->Remove(old);
      return false;
   }
   else {
      // if previous status was on and now is off
      TF1 * old = dynamic_cast<TF1*>( gROOT->GetListOfFunctions()->FindObject(this) );
      if (!old) {
         //Warning("AddToGlobalList","Function is supposed to be in the global list but it is not present");
         return false;
      }
      gROOT->GetListOfFunctions()->Remove(this);
      return true;
   }
   return true;
}

Bool_t TGRSIFit::AddToGlobalList(TF1* func, Bool_t on){
   // Add to global list of functions (gROOT->GetListOfFunctions() )
   // return previous status (true of functions was already in the list false if not)
   
   if (!gROOT) return false;

   if (on )  {
      if(gROOT->GetListOfFunctions()->FindObject(func) != nullptr){
         return true;
      }
      gROOT->GetListOfFunctions()->Add(func);
      // do I need to delete previous one with the same name ???
      //TF1 * old = dynamic_cast<TF1*>( gROOT->GetListOfFunctions()->FindObject(GetName()) );
      //if (old) gROOT->GetListOfFunctions()->Remove(old);
      return false;
   }
   else {
      // if previous status was on and now is off
      TF1 * old = dynamic_cast<TF1*>( gROOT->GetListOfFunctions()->FindObject(func) );
      if (!old) {
         //func->Warning("AddToGlobalList","Function is supposed to be in the global list but it is not present");
         return false;
      }
      gROOT->GetListOfFunctions()->Remove(func);
      return true;
   }
   return true;
}


