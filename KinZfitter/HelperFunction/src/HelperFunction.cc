// -*- C++ -*-
//
// Package:     Subsystem/Package
// Class  :     HelperFunction
// 
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Tongguang Cheng
//         Created:  Mon, 21 Dec 2015 12:47:33 GMT
//

#ifndef HelperFunction_CC
#define HelperFunction_CC
// system include files

// user include files
#include "KinZfitter/HelperFunction/interface/HelperFunction.h"

// fileinPath
#include "FWCore/ParameterSet/interface/FileInPath.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HelperFunction::HelperFunction(bool isData)
{

        //declarations
        debug_ = 0;

        /* Legacy 7+8 TeV
        TString fmu_s = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/ebeOverallCorrections.Legacy2013.v0.root" ).fullPath());
        TString fel_s = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/ebeOverallCorrections.Legacy2013.v0.root" ).fullPath());

        fel = boost::shared_ptr<TFile>( new TFile(fel_s));
        fmu = boost::shared_ptr<TFile>( new TFile(fmu_s));

        muon_corr_data = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fmu->Get( "mu_reco53x" )->Clone() )) );
        muon_corr_mc = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fmu->Get( "mu_mc53x" )->Clone() )) );
                
        electron_corr_data = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fel->Get( "el_reco53x" )->Clone() )) );
        electron_corr_mc = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fel->Get( "el_mc53x" )->Clone() )) );
        */

        // ICHEP 2016 13 TeV
        /*
        TString fel_s_mc = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e.root" ).fullPath());
        TString fmu_s_mc = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2muLUT_m2mu.root" ).fullPath());
        TString fel_s_data = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/DoubleLepton_m2eLUT_m2e.root" ).fullPath());
        TString fmu_s_data = TString(edm::FileInPath ( "KinZfitter/HelperFunction/hists/DoubleLepton_m2muLUT_m2mu.root" ).fullPath());

        fel_mc = boost::shared_ptr<TFile>( new TFile(fel_s_mc));
        fmu_mc = boost::shared_ptr<TFile>( new TFile(fmu_s_mc));
        fel_data = boost::shared_ptr<TFile>( new TFile(fel_s_data));
        fmu_data = boost::shared_ptr<TFile>( new TFile(fmu_s_data));
                
        muon_corr_data = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fmu_data->Get( "2mu" )->Clone() )) );
        muon_corr_mc = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fmu_mc->Get( "2mu" )->Clone() )) );
                
        electron_corr_data = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fel_data->Get( "2e" )->Clone() )) );
        electron_corr_mc = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(fel_mc->Get( "2e" )->Clone() )) );
        */

        // MORIOND 17
        TString s_corr_e_1_mc = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_1.root" ).fullPath());
        TString s_corr_e_2_mc = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_2.root" ).fullPath());
        TString s_corr_e_3_mc = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2eLUT_m2e_3.root" ).fullPath());

        TString s_corr_mu_mc = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DYJetsToLL_M-50_m2muLUT_m2mu.root" ).fullPath());

        TString s_corr_e_1_data = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DoubleLepton_m2eLUT_m2e_1.root" ).fullPath());
        TString s_corr_e_2_data = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DoubleLepton_m2eLUT_m2e_2.root" ).fullPath());
        TString s_corr_e_3_data = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DoubleLepton_m2eLUT_m2e_3.root" ).fullPath());

        TString s_corr_mu_data = TString(edm::FileInPath ("KinZfitter/HelperFunction/hists/DoubleLepton_m2muLUT_m2mu.root" ).fullPath());

        TString s_corr_e_1, s_corr_e_2, s_corr_e_3, s_corr_mu;

        if (isData) {

           s_corr_e_1 = s_corr_e_1_data;
           s_corr_e_2 = s_corr_e_2_data;
           s_corr_e_3 = s_corr_e_3_data;
           s_corr_mu = s_corr_mu_data;

           } else {

                  s_corr_e_1 = s_corr_e_1_mc;
                  s_corr_e_2 = s_corr_e_2_mc;
                  s_corr_e_3 = s_corr_e_3_mc;
                  s_corr_mu = s_corr_mu_mc;

                  }

        f_corr_e_1 = boost::shared_ptr<TFile>( new TFile(s_corr_e_1)); 
        f_corr_e_2 = boost::shared_ptr<TFile>( new TFile(s_corr_e_2)); 
        f_corr_e_3 = boost::shared_ptr<TFile>( new TFile(s_corr_e_3)); 
        f_corr_mu = boost::shared_ptr<TFile>( new TFile(s_corr_mu));
        
        el_corr_1 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_1->Get("2e")->Clone() )) );
        el_corr_2 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_2->Get("2e")->Clone() )) );
        el_corr_3 = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_e_3->Get("2e")->Clone() )) );

        mu_corr = boost::shared_ptr<TH2F>(  (static_cast<TH2F*>(f_corr_mu->Get("2mu")->Clone() )) );
        
        x_elpTaxis_1 = el_corr_1->GetXaxis(); y_eletaaxis_1 = el_corr_1->GetYaxis();
        maxPtEl_1 = x_elpTaxis_1->GetXmax(); minPtEl_1 = x_elpTaxis_1->GetXmin();

        x_elpTaxis_2 = el_corr_2->GetXaxis(); y_eletaaxis_2 = el_corr_2->GetYaxis();
        maxPtEl_2 = x_elpTaxis_2->GetXmax(); minPtEl_2 = x_elpTaxis_2->GetXmin();

        x_elpTaxis_3 = el_corr_3->GetXaxis(); y_eletaaxis_3 = el_corr_3->GetYaxis();
        maxPtEl_3 = x_elpTaxis_3->GetXmax(); minPtEl_3 = x_elpTaxis_3->GetXmin();

        x_mupTaxis = mu_corr->GetXaxis(); y_muetaaxis = mu_corr->GetYaxis();
        maxPtMu = x_mupTaxis->GetXmax(); minPtMu = x_mupTaxis->GetXmin();

}

// HelperFunction::HelperFunction(const HelperFunction& rhs)
// {
//    // do actual copying here;
// }

HelperFunction::~HelperFunction()
{
}

//
// assignment operators
//
// const HelperFunction& HelperFunction::operator=(const HelperFunction& rhs)
// {
//   //An exception safe implementation is
//   HelperFunction temp(rhs);
//   swap(rhs);
//
//   return *this;
// }

//
// member functions
//

double HelperFunction:: masserrorFullCov(std::vector<TLorentzVector> p4s, TMatrixDSym covMatrix){

        int ndim = 3*p4s.size();
        if(debug_) cout<<""<<endl;

        TMatrixD jacobian(1,ndim);

        double e = 0; double mass = 0;
        double px = 0; double py = 0; double pz = 0;
        for (unsigned int ip = 0; ip < p4s.size(); ip++) {
         
            e = e + p4s[ip].E();
            px = px + p4s[ip].Px();
            py = py + p4s[ip].Py();
            pz = pz + p4s[ip].Pz();
        }

        mass = TMath::Sqrt(e*e-px*px-py*py-pz*pz);

        for (unsigned int i = 0, o = 0; i < p4s.size(); i++, o += 3) {

                double pxi = p4s[i].Px();
                double pyi = p4s[i].Py();
                double pzi = p4s[i].Pz();
                double ei = p4s[i].E();

                jacobian(0, o+0) = (e*(pxi/ei) - px)/mass;
                jacobian(0, o+1) = (e*(pyi/ei) - py)/mass;
                jacobian(0, o+2) = (e*(pzi/ei) - pz)/mass;
        }

        TMatrixDSym massCov = covMatrix.Similarity(jacobian);

        double dm2 = massCov(0,0);
        return (dm2 > 0 ? std::sqrt(dm2) : 0.0);

}


double HelperFunction::masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr){
        // if(Lep.size()!= pterr.size()!=4) {std::cout<<" Lepsize="<<Lep.size()<<", "<<pterr.size()<<std::endl;}
        TLorentzVector compositeParticle ;
        for(unsigned int i=0; i<Lep.size(); i++){
                compositeParticle+=Lep[i];
        }
        double mass  =  compositeParticle.M();

        double masserr = 0;

        for(unsigned int i=0; i<Lep.size(); i++){
                TLorentzVector variedLep; // = Lep[i];

                variedLep.SetPtEtaPhiM(Lep[i].Pt()+ pterr[i], Lep[i].Eta(), Lep[i].Phi(), Lep[i].M());
                TLorentzVector compositeParticleVariation ;
                for(unsigned int j=0; j<Lep.size(); j++){
                        if(i!=j)compositeParticleVariation+=Lep[j];
                        else compositeParticleVariation+=variedLep;
                }

                masserr += (compositeParticleVariation.M()-mass)*(compositeParticleVariation.M()-mass);
        }

        return sqrt(masserr);
}


double HelperFunction::pterr( reco::Candidate *c, bool isData){

  reco::GsfElectron *gsf; reco::Muon *mu;
  reco::PFCandidate *pf;

  pat::Muon *patmu;

  double pterrLep = 0.0;

  if ((gsf = dynamic_cast<reco::GsfElectron *> (&(*c)) ) != 0)
  {

    pterrLep=pterr(gsf, isData);

    double pT_e = gsf->pt();
    double eta_e = gsf->eta();
    
    if (gsf->ecalDriven()) {

        if (fabs(eta_e) < 1) {
            if (pterrLep/pT_e < 0.03) { // hardcode 1
                int xbin = x_elpTaxis_1->FindBin(pT_e);
                int ybin = y_eletaaxis_1->FindBin(fabs(eta_e));
                if(pT_e > minPtEl_1 && pT_e < maxPtEl_1 ){
                    pterrLep*=el_corr_1->GetBinContent(xbin,ybin);
                } else {
                    pterrLep*=1.0;
                }
            } else {
                   if (isData){pterrLep*=1.187;}
                      else{pterrLep*=1.224;}
                   } // hardcode 2        
        } else if (fabs(eta_e) > 1 && fabs(eta_e) < 2.5) {
            if (pterrLep/pT_e < 0.07) { // hardcode 3
                int xbin = x_elpTaxis_2->FindBin(pT_e);
                int ybin = y_eletaaxis_2->FindBin(fabs(eta_e));
                if(pT_e > minPtEl_2 && pT_e < maxPtEl_2 ){
                    pterrLep*=el_corr_2->GetBinContent(xbin,ybin);
                } else {
                    pterrLep*=1.0;
                }
            } else {
                   if (isData){pterrLep*=0.815;}
                      else{pterrLep*=0.786;} // hardcode 4
                   }
        } // 1 < |eta| < 2.5      
    } else {
      
        int xbin = x_elpTaxis_3->FindBin(pT_e);
        int ybin = y_eletaaxis_3->FindBin(fabs(eta_e));
        if(pT_e > minPtEl_3 && pT_e < maxPtEl_3 ){
            pterrLep*=el_corr_3->GetBinContent(xbin,ybin);
        } else {
            pterrLep*=1.0;
        }
    }
    

  }
  else if ((mu = dynamic_cast<reco::Muon *> (&(*c)) ) != 0)
  {
    pterrLep=pterr(mu, isData);
    if(debug_)cout<<"reco pt err is "<<pterrLep<<endl;

    if( (patmu = dynamic_cast<pat::Muon *> (&(*c)) )!=0){

     if ( patmu->hasUserFloat("correctedPtError") == true ) {
       if(debug_) cout<<"use userFloat for muon pt err"<<endl;
       pterrLep = patmu->userFloat("correctedPtError");
       if(debug_) cout<<"calib pt err is "<<pterrLep<<endl;
     }
 
    }

    int xbin = x_mupTaxis->FindBin(mu->pt());
    int ybin = y_muetaaxis->FindBin(fabs(mu->eta()));
    if(mu->pt()>minPtMu && mu->pt()<maxPtMu ){
        pterrLep*=mu_corr->GetBinContent(xbin,ybin);
    } else {
        pterrLep*=1.0;
    }

  }
  else if ((pf = dynamic_cast<reco::PFCandidate *> (&(*c)) ) != 0)
  { 
    pterrLep=pterr(c, isData);
  }

  return pterrLep;

}

double HelperFunction::pterr( reco::Muon* mu, bool isData){

        double pterr = mu->muonBestTrack()->ptError();

        return pterr;
}

double HelperFunction::pterr( reco::GsfElectron * elec, bool isData ){

        if(debug_) cout<<"reco:gsfelectron pt err"<<endl; 

        double perr = elec->p();

        if (elec->ecalDriven()){
           perr = elec->p4Error(reco::GsfElectron::P4_COMBINATION);         
        }
        else{

                 double ecalEnergy = elec->correctedEcalEnergy() ;

                 if(debug_)cout<<"ecalEnergy "<<ecalEnergy<<endl;
                 double err2 = 0.0;
                 if (elec->isEB()) {
                        err2 += (5.24e-02*5.24e-02)/ecalEnergy;
                        err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
                        err2 += 1.00e-02*1.00e-02;
                 } else if (elec->isEE()) {
                        err2 += (1.46e-01*1.46e-01)/ecalEnergy;
                        err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
                        err2 += 1.94e-03*1.94e-03;
                 }
                 perr = ecalEnergy * sqrt(err2);

        }
         
        double pterr = perr*elec->pt()/elec->p();

        return pterr;

}


double HelperFunction::pterr(TLorentzVector ph){

         if(debug_) cout<<"perr for pf photon"<<endl;

         double perr = PFEnergyResolution().getEnergyResolutionEm(ph.E(), ph.Eta());

         double pterr = perr*ph.Pt()/ph.P();

         return pterr;
}

//
// const member functions
//

//
// static member functions
//
#endif
