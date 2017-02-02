/** \file FissionStudyProcessor.cpp
 * \brief A class to process data from Fission study experiments using
 * VANDLE.
 *
 *\author W. A. Peters
 *\date Oct 3, 2016
 */
#include <fstream>
#include <iostream>

#include <cmath>

#include "BarBuilder.hpp"
#include "BarDetector.hpp"
#include "DammPlotIds.hpp"
#include "DoubleBetaProcessor.hpp"
#include "DetectorDriver.hpp"
//#include "GeProcessor.hpp"
#include "GetArguments.hpp"
#include "Globals.hpp"
#include "FissionStudyProcessor.hpp"
#include "RawEvent.hpp"
#include "TimingMapBuilder.hpp"
#include "VandleProcessor.hpp"

#ifdef useroot
//static double tof_;
//static double qdc_;
#endif

namespace dammIds {
    namespace experiment {
        const int DD_FS0  = 0; //!<Si_E vs TOF
        const int DD_FS1  = 1; //!<QDC ToF Ungated
        const int DD_FS2  = 2; //!<Cor ToF vs. Gamma E
        const int D_FS3  = 3; //!<Vandle Multiplicity
        const int DD_FS4  = 4; //!<QDC vs Cor Tof Mult1
        const int D_FS5  = 5; //!<Mult2 Sym Plot Tof 
        const int DD_FS6  = 6; //!<Si-Gated Bar vs. TOF
        const int D_FS7  = 7; //!<TOF mass lo
        const int D_FS8  = 8; //!<TOF mass hi
        const int D_FS9  = 9; //!<TOF mass lo1
        const int D_FS10  = 10; //!<TOF mass lo2
        const int DD_FS11 = 11; //!< mass v TOF forward
        const int DD_FS13 = 13; //!< mass v TOF all
        const int D_ENERGY = 14; //!< Si singles ungated
        const int D_HPGE_E = 15; //Si:hpge energy
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

void FissionStudyProcessor::DeclarePlots(void) {
       DeclareHistogram2D(DD_FS0, SC, SC, "Si_raw/4 vs TOF");
       DeclareHistogram2D(DD_FS1, SC, SC, "QDC/4 ToF Ungated");
    DeclareHistogram2D(DD_FS2, SC, SC, "Beta vs mass");
    DeclareHistogram1D(D_FS3, SC, "Si_E *1000");
    DeclareHistogram2D(DD_FS4, SC, SC, "Si_E vs mass");
    DeclareHistogram1D(D_FS5, SC, "Beta frag *1000");
    DeclareHistogram2D(DD_FS6, SC, 256, "Si-Gated Bar vs. TOF");
    DeclareHistogram1D(D_FS7, SC, "CorTOF lo mass");
    DeclareHistogram1D(D_FS8, SC, "CorTOF hi mass");
    DeclareHistogram1D(D_FS9, SC, "CorTOF lo1 mass");
    DeclareHistogram1D(D_FS10, SC, "CorTOF lo2 mass");
    DeclareHistogram2D(DD_FS11, SC, SA, "mass vs. forward CorTOF");
    DeclareHistogram2D(DD_FS13, SC, SA, "mass vs. CorTOF");
    DeclareHistogram1D(D_ENERGY, SD, "Si_raw");
    DeclareHistogram1D(D_HPGE_E, SD, "HPGe_raw / 2");
}

FissionStudyProcessor::FissionStudyProcessor() : EventProcessor(OFFSET, RANGE, "FissionStudyProcessor") {
    associatedTypes.insert("vandle");
    associatedTypes.insert("si");
    associatedTypes.insert("beta");

    char hisFileName[32];
    GetArgument(1, hisFileName, 32);
    string temp = hisFileName;
    temp = temp.substr(0, temp.find_first_of(" "));
    stringstream name;
    name << temp << ".dat";
    outstream = new ofstream(name.str().c_str());
#ifdef useroot
    stringstream rootname;
     rootname << temp << ".root";
    rootfile_ = new TFile(rootname.str().c_str(),"RECREATE");
///    roottree_ = new TTree("vandle","");
///    roottree_->Branch("tof",&tof_,"tof/D");
///    roottree_->Branch("qdc",&qdc_,"qdc/D");
///    qdctof_ = new TH2D("qdctof","",1000,-100,900,16000,0,16000);
///    vsize_ = new TH1D("vsize","",40,0,40);
#endif
}

FissionStudyProcessor::~FissionStudyProcessor() {
    outstream->close();
    delete(outstream);
#ifdef useroot
    rootfile_->Write();
    rootfile_->Close();
    delete(rootfile_);
#endif
}

///We do nothing here since we're completely dependent on the resutls of others
bool FissionStudyProcessor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
        return(false);
    return(true);
}

bool FissionStudyProcessor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return(false);
    double plotMult_ = 2;
    double plotOffset_ = 1000;

    BarMap vbars, betas;
    map<unsigned int, pair<double,double> > lrtBetas;
    vector<ChanEvent*> geEvts;

    if(event.GetSummary("vandle")->GetList().size() != 0)
        vbars = ((VandleProcessor*)DetectorDriver::get()->
            GetProcessor("VandleProcessor"))->GetBars();
    if(event.GetSummary("beta:double")->GetList().size() != 0) {
        betas = ((DoubleBetaProcessor*)DetectorDriver::get()->
            GetProcessor("DoubleBetaProcessor"))->GetBars();
        lrtBetas = ((DoubleBetaProcessor*)DetectorDriver::get()->
            GetProcessor("DoubleBetaProcessor"))->GetLowResBars();
       }
    
// Get si events here 
//But no Si processor
// need raw si trap energy from module 0 channel 4 and module 5 ch 0   
    static const vector<ChanEvent*> &SiPinEvts =
	event.GetSummary("si:pin")->GetList();
 
    static const vector<ChanEvent*> &SiHpgeEvts =
    event.GetSummary("si:hpge")->GetList();
    
//    static const vector<ChanEvent*> &labr3Evts =
//	event.GetSummary("labr3:mrbig")->GetList();

#ifdef useroot
///    vsize_->Fill(vbars.size());
#endif

    //Obtain some useful logic statuses
//    double lastProtonTime =
//	TreeCorrelator::get()->place("logic_t1_0")->last().time;
//    bool isTapeMoving = TreeCorrelator::get()->place("TapeMove")->status();

//    int bananaNum = 2;
//    bool hasMultOne = vbars.size() == 1;
//    bool hasMultTwo = vbars.size() == 2;
    //bool isFirst = true;
//MVT settings:
//	   double E_offset = 1200.0; // arb units 700 1200
//	   double E_scale = (1.0/20000.0); //abr units 1/20000
//	   double V_offset = 0.09; // or 0.0 with 700 above 0.09 with 1200

////// Variables-------------------------------    
	   double E_offset = 1500.0; // arb units 700 1200
	   double E_scale = (1.0/30000.0); //abr units 1/20000
	   double V_offset = 0.0; // or 0.0 with 700 above 0.09 with 1200
	   double si_energy = 0.0;
	   double hpge_energy = 0.0;
	   double tof = 0.0;
	   double CorTof =0.0;
	   double frag_mass = 0.0;
	   double Si_E = 0.0 ;
	   double v_c = 29.98; // V of light cm/ns
	   double v_frag = 0.0;
	   double beta_frag = 0.0;
	    
	//-------------- SiPin Processing ---------------
        for(vector<ChanEvent*>::const_iterator it = SiPinEvts.begin();
	    it != SiPinEvts.end(); it++) {
	    // need to identify which channel it was
	       si_energy = (*it)->GetEnergy();
	       plot(D_ENERGY, si_energy);
	    }
	    
	    for(vector<ChanEvent*>::const_iterator it = SiHpgeEvts.begin();
	    it != SiHpgeEvts.end(); it++) {
            hpge_energy = (*it)->GetEnergy();
            plot(D_HPGE_E, (hpge_energy/2.0));
	    }
	    
///    //Begin processing for VANDLE bars
   for (BarMap::iterator it = vbars.begin(); it !=  vbars.end(); it++) {
        TimingDefs::TimingIdentifier barId = (*it).first;
        BarDetector bar = (*it).second;

        if(!bar.GetHasEvent())
            continue;

        unsigned int barLoc = barId.first;
        TimingCalibration cal = bar.GetCalibration();

        for(BarMap::iterator itStart = betas.begin();
	    itStart != betas.end(); itStart++) {
	    unsigned int startLoc = (*itStart).first.first;
            BarDetector start = (*itStart).second;
            if(!start.GetHasEvent())
                continue;
///         unsigned int barPlusStartLoc = barLoc + startLoc;
                            
            v_frag = start.GetBetaVel();
            beta_frag = (v_frag + V_offset)/v_c;

            tof = bar.GetTimeAverage() - start.GetTimeZero() + 
                    cal.GetTofOffset(startLoc);

            CorTof =
                ((VandleProcessor*)DetectorDriver::get()->
        		 GetProcessor("VandleProcessor"))->
        		CorrectTOF(tof, bar.GetFlightPath(), cal.GetZ0());
                    
            Si_E = ((si_energy + E_offset) * E_scale);       
                           
            frag_mass = ( (2.0*Si_E) / (beta_frag * beta_frag));        
            	    
	            plot(D_FS3, Si_E * 1000.0);
              	    
	            plot(D_FS5, beta_frag *1000.0);
      
            if( (si_energy > 651) && (v_frag > 0.1) ){                
 
	            plot(DD_FS0, CorTof*plotMult_ + plotOffset_, (si_energy/4.0));
	    
	            plot(DD_FS1, CorTof*plotMult_ + plotOffset_, (bar.GetQdc()/4.0));
	    
	            plot(DD_FS2, frag_mass, beta_frag*1000.0);
	            
	            plot(DD_FS4, frag_mass, Si_E*1000.0);
	            
                plot(DD_FS6, tof*plotMult_+plotOffset_ , barLoc);
                
	            if(frag_mass > 280 && frag_mass < 371) plot(D_FS7, 
	                CorTof*plotMult_ + plotOffset_);
	            if(frag_mass > 371 && frag_mass < 483) plot(D_FS8, 
	                CorTof*plotMult_ + plotOffset_);	            
	            if(frag_mass > 280 && frag_mass < 326) plot(D_FS9, 
	                CorTof*plotMult_ + plotOffset_);
	            if(frag_mass > 326 && frag_mass < 371) plot(D_FS10, 
	                CorTof*plotMult_ + plotOffset_);
 
	            plot(DD_FS13, CorTof*plotMult_ + plotOffset_, frag_mass);
	            
	            if(barLoc >= 21) {
	               plot(DD_FS11, CorTof*plotMult_ + plotOffset_, frag_mass);
	            
	           
	            }// forward bars
        	}// si_energy and velocity cuts
	    }
	} //End processing for VANDLE bars
	
	


//	    *outstream << tof << " " << bar.GetQdc() << endl;
#ifdef useroot
//        qdctof_->Fill(tof,bar.GetQdc());
//        qdc_ = bar.GetQdc();
///        tof_ = tof;
///        roottree_->Fill();
//        qdc_ = tof_ = -9999;
#endif


	    ///Starting to look for 2n coincidences in VANDLE
//	    BarMap::iterator itTemp = it;
//	    itTemp++;
//	    for (BarMap::iterator it2 = itTemp; it2 !=  vbars.end(); it2++) {
//		TimingDefs::TimingIdentifier barId2 = (*it2).first;
//		BarDetector bar2 = (*it2).second;
//		if(!bar2.GetHasEvent())
//		    continue;
//		unsigned int barLoc2 = barId2.first;
//		bool isAdjacent = abs((int)barLoc2 - (int)barLoc) < 1;
//		TimingCalibration cal2 = bar2.GetCalibration();
//		double tofOffset2 = cal2.GetTofOffset(startLoc);
//		double tof2 = bar2.GetCorTimeAve() -
//		    start.GetCorTimeAve() + tofOffset2;
//		double corTof2 =
//		    ((VandleProcessor*)DetectorDriver::get()->
//		     GetProcessor("VandleProcessor"))->
//		    CorrectTOF(tof2, bar2.GetFlightPath(), cal2.GetZ0());
//		bool inPeel2 = histo.BananaTest(bananaNum,
//						corTof2*plotMult_+plotOffset_,
//						bar2.GetQdc());
//	    }//End 2n coincidence routine
	    
    //-------------- LaBr3 Processing ---------------
//    for(vector<ChanEvent*>::const_iterator it = labr3Evts.begin();
//	it != labr3Evts.end(); it++)
//	plot(DD_FS6, (*it)->GetEnergy());


    //----------------- GE Processing -------------------
//    bool hasBeta = TreeCorrelator::get()->place("Beta")->status();
//    double clockInSeconds = Globals::get()->clockInSeconds();
    // plot with 10 ms bins
//    const double plotResolution = 10e-3 / clockInSeconds;
//    for (vector<ChanEvent*>::iterator it1 = geEvts.begin();
//	 it1 != geEvts.end(); ++it1) {
//       ChanEvent *chan = *it1;



    EndProcess();
    return(true);
 }
