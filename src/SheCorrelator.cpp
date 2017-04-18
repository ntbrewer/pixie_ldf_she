/*! \file SheCorrelator.cpp
 *
 * The SheCorrelator is designed to recontruct dssd events in a SHE experiment
 * and to correlate chains of alphas in dssd pixels
 */

#include <ctime>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <signal.h>

#include "DammPlotIds.hpp"
#include "Dssd4SHEProcessor.hpp"
#include "SheCorrelator.hpp"
#include "DetectorDriver.hpp"
#include "Exceptions.hpp"
#include "Notebook.hpp"


using namespace std;

SheEvent::SheEvent() {
    energy_ = -1.0;
    ratio_ = -1.0;
    time_= -1.0;
    mwpc_= -1;
    mwpcTime_= -1;
    mwpcEnergy_ = -1;
    has_beam_= false;
    has_veto_= false;
    has_escape_ = 0.0;
    type_= unknown;
}

SheEvent::SheEvent(double energy, double ratio, double time, int mwpc, double mwpcTime, double mwpcEnergy,
                   bool has_beam, bool has_veto, double has_escape, int pixel[3], 
 		   SheEventType type /* = unknown*/)
{
    set_energy(energy);
    set_ratio(ratio);
    set_time(time);
    set_mwpc(mwpc);
    set_mwpcTime(mwpcTime);
    set_mwpcEnergy(mwpcEnergy);
    set_beam(has_beam);
    set_veto(has_veto);
    set_escape(has_escape);
    set_Xpixel(pixel);
    set_Ypixel(pixel);
    set_Epixel(pixel);
    set_type(type);
}


SheCorrelator::SheCorrelator(int size_x, int size_y) 
{
    size_x_ = size_x + 1;
    size_y_ = size_y + 1;
    pixels_ = new deque<SheEvent>*[size_x_];
    for(int i = 0; i < size_x_; ++i)
        pixels_[i] = new deque<SheEvent>[size_y_];

}


SheCorrelator::~SheCorrelator() {
    for(int i = 0; i < size_y_; ++i) {
        delete[] pixels_[i];
    }
    delete[] pixels_;
}


bool SheCorrelator::add_event(SheEvent& event, int x, int y, Plots& histo) {

    if (x < 0 || x >= size_x_) {
        stringstream ss;
        ss << "Requested event at non-existing X strip " << x << endl;
        throw GeneralWarning(ss.str());
    }
    if (y < 0 || y >= size_y_) {
        stringstream ss;
        ss << "Requested event at non-existing Y strip " << y << endl;
        throw GeneralWarning(ss.str());
    }
    
    if (event.get_type() == heavyIon) //test for looking at R - (R-False) - F
//      if (event.get_type() == unknown)
          flush_chain(x, y, histo);

    pixels_[x][y].push_back(event);

    if (event.get_type() == fission)
        flush_chain(x, y, histo);
    
    return true;
}


bool SheCorrelator::flush_chain(int x, int y, Plots& histo){

    unsigned chain_size = pixels_[x][y].size();

    /** If chain too short just clear it */
    if (chain_size < 2) {
        pixels_[x][y].clear();
        return false;
    }

    SheEvent first = pixels_[x][y].front();

    /** Conditions for interesing chain:
     *      * starts with heavy ion implantation
     *      * has ion + x number of decays (0 <= x) + fission
     *      * or includes at least two alphas
     */

    /** If it doesn't start with heavyIon, clear and exit**/
    if (first.get_type() != heavyIon) {
        pixels_[x][y].clear();
        return false;
    } 
    /* If it is greater than 1 element long, check if the last is fission,
     *  if not - clear and exit**/
    if (chain_size <= 2 && pixels_[x][y].back().get_type() != fission) {
        pixels_[x][y].clear();
        return false;
    }
    int alphas = 0;
    double alphaE[8]={0};
    double alphaTime[8]={0};
    double fissionE =0, fissionTime =0;
    double VRecoilTime=0, VRecoilE=0;
    static double mwpcTime, mwpcEnergy;//mwpcTime -> TOF
    int interest=0;
    static int ctra=0 , ctrsf=0, ctro=0;
    int l=0; //ittr;
    stringstream ss;

    time_t wallTime = DetectorDriver::get()->GetWallTime(first.get_time());
    string humanTime = ctime(&wallTime);
    humanTime.erase(humanTime.find('\n', 0), 1);
    ss << humanTime << "\t X = " << x <<  " Y = " << y << " num " << ctra << endl; //" time: " << wallTime <<  endl;


    for (deque<SheEvent>::iterator it = pixels_[x][y].begin();
         it != pixels_[x][y].end();
         ++it)
    { 	// Process Correlations in the SHE Event. Make a logic for a good event and pass it to the plot routines. 
	int type = (*it).get_type();
	int pixel[2] = {(*it).get_Xpixel(),(*it).get_Ypixel()}; 
	double energy = (*it).get_energy();
	double time = (*it).get_time()*(Globals::get()->clockInSeconds());//units in s *1e-3 gives units in ms
	
	if( type == heavyIon) {
	    mwpcTime = (*it).get_mwpcTime();
	    mwpcEnergy=(*it).get_mwpcEnergy();
   	    VRecoilE=energy;
	    VRecoilTime=time;
            
	}    	
	if(type == fission) {
	    fissionE = energy;
	    fissionTime = time;
	}
	if (type == alpha) {
            alphas += 1;
	    if(alphas <= 8) {
		alphaE[alphas]=energy;
		alphaTime[alphas] = time; 
	    }
        }

	
	if (type != unknown && type !=check && type != lightIon) {
	   human_event_info((*it), ss, first.get_time());
	   if (pixel[0] >0) {
	      ss << " in X: " << int(pixel[0]/1000) << " and " << int(pixel[0] - round(pixel[0]/1000)*1000);
	   } else if (pixel[1]>0 ) {
	      ss << " in Y: " << int(pixel[1]/1000) << " and " << int(pixel[1] - round(pixel[1]/1000)*1000);
	   }
	    if (type == heavyIon ) {
	      if ( mwpcEnergy < 500) {
    	        ss << " MWE: " <<  mwpcEnergy << " " ;
	      }    
	      if (abs(mwpcTime - 2000)>10 ) {
		ss << " TOF: " <<  mwpcTime << " " ;
	      }
	    }
           ss << endl;
	}
	
	if (alphas >= 2 && alphaE[1] > 9000 && l==0 && (alphaTime[1]-VRecoilTime) > 0 && (alphaTime[1]-VRecoilTime) < 0.1 && fissionE==0){//In units of s    
	   interest=1;
	   l++;
	} else if (fissionE>0 && (fissionTime - VRecoilTime) <270000 && alphas >1 && (alphaTime[1]-VRecoilTime) < 0.5 ) {
	   interest=2;
	} else if (fissionE > 0 && alphas==0 && (fissionTime - VRecoilTime)<1)  { 
	   interest=3;
	}  else if (alphas > 0)  { 
	   interest=4;
	}
	//tbd process out the unknown events in correlation with the good chains.
	
    }

    // Plot and report characteristics of interesting events. 
    if (interest>0 && interest <3) {
	Notebook::get()->report(ss.str());
    }
    if (interest==3 && fissionE > 50000 ) {
        Notebook::get()->report(ss.str());
    }
    if (interest==1){
        //cout << "INTERESTING alpha!!" << endl;
	/* Event vs. Number for Recoil */
        /*ittr =0;
        while(ittr < 2 ) { 
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_ALPHA,VRecoilE/10,ctra);
	    ittr++; 
	}*/
	histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_TIME,10,VRecoilE/10);
	/*Event in chain plots*/
	for (int i=1; i!=alphas; ++i) {//((alphaE[i]>0) & ((alphaTime[i]-VRecoilTime)<=100)){
	 
	    
     	    histo.Plot(dammIds::dssd4she::DD_CHAIN_ALPHA_V_ALPHA, alphaE[1]/10 , alphaE[i+1]/10);
 	    histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_TIME,(log10(alphaTime[i]-VRecoilTime)+10)*100,alphaE[i]/10);
	   

	   /*Event in chain vs. Number plots*/
	   /*ittr =0;
	   while (ittr <= (alphaTime[i]-VRecoilTime)*1e6 && ittr <8000) {
		histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_ALPHA,alphaE[i]/10,ctra);
		ittr++;
	   }*/
	}

	/*MWPC Plots*/
        histo.Plot(dammIds::dssd4she::DD_TOF_A_EVENT,mwpcTime,ctra);
	histo.Plot(dammIds::dssd4she::DD_MWPC_ENERGY_A_EVENT,mwpcEnergy,ctra);

	/*For Calibration runs*/	
	histo.Plot(dammIds::dssd4she::DD_FRONT_V_ALPHA1,alphaE[1],x); //For gainmatching from traces
	histo.Plot(dammIds::dssd4she::DD_BACK_V_ALPHA1,alphaE[1],y);
	
	ctra++;  
    }

    if (interest==2 ) {
	//cout << "Fission after alpha chain!" << endl;
	/* Event vs. Number plots */
	/*Recoil*/
	histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_TIME,10,VRecoilE/10);
	/*ittr=0;
        while(ittr < 2 ) { 
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_ALPHA,VRecoilE/10,ctra);
	    ittr++; 
	}*/
	/*Fission*/
	/*ittr=0;
	while (ittr < (fissionTime - VRecoilTime) ) {
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_ALPHA,fissionE/100+1500,ctra);
	    ittr++;
	}*/
	histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_TIME,(log10(fissionTime-VRecoilTime)+10)*100,fissionE/100+1500);
	/*Event in chain plots*/
	for (int j=1; j!=alphas; ++j) {
    	    histo.Plot(dammIds::dssd4she::DD_CHAIN_ALPHA_V_ALPHA, alphaE[1]/10 , alphaE[j+1]/10);
	    histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_TIME,(log10(alphaTime[j]-VRecoilTime)+10)*100,alphaE[j]/10);
	   /*Event in chain vs. Number plots*/
	   /* ittr=0;
	    while (ittr <= (alphaTime[i]-VRecoilTime)*1e6 && ittr <8000) {
		histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_ALPHA,alphaE[i]/10,ctra);
		ittr++;
	    }*/

	}

	/*MWPC Plots*/
	histo.Plot(dammIds::dssd4she::DD_TOF_A_EVENT,mwpcTime,ctra);
	histo.Plot(dammIds::dssd4she::DD_MWPC_ENERGY_A_EVENT,mwpcEnergy,ctra);
	
        ctra++;
    }

    if (interest==3)  {
	/* Event vs. Number plots */
	/*Recoil*/ 
	/*ittr =0;
        while(ittr < 2 ) { 
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_FISSION,VRecoilE/10,ctr);
	    ittr++; 
	}*/
	/*Fission*/
	/*ittr =0;
	while (ittr < (fissionTime - VRecoilTime)*1e6 ) {
	    histo.Plot(dammIds::dssd4she::DD_CHAIN_NUM_FISSION,fissionE/200,ctr);
	    ittr++;
	}*/
	//plot SF context plots here. 
        // I-SF dts for <1s in .1 ms (10k) and <10ms in 1 us units (10k).
        //SF E in E/100 VRecoilE in E/10
        histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_DTIME_1S,fissionE/100.,(fissionTime-VRecoilTime)*1.0e4);
        histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_DTIME_10MS,fissionE/100.,(fissionTime-VRecoilTime)*1.0e7);
	histo.Plot(dammIds::dssd4she::DD_RECOIL_ENERGY_V_SFE,fissionE/100.,VRecoilE/10.);        
        //cout << fissionTime-VRecoilTime << endl;
//	histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_TIME,10,VRecoilE/10);
        //cout << mwpcTime <<endl;
	/*Event in sf plot*/
//	histo.Plot(dammIds::dssd4she::DD_CHAINS_ENERGY_V_TIME,(log10(fissionTime-VRecoilTime)+10)*100,fissionE/100+1500);
	/*MWPC Plots*/
//        histo.Plot(dammIds::dssd4she::DD_TOF_SF_EVENT,mwpcTime,ctrsf);       
//	histo.Plot(dammIds::dssd4she::DD_MWPC_ENERGY_SF_EVENT,mwpcEnergy,ctrsf);
//	histo.Plot(dammIds::dssd4she::DD_MWPC_ENERGY_V_SFE,fissionE/100,mwpcEnergy/10);
    //cout << "SFE " << fissionE/100 << " MWPCE " << mwpcEnergy << endl;
	//ctrsf++;
    }

    if (interest==4 ) {
	//cout << "other alpha chains" << endl;
	/*Fission*/
	histo.Plot(dammIds::dssd4she::DD_OTHER_ENERGY_V_TIME,10,VRecoilE/10);
	if (fissionE>0) {
	    histo.Plot(dammIds::dssd4she::DD_OTHER_ENERGY_V_TIME,(log10(fissionTime-VRecoilTime)+10)*100,fissionE/100+1500);
	}
	/*Event in chain plots*/
	for (int j=1; j!=alphas; ++j) {
    	    histo.Plot(dammIds::dssd4she::DD_OTHER_ALPHA_V_ALPHA, alphaE[1]/10 , alphaE[j+1]/10);
	    histo.Plot(dammIds::dssd4she::DD_OTHER_ENERGY_V_TIME,(
		log10(alphaTime[j]-VRecoilTime)+10)*100,alphaE[j]/10);
	}

	/*MWPC Plots*/
	histo.Plot(dammIds::dssd4she::DD_TOF_O_EVENT,mwpcTime,ctro);
	histo.Plot(dammIds::dssd4she::DD_MWPC_ENERGY_O_EVENT,mwpcEnergy,ctro);
	ctro++;
    }

    pixels_[x][y].clear();
    
    return true;
}

// Save event to file
void SheCorrelator::human_event_info(SheEvent& event, stringstream& ss,
                                     double clockStart) {
    string humanType;
    switch (event.get_type()) {
        case alpha:   humanType = "A";
                        break;
        case heavyIon: humanType = "I";
                        break;
        case fission: humanType = "F";
                        break;
	case check: humanType = "C";
                        break;
        case lightIon: humanType = "L";
                        break;
        case unknown: humanType = "U";
                        break;
    }


    ss << fixed 
       << humanType 
       << " " 
       << setprecision(0) << setw(12) << event.get_energy();

	if (event.get_ratio()==0.0) {
	     ss  << " " << setprecision(3) << setw(13) ;
	} else {
	     ss  << " " << setprecision(2) << event.get_ratio() << " " << setprecision(3) << setw(8);
	}

    ss << (event.get_time() - clockStart) * 1.0e-5  //ms
       << " M" << event.get_mwpc() 
       << "B" << event.get_beam() 
       << "V" << event.get_veto();
 	
	if (event.get_escape()==0.0) {
	     ss  << "E0";
	} else {
	     ss << "E"<< event.get_Epixel() << " " << event.get_escape();
	}
}
