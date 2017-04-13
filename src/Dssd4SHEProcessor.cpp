/*! \file DssdProcessor.cpp
 *
 * The DSSD processor handles detectors of type dssd_front and dssd_back and
 *   determines whether the events are implants or decays and informs the
 *   correlator accordingly
 */

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <signal.h>

#include "Dssd4SHEProcessor.hpp"
#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "Messenger.hpp"
#include "Notebook.hpp"
#include "RawEvent.hpp"

#include "CfdAnalyzer.hpp"
#include "DoubleTraceAnalyzer.hpp"
#include "FittingAnalyzer.hpp"
#include "TauAnalyzer.hpp"
#include "TraceAnalyzer.hpp"
#include "TraceExtracter.hpp"
#include "WaveformAnalyzer.hpp"
#include "DetectorDriver.hpp"
#include "Notebook.hpp"
using namespace dammIds::dssd4she;
using namespace std;

Dssd4SHEProcessor::Dssd4SHEProcessor(double timeWindowFB,
				     double timeWindowDS,
                                     double deltaEnergy,
				     double recoilEnergyCut,
                                     double highEnergyCut,
                                     double lowEnergyCut,
                                     double fissionEnergyCut,
                      		     int XToffset1,
		                     int YToffset1,
		                     int YToffset2,
		                     double TempCheck,
                                     int numBackStrips,
                                     int numFrontStrips) :
    EventProcessor(OFFSET, RANGE, "dssd4she"),

    correlator_(numBackStrips, numFrontStrips)
{
    timeWindowFB_ = timeWindowFB;
    timeWindowDS_ = timeWindowDS;
    deltaEnergy_ = deltaEnergy;
    recoilEnergyCut_ = recoilEnergyCut;
    highEnergyCut_ = highEnergyCut;
    lowEnergyCut_ = lowEnergyCut;
    fissionEnergyCut_ = fissionEnergyCut;
    XToffset1_ = XToffset1;
    YToffset1_ = YToffset1;
    YToffset2_ = YToffset2;
    TempCheck_ = TempCheck;
    name = "dssd";
    associatedTypes.insert("dssd_front");
    associatedTypes.insert("dssd_back");

    stringstream ss;
    ss << fixed 
       << "#T " 
       << setw(12) << "E (keV) " 
       << "XY Ratio "
       << setw(12) << "t (ms) "  
       << " M  B  V  E" 
       << endl;
    Notebook::get()->report(ss.str());
}


void Dssd4SHEProcessor::DeclarePlots(void)
{
    
    const int energyBins = SE;
    const int energyBins2 = SB;
    const int xBins = S7;
    const int yBins = S6;
    const int timeBins = S8;
    const int traceBins = SA;
    unsigned short numTraces = Globals::get()->numTraces();

    
    DeclareHistogram2D(DD_TRACES_DECAY_CORR_X, traceBins, numTraces, "Traces correlated as a decay X");
    DeclareHistogram2D(DD_TRACES_DECAY_CORR_Y, traceBins, numTraces, "Traces correlated as a decay Y");
    DeclareHistogram2D(DD_TRACES_DECAY_CORR_E_V_NX, SD, numTraces, "Traces correlated as a decay X E vs N ");
    DeclareHistogram2D(DD_TRACES_DECAY_CORR_POS_V_NX, xBins, numTraces, "Traces correlated as a decay X POS vs N ");
    DeclareHistogram2D(DD_TRACES_DECAY_CORR_E_V_NY, SD, numTraces, "Traces correlated as a decay Y E vs N ");
    DeclareHistogram2D(DD_TRACES_DECAY_CORR_POS_V_NY, yBins, numTraces, "Traces correlated as a decay Y POS vs N ");

    /*DeclareHistogram1D(D_TRACES_DECAY_CORR_GATE1_X, numTraces, "Traces correlated in X Gate 1");
    DeclareHistogram1D(D_TRACES_DECAY_CORR_GATE2_X, numTraces, "Traces correlated in X Gate 2");
    DeclareHistogram1D(D_TRACES_DECAY_CORR_GATE3_X, numTraces, "Traces correlated in X Gate 3");
    DeclareHistogram1D(D_TRACES_DECAY_CORR_GATE1_Y, numTraces, "Traces correlated in Y Gate 1");
    DeclareHistogram1D(D_TRACES_DECAY_CORR_GATE2_Y, numTraces, "Traces correlated in Y Gate 2");
    DeclareHistogram1D(D_TRACES_DECAY_CORR_GATE3_Y, numTraces, "Traces correlated in Y Gate 3");
  */
    DeclareHistogram1D(D_ENERGY_X, energyBins, "Energy/10 dssd X strips");
    DeclareHistogram1D(D_ENERGY_Y, energyBins, "Energy/10 dssd Y strips");

    DeclareHistogram1D(D_DTIME, S8, "Pairs time diff in 10 ns (+ 1 bin)");

    DeclareHistogram1D(D_MWPC_MULTI, S5, "MWPC multiplicity");
    DeclareHistogram1D(D_ENERGY_CORRELATED_SIDE, energyBins, 
                       "Energy Side corr. with DSSD");
    DeclareHistogram2D(DD_ENERGY_CORRELATED_SIDE, energyBins, S3,
                       "Side Det. Num vs Energy Side corr. with DSSD");
    DeclareHistogram1D(D_DTIME_SIDE, S8, 
                        "Side det. time diff in 10 ns (+ 1 bin)");
    DeclareHistogram2D(DD_DTIME_SIDE, S7, S3, 
                        "Side det. Num vs time diff in 10 ns (+ 1 bin) ");

    DeclareHistogram2D(DD_ENERGY_DT__DSSD_MWPC, 
		       SB, S8, "DSSD energy/100 vs DT (10 ns) to MWPC");

    DeclareHistogram2D(DD_DE_E__DSSD_VETO, 
		       SB, SB, "DSSD energy/100 vs veto/100");

    /**Declare Saturation diagnostics **/
    DeclareHistogram1D(D_BACKE_W_FRONT_SAT, energyBins, "Back Energy/20 with Front Saturation");
    DeclareHistogram1D(D_FRONTE_W_BACK_SAT, energyBins, "Front Energy/20 with Back Saturation");

    DeclareHistogram1D(D_ENERGY_ALL,
		       energyBins2, "Event/10 All");
    DeclareHistogram1D(D_ENERGY_IMPLANT,
		       energyBins2, "Event/10 implant");
    DeclareHistogram1D(D_ENERGY_DECAY,
		       energyBins2, "Event/10 decay");
    DeclareHistogram1D(D_ENERGY_LIGHT,
		       energyBins2, "Event/10 light ion");
    DeclareHistogram1D(D_ENERGY_UNKNOWN,
		       energyBins2, "Event/10 unknown");
    DeclareHistogram1D(D_ENERGY_FISSION,
		       energyBins2, "Event/100 fission");
    DeclareHistogram1D(D_ENERGY_DECAY_BEAMSTOP,
		       energyBins, "Event/10 alpha beam stopped");
//    DeclareHistogram1D(D_ENERGY_DECAY_NO_SI,
//		       energyBins2, "Event/10 Decay, No Si");
    DeclareHistogram1D(D_ENERGY_ALL_BEAMSTOP,
               energyBins, "Event/10 All with beam off");

    DeclareHistogram1D(D_ENERGY_WITH_VETO, energyBins, 
                      "Fr.E/10 veto");
    DeclareHistogram1D(D_ENERGY_WITH_MWPC, energyBins, 
                      "Fr.E/10 mwpc, no veto");
    DeclareHistogram1D(D_ENERGY_WITH_VETO_MWPC, energyBins, 
                      "Fr.E/10 veto and mwpc");
    DeclareHistogram1D(D_ENERGY_NO_VETO_NO_MWPC, energyBins, 
                      "Fr.E/10 no veto/si, no mwpc");
    DeclareHistogram1D(D_ENERGY_NO_SI_NO_VETO, energyBins, 
                      "Fr.E/10 no veto/si");

    DeclareHistogram2D(DD_EVENT_POSITION, 
		       xBins, yBins, "DSSD all events positions");
    DeclareHistogram2D(DD_EVENT_POSITION_FROM_E, 
		       xBins, yBins, "DSSD position all max event");
    DeclareHistogram2D(DD_LIGHT_POSITION, 
		       xBins, yBins, "DSSD position light ion");
    DeclareHistogram2D(DD_UNKNOWN_POSITION, 
		       xBins, yBins, "DSSD position unknown");
    DeclareHistogram2D(DD_FISSION_POSITION, 
		       xBins, yBins, "DSSD position fission");
    DeclareHistogram2D(DD_IMPLANT_POSITION, 
		       xBins, yBins, "DSSD position implant");
    DeclareHistogram2D(DD_DECAY_POSITION, 
		       xBins, yBins, "DSSD position decay");

    /** Trace Information **/
    DeclareHistogram2D(DD_ENERGY_POS_X_TRACE,
		       energyBins, xBins, "DSSD E vs X from Traces");
    DeclareHistogram2D(DD_ENERGY_POS_Y_TRACE,
		       energyBins, yBins, "DSSD E vs Y from Traces");
    DeclareHistogram2D(DD_FRONTE__BACKE, energyBins2, energyBins2,
            "Front vs Back energy (calib / 100)");
    DeclareHistogram2D(DD_POS_X_VS_TRACE,
		       traceBins, xBins, "DSSD Position vs X from Traces");
    DeclareHistogram2D(DD_POS_Y_VS_TRACE,
		       traceBins, xBins, "DSSD Position vs Y from Traces");//Wrong Labels
    DeclareHistogram2D(DD_POS_SI_VS_TRACE,
		       traceBins, numTraces, "Si Box Traces");
    /** Check Gain match alpha region events **/
    DeclareHistogram2D(DD_EVENT_ENERGY__X_POSITION_IMP,
		       energyBins, xBins, "DSSD X strips E vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY__Y_POSITION_IMP,
		       energyBins, yBins, "DSSD Y strips E vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY__X_POSITION_DEC,
		       energyBins, xBins, "DSSD X strips E vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY__Y_POSITION_DEC,
		       energyBins, yBins, "DSSD Y strips E vs. position");
    /** Check Gain match for HE events **/
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_X_POSITION_IMP,
		       energyBins2, xBins, "DSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_Y_POSITION_IMP,
		       energyBins2, yBins, "DSSD Y strips E/100 vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_X_POSITION_DEC,
		       energyBins2, xBins, "DSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_EVENT_ENERGY_COMP_Y_POSITION_DEC,
		       energyBins2, yBins, "DSSD Y strips E/100 vs. position");

    /** Check Gain match spectra via MaxEvents routine **/

    DeclareHistogram2D(DD_MAXEVENT_ENERGY__X_POSITION_IMP,
		       energyBins, xBins, "MAXDSSD X strips E vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY__Y_POSITION_IMP,
		       energyBins, yBins, "MAXDSSD Y strips E vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY__X_POSITION_DEC,
		       energyBins, xBins, "MAXDSSD X strips E vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY__Y_POSITION_DEC,
		       energyBins, yBins, "MAXDSSD Y strips E vs. position");

    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_X_POSITION_IMP,
		       energyBins2, xBins, "MAXDSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_IMP,
		       energyBins2, yBins, "MAXDSSD Y strips E/100 vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_X_POSITION_DEC,
		       energyBins2, xBins, "MAXDSSD X strips E/100 vs. position");
    DeclareHistogram2D(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_DEC,
		       energyBins2, yBins, "MAXDSSD Y strips E/100 vs. position");
    DeclareHistogram2D(DD_FRONT_BACK_DTIME_VS_POSITION,SA,S5,"Time difference between preamp groups of 16, total 24 (8x3).");

    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 0, energyBins, timeBins,
		       "DSSD Ty,Ex (10ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 1, energyBins, timeBins,
		       "DSSD Ty,Ex (100ns/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 2, energyBins, timeBins,
		       "DSSD Ty,Ex (1us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 3, energyBins, timeBins,
		       "DSSD Ty,Ex (10us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 4, energyBins, timeBins,
		       "DSSD Ty,Ex (100us/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 5, energyBins, timeBins,
		       "DSSD Ty,Ex (1ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 6, energyBins, timeBins,
		       "DSSD Ty,Ex (10ms/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 7, energyBins, timeBins,
		       "DSSD Ty,Ex (1s/ch)(xkeV)");
    DeclareHistogram2D(DD_ENERGY_DECAY_TIME_GRANX + 8, energyBins, timeBins,
		       "DSSD Ty,Ex (10s/ch)(xkeV)");
    /** Diagnostics **/ 
    DeclareHistogram2D(DD_ENERGY__POSX_T_MISSING,
		       energyBins, xBins, "DSSD T missing X strips E vs. position");
    DeclareHistogram2D(DD_ENERGY__POSY_T_MISSING,
		       energyBins, yBins, "DSSD T missing Y strips E vs. position");
    if (TempCheck_>=0) {
	DeclareHistogram2D(DD_ALPHA_XENERGY_VS_TIME,energyBins,SC ,"XEnergy between 6 and 10 MeV vs. Time." );
	DeclareHistogram2D(DD_ALPHA_YENERGY_VS_TIME,energyBins,SC ,"YEnergy between 6 and 10 MeV vs. Time." );
    }
    /** Check how many strips and how far fired **/
    /**DeclareHistogram2D(DD_DENERGY__DPOS_X_CORRELATED,
		       energyBins, xBins, "DSSD dE dX correlated events");
    DeclareHistogram2D(DD_DENERGY__DPOS_Y_CORRELATED,
		       energyBins, yBins, "DSSD dE dY correlated events");**/
    /** Pixel Correlated Events **/    
    DeclareHistogram2D(DD_CHAINS_ENERGY_V_TIME,SB,SC,
	    "Event Chains vs. Log of dT between Events");
    DeclareHistogram2D(DD_OTHER_ENERGY_V_TIME,SB,SC,
	    "Other Event Chains vs. Log of dT between");
/*    DeclareHistogram2D(DD_CHAIN_NUM_FISSION, SC, S9, 
            "Event Number vs. Fission Chain Energy/20");
    DeclareHistogram2D(DD_CHAIN_NUM_ALPHA, SC, S9, 
            "Event Number vs. Alpha Chain Energy/10");*/
    DeclareHistogram2D(DD_CHAIN_ALPHA_V_ALPHA, SC, SC,
	    "First Alpha vs. Alpha Chain Energy/10");
    DeclareHistogram2D(DD_OTHER_ALPHA_V_ALPHA, SC, SC,
	    "First Alpha vs. Other Alpha Chain Energy/10");
    DeclareHistogram2D(DD_TOF_A_EVENT, SC, S9,
	    "TOF vs. Event Number Alpha");
    DeclareHistogram2D(DD_TOF_SF_EVENT, SC, S9,
	    "TOF vs. Event Number R-SF");
    DeclareHistogram2D(DD_TOF_O_EVENT, SC, SA,
	    "TOF vs. Event Number R-Other");
    DeclareHistogram2D(DD_MWPC_ENERGY_A_EVENT, SD, S9,
	    "MWPC Energy vs. Event Number Alpha");
    DeclareHistogram2D(DD_MWPC_ENERGY_SF_EVENT, SD, S9,
	    "MWPC Energy vs. Event Number R-SF"); 
    DeclareHistogram2D(DD_MWPC_ENERGY_V_SFE, SB, SB,
	    "MWPC Energy/10 vs. SF Energy/100 R-SF"); 
    DeclareHistogram2D(DD_MWPC_ENERGY_O_EVENT, SD, SA,
	    "MWPC Energy vs. Event Number R-Other");    
    DeclareHistogram2D(DD_FRONT_V_ALPHA1, SE, xBins, 
  	    "Front Position vs. Alpha Energy 1 from interesting events in SheCorrelator chain");
    DeclareHistogram2D(DD_BACK_V_ALPHA1, SE, yBins, 
	    "Back Position vs. Alpha Energy 1 from interesting events in SheCorrelator chain");	
    DeclareHistogram1D(D_TOF_HEAVY, SE,"Avg. TOF of Heavy Implants");
    DeclareHistogram1D(D_TOF_LIGHT, SE,"Avg. TOF of Light Particles"); 
    DeclareHistogram1D(D_TOF_UNK, SE,"Avg. TOF of 'Unknown' Particles"); 
    DeclareHistogram1D(D_MWPC_ENERGY_HEAVY, SE,"Avg. MWPC Energy of Heavy Implants");
    DeclareHistogram1D(D_MWPC_ENERGY_LIGHT, SE,"Avg. MWPC Energy of Light Particles"); 
    DeclareHistogram1D(D_MWPC_ENERGY_UNK, SE,"Avg. MWPC Energy of 'Unknown' Particles");


    DeclareHistogram2D(DD_ENERGY_YVX_BAD_MATCH,energyBins2,energyBins2, "Energy/100 Y vs X if x or y is poorly matched in energy"); 
    DeclareHistogram2D(DD_EVT_Y_BAD_MATCH,energyBins2 ,SB,"Energy/10 V FrontBack time if X is missing energy");
    DeclareHistogram2D(DD_EVT_X_BAD_MATCH,energyBins2 ,SB,"Energy/10 V FrontBack time if Y is missing energy");

    DeclareHistogram2D(DD_SI_RAW_EDSSD+0,SD,SC,"Raw Si vs Cal DSSD 1");
    DeclareHistogram2D(DD_SI_RAW_EDSSD+1,SD,SB,"Raw Si vs Cal DSSD 2");
    DeclareHistogram2D(DD_SI_RAW_EDSSD+2,SD,SB,"Raw Si vs Cal DSSD 3");
    DeclareHistogram2D(DD_SI_RAW_EDSSD+3,SD,SB,"Raw Si vs Cal DSSD 4");
    DeclareHistogram2D(DD_SI_RAW_EDSSD+4,SD,SB,"Raw Si vs Cal DSSD 5");
    DeclareHistogram2D(DD_SI_RAW_EDSSD+5,SD,SB,"Raw Si vs Cal DSSD 6");
    DeclareHistogram2D(DD_SI_RAW_EDSSD+6,SD,SC,"Raw Si vs Cal DSSD 0");
    
    
    
}


bool Dssd4SHEProcessor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return false;

    xyEventsTMatch_.clear();
    xyEventsEMatch_.clear();

    vector<ChanEvent*> xEvents = 
        event.GetSummary("dssd_back:dssd_back", true)->GetList();
    vector<ChanEvent*> yEvents = 
        event.GetSummary("dssd_front:dssd_front", true)->GetList();

    /**
     * Matching the front-back by the time correlations
     */
    vector< pair<StripEvent, bool> > xEventsTMatch;
    vector< pair<StripEvent, bool> > yEventsTMatch;
    StripEvent ev2x;
    StripEvent ev2y;
    int numxTrace=0;

    for (vector<ChanEvent*>::iterator itx = xEvents.begin();
         itx != xEvents.end();
         ++itx) { //PreProcess traces in X
	int xid = (*itx)->GetChanID().GetLocation();
	double xtime = (*itx)->GetTime();
        if (xid >48 && xid<65) { //Correct for delays by an offset 
	    (*itx)->SetTime(xtime+XToffset1_);
	}
	StripEvent ev((*itx)->GetCalEnergy(), 
                      (*itx)->GetTime(),
                      (*itx)->GetChanID().GetLocation(),
                      (*itx)->IsSaturated(),
		      0,	      
		      (*itx)->GetTrace());
        pair<StripEvent, bool> match(ev, false);
        xEventsTMatch.push_back(match);
	
        Trace& trace = (*itx)->GetTrace();
        int pulses = trace.GetValue("numPulses");
	if (pulses >= 0) {

	    trace.Plot(DD_POS_X_VS_TRACE, numxTrace);
	    numxTrace++;
	}
        /** Handle additional pulses (no. 2, 3, ...) */
        
        for (int i = 1; i < pulses; ++i) {
            stringstream energyCalName;
            energyCalName << "filterEnergy" << i + 1 << "Cal";
            stringstream timeName;
            timeName << "filterTime" << i + 1;

            ev.pileup = true;

            
            ev2x.E = trace.GetValue(energyCalName.str());
            ev2x.t = (trace.GetValue(timeName.str()) - 
                     trace.GetValue("filterTime") + ev.t);
            ev2x.pos = ev.pos;
            ev2x.sat = false;
            ev2x.pileup = true;
	    
            pair<StripEvent, bool> match2(ev2x, false);
            xEventsTMatch.push_back(match2);
	    
    
        }

        /*theTime = driver -> GetWallTime(ev.t);
        if ( abs(theTime-1484649104) < 1) {
                 sst << 0 << " " << ev.E << " " << ev.t << " " << ev.pos << " " << ev.sat << endl;
         	 Notebook::get()->report(sst.str());
	}*/
    }
    int numyTrace=0;
    for (vector<ChanEvent*>::iterator ity = yEvents.begin();
         ity != yEvents.end();
         ++ity) { //PreProcess traces in Y.
	
	int yid = (*ity)->GetChanID().GetLocation(); 
	double ytime = (*ity)->GetTime();
        if (yid<=16) { //Correct for delays by an offset ?!
	    (*ity)->SetTime(ytime+YToffset1_);
	} else if (yid>16 ) { 
	    (*ity)->SetTime(ytime+YToffset2_);
	} 
        StripEvent ev((*ity)->GetCalEnergy(), 
                      (*ity)->GetTime(),
                      (*ity)->GetChanID().GetLocation(),
                      (*ity)->IsSaturated(),
		      0,
		      (*ity)->GetTrace());
        pair<StripEvent, bool> match(ev, false);
        yEventsTMatch.push_back(match);

        const Trace& trace = (*ity)->GetTrace();
        int pulses = trace.GetValue("numPulses");
	if (pulses >=1) {
	    trace.plot(DD_POS_Y_VS_TRACE, numyTrace);
	    numyTrace++;
	}
        for (int i = 1; i < pulses; ++i) {
            stringstream energyCalName;
            energyCalName << "filterEnergy" << i + 1 << "Cal";
            stringstream timeName;
            timeName << "filterTime" << i + 1;

            ev.pileup = true;


            ev2y.E = trace.GetValue(energyCalName.str());
            ev2y.t = (trace.GetValue(timeName.str()) - 
                     trace.GetValue("filterTime") + ev.t);
            ev2y.pos = ev.pos;
            ev2y.sat = false;
            ev2y.pileup = true;
            pair<StripEvent, bool> match2(ev2y, false);
            yEventsTMatch.push_back(match2);
	    
            if (i > 1 && abs(1-ev2x.E/ev2y.E) < 0.3) {

                plot(DD_ENERGY_POS_Y_TRACE, ev2y.E, ev2y.pos);
		plot(DD_ENERGY_POS_X_TRACE, ev2x.E, ev2x.pos);


            }
        }

        /* theTime = driver -> GetWallTime(ev.t);
        if ( abs(theTime-1484649104) < 1) {
                 sst << 1 << " " << ev.E << " " << ev.t << " " << ev.pos << " " << ev.sat << endl;
         	 Notebook::get()->report(sst.str());
	}*/
    } // It may be possible to move these into the next loop to avoid the repetition.

    for (vector< pair<StripEvent, bool> >::iterator itx = xEventsTMatch.begin();
         itx != xEventsTMatch.end();
         ++itx) {
        double bestDtime = numeric_limits<double>::max();
        vector< pair<StripEvent, bool> >::iterator bestMatch =
            yEventsTMatch.end();
        for (vector< pair<StripEvent, bool> >::iterator ity = 
                                                     yEventsTMatch.begin();
            ity != yEventsTMatch.end();
            ++ity) 
        {  
            
            double energyX = (*itx).first.E;
            double energyY = (*ity).first.E;
	    //double energyXb, energyYb;
	    double xtime=(*itx).first.t;
	    double ytime=(*ity).first.t;
	    //double Xbtime, Ybtime;
	    double xpos=(*itx).first.pos;
	    double ypos=(*ity).first.pos;
	    //double Xbpos, Ybpos;
	    int preampNum=floor(xpos/16)+8*floor(ypos/16);
	    plot(DD_FRONT_BACK_DTIME_VS_POSITION,xtime-ytime+500,preampNum);
            
	    //Make condition if more than one pixel hits. 
	    /*if ( (xtime-Xbtime==0 || ytime-Ybtime==0 )&& abs(xtime-ytime) < timeWindowFB_*1.0e8) {
		if ( ( Ybpos==ypos && abs(Xbpos-xpos)==1 && 
		    (abs(energyYb-(energyXb+energyX))<0.3*energyYb) && energyX>0 && energyXb>0) ) {
		   (*itx).first.pxl =Xbpos*1000 +xpos;
		   (*itx).first.E = energyX + energyXb;
                   if ( energyXb > energyX ) {
                       (*itx).first.t = Xbtime;
                   }
		} else if ( ( Xbpos==xpos && abs(Ybpos-ypos)==1 && 
		    (abs(energyXb-(energyYb+energyY))<0.3*energyXb) && energyY>0 && energyYb>0 ) ){
		   (*ity).first.pxl = Ybpos*1000 +ypos;
		   (*ity).first.E = energyY + energyYb;
                   if ( energyYb > energyY ) {
                       (*ity).first.t = Ybtime;
                   }
		}   
 				
     	    }*/

	    /*theTime = driver -> GetWallTime((*itx).first.t);
	    if ( abs(theTime-1484649104) < 1) {
                 sst << 0 << " " << (*itx).first.E << " " << (*itx).first.t << " " << (*itx).first.pos << " " << (*itx).first.sat << endl;
                 sst << 1 << " " << (*ity).first.E << " " << (*ity).first.t << " " << (*ity).first.pos << " " << (*ity).first.sat << endl;
         	 Notebook::get()->report(sst.str());
	    }
	    energyXb=energyX;
	    energyYb=energyY;
	    Xbtime=xtime;
	    Ybtime=ytime;
	    Xbpos=xpos;
	    Ybpos=ypos;
	    energyX=(*itx).first.E;
	    energyY=(*ity).first.E;
	    */


	    //Make matched pairs.
            // If already matched, skip
            if ((*ity).second)
                continue;

            /** If energies are in lower range and/or not satured
             *  check if delta energy condition is not met, 
             *  if not, skip this event
             *
             *  For saturated events and set 20 MeV energy for difference check. 
	     *  The calibration in this is imprecise, so one cannot correlate
             *  by energy difference and a wider tolerance is given. 
             *  
             **/

            if (!( (*itx).first.sat || (*ity).first.sat ) ) {

	        if( energyX > highEnergyCut_ || energyY > highEnergyCut_ ) {
		    if (abs(energyX - energyY)*2/(energyX+energyY) > 10*deltaEnergy_   ){
                       continue;
	            }
	        } else if (abs(energyX - energyY)*2/(energyX+energyY) > deltaEnergy_   ){
                    continue;
    	        }

	    }
            double dTime = abs((*itx).first.t - (*ity).first.t) *
                                Globals::get()->clockInSeconds();
            if (dTime < bestDtime){
                bestDtime = dTime;
                bestMatch = ity;
            }
        }
	
        if (bestDtime < timeWindowFB_) {
            (*itx).second = true;
            (*bestMatch).second = true;
            xyEventsTMatch_.push_back(
                pair<StripEvent, StripEvent>((*itx).first, (*bestMatch).first));
            plot(D_DTIME, int(bestDtime / 1.0e-8) + 2);
        } else {
            bestDtime = int(bestDtime / 1.0e-8);
            if (bestDtime > S8)
                bestDtime = S8 - 1;
            else if (bestDtime < 0)
                bestDtime = 0;
            plot(D_DTIME, bestDtime);
        }
    }

	int ev2xpos;
	double ev2xE;

           
    for (vector< pair<StripEvent, bool> >::iterator itx = xEventsTMatch.begin();
         itx != xEventsTMatch.end();
         ++itx) {

        ev2xpos = (*itx).first.pos;
        ev2xE = (*itx).first.E;        

        if ((*itx).second) { //if matched search for neighbors
            continue; //to plot missing
        }
        
        plot(DD_ENERGY__POSX_T_MISSING, ev2xE, ev2xpos);
  
    }
	int ev2ypos;

	double ev2yE;

    for (vector< pair<StripEvent, bool> >::iterator ity = yEventsTMatch.begin();
         ity != yEventsTMatch.end();
         ++ity) {
	
        if ((*ity).second)
            continue;

        ev2ypos = (*ity).first.pos;
        ev2yE = (*ity).first.E;

        plot(DD_ENERGY__POSY_T_MISSING, ev2yE, ev2ypos);
    }
    
    /*DetectorDriver* driver = DetectorDriver::get();
    time_t theTime;
    stringstream sst;
    sst.precision(15);
    */
    for (vector< pair<StripEvent, StripEvent> >::iterator it =
        xyEventsTMatch_.begin();it != xyEventsTMatch_.end(); ++it) {

        int xPosition = (*it).first.pos;
        int yPosition = (*it).second.pos;
        double xEnergy = (*it).first.E;
        double yEnergy = (*it).second.E;
	double xtime = (*it).first.t;
	double ytime = (*it).second.t;
        double newXEnergy = 0 , newYEnergy = 0;

        for (vector< pair<StripEvent, bool> >::iterator itx = xEventsTMatch.begin();
            itx != xEventsTMatch.end(); ++itx) {
	    if (abs(xPosition-(*itx).first.pos)==1 && abs(xtime-(*itx).first.t) < 2*timeWindowFB_*1.0e8 ){
               if ( (*it).first.pxl == 0 ) {
                   (*it).first.pxl = xPosition;
                   newXEnergy = xEnergy;
	       }
	       (*it).first.pxl += (*itx).first.pos*1000;
               newXEnergy += (*itx).first.E;
	    }
	    /*theTime = driver -> GetWallTime(xtime);
	    if ( abs(theTime-1484649104) < 1) {
                 sst << 0 << " " << xEnergy << " " << xtime << " " << xPosition << " " << endl;
                 sst << (*itx).first.E << " " << (*itx).first.t << " " << (*itx).first.pos << " " << newXEnergy << endl;
	    }*/


	}

        for (vector< pair<StripEvent, bool> >::iterator ity = yEventsTMatch.begin();
            ity != yEventsTMatch.end(); ++ity) {
 	    if (abs(yPosition-(*ity).first.pos)==1 && abs(ytime-(*ity).first.t) < 2*timeWindowFB_*1.0e8 ){
               if ( (*it).second.pxl == 0 ) {
                   (*it).second.pxl = yPosition;
                   newYEnergy = yEnergy;
	       }
	       (*it).second.pxl += (*ity).first.pos*1000;
               newYEnergy += (*ity).first.E;
	    }
	    /*theTime = driver -> GetWallTime(ytime);
	    if ( abs(theTime-1484649104) < 1) {
                 sst << 1 << " " << yEnergy << " " << ytime << " " << yPosition << " " << endl;
                 sst << (*ity).first.E << " " << (*ity).first.t << " " << (*ity).first.pos << " " << newYEnergy << endl;
         	 Notebook::get()->report(sst.str());
	    }*/

	}

        if ( newXEnergy > highEnergyCut_ || newYEnergy > highEnergyCut_ ) {
            if ( abs(newXEnergy - newYEnergy)*2/(newYEnergy+newYEnergy) < 10*deltaEnergy_ ){
                if (newXEnergy !=0 )
	           (*it).first.E = newXEnergy;
                if (newYEnergy !=0 )
                   (*it).second.E = newYEnergy;
	    } else {
		(*it).first.pxl = -1;
		(*it).second.pxl = -1;		
            }
       } else {
           if ( abs(newXEnergy - newYEnergy)*2/(newXEnergy+newYEnergy) > deltaEnergy_ ){
                if (newXEnergy !=0 )
	           (*it).first.E = newXEnergy;
                if (newYEnergy !=0 )
                   (*it).second.E = newYEnergy;
           } else {
		(*it).first.pxl = -1;
		(*it).second.pxl = -1;		
           }
       }



    }



    /**
     * Matching the front-back by the Energy of the event
     * Using the old style GetMaxEvent for comparison
     */
    if (xEvents.size() > 0 && yEvents.size() > 0) {
            ChanEvent* maxFront =
                event.GetSummary("dssd_back:dssd_back")->GetMaxEvent(true);
            ChanEvent* maxBack = 
                event.GetSummary("dssd_front:dssd_front")->GetMaxEvent(true);
            StripEvent evf(maxFront->GetCalEnergy(), 
                           maxFront->GetTime(),
                           maxFront->GetChanID().GetLocation(),
                           maxFront->IsSaturated(),
                           0,
			   maxFront->GetTrace());
            StripEvent evb(maxBack->GetCalEnergy(), 
                           maxBack->GetTime(),
                           maxBack->GetChanID().GetLocation(),
                           maxBack->IsSaturated(),
                           0,
			   maxBack->GetTrace());
        xyEventsEMatch_.push_back(pair<StripEvent, StripEvent>(evf, evb));
    }

    return true; 
}


bool Dssd4SHEProcessor::Process(RawEvent &event)
{
    using namespace dammIds::dssd4she;

    if (!EventProcessor::Process(event))
        return false;

    vector<ChanEvent*> vetoEvents = 
        event.GetSummary("si:veto", true)->GetList();
    vector<ChanEvent*> sideEvents = 
        event.GetSummary("si:si", true)->GetList();
    vector<ChanEvent*> mwpcEvents = 
        event.GetSummary("mcp:mcp", true)->GetList();
    int mwpc = event.GetSummary("mcp", true)->GetMult();
    bool hasBeam = TreeCorrelator::get()->place("Beam")->status();

   plot(D_MWPC_MULTI, mwpc);

   
	
    for (vector< pair<StripEvent, StripEvent> >::iterator it =
                                                 xyEventsTMatch_.begin();
         it != xyEventsTMatch_.end(); ++it)
    {
        int xPosition = (*it).first.pos;
        int yPosition = (*it).second.pos;
        double xEnergy = (*it).first.E;
        double yEnergy = (*it).second.E;
	int pixel[3] = {(*it).first.pxl,(*it).second.pxl,0};

	int hasPileup = (*it).first.pileup;
	int hasSat =0;
	
        /** If both saturated set to 250 MeV **/
        if ((*it).first.sat && !(*it).second.sat) {
	    plot(D_FRONTE_W_BACK_SAT, yEnergy/20);
            xEnergy = yEnergy;
	    hasSat +=1;
	} else if (!(*it).first.sat && (*it).second.sat) {
	    plot(D_BACKE_W_FRONT_SAT, xEnergy/20);
            yEnergy = xEnergy;
	    hasSat +=1;
	} else if ((*it).first.sat && (*it).second.sat) {
            xEnergy = 250000.0;
            yEnergy = 250000.0;
	    hasSat+=2;
        }
	
	
	



        double time = min((*it).first.t, (*it).second.t);
	
	plot(D_ENERGY_X, xEnergy / 10.0); 
        plot(D_ENERGY_Y, yEnergy / 10.0);
        plot(DD_FRONTE__BACKE, xEnergy / 100.0, yEnergy / 100.0);

	static double firstTempCheck;	
	if (TempCheck_==0 && firstTempCheck==0) {
	    firstTempCheck=time;
	    cout << "using time " << firstTempCheck << " as time=0 for temperature checking." << endl;
	} else if (TempCheck_>0 && firstTempCheck==0) {
	    firstTempCheck=TempCheck_;
	    cout << "using time " << firstTempCheck << " as time=0 for temperature checking." << endl;
	}
        if (firstTempCheck>0 && (xEnergy>4000||yEnergy>4000)) {
	    plot(DD_ALPHA_XENERGY_VS_TIME,(time-firstTempCheck)*1e-9,xEnergy-4000);
	    plot(DD_ALPHA_YENERGY_VS_TIME,(time-firstTempCheck)*1e-9 ,yEnergy-4000);
	    
	}
	/*Needs redone*/
        /*double static mwpcTime;
        double static mwpcdT;
        double static MwpcE;

	int k=1;

        for (vector<ChanEvent*>::iterator itm = mwpcEvents.begin();
             itm != mwpcEvents.end(); ++itm) {
            double static dt=(*mwpcEvents.begin())->GetTime();
            mwpcdT = (*itm)->GetTime()-dt;
	    mwpcTime = (*itm)->GetTime();
            MwpcE = (*itm) -> GetCalEnergy();
       	    if ( hasPileup) {
	        k++;
            }
	    
	}*/
	/*consider first / last method)*/
	double static mwpcTime;
	double tof, MwpcE;
	int k=1;
	for (vector<ChanEvent*>::iterator itm = mwpcEvents.begin();
             itm != mwpcEvents.end(); ++itm) {
	    ChanEvent *chan = *itm;

	    int number 	= chan->GetChanID().GetLocation();
	    MwpcE = chan->GetCalEnergy();
	    mwpcTime   = chan->GetTime();
	    double mwpcT[2];

	    if (number==0) {
		mwpcT[0]=mwpcTime;
	    } else if (number==1) {
		mwpcT[1]=mwpcTime;
	    }

	    if (mwpc >=1) {
		tof=  mwpcT[1] - mwpcT[0]+2000;
 	    }

	    if ( hasPileup) {
	        k++;
            }
	}
        if (vetoEvents.size() > 0) {
            for (vector<ChanEvent*>::iterator itv = vetoEvents.begin();
                itv != vetoEvents.end();
                ++itv) {
                double vetoEnergy = (*itv)->GetCalEnergy();
                plot(DD_DE_E__DSSD_VETO, (vetoEnergy + xEnergy) / 100.0, xEnergy / 100.0);
            }
        }

        double bestSiTime = numeric_limits<double>::max();
        ChanEvent* correlatedSide = 0;
        bool hasEscape = false;
        double escapeEnergy = 0.0;
        double escapeRaw = 0.0;
        int escapePos = 0;
        static int numsTrace = 0;
        for (vector<ChanEvent*>::iterator its = sideEvents.begin();
            its != sideEvents.end();
            ++its) {
            double dt = abs(time - (*its)->GetTime()) *
                        Globals::get()->clockInSeconds();
            if (dt < bestSiTime) {
                bestSiTime = dt;
                correlatedSide = *its;
            }
            Trace& traceSi = (*its)->GetTrace();
            if (!traceSi.empty()) {
                for (unsigned int i=0; i < traceSi.size(); i++) {
                   histo.Plot(DD_POS_SI_VS_TRACE, i, numsTrace, traceSi.at(i));
                }
   	        numsTrace++;
   	    }
        }

        if (correlatedSide != 0) {
            int siTime = int(bestSiTime / 1.0e-8) + 1;
            if (siTime > S8)
                siTime = S8 - 1;
            else if (siTime < 0)
                siTime = 0;
            plot(D_DTIME_SIDE, siTime );
	    escapePos = (correlatedSide->GetChanID().GetLocation() -1);
            escapeEnergy = correlatedSide->GetCalEnergy();
            plot(DD_DTIME_SIDE, siTime, escapePos);
            plot(DD_ENERGY_CORRELATED_SIDE, escapeEnergy,escapePos);
          
	   
	    if (bestSiTime < timeWindowDS_) {

		escapeRaw = correlatedSide->GetEnergy();
		escapePos = (correlatedSide->GetChanID().GetLocation() -1);
		pixel[2] = escapePos;
                plot(D_ENERGY_CORRELATED_SIDE, escapeEnergy);
                hasEscape = true;

            }
        }

	/** if pileup second signal has no mwpc **/
	if (hasPileup && mwpc == 2 && k > 1) {
	    mwpc=0; 
        }
	
	/*if (mwpc >=2) {
	    
	}*/

        bool hasVeto = false;
        static int traceNum1=0;
        static int traceNum2=0;
        if (vetoEvents.size() > 0)
            hasVeto = true;

        if (hasVeto) {
            plot(D_ENERGY_WITH_VETO, yEnergy / 10.0);
    	}
        if (mwpc > 0 && !hasVeto) {
            plot(D_ENERGY_WITH_MWPC, yEnergy / 10.0);
	}
        if (hasVeto && mwpc > 0) {
	    plot(D_ENERGY_WITH_VETO_MWPC, yEnergy / 10.0);
	}
        if (!hasVeto && !hasEscape && mwpc == 0) {
            plot(D_ENERGY_NO_VETO_NO_MWPC, yEnergy / 10.0);
            if (!(*it).first.tr.empty() ) {
	        for (unsigned int j=0; j<(*it).first.tr.size(); j++) {
                    histo.Plot(DD_TRACES_DECAY_CORR_X, j, traceNum1, (*it).first.tr.at(j));
		}
		plot(DD_TRACES_DECAY_CORR_E_V_NX, xEnergy, traceNum1);
		plot(DD_TRACES_DECAY_CORR_POS_V_NX, xPosition, traceNum1);
	      	/*if (xEnergy <= 6350 && xEnergy >= 5800) {
		   for (int m=0; m<xPosition+1; m++) {
		      if (xEnergy <= 5800 && xEnergy >= 5700 ) {
			   plot(D_TRACES_DECAY_CORR_GATE1_X,traceNum1);
		      } else if (xEnergy <= 6100 && xEnergy >= 6000) {
			   plot(D_TRACES_DECAY_CORR_GATE2_X,traceNum1);
		      } else if (xEnergy <= 6350 && xEnergy >= 6200) {
			   plot(D_TRACES_DECAY_CORR_GATE3_X,traceNum1);
		      }
		   }
		}*/
	        
                traceNum1++;
            } else if (!(*it).second.tr.empty()) {
	        for (unsigned int j=0; j<(*it).second.tr.size(); j++) {
	            histo.Plot(DD_TRACES_DECAY_CORR_Y, j, traceNum2, (*it).second.tr.at(j));
		}
		plot(DD_TRACES_DECAY_CORR_E_V_NY, yEnergy, traceNum2);
		plot(DD_TRACES_DECAY_CORR_POS_V_NY, yPosition, traceNum2);
		/*if (yEnergy <= 6350 && yEnergy >= 5800) {
		   for (int m=0; m<yPosition+1; m++) {
	              if (yEnergy <= 5800 && yEnergy >= 5700 ) {
			   plot(D_TRACES_DECAY_CORR_GATE1_Y,traceNum2);
		      } else if (yEnergy <= 6100 && yEnergy >= 6000) {
			   plot(D_TRACES_DECAY_CORR_GATE2_Y,traceNum2);
		      } else if (yEnergy <= 6350 && yEnergy >= 6200) {
			   plot(D_TRACES_DECAY_CORR_GATE3_Y,traceNum2);
		      }
		   }
		}*/
	        traceNum2++;                
            }
	}
        if (!hasVeto && !hasEscape) {
            plot(D_ENERGY_NO_SI_NO_VETO, yEnergy / 10.0);
	}

        double XYRatio=xEnergy/yEnergy;
	if (XYRatio>=0 && abs((XYRatio-1)/(XYRatio+1))>=0.01){ //Bad Match 
	    double dtxy = (*it).first.t-(*it).second.t;
	    plot(DD_ENERGY_YVX_BAD_MATCH,yEnergy/100,xEnergy/100);
	    plot(DD_EVT_Y_BAD_MATCH,yEnergy/10,dtxy+500);
    	plot(DD_EVT_X_BAD_MATCH,xEnergy/10,dtxy+500);
	    if ( hasSat > 0)  {
	        XYRatio=-hasSat;
	    }
	} /*else {
	    XYRatio=0.0;
	}
	if ( abs((XYRatio-1)/(XYRatio+1))<=0.01 ) {
	    XYRatio=0.0;
	} else if ( hasSat > 0)  {
	    XYRatio=-hasSat;
	}*/
	
 
	SheEvent event = SheEvent((xEnergy +yEnergy)/2 + escapeEnergy, XYRatio,time, mwpc, tof ,MwpcE, hasBeam, hasVeto, escapeEnergy*2/(xEnergy +yEnergy), pixel, unknown);
	
        pickEventType(event);
	
	double eventE=event.get_energy();
        plot(DD_EVENT_POSITION, xPosition, yPosition);    
    if (!event.get_beam() ) {
        plot(D_ENERGY_ALL_BEAMSTOP, eventE/10);
        if (event.get_type() == alpha)
            plot(D_ENERGY_DECAY_BEAMSTOP, eventE/10);
    }
	if (eventE>=0) {
	    plot(D_ENERGY_ALL, eventE/10);
	}
	if (event.get_type() == heavyIon) {
	    plot(D_TOF_HEAVY,event.get_mwpcTime());//+2000);
	    plot(D_MWPC_ENERGY_HEAVY, event.get_mwpcEnergy());
            plot(DD_IMPLANT_POSITION, xPosition, yPosition);
            plot(D_ENERGY_IMPLANT, eventE/10);
      	    plot(DD_EVENT_ENERGY__X_POSITION_IMP, eventE, xPosition);
            plot(DD_EVENT_ENERGY__Y_POSITION_IMP, eventE, yPosition);
            plot(DD_EVENT_ENERGY_COMP_X_POSITION_IMP, eventE/100, xPosition);
            plot(DD_EVENT_ENERGY_COMP_Y_POSITION_IMP, eventE/100, yPosition);
		}
        else if (event.get_type() == alpha) {
	    plot(D_ENERGY_DECAY, eventE / 10.0);
            plot(DD_DECAY_POSITION, xPosition, yPosition);
            plot(DD_EVENT_ENERGY__X_POSITION_DEC,eventE, xPosition);
            plot(DD_EVENT_ENERGY__Y_POSITION_DEC, eventE, yPosition);
            plot(DD_EVENT_ENERGY_COMP_X_POSITION_DEC, eventE/100, xPosition);
            plot(DD_EVENT_ENERGY_COMP_Y_POSITION_DEC, eventE/100, yPosition);
	   /* if (!hasEscape) {
		plot(D_ENERGY_DECAY_NO_SI,eventE/10);
	    }*/
	 }
        else if (event.get_type() == lightIon) {
	    plot(D_TOF_LIGHT,event.get_mwpcTime());//+2000);
	    plot(D_MWPC_ENERGY_LIGHT, event.get_mwpcEnergy());
            plot(DD_LIGHT_POSITION, xPosition, yPosition);
            plot(D_ENERGY_LIGHT, eventE / 10.0);
       	    plot(DD_EVENT_ENERGY__X_POSITION_IMP, eventE, xPosition);
            plot(DD_EVENT_ENERGY__Y_POSITION_IMP, eventE, yPosition);
            plot(DD_EVENT_ENERGY_COMP_X_POSITION_IMP, eventE/100, xPosition);
            plot(DD_EVENT_ENERGY_COMP_Y_POSITION_IMP, eventE/100, yPosition);
	 }
        else if (event.get_type() == unknown) {
	    plot(D_TOF_UNK,event.get_mwpcTime());//+2000);
	    plot(D_MWPC_ENERGY_UNK, event.get_mwpcEnergy());
            plot(DD_UNKNOWN_POSITION, xPosition, yPosition);
            plot(D_ENERGY_UNKNOWN, eventE / 10.0);
        }
        else if (event.get_type() == fission) {
            plot(DD_FISSION_POSITION, xPosition, yPosition);
            plot(D_ENERGY_FISSION, eventE / 100.0);
            plot(DD_EVENT_ENERGY__X_POSITION_DEC, eventE, xPosition);
            plot(DD_EVENT_ENERGY__Y_POSITION_DEC, eventE, yPosition);
            plot(DD_EVENT_ENERGY_COMP_X_POSITION_DEC, eventE/100, xPosition);
            plot(DD_EVENT_ENERGY_COMP_Y_POSITION_DEC, eventE/100, yPosition);
	   }

	if (event.get_type()== alpha) {
	    const unsigned int NumGranularities = 9;
	    //time resolution in seconds per bin
	    const double timeResolution[NumGranularities] = 
		{10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-1, 10};
	   
            for (unsigned int i = 0; i < NumGranularities; i++) {
		double timeBin = ((event.get_time()-mwpcTime)
			*Globals::get()->clockInSeconds())/timeResolution[i];
		double energyBin = eventE;
//		cout << energyBin << " " << timeBin << endl; 
//' ' << ( event.get_time()*Globals::get()->clockInSeconds() )  << " " << ( mwpcTime*Globals::get()->clockInSeconds() ) <<  endl;
		plot(DD_ENERGY_DECAY_TIME_GRANX + i,energyBin,timeBin);
            }
	}


	if (hasEscape && event.get_type() == alpha) {
	    plot(DD_SI_RAW_EDSSD+escapePos, yEnergy, escapeRaw);
	}
	//if (  event.get_type()!=heavyIon || ( abs(tof-2000)<1000 && MwpcE > 10 )  ) {
	    correlator_.add_event(event, xPosition, yPosition, histo);
	//}
    }
	
    /** Old style max evefnt for comparison */
    for (vector< pair<StripEvent, StripEvent> >::iterator it =
                                                 xyEventsEMatch_.begin();
         it != xyEventsEMatch_.end();
         ++it) {
	//double time = min((*it).first.t, (*it).second.t);
        double xEnergy = (*it).first.E;
        double yEnergy = (*it).second.E;

        int xPosition = (*it).first.pos;
        int yPosition = (*it).second.pos;
	int mwpc = event.GetSummary("mcp", true)->GetMult();
	vector<ChanEvent*> vetoEvents = event.GetSummary("si:veto", true)->GetList();
        bool hasVeto = false;

        if (vetoEvents.size() > 0){
            hasVeto = true;
	}

        plot(DD_EVENT_POSITION_FROM_E, xPosition, yPosition);

        if (mwpc > 0 && !hasVeto) {
            plot(DD_MAXEVENT_ENERGY__X_POSITION_IMP, xEnergy, xPosition);
            plot(DD_MAXEVENT_ENERGY__Y_POSITION_IMP, yEnergy, yPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_X_POSITION_IMP, xEnergy/100, xPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_IMP, yEnergy/100, yPosition);
	}
        if (!hasVeto && mwpc == 0) {
            plot(DD_MAXEVENT_ENERGY__X_POSITION_DEC, xEnergy, xPosition);
            plot(DD_MAXEVENT_ENERGY__Y_POSITION_DEC, yEnergy, yPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_X_POSITION_DEC, xEnergy/100, xPosition);
            plot(DD_MAXEVENT_ENERGY_COMP_Y_POSITION_DEC, yEnergy/100, yPosition);
	}

   }        

    EndProcess();
    return true;
}


bool Dssd4SHEProcessor::pickEventType(SheEvent& event) {
    /**
     * Logic table (V - veto, M - mwpc, B - beam )
     * Logic state is converted into a numerical value N
     * like a binary number:
     *
     * V M B | N | decision
     * --------------------
     * 0 0 0 | 0 | check / alpha / fission (depending on energy)
     * 0 0 1 | 1 | ""
     * 0 1 0 | 2 | unknown (likely scattered)
     * 0 1 1 | 3 | heavyIon / check (depending on energy)
     * 1 0 0 | 4 | unknown  (likely scattered)
     * 1 0 1 | 5 | lightIon
     * 1 1 0 | 6 | unknown (likely scattered)
     * 1 1 1 | 7 | lightIon
     *
     **/
    int condition = 0;
    if (event.get_beam()) 
        condition += 1;
    if (event.get_mwpc() > 0) 
        condition += 2;
    if (event.get_veto()) 
        condition += 4;

    if (condition == 0) {
        double energy = event.get_energy();
        if (energy < lowEnergyCut_)
            event.set_type(check);
        else if (energy < highEnergyCut_)
            event.set_type(alpha);
        else if (energy < fissionEnergyCut_)
            event.set_type(unknown); //Check Energy Conditions!
        else 
            event.set_type(fission);
    } 
    else if (condition == 1) {
        double energy = event.get_energy();
        if (energy < lowEnergyCut_)
            event.set_type(check);
        else if (energy < highEnergyCut_)
            event.set_type(alpha);
        else if (energy < fissionEnergyCut_)
            event.set_type(unknown); //Check Energy Conditions!
        else 
            event.set_type(fission);
    } 
    else if (condition == 2 || 
             condition == 4 ||
             condition == 6) {
        event.set_type(unknown);
    } 
    else if (condition == 3 ) {
	if (event.get_energy() > recoilEnergyCut_ && event.get_energy() < highEnergyCut_) {
            event.set_type(heavyIon);
	} else {
	    event.set_type(check);
	}
    }
    else if (condition == 5 || condition == 7) {
        event.set_type(lightIon);
    }
    else
        event.set_type(unknown);

    return true;
}


