/**    \file DoubleTraceAnalyzer.cpp
 *     \brief Identifies double traces.
 *
 *     Implements a quick online trapezoidal filtering mechanism
 *     for the identification of double pulses
 *
 *     - SNL - 7-2-07 - created
 *     - SNL - 2-4-08 - Add plotting spectra
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <numeric>

#include <cstdlib>

#include "DammPlotIds.hpp"
#include "RandomPool.hpp"
#include "Trace.hpp"
#include "DoubleTraceAnalyzer.hpp"
#include "Messenger.hpp"
#include "Globals.hpp"

using namespace std;

int DoubleTraceAnalyzer::numDoubleTraces = 0;
int DoubleTraceAnalyzer::numDoubleTracesCut = 0;
int DoubleTraceAnalyzer::numRejectTraces =0;
/**
 * Set default values for time and energy
 */
DoubleTraceAnalyzer::DoubleTraceAnalyzer(double energyScaleFactor,
                                         short fast_rise, short fast_gap,
                                         short fast_threshold,
                                         short energy_rise, short energy_gap,
                                         short slow_rise, short slow_gap,
                                         short slow_threshold) :
    TraceFilterer(energyScaleFactor,
                  fast_rise, fast_gap, fast_threshold,
                  energy_rise, energy_gap,
                  slow_rise, slow_gap, slow_threshold)
{
}


DoubleTraceAnalyzer::~DoubleTraceAnalyzer() 
{
    // do nothing
}

void DoubleTraceAnalyzer::DeclarePlots()
{
    using namespace dammIds::trace::doubletraceanalyzer;

    TraceFilterer::DeclarePlots();

    const int energyBins = SE;
    const int energyBins2 = SB;
    const int timeBins = SA;
    const int traceBins = dammIds::trace::traceBins;

    Trace sample_trace = Trace();
    unsigned short numTraces = Globals::get()->numTraces();

    sample_trace.DeclareHistogram1D(D_ENERGY2, energyBins, "E2 from traces");
    sample_trace.DeclareHistogram1D(D_ENERGY2_E1CUT, energyBins, "E2 from traces, when E1 > 5 Mev");
    sample_trace.DeclareHistogram2D(DD_DOUBLE_TRACE, traceBins, numTraces,
                                    "Double traces");
    sample_trace.DeclareHistogram2D(DD_REJECT_DOUBLE_TRACE, traceBins, numTraces,
                                    "Rejected traces based on num pulse");
    sample_trace.DeclareHistogram2D(DD_DOUBLE_TRACE_E1CUT, traceBins, numTraces,
                                    "Double traces, when E1 > 5MeV");
    sample_trace.DeclareHistogram2D(DD_ENERGY2__TDIFF, energyBins2, timeBins,
                                    "E2 vs DT");
    sample_trace.DeclareHistogram2D(DD_ENERGY1__TDIFF, energyBins2, timeBins,
                                    "E1 vs DT");
    sample_trace.DeclareHistogram2D(DD_ENERGY2__ENERGY1, energyBins2,
                                    energyBins2, "E2/1000 vs E1/1000");

    sample_trace.DeclareHistogram2D(DD_TRIPLE_TRACE, traceBins, 
                                    numTraces, "Interesting triple traces");
    sample_trace.DeclareHistogram2D(DD_TRIPLE_TRACE_FILTER1, traceBins,
                                numTraces, "Interesting traces (fast filter)");
    sample_trace.DeclareHistogram2D(DD_TRIPLE_TRACE_FILTER2, traceBins,
                            numTraces, "Interesting traces (energy filter)");
    if (useThirdFilter) 
    {                           
        sample_trace.DeclareHistogram2D(DD_TRIPLE_TRACE_FILTER3, traceBins,
                                numTraces, "Interesting traces (3rd filter)");
    }
}

/**
 *   Detect a second crossing of the fast filter corresponding to a piled-up
 *     trace and deduce its energy
 */
void DoubleTraceAnalyzer::Analyze(Trace &trace, 
				  const string &type, const string &subtype)
{    
    if (subtype == "top" || subtype == "bottom")
        return;

    Messenger m;

    TraceFilterer::Analyze(trace, type, subtype);
    // class to see when the fast filter falls below threshold
    static binder2nd< less<Trace::value_type> > recrossesThreshold
	(less<Trace::value_type>(), fastThreshold);

    if ( pulse.isFound && level >= 10 ) {
        /*
         * Show number of traces in messenger
        stringstream ss;
        ss << "Double trace #" << numDoubleTraces << " for type " 
           << type << ":" << subtype;
        m.run_message(ss.str());
        */

        // trace filterer found a first pulse

        Trace::iterator iThr = fastFilter.begin() + pulse.time;
        Trace::iterator iHigh = fastFilter.end();

        vector<PulseInfo> pulseVec;
        // put the original pulse in the vector
        pulseVec.push_back(pulse);
        

        while (iThr < iHigh) {
            // find the trailing edge (use rise samples?)
            advance(iThr, fastParms.GetGapSamples());
            iThr = find_if(iThr, iHigh, recrossesThreshold);					
            // advance(iThr, fastParms.GetSize());
            advance(iThr, fastParms.GetRiseSamples());

            FindPulse(iThr, iHigh);
            if (pulse.isFound) {
                pulseVec.push_back(pulse);
                iThr = fastFilter.begin() + pulse.time;
            } else break;

        } // while searching for multiple traces
        
        trace.SetValue("numPulses", (int)pulseVec.size());
        // now plot stuff        
        if ( pulseVec.size() > 1 ) {
            using namespace dammIds::trace::doubletraceanalyzer;

            const size_t pulseLimit = 5; // maximum number of pulses to find
            if (pulseVec.size() > pulseLimit) {

                /*stringstream ss;
                ss << "Too many pulses, limit = " 
                   << pulseLimit << ", breaking out.";
                m.warning(ss.str());*/
                trace.Plot(DD_REJECT_DOUBLE_TRACE,numRejectTraces);
                numRejectTraces++;
        
                EndAnalyze(); // update timing
                return;
            }
            if (pulseVec.size() > 1) {
                /*cout << pulseVec[0].time << " t1 " << pulseVec[1].time << " t2 " ;
                cout << pulseVec[2].time << " t3 " << pulseVec.size() << " s " << numRejectTraces << endl;*/
                double dt = (pulseVec[1].time - pulseVec[0].time);
                if (abs(dt) < 10.)
                {
                    trace.Plot(DD_REJECT_DOUBLE_TRACE,numRejectTraces);
                    numRejectTraces++;
        
                    EndAnalyze(); // update timing
                    return;
                }
            }
            // fill the trace info
            // first pulse info is set in TraceFilterer
            for (Trace::size_type i=1; i < pulseVec.size(); i++) {
                stringstream str;
                // the first pulse in the vector is the SECOND pulse in the trace
                str << "filterEnergy" << i+1;
                trace.SetValue(str.str(), pulseVec[i].energy);
                str.str(""); // clear the string
                str << "filterTime" << i+1;
                trace.SetValue(str.str(), (int)pulseVec[i].time);
            }
            
            // plot the double pulse stuff
            trace.Plot(DD_DOUBLE_TRACE,numDoubleTraces);

	    if (pulseVec[0].energy > 5000. || pulseVec[1].energy > 50000.) {
                trace.plot(D_ENERGY2_E1CUT, pulseVec[1].energy/100.);
		trace.Plot(DD_DOUBLE_TRACE_E1CUT,numDoubleTracesCut);
		numDoubleTracesCut++;
	    } 

            if (pulseVec.size() > 2  && pulseVec.size()<10) { //&& pulseVec.size() < pulseLimit?
                static int numTripleTraces = 0;
                if (numTripleTraces == 82) 
                {

                stringstream ss;
                ss << "Found triple trace " << numTripleTraces 
                   << ", num pulses = " << pulseVec.size();
                   for (int i =0; i<= pulseVec.size(); i++) 
                   {
                     ss << pulseVec[i].energy << " ";
                   }
                   ss << ", sigma baseline = " << trace.GetValue("sigmaBaseline");
                m.run_message(ss.str());
                }
                //if (numTripleTraces <20) cout << numTripleTraces << " " << pulseVec.size() << endl;
                trace.Plot(DD_TRIPLE_TRACE, numTripleTraces);
                fastFilter.ScalePlot(DD_TRIPLE_TRACE_FILTER1, numTripleTraces,
                                    fastParms.GetRiseSamples());
                energyFilter.ScalePlot(DD_TRIPLE_TRACE_FILTER2, numTripleTraces,
                                    energyParms.GetRiseSamples());
                if (useThirdFilter)
                    thirdFilter.ScalePlot(DD_TRIPLE_TRACE_FILTER3,
                                numTripleTraces, thirdParms.GetRiseSamples());
                
                numTripleTraces++;
            }

            trace.plot(D_ENERGY2, pulseVec[1].energy);
            trace.plot(DD_ENERGY1__TDIFF, 
                pulseVec[0].energy/100, pulseVec[1].time - pulseVec[0].time);
            trace.plot(DD_ENERGY2__TDIFF, 
                pulseVec[1].energy/100, pulseVec[1].time - pulseVec[0].time);
            trace.plot(DD_ENERGY2__ENERGY1, 
                pulseVec[1].energy/100, pulseVec[0].energy/100);

            numDoubleTraces++;
        } // if found double trace
    } // sufficient analysis level

    EndAnalyze(trace);
}
