/** \file McpProcessor.cpp
 * Mcp class takes the signals from the mcp detector
 *   and calculates a 2D position based on the readout
 * SNL 2-2-08, Modified DTM 9-09
 */
#ifdef useroot
#include <TTree.h>
#endif

#include "DammPlotIds.hpp"
#include "McpProcessor.hpp"
#include "RawEvent.hpp"

using std::string;
using std::vector;
using namespace dammIds::mcp;

namespace dammIds {
    namespace mcp {	
        const int D_POSX   = 1;
        const int D_POSY   = 2;
        const int DD_POSXY = 3;
        const int D_AMP1 = 4;
        const int D_AMP2 = 5;
        const int DD_XE_YE = 6;
    }
}


void McpProcessor::McpData::Clear(void)
{
  for (size_t i=0; i < nPos; i++)
    raw[i] = 0;
  xpos = ypos = 0.;

  mult = 0;
}

McpProcessor::McpProcessor(void) : EventProcessor(OFFSET, RANGE, "mcp")
{
  associatedTypes.insert("mcp");
}

void McpProcessor::DeclarePlots(void)
{

  const int posBins   = SE;
  const int posBins2D = S9;

  DeclareHistogram1D(D_POSX, posBins, "MWPC Time Difference");
  DeclareHistogram1D(D_POSY, posBins, "Cathode");
  DeclareHistogram2D(DD_POSXY, posBins2D, SA, "(Energy MPC2)/100 vs. Time Diff (dT-1500)/2");
  DeclareHistogram2D(DD_XE_YE, posBins2D, posBins2D, "Energy_1/5 vs. Energy_2/5");
  DeclareHistogram1D(D_AMP1, posBins, "MWPC 1 Amplitude");
  DeclareHistogram1D(D_AMP2, posBins, "MWPC 2 Amplitude");
}


bool McpProcessor::Process(RawEvent &event)
{
  if (!EventProcessor::Process(event))
    return false;

  double qTotal, qRight, qTop, Cathode;
  
// static const vector<ChanEvent*> &mcpEvents = sumMap["mcp"]->GetList();
     static const vector<ChanEvent*> &mcpEvents = 
        event.GetSummary("mcp:mcp", true)->GetList();
  data.Clear();
  for (vector<ChanEvent*>::const_iterator it = mcpEvents.begin();
       it != mcpEvents.end(); it++) {
      ChanEvent *chan = *it;

      string subtype   = chan->GetChanID().GetSubtype();
      int number 	= chan->GetChanID().GetLocation();
      double calEnergy = chan->GetCalEnergy();
      double mcpTime   = chan->GetTime();
      int hasPileup = chan -> IsPileup();

      using std::cout;
      using std:: endl; 
      /*static double cnt1,cnt2;
      double Ratio;
	Ratio=cnt1/cnt2;
        if (hasPileup) {
	    cout << calEnergy << " mwpc Energy" << endl;        
       	} */
	

      if(subtype == "1time") { 
	  // do nothing
        Cathode=calEnergy;
      } else if (subtype == "mcp" && number==0){
	//  cout << number << endl;
	  data.raw[0] = calEnergy;
	  data.time[0] = mcpTime;
	  data.mult++;
      } else if (subtype == "mcp" && number==1){
	 // cout << number << endl;
	  data.raw[1] = calEnergy;
	  data.mult++;
	  data.time[1] = mcpTime;
      } else if (subtype == "mcp" && number==2){
	  //cout << number << endl;
	  data.raw[2] = calEnergy;
	  data.mult++;
	  data.time[2] = mcpTime;
      } else if (subtype == "mcp" && number==3){
	  data.raw[3] = calEnergy;
	  data.time[3] = mcpTime;
	  data.mult++;
      }
  }
    
  // calculation of position from charge collected on the four corners
    
  // magic numbers here
  //   data.raw[0] *= 1.3;
    
     qTotal = data.raw[0] + data.raw[1] + data.raw[2] + data.raw[3];
    
     qRight = data.raw[3] + data.raw[0];
     //qLeft   = data.raw[2] + data.raw[1];  
     qTop   = data.raw[0] + data.raw[1];
  // qBottom = data.raw[2] + data.raw[3];
    
  data.xpos = (qRight / qTotal) * 512. - 75; //horizontal MCP pos
  data.ypos = (qTop   / qTotal) * 512. - 75; //vertical MCP pos
  // qLeft, qBottom not used
  
  if (data.mult >= 1) {
    using namespace dammIds::mcp;
    double tDiff=  data.time[1] - data.time[0]+2000;
//    data.raw[1]=data.raw[1]/5.;
    plot(D_POSX, tDiff);
    plot(D_POSY, Cathode);
    plot(DD_POSXY,data.raw[1]/100,(tDiff-1500)/2);
    plot(D_AMP1, data.raw[0]/2);
    plot(D_AMP2, data.raw[1]/2);
    plot(DD_XE_YE, data.raw[0]/5,data.raw[1]/5);

    // plot(D_POSX, data.xpos);
    //plot(D_POSY, data.ypos);
    // plot(DD_POSXY, data.xpos, data.ypos);
  }
  
  EndProcess();
  return (data.mult == 3);
}

#ifdef useroot
bool McpProcessor::AddBranch(TTree *tree)
{
  if (tree) {
    TBranch *mcpBranch = 
      tree->Branch(name.c_str(), &data, "raw[4]/D:xpos:ypos:mult/I");

    return (mcpBranch != NULL);
  } 

  return false;
}

void McpProcessor::FillBranch(void)
{
  if (!HasEvent())
    data.Clear();  
}
#endif //useroot
