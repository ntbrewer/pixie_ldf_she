/**
 *   \file Trace.cpp
 *   
 *   Implement how to do our usual tricks with traces
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <vector>

#include "Trace.hpp"

using namespace std;
using namespace dammIds::trace;

namespace dammIds {
    namespace trace {
    }
} // trace namespace

const Trace emptyTrace; ///< an empty trace for const references to point to

/*
 * Plots are static, class-wide variable, so every trace instance has
 * an access to the same histogram range
 */
Plots Trace::histo(OFFSET, RANGE, "traces");

/**
 * Defines how to implement a trapezoidal filter characterized by two
 * moving sum windows of width risetime separated by a length gaptime.
 * Filter is calculated from channels lo to hi.
 */
void Trace::TrapezoidalFilter(Trace &filter, 
			      const TrapezoidalFilterParameters &parms,
			      unsigned int lo, unsigned int hi) const
{
    // don't let the filter work outside of its reasonable range
    lo = max(lo, (unsigned int)parms.GetSize());

    filter.assign(lo, 0);
    
    //! check if we're going to do something bad here
    for (unsigned int i = lo; i < hi; i++) {
        int leftSum = accumulate(begin() + i - parms.GetSize(),
                                 begin() + i - parms.GetRiseSamples() 
                                 - parms.GetGapSamples(), 0);
        int rightSum = accumulate(begin() + i - parms.GetRiseSamples(),
                                  begin() + i, 0);
        filter.push_back(rightSum - leftSum);
    }
}


double Trace::DoBaseline(unsigned int lo, unsigned int numBins)
{
// fits a line to the baseline sets the mean and sigma and returns the min.
    if (size() < lo + numBins) {
        cerr << "Bad range in baseline calculation." << endl;
        return NAN;
    }
    lo = double(lo);
    numBins = double(numBins);
    double hi = lo + numBins;
    
    if (baselineLow == lo && baselineHigh == hi)
        return GetValue("baseline");
    std::vector<double> xVector;
    for (int i = lo; i<= hi; i++) xVector.push_back(i);
    
    double sumX = accumulate(xVector.begin(),xVector.end(),0.0);
    double meanX = sumX/numBins;
    double sumY = accumulate(begin() + lo, begin() + hi, 0.0);
    //double sumYoff = accumulate(begin()+lo+1,begin() + hi,0.0);
    double meanY = sumY / numBins;
    double minY = *min_element(begin()+lo,begin()+hi); 
    double upper= 0, lower=0;
    int numNegSpike =0;

    for (int i=lo; i<=hi; i++) 
    {
        double dx = *(xVector.begin()+i-lo)-meanX;
        double dy = *(begin()+i)-meanY;
        upper += dx*dy;
        lower += pow(dx,2.);
    }
    for (int i=lo; i<=size()/4; i++) 
    {
        double diff = *(begin()+i+1) - *(begin()+i); 
	//cout << diff << " " ; 
	if (diff < -5 ) 
	{
	  numNegSpike++;
	}
    } 
    double slope = upper/lower;
    double offset = meanY - slope*meanX;
    double stdDev = abs(slope) * sqrt(lower/numBins);
    //cout << upper << " "<< lower << " "<< slope << " " << offset << " " << stdDev << " " << minY << endl;
    SetValue("baseline", meanY);
    SetValue("sigmaBaseline", stdDev);
    SetValue("numNegSpike",numNegSpike);
    //cout << "nns " << numNegSpike << endl; 
    baselineLow  = lo;
    baselineHigh = hi;

    return minY;
}
/*double Trace::DoBaseline(unsigned int lo, unsigned int numBins)
{
    if (size() < lo + numBins) {
        cerr << "Bad range in baseline calculation." << endl;
        return NAN;
    }

    unsigned int hi = lo + numBins;

    if (baselineLow == lo && baselineHigh == hi)
        return GetValue("baseline");

    double sum = accumulate(begin() + lo, begin() + hi, 0.0);
    double mean = sum / numBins;
    double sq_sum = inner_product(begin() + lo, begin() + hi,
                                  begin() + lo, 0.0);
    double std_dev = sqrt(sq_sum / numBins - mean * mean);

    SetValue("baseline", mean);
    SetValue("sigmaBaseline", std_dev);

    baselineLow  = lo;
    baselineHigh = hi;
    //cout << mean << " " << std_dev << endl;
    return mean;
}*/
/* //backup notes tb deleted. 
    double sum_x = (lo+hi)*(hi-lo+1.)/2.;
    double sq_sum_x = (hi-lo+1.)*(2.*(pow(lo,2.) + pow(hi,2.) + lo*hi) -lo + hi)/6.;
    double sqSumX = inner_product(xVector.begin(), xVector.end(),
                                  xVector.begin(), 0.0);
    
    double sqSumY = inner_product(begin() + lo, begin() + hi,
                                  begin() + lo, 0.0);
    double sumXY = inner_product(begin() + lo, begin() + hi,
                                  xVector.begin(), 0.0);
 
    double denom = numBins*sqSumX - pow(sumX,2.);
    double n1 = numBins*sumXY;
    double n2 = sumX*sumY;
    double slope  =   abs(n1/denom - n2/denom);
    double offset =  ( sqSumX*sumY - sumX*sumXY ) / ( numBins*sqSumX - pow(sumX,2.) );
       
    double stdDev = sqrt( ( sqSumY +pow(slope,2.)*sqSumX+numBins*pow(offset,2.)+
                     2*slope*offset*sumX - 2*slope*sumXY-2*offset*sumY ) / ( numBins-2 ) );*/
unsigned int Trace::DoDiscrimination(unsigned int lo, unsigned int numBins)
{
    unsigned int high = lo+numBins;

    if(size() < high)
        return pixie::U_DELIMITER;
    
    int discrim = 0, max = GetValue("maxpos");
    double baseline = GetValue("baseline");

    for(unsigned int i = max+lo; i <= max+high; i++)
	discrim += at(i)-baseline;
    
    InsertValue("discrim", discrim);
    
    return(discrim);
}

unsigned int Trace::DoQDC(unsigned int lo, unsigned int numBins) 
{
    unsigned int high = lo+numBins;

    if(size() < high)
	return pixie::U_DELIMITER;

    double baseline = GetValue("baseline");
    double qdc = 0;

    for(unsigned int i = lo; i < high; i++)
	qdc += at(i)-baseline;

    InsertValue("tqdc", qdc);

    return(qdc);
}

unsigned int Trace::FindMaxInfo(unsigned int lo, unsigned int numBins)
{
    lo = constants.GetConstant("waveformLow");
    unsigned int hi = constants.GetConstant("waveformHigh");
    numBins = lo + hi;
    
    if(size() < lo + numBins)
       return pixie::U_DELIMITER;
    
    Trace::const_iterator itTrace = max_element(begin()+lo, end()-lo);
    
    int maxPos = int(itTrace-begin());

    DoBaseline(0,maxPos-constants.GetConstant("waveformLow"));

    InsertValue("maxpos", int(itTrace-begin()));
    InsertValue("maxval", int(*itTrace)-GetValue("baseline"));

    return (itTrace-begin());
}

void Trace::Plot(int id)
{
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, 1, at(i));
    }
}

void Trace::Plot(int id, int row)
{
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, row, at(i));
    }
}

void Trace::ScalePlot(int id, double scale)
{
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, 1, abs(at(i)) / scale);
    }
}

void Trace::ScalePlot(int id, int row, double scale)
{
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, row, abs(at(i)) / scale);
    }
}

void Trace::OffsetPlot(int id, double offset)
{
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, 1, max(0., at(i) - offset));
    }
}

void Trace::OffsetPlot(int id, int row, double offset)
{
    for (size_type i=0; i < size(); i++) {
        histo.Plot(id, i, row, max(0., at(i) - offset));
    }
}
