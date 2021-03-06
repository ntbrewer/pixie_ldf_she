#include <cmath>
//TEMP
#include <iostream>

#include "Calibrator.hpp"
#include "Exceptions.hpp"

using namespace std;


void Calibrator::AddChannel(const Identifier& chanID, const std::string model,
                            double min, double max, 
                            const std::vector<double>& par) {
    CalibrationParams cf;
    unsigned required_parameters = 0;
    if (model == "raw") {
        cf.model = cal_raw;
    } else if (model == "off") {
        cf.model = cal_off;
    } else if (model == "linear") {
        cf.model = cal_linear;
        required_parameters = 2;
    } else if (model == "quadratic") {
        cf.model = cal_quadratic;
        required_parameters = 3;
    } else if (model == "polynomial") {
        cf.model = cal_polynomial;
        required_parameters = 1;
    } else if (model == "hyplin") {
        cf.model = cal_hyplin;
        required_parameters = 3;
    } else if (model == "linlin") {
	cf.model = cal_linlin;
        required_parameters = 4;
    } else if (model == "linlog") {
	cf.model = cal_linlog;
        required_parameters = 4;
    } else {
        stringstream ss;
        ss << "Calibrator: unknown calibration model " << model;
        throw GeneralException(ss.str());
    }

    if (min > max) {
        stringstream ss;
        ss << "Calibrator: wrong calibration range, channels "
           << min << " to " << max;
        throw GeneralException(ss.str());
    }

    cf.min = min;
    cf.max = max;
    
    for (vector<double>::const_iterator it = par.begin(); it != par.end();
        ++it) {
        cf.parameters.push_back(*it);
    }

    if (cf.parameters.size() < required_parameters) {
        stringstream ss;
        ss << "Calibrator: selected model needs at least " 
           << required_parameters
           << " but only " << cf.parameters.size() << " where found";
        throw GeneralException(ss.str());
    }

    if (channels_.find(chanID) != channels_.end()) {
        channels_[chanID].push_back(cf);
    } else {
        vector<CalibrationParams> vcf;
        vcf.push_back(cf);
        channels_[chanID] = vcf;
    }
}

double Calibrator::GetCalEnergy(const Identifier& chanID, double raw) const {
    map<Identifier, vector<CalibrationParams> >::const_iterator itch =
        channels_.find(chanID);
    if (itch != channels_.end()) {
        vector<CalibrationParams>::const_iterator itf;
        for (itf = itch->second.begin(); itf != itch->second.end(); ++itf) {
            if (itf->min <= raw && raw <= itf->max)
                break;
        }
        // Parts of spectrum that are not within some min-max range are
        // zeroed
        if (itf == itch->second.end()) {
            return 0;
        }
        switch(itf->model) {
            case cal_raw:
                return ModelRaw(raw);
                break;
            case cal_off: 
                return ModelOff();
                break;
            case cal_linear: 
                return ModelLinear(itf->parameters, raw);
                break;
            case cal_quadratic: 
                return ModelQuadratic(itf->parameters, raw);
                break;
            case cal_polynomial: 
                return ModelPolynomial(itf->parameters, raw);
                break;
            case cal_hyplin: 
                return ModelHypLin(itf->parameters, raw);
                break;
	    case cal_linlin: 
                return ModelLinLin(itf->parameters, raw);
                break;
	    case cal_linlog: 
                return ModelLinLog(itf->parameters, raw);
                break;
            default: 
                break;
        }
    }

    // If no calibration found, return raw channel
    return raw;
}

double Calibrator::ModelRaw(double raw) const {
    return raw;
}

double Calibrator::ModelOff() const {
    return 0;
}

double Calibrator::ModelLinear(const std::vector<double>& par,
                                    double raw) const {
    return par[0] + par[1] * raw;
}

double Calibrator::ModelQuadratic(const std::vector<double>& par,
                                    double raw) const {
    return par[0] + par[1] * raw + par[2] * raw * raw;
}

double Calibrator::ModelPolynomial(const std::vector<double>& par,
                                    double raw) const {
    int p = 0;
    double r = 0;
    for (vector<double>::const_iterator it = par.begin(); it != par.end();
        ++it) {
        r += (*it) * pow(raw, p);
        ++p;
    }
    return r;
}

double Calibrator::ModelHypLin(const std::vector<double>& par,
                               double raw) const {
    if (raw > 0)
        return par[0] / raw + par[1] + par[2] * raw;
    else
        return 0;
}

double Calibrator::ModelLinLin(const std::vector<double>& par,
                                    double raw) const {
    if (raw*par[1] < 30000) {
    	return (par[0] + par[1] * raw);
    } else {
	return (par[2] + par[3] *raw);
    }
    
	
}

double Calibrator::ModelLinLog(const std::vector<double>& par,
                                    double raw) const {
   
//    cout << raw <<" raw" << endl;
    
//	return(raw);
    if (raw*par[0] < 30000) {
//	cout << par[0] <<" " << par[1] << " " << par[0]*raw+par[1] <<endl;
    	return (par[0] * raw + par[1]);
    } else {
//	cout << par[2] <<" " << par[3] << " " << exp(par[2] + par[3] *raw) <<endl;
	return (exp(par[2] + par[3] *raw));
    }


    
	
}
