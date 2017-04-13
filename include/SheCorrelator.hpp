/** \file SHECorrelator.hpp
 *
 */

#ifndef __SHECORRELATOR_HPP_
#define __SHECORRELATOR_HPP_

#include <vector>
#include <deque> 
#include <sstream>

#include "Plots.hpp"
#include "Globals.hpp"
#include "ChanEvent.hpp"
#include "Messenger.hpp"
#include "WalkCorrector.hpp"
#include "Calibrator.hpp"
#include "DammPlotIds.hpp"


enum SheEventType {
    alpha,
    heavyIon,
    fission,
    check,
    lightIon,
    unknown
};

class SheEvent {
    public:
        SheEvent();
        SheEvent(double energy, double ratio, double time, int mwpc, double mwpcTime, double mwpcEnergy, 
                 bool has_beam, bool has_veto, double has_escape, int pixel[3],
                 SheEventType type=unknown);

        ~SheEvent() {}

        double get_energy() const {
            return energy_;
        }
        double get_ratio() const {
            return ratio_;
        }
        double get_time() const {
            return time_;
        }
        int get_mwpc() const {
            return mwpc_;
        }
	double get_mwpcTime() const {
            return mwpcTime_;
        }
	double get_mwpcEnergy() const {
            return mwpcEnergy_;
        }
        bool get_beam() const {
            return has_beam_;
        }
        bool get_veto() const {
            return has_veto_;
        }
        double get_escape() const {
            return has_escape_;
        }
	int get_Xpixel() const {
            return Xpixel_;
        }
	int get_Ypixel() const {
            return Ypixel_;
        }
	int get_Epixel() const {
	    return Epixel_;
        }
        SheEventType get_type() const {
            return type_;
        }

        void set_energy(double energy) {
            energy_  = energy;
        }
        void set_ratio(double ratio) {
            ratio_  = ratio;
        }
        void set_time(double time) {
            time_ = time;
        }
        void set_mwpc(int mwpc) {
            mwpc_ = mwpc;
        }
        void set_mwpcTime(double mwpcTime) {
            mwpcTime_ = mwpcTime;
        }
        void set_mwpcEnergy(double mwpcEnergy) {
            mwpcEnergy_ = mwpcEnergy;
        }
        void set_beam(bool has_beam) {
            has_beam_ = has_beam;
        }
        void set_veto(bool has_veto) {
            has_veto_ = has_veto;
        }
        void set_escape(double has_escape) {
            has_escape_ = has_escape;
        }
	void set_Xpixel(int pixel[3]) {
            Xpixel_ = pixel[0];
        }
	void set_Ypixel(int pixel[3]) {
            Ypixel_ = pixel[1];
        }
	void set_Epixel(int pixel[3]) {
            Epixel_ = pixel[2];
        }
        void set_type(SheEventType type) {
            type_ = type;
        }

    private:
        /** Total (reconstructed) energy, may include escape **/
        double energy_;
        /** Ratio of the energy front/back reported if outside a delta_energy defined in the Config.xml **/
        double ratio_;
        /** Shortest time of all subevents (e.g. back and front) */
        double time_;
        /** Number of MWPC chambers hits */
        int mwpc_;
	/** Time MWPC of hits */
        double mwpcTime_;
	/** Energy MWPC of hits */
        double mwpcEnergy_;
        /** Veto hit flag **/
        bool has_veto_;
        /** Beam on flag **/
        bool has_beam_;
        /** If reconstructed energy includes escape **/
        double has_escape_;
        /** If chain has a multiplicity in the Dssd4SHEProcessor the other pixels are recorded**/
        int Xpixel_;
	int Ypixel_;
	int Epixel_;
        /** Type of event decided by Correlator **/
        SheEventType type_;
};


class SheCorrelator {
    public:
        SheCorrelator(int size_x, int size_y);
        ~SheCorrelator();
        bool add_event(SheEvent& event, int x, int y, Plots& histo);
        void human_event_info(SheEvent& event, std::stringstream& ss, double clockStart);
  

    private:
        int size_x_;
        int size_y_;
        std::deque<SheEvent>** pixels_;

        bool flush_chain(int x, int y, Plots& histo);
};


#endif
