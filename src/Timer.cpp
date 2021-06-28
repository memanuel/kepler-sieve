/** @file Timer.cpp
 *  @brief Implmentation of Timer class.
 *
 *  @author Michael S. Emanuel
 *  @date 2019-02-04
 */

// *********************************************************************************************************************
#include "Timer.hpp"

// *********************************************************************************************************************
// Constants used in converting units
constexpr long aBillion {static_cast<long>(1.0e9)};
constexpr long aMillion {static_cast<long>(1.0e6)};
	
// *********************************************************************************************************************
// Constructor
Timer::Timer() {
    tick();
}

// *********************************************************************************************************************
// Start the timer.
void Timer::tick() {
    tp0 = high_resolution_clock::now();
}

// *********************************************************************************************************************
double Timer::tock() {
	// Time point when tock() is called
	highResTimePoint tp1 = high_resolution_clock::now();

	// The elapsed time in nanoseconds
	time_unit_t t = duration_cast<nanoseconds>(tp1 - tp0).count();

	// Compute the elapsed time in seconds.
	double tSeconds = static_cast<double>(t) / aBillion;

	// Return the elapsed time in seconds
	return tSeconds;
}

// *********************************************************************************************************************
double Timer::tock_msg(string blurb, ostream &os) {
	// Time point when tock() is called
	highResTimePoint tp1 = high_resolution_clock::now();

	// The elapsed time in nanoseconds
	time_unit_t t = duration_cast<nanoseconds>(tp1 - tp0).count();

	// Create a formatted output with an appropriate amount of resolution for legibility.
	// This message is either:
	// Elapsed time <blurb>: nn.ddd <TimeUnits>.
	// Elapsed time: nn.ddd <TimeUnits>.
	// TimeUnits is one of seconds, milliseoncds, microseconds, or nanoseconds.
	format msg("Elapsed time: %.3f %s.\n");
	if (blurb.length() > 0) {
		msg = format("Elapsed time %s: %.3f %s.\n");
		msg % blurb;
	}

	// Compute the elapsed time in seconds.
	double tSeconds = static_cast<double>(t) / aBillion;

	// Send a message to cout stating the elapsed time.
	if (t > aBillion) {
		cout << msg % tSeconds % "seconds";
	}
	else if (t > aMillion) {
		double tMilliSeconds = tSeconds * 1000;
		cout << msg % tMilliSeconds % "milliseconds";
	}
	else if (t > 1000) {
		double tMicroSeconds = tSeconds * aMillion;
		cout << msg % tMicroSeconds% "microseconds";
	}
	else {
		double tNanoSeconds = static_cast<double>(t);
		cout << msg % tNanoSeconds % "microseconds";
	}

	// Return the elapsed time in seconds
	return tSeconds;
}

// *********************************************************************************************************************
void TimerTest() {
	Timer t;

	t.tick();
	t.tock_msg("doing nothing");
}
