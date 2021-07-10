/** @file Timer.cpp
 *  @brief Implmentation of Timer class.
 *
 *  @author Michael S. Emanuel
 *  @date 2019-02-04
 */

// *********************************************************************************************************************
#include "Timer.hpp"
	using ks::Timer;

// *********************************************************************************************************************
// Constants used in converting units
constexpr long aBillion {static_cast<long>(1.0e9)};
constexpr long aMillion {static_cast<long>(1.0e6)};
	
// *********************************************************************************************************************
// Constructor
Timer::Timer() 
	{tick();}

// *********************************************************************************************************************
// Start the timer.
void Timer::tick() 
	{tp0 = high_resolution_clock::now();}

// *********************************************************************************************************************
double Timer::tock() 
{
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
double Timer::tock_msg(const string blurb) 
{
	// Time point when tock() is called
	highResTimePoint tp1 = high_resolution_clock::now();

	// The elapsed time in nanoseconds
	time_unit_t t = duration_cast<nanoseconds>(tp1 - tp0).count();

	// Create a formatted output with an appropriate amount of resolution for legibility.
	// This message is either:
	// Elapsed time <blurb>: nn.ddd <TimeUnits>.
	// Elapsed time: nn.ddd <TimeUnits>.
	// TimeUnits is one of seconds, milliseoncds, microseconds, or nanoseconds.

	// The template message depends on whether a blurb was provided or not
	string msg = (blurb.length() > 0) ? 
		"Elapsed time {:s}: {:.3f} {:s}.\n" : 
		"Elapsed time: {:.3f} {:s}.\n";

	// Compute the elapsed time in seconds.
	double tSeconds = static_cast<double>(t) / aBillion;

	// Print a message stating the elapsed time.
	if (t > aBillion) 
	{
		print(msg, tSeconds, "seconds");
	}
	else if (t > aMillion) 
	{
		double tMilliSeconds = tSeconds * 1000;
		print(msg, tMilliSeconds, "milliseconds");
	}
	else if (t > 1000) 
	{
		double tMicroSeconds = tSeconds * aMillion;
		print(msg, tMicroSeconds, "microseconds");
	}
	else 
	{
		double tNanoSeconds = static_cast<double>(t);
		print(msg, tNanoSeconds, "nanoseconds");
	}

	// Return the elapsed time in seconds
	return tSeconds;
}

// *********************************************************************************************************************
void TimerTest() 
{
	Timer t;

	t.tick();
	t.tock_msg("doing nothing");
}
