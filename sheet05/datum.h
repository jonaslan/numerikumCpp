#ifndef _DATUM_H_
#define _DATUM_H_

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>

// Copy of datum.h with modified constructor

// time(0) gets current system time
static time_t t = time(0);
// the time struct contains (among other data members):
// current day of the month in tm_mday (1-31)
// current month since January in tm_mon (0-11)
// current year since 1900 in tm_year
static struct tm* curr_date = localtime(&t);

class Datum
{
	int day, month, year;
	
public:
	
	// constructor now has all default parameters, through the struct time the default 
	// parameters are those of the current date
	Datum(const int& dayIn=curr_date->tm_mday, const int& monthIn=curr_date->tm_mon + 1, const int& yearIn=curr_date->tm_year+1900)
	{
		day = dayIn;
		month = monthIn;
		year = yearIn;
	}
  
	friend std::ostream& operator<< (std::ostream& out, const Datum& d)
	{
		// Print the American date representation mm/dd/yyyy with leading zeros
		out << std::setfill('0') << std::setw(2) << d.month << "/";
		out << std::setfill('0') << std::setw(2) << d.day << "/";
		out << std::setfill('0') << std::setw(2) << d.year;
		return out;
	}


	const std::string& weekday()
	{
		if (year<= 1582)
			throw ("Input parameter 'year' must be greater than 1582!");

		const int m = month<3 ? month+10 : month-2;

		int c = year/100;
		if (!(c%4) && m>10)
			--c;

		int y = year-100*c;
		if (m>10)
			--y;

		const int d = day;
		const int w = int(d + int(2.6*m-0.2) + y + int(y/4) + int(c/4) - 2*c ) % 7;

		static const std::string name[7] = {"Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"};
		return name[w];
	}

};

#endif