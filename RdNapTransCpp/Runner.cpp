// ***********************************************************************
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-06-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-07-2019
// ***********************************************************************
// C++ PORT from C version of RDNAPTRANS
// ***********************************************************************
/*
**--------------------------------------------------------------
**    Function name: main
**    Description:   very simple interface
**
**    Parameter      Type        In/Out Req/Opt Default
**    none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
#include <iostream>
#include "Transformer.h"
#include <iomanip>
using namespace  std;

int main()
{
	auto choice = 9;
	double x_rd, y_rd, nap, phi, lambda, h;

	cout << endl;
	cout << "RDNAPTRANS(TM)2008" << endl;

	while (choice != 0)
	{
		cout << endl;
		cout << "Make your choice:" << endl;
		cout << "[1] ETRS89 to RD and NAP" << endl;
		cout << "[2] RD and NAP to ETRS89" << endl;
		cout << "[0] Close program" << endl << " ";
		cin >> choice;
		if (choice == 1)
		{
			cout << endl;
			cout << "Enter latitude (degrees): ";
			cin >> phi;
			cout << "Enter longitude (degrees): ";
			cin >> lambda;
			cout << "Enter ellipsoidal height (meters, enter 43 if unknown): ";
			cin >> h;
			cout << endl;

			/*
			**--------------------------------------------------------------
			**    Calculation ETRS89 to RD and NAP
			**--------------------------------------------------------------
			*/
			Transformer::etrs2rdnap(phi, lambda, h, x_rd, y_rd, nap);
			{
				std::cout << fixed << setprecision(4);
				cout << "RD x = " << x_rd << endl;
				cout << "   y = " << y_rd << endl;
				cout << "NAP  = " << nap << endl;
			}
		}
		else if (choice == 2)
		{
			cout << endl;
			cout << "Enter RD x (meters): ";
			cin >> x_rd;
			cout << "Enter RD y (meters): ";
			cin >> y_rd;
			cout << "Enter NAP  (meters): ";
			cin >> nap;
			cout << endl;

			/*
			**--------------------------------------------------------------
			**    Calculation RD and NAP to ETRS89
			**--------------------------------------------------------------
			*/
			Transformer::rdnap2etrs(x_rd, y_rd, nap, phi, lambda, h);
			{
				cout << fixed << setprecision(9);
				cout << "phi                 " << phi << endl;
				cout << "lambda             = " << lambda << endl;
				cout << fixed << setprecision(4);
				cout << "ellipsoidal height = " << h << endl;
			}
		}
		else if (choice == 0)
		{
			cout << endl;
			cout << "Regular end of program" << endl;
			cout << endl;
		}
	}
}
