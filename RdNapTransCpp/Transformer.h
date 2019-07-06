#pragma once
#include "Cartesian.h"
#include "Geographic.h"

class Transformer
{
public:
	static int etrs2rd(double phi_etrs, double lambda_etrs, double h_etrs,
		double& x_rd, double& y_rd, double& h_bessel);
	static int rd2etrs(double x_rd, double y_rd, double nap,
		double& phi_etrs, double& lambda_etrs, double& h_etrs);
	static int etrs2nap(double phi, double lambda, double h,
		double& nap);
	static int nap2etrs(double phi, double lambda, double nap,
		double& h);
	static int etrs2rdnap(double phi, double lambda, double h,
		double& x_rd, double& y_rd, double& nap);
	static int rdnap2etrs(double x_rd, double y_rd, double nap,
		double& phi, double& lambda, double& h);

	static Cartesian etrs2rd(Geographic etrs);
	static Geographic rd2etrs(Cartesian rd);
	static Cartesian etrs2rdnap(Geographic etrs);
	static Geographic rdnap2etrs(Cartesian rd);
};