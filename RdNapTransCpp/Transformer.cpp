// ***********************************************************************
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-06-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-07-2019
// ***********************************************************************
// C++ PORT from C version of RDNAPTRANS
// ***********************************************************************
#include "Transformer.h"
#include "Constants.h"
#include "Helpers.h"
#include "GrdFile.h"

const std::string GRID_FILE_GEOID = "nlgeo04.grd";

/*
**--------------------------------------------------------------
**    Function name: etrs2rd
**    Description:   convert ETRS89 coordinates to RD coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    etrs           Geographic  in     req     none
**    -              Cartesian   out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi_etrs, lambda_etrs, h_etrs  input ETRS89 coordinates
**    x_rd, y_rd                     output RD coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
Cartesian Transformer::etrs2rd(const Geographic& etrs)
{
	double x_rd;
	double y_rd; 
	double h_bessel;
	etrs2rd(etrs.get_phi(), etrs.get_lambda(), etrs.get_h(), x_rd, y_rd, h_bessel);
	return Cartesian(x_rd, y_rd, h_bessel);
}

/*
**--------------------------------------------------------------
**    Function name: rd2etrs
**    Description:   convert RD coordinates to ETRS89 coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    rd             Cartesian    in     req     none
**    -              Geographic   out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_rd, y_rd                    input RD coordinates
**    phi_etrs, lambda_etrs, h_etrs output ETRS89 coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
Geographic Transformer::rd2etrs(const Cartesian& rd)
{
	double phi_etrs;
	double lambda_etrs;
	double h_etrs;
	rd2etrs(rd.get_x(), rd.get_y(), rd.get_z(), phi_etrs, lambda_etrs, h_etrs);
	return Geographic(phi_etrs, lambda_etrs, h_etrs);
}

/*
**--------------------------------------------------------------
**    Function name: etrs2rdnap
**    Description:   convert ETRS89 coordinates to RD coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    etrs           Geographic  in     req     none
**    -              Cartesian   out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi_etrs, lambda_etrs, h_etrs  input ETRS89 coordinates
**    x_rd, y_rd, z_rd               output RD coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
Cartesian Transformer::etrs2rdnap(const Geographic& etrs)
{
	double x_rd;
	double y_rd;
	double h_bessel;
	etrs2rdnap(etrs.get_phi(), etrs.get_lambda(), etrs.get_h(), x_rd, y_rd, h_bessel);
	return Cartesian(x_rd, y_rd, h_bessel);
}

/*
**--------------------------------------------------------------
**    Function name: rdnap2etrs
**    Description:   convert RD coordinates to ETRS89 coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    rd             Cartesian    in     req     none
**    -              Geographic   out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_rd, y_rd, z_rd                    input RD coordinates
**    phi_etrs, lambda_etrs, h_etrs output ETRS89 coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
Geographic Transformer::rdnap2etrs(const Cartesian& rd)
{
	double phi_etrs;
	double lambda_etrs;
	double h_etrs;
	rdnap2etrs(rd.get_x(), rd.get_y(), rd.get_z(), phi_etrs, lambda_etrs, h_etrs);
	return Geographic(phi_etrs, lambda_etrs, h_etrs);
}

/*
**--------------------------------------------------------------
**    Function name: etrs2rd
**    Description:   convert ETRS89 coordinates to RD coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    phi_etrs       double      in     req     none
**    lambda_etrs    double      in     req     none
**    h_etrs         double      in     req     none
**    x_rd           double      out    -       none
**    y_rd           double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi_etrs, lambda_etrs, h_etrs  input ETRS89 coordinates
**    x_rd, y_rd                     output RD coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
int Transformer::etrs2rd(double phi_etrs, double lambda_etrs, double h_etrs,
	double& x_rd, double& y_rd, double& h_bessel)
{
	double x_etrs, y_etrs, z_etrs;
	double x_bessel, y_bessel, z_bessel;
	double phi_bessel, lambda_bessel;
	double x_pseudo_rd, y_pseudo_rd;

	double x_amersfoort_bessel;
	double y_amersfoort_bessel;
	double z_amersfoort_bessel;

	/*
	**--------------------------------------------------------------
	**    Calculate the cartesian ETRS89 coordinates of the pivot point Amersfoort
	**--------------------------------------------------------------
	*/
	Helpers::geographic2cartesian(Constants::PHI_AMERSFOORT_BESSEL, Constants::LAMBDA_AMERSFOORT_BESSEL, Constants::H_AMERSFOORT_BESSEL,
	                     Constants::A_BESSEL, Constants::INV_F_BESSEL,
		x_amersfoort_bessel, y_amersfoort_bessel, z_amersfoort_bessel);
	const auto x_amersfoort_etrs = x_amersfoort_bessel + Constants::TX_BESSEL_ETRS;
	const auto y_amersfoort_etrs = y_amersfoort_bessel + Constants::TY_BESSEL_ETRS;
	const auto z_amersfoort_etrs = z_amersfoort_bessel + Constants::TZ_BESSEL_ETRS;

	/*
	**--------------------------------------------------------------
	**    Convert ETRS89 coordinates to RD coordinates
	**    (To convert from degrees, minutes and seconds use the function deg_min_sec2decimal() here)
	**--------------------------------------------------------------
	*/
	Helpers::geographic2cartesian(phi_etrs, lambda_etrs, h_etrs,
	                              Constants::A_ETRS, Constants::INV_F_ETRS,
		x_etrs, y_etrs, z_etrs);
	Helpers::sim_trans(x_etrs, y_etrs, z_etrs,
	          Constants::TX_ETRS_BESSEL, Constants::TY_ETRS_BESSEL, Constants::TZ_ETRS_BESSEL,
	          Constants::ALPHA_ETRS_BESSEL, Constants::BETA_ETRS_BESSEL, Constants::GAMMA_ETRS_BESSEL,
	          Constants::DELTA_ETRS_BESSEL,
		x_amersfoort_etrs, y_amersfoort_etrs, z_amersfoort_etrs,
		x_bessel, y_bessel, z_bessel);
	Helpers::cartesian2geographic(x_bessel, y_bessel, z_bessel,
	                     Constants::A_BESSEL, Constants::INV_F_BESSEL,
		phi_bessel, lambda_bessel, h_bessel);
	Helpers::rd_projection(phi_bessel, lambda_bessel, x_pseudo_rd, y_pseudo_rd);
	const auto error = Helpers::rd_correction(x_pseudo_rd, y_pseudo_rd, x_rd, y_rd);
	return error;
}

/*
**--------------------------------------------------------------
**    Function name: rd2etrs
**    Description:   convert RD coordinates to ETRS89 coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    x_rd           double      in     req     none
**    y_rd           double      in     req     none
**    nap            double      in     req     none
**    phi_etrs       double      out    -       none
**    lambda_etrs    double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_rd, y_rd, nap        input RD and NAP coordinates
**    phi_etrs, lambda_etrs  output ETRS89 coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
int Transformer::rd2etrs(double x_rd, double y_rd, double nap,
	double& phi_etrs, double& lambda_etrs, double& h_etrs)
{
	double x_pseudo_rd, y_pseudo_rd;
	double phi_bessel, lambda_bessel;
	double x_bessel, y_bessel, z_bessel;
	double x_etrs, y_etrs, z_etrs;
	double x_amersfoort_bessel;
	double y_amersfoort_bessel;
	double z_amersfoort_bessel;

	/*
	**--------------------------------------------------------------
	**    Calculate the cartesian Bessel coordinates of the pivot point Amersfoort
	**--------------------------------------------------------------
	*/
	Helpers::geographic2cartesian(Constants::PHI_AMERSFOORT_BESSEL, Constants::LAMBDA_AMERSFOORT_BESSEL, Constants::H_AMERSFOORT_BESSEL,
	                     Constants::A_BESSEL, Constants::INV_F_BESSEL,
		x_amersfoort_bessel, y_amersfoort_bessel, z_amersfoort_bessel);

	/*
	**--------------------------------------------------------------
	**    Calculate approximated value of ellipsoidal Bessel height
	**    The error made by using a constant for de Bessel geoid height is max. circa 1 meter in the 
	**    ellipsoidal height (for the NLGEO2004 geoid model). This introduces an error in the phi, lambda position too,
	**    this error is nevertheless certainly smaller than 0.0001 m.
	**--------------------------------------------------------------
	*/
	const auto h_bessel = nap + Constants::MEAN_GEOID_HEIGHT_BESSEL;

	/*
	**--------------------------------------------------------------
	**    Convert RD coordinates to ETRS89 coordinates
	**--------------------------------------------------------------
	*/
	const auto error = Helpers::inv_rd_correction(x_rd, y_rd,
	                                       x_pseudo_rd, y_pseudo_rd);
	Helpers::inv_rd_projection(x_pseudo_rd, y_pseudo_rd,
		phi_bessel, lambda_bessel);
	Helpers::geographic2cartesian(phi_bessel, lambda_bessel, h_bessel,
	                              Constants::A_BESSEL, Constants::INV_F_BESSEL,
		x_bessel, y_bessel, z_bessel);
	Helpers::sim_trans(x_bessel, y_bessel, z_bessel,
	          Constants::TX_BESSEL_ETRS, Constants::TY_BESSEL_ETRS, Constants::TZ_BESSEL_ETRS,
	          Constants::ALPHA_BESSEL_ETRS, Constants::BETA_BESSEL_ETRS, Constants::GAMMA_BESSEL_ETRS,
	          Constants::DELTA_BESSEL_ETRS,
		x_amersfoort_bessel, y_amersfoort_bessel, z_amersfoort_bessel,
		x_etrs, y_etrs, z_etrs);
	Helpers::cartesian2geographic(x_etrs, y_etrs, z_etrs,
	                     Constants::A_ETRS, Constants::INV_F_ETRS,
		phi_etrs, lambda_etrs, h_etrs);
	/*
	**--------------------------------------------------------------
	**    To convert to degrees, minutes and seconds use the function decimal2deg_min_sec() here
	**--------------------------------------------------------------
	*/
	return error;
}

/*
**--------------------------------------------------------------
**    Function name: etrs2nap
**    Description:   convert ellipsoidal ETRS89 height to NAP height
**
**    Parameter      Type        In/Out Req/Opt Default
**    phi            double      in     req     none
**    lambda         double      in     req     none
**    h              double      in     req     none
**    nap            double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi, lambda, h  input ETRS89 coordinates
**    nap             output NAP height
**
**    Return value: (besides the standard return values) none
**    on error (outside geoid grid) nap is not compted here
**    instead in etrs2rdnap nap=h_bessel
**--------------------------------------------------------------
*/
int Transformer::etrs2nap(double phi, double lambda, double h,
	double& nap)
{
	double n;
	const auto error = GrdFile::grid_interpolation(lambda, phi, GRID_FILE_GEOID, n);
	if (error != 0)
	{
		return error;
	}
	nap = h - n + 0.0088;
	return 0;
}

/*
**--------------------------------------------------------------
**    Function name: nap2etrs
**    Description:   convert NAP height to ellipsoidal ETRS89 height
**
**    Parameter      Type        In/Out Req/Opt Default
**    phi            double      in     req     none
**    lambda         double      in     req     none
**    nap            double      in     req     none
**    h              double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi, lambda  input ETRS89 position
**    nap          input NAP height at position phi, lambda
**    h            output ellipsoidal ETRS89 height
**
**    Return value: (besides the standard return values)
**    none
**    on error (outside geoid grid) h is not compted here
**    instead in rdnap2etrs h=h_etrs_sim (from similarity transformation)
**--------------------------------------------------------------
*/
int Transformer::nap2etrs(double phi, double lambda, double nap,
	double& h)
{
	double n;
	const auto error = GrdFile::grid_interpolation(lambda, phi, GRID_FILE_GEOID, n);
	if (error != 0)
	{
		return error;
	}
	h = nap + n - 0.0088;
	return 0;
}

/*
**--------------------------------------------------------------
**    Function name: etrs2rdnap
**    Description:   convert ETRS89 coordinates to RD and NAP coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    phi            double      in     req     none
**    lambda         double      in     req     none
**    h              double      in     req     none
**    x_rd           double      out    -       none
**    y_rd           double      out    -       none
**    nap            double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi, lambda, h   input ETRS89 coordinates
**    x_rd, y_rd, nap  output RD and NAP coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
int Transformer::etrs2rdnap(double phi, double lambda, double h,
	double& x_rd, double& y_rd, double& nap)
{
	double h_bessel, h_geoid;
	auto error = etrs2rd(phi, lambda, h, x_rd, y_rd, h_bessel);
	if (error != 0)
	{
		return error;
	}
	error = etrs2nap(phi, lambda, h, h_geoid);
	if (error == 3)
	{
		nap = h_bessel;
	}
	else nap = h_geoid;
	return error;
}

/*
**--------------------------------------------------------------
**    Function name: rdnap2etrs
**    Description:   convert RD and NAP coordinates to ETRS89 coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    x_rd           double      in     req     none
**    y_rd           double      in     req     none
**    nap            double      in     req     none
**    phi            double      out    -       none
**    lambda         double      out    -       none
**    h              double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_rd, y_rd, nap  input RD and NAP coordinates
**    phi, lambda, h   output ETRS89 coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
int Transformer::rdnap2etrs(double x_rd, double y_rd, double nap,
	double& phi, double& lambda, double& h)
{
	double h_etrs_sim, h_etrs_geoid;
	auto error = rd2etrs(x_rd, y_rd, nap, phi, lambda, h_etrs_sim);
	if (error != 0)
	{
		return error;
	}
	error = nap2etrs(phi, lambda, nap, h_etrs_geoid);
	if (error == 3)
	{
		h = h_etrs_sim;
	}
	else
	{
		h = h_etrs_geoid;
	}
	return error;
}