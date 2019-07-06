#include <cmath>
#include "Helpers.h"
#include "Constants.h"
#include  "GrdFile.h"


const std::string GRID_FILE_DX = "x2c.grd";
const std::string GRID_FILE_DY = "y2c.grd";
const std::string GRID_FILE_GEOID = "nlgeo04.grd";
/*
**--------------------------------------------------------------
**    Functions
**--------------------------------------------------------------
*/

/*
**--------------------------------------------------------------
**    Function name: deg_sin
**    Description:   sine for angles in degrees
**
**    Parameter      Type        In/Out Req/Opt Default
**    alpha          double      in     req     none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    sin(alpha)
**--------------------------------------------------------------
*/
double Helpers::deg_sin(double alpha)
{
	return sin(alpha / 180.0 * Constants::PI);
}

/*
**--------------------------------------------------------------
**    Function name: deg_cos
**    Description:   cosine for angles in degrees
**
**    Parameter      Type        In/Out Req/Opt Default
**    alpha          double      in     req     none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    cos(alpha)
**--------------------------------------------------------------
*/
double Helpers::deg_cos(double alpha)
{
	return cos(alpha / 180.0 * Constants::PI);
}

/*
**--------------------------------------------------------------
**    Function name: deg_tan
**    Description:   tangent for angles in degrees
**
**    Parameter      Type        In/Out Req/Opt Default
**    alpha          double      in     req     none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    tan(alpha)
**--------------------------------------------------------------
*/
double Helpers::deg_tan(double alpha)
{
	return tan(alpha / 180.0 * Constants::PI);
}

/*
**--------------------------------------------------------------
**    Function name: deg_asin
**    Description:   inverse sine for angles in degrees
**
**    Parameter      Type        In/Out Req/Opt Default
**    a              double      in     req     none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    asin(a)
**--------------------------------------------------------------
*/
double Helpers::deg_asin(double a)
{
	return (asin(a) * 180.0 / Constants::PI);
}

/*
**--------------------------------------------------------------
**    Function name: deg_atan
**    Description:   inverse tangent for angles in degrees
**
**    Parameter      Type        In/Out Req/Opt Default
**    a              double in     req     none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    atan(a)
**--------------------------------------------------------------
*/
double Helpers::deg_atan(double a)
{
	return (atan(a) * 180.0 / Constants::PI);
}

/*
**--------------------------------------------------------------
**    Function name: atanh
**    Description:   inverse hyperbolic tangent
**
**    Parameter      Type        In/Out Req/Opt Default
**    a              double      in     req     none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    atanh(a)
**--------------------------------------------------------------
*/
double Helpers::atanh(double a)
{
	return (0.5 * log((1.0 + a) / (1.0 - a)));
}

/*
**--------------------------------------------------------------
**    Function name: deg_min_sec2decimal
**    Description:   converts from degrees, minutes and seconds to decimal degrees
**
**    Parameter      Type        In/Out Req/Opt Default
**    deg            double      in     req     none
**    min            double      in     req     none
**    sec            double      in     req     none
**    dec_deg        double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    All parameters are doubles, so one can also enter decimal minutes or degrees.
**    Note: Nonsense input is accepted too.
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void Helpers::deg_min_sec2decimal(double deg, double min, double sec, double& dec_deg)
{
	dec_deg = (deg + min / 60.0 + sec / 3600.0);
}

/*
**--------------------------------------------------------------
**    Function name: decimal2deg_min_sec
**    Description:   converts from decimal degrees to degrees, minutes and seconds
**
**    Parameter      Type        In/Out Req/Opt Default
**    dec_deg        double      in     req     none
**    deg            int         out    -       none
**    min            int         out    -       none
**    sec            double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    none
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void Helpers::decimal2deg_min_sec(double dec_deg, int& deg, int& min, double& sec)
{
	deg = int(dec_deg);
	min = int((dec_deg - deg) * 60.0);
	sec = ((dec_deg - deg) * 60.0 - min) * 60.0;
}

/*
**--------------------------------------------------------------
**    Function name: geographic2cartesian
**    Description:   from geographic coordinates to cartesian coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    phi            double      in     req     none
**    lambda         double      in     req     none
**    h              double      in     req     none
**    a              double      in     req     none
**    inv_f          double      in     req     none
**    x              double      out    -       none
**    y              double      out    -       none
**    z              double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi      latitude in degrees
**    lambda   longitude in degrees
**    h        ellipsoidal height
**    a        half major axis of the ellisoid
**    inv_f    inverse flattening of the ellipsoid
**    x, y, z  output of cartesian coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void Helpers::geographic2cartesian(double phi, double lambda, double h,
	double a, double inv_f,
	double& x, double& y, double& z)
{
	/*
	**--------------------------------------------------------------
	**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
	**--------------------------------------------------------------
	*/

	/*
	**--------------------------------------------------------------
	**    Explanation of the meaning of variables:
	**        f    flattening of the ellipsoid
	**        ee   first eccentricity squared (e squared in some notations)
	**        n    second (East West) principal radius of curvature (N in some notations)
	**--------------------------------------------------------------
	*/
	const auto f = 1.0 / inv_f;
	const auto ee = f * (2.0 - f);
	const auto n = a / sqrt(1.0 - ee * pow(deg_sin(phi), 2));
	x = (n + h) * deg_cos(phi) * deg_cos(lambda);
	y = (n + h) * deg_cos(phi) * deg_sin(lambda);
	z = (n * (1.0 - ee) + h) * deg_sin(phi);
}

/*
**--------------------------------------------------------------
**    Function name: cartesian2geographic
**    Description:   from cartesian coordinates to geographic coordinates
**
**    Parameter      Type        In/Out Req/Opt Default
**    x              double      in     req     none
**    y              double      in     req     none
**    z              double      in     req     none
**    a              double      in     req     none
**    inv_f          double      in     req     none
**    phi            double      out    -       none
**    lambda         double      out    -       none
**    h              double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x, y, z  input of cartesian coordinates
**    a        half major axis of the ellisoid
**    inv_f    inverse flattening of the ellipsoid
**    phi      output latitude in degrees
**    lambda   output longitude in degrees
**    h        output ellipsoidal height
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void Helpers::cartesian2geographic(double x, double y, double z,
	double a, double inv_f,
	double& phi, double& lambda, double& h)
{
	/*
	**--------------------------------------------------------------
	**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
	**--------------------------------------------------------------
	*/

	/*
	**--------------------------------------------------------------
	**    Explanation of the meaning of variables:
	**        f    flattening of the ellipsoid
	**        ee   first eccentricity squared (e squared in some notations)
	**        rho  distance to minor axis
	**        n    second (East West) principal radius of curvature (N in some notations)
	**--------------------------------------------------------------
	*/
	const double f = 1.0 / inv_f;
	const double ee = f * (2.0 - f);
	const double rho = sqrt(x * x + y * y);
	double n = 0;

	/*
	**--------------------------------------------------------------
	**    Iterative calculation of phi
	**--------------------------------------------------------------
	*/
	phi = 0;
	double diff = 90;
	while (diff > Constants::DEG_PRECISION)
	{
		const double previous = phi;
		n = a / sqrt(1.0 - ee * pow(deg_sin(phi), 2));
		phi = deg_atan(z / rho + n * ee * deg_sin(phi) / rho);
		diff = fabs(phi - previous);
	}

	/*
	**--------------------------------------------------------------
	**     Calculation of lambda and h
	**--------------------------------------------------------------
	*/
	lambda = deg_atan(y / x);
	h = rho * deg_cos(phi) + z * deg_sin(phi) - n * (1.0 - ee * pow(deg_sin(phi), 2));
}

/*
**--------------------------------------------------------------
**    Function name: sim_trans
**    Description:   3 dimensional similarity transformation (7 parameters) around another pivot point "a" than the origin
**
**    Parameter      Type        In/Out Req/Opt Default
**    x_in           double      in     req     none
**    y_in           double      in     req     none
**    z_in           double      in     req     none
**    tx             double      in     req     none
**    ty             double      in     req     none
**    tz             double      in     req     none
**    alpha          double      in     req     none
**    beta           double      in     req     none
**    gamma          double      in     req     none
**    delta          double      in     req     none
**    xa             double      in     req     none
**    ya             double      in     req     none
**    za             double      in     req     none
**    x_out          double      out    -       none
**    y_out          double      out    -       none
**    z_out          double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_in, y_in, z_in     input coordinates
**    tx                   translation in direction of x axis
**    ty                   translation in direction of y axis
**    tz                   translation in direction of z axis
**    alpha                rotation around x axis in radials
**    beta                 rotation around y axis in radials
**    gamma                rotation around z axis in radials
**    delta                scale parameter (scale = 1 + delta)
**    xa, ya, za           coordinates of pivot point a (in case of rotation around the center of the ellipsoid these parameters are zero)
**    x_out, y_out, z_out  output coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void Helpers::sim_trans(const double x_in, const double y_in, const double z_in,
                        const double tx, const double ty, const double tz,
                        const double alpha, const double beta, const double gamma,
                        const double delta,
                        const double xa, double ya, const double za,
	double& x_out, double& y_out, double& z_out)

{
	/*
	**--------------------------------------------------------------
	**    Source: HTW, "Handleiding voor de Technische Werkzaamheden van het Kadaster". Apeldoorn: Kadaster, 1996.
	**--------------------------------------------------------------
	*/

	/*
	**--------------------------------------------------------------
	**    Calculate the elements of the rotation_matrix:
	**
	**    a b c
	**    d e f
	**    g h i
	**
	**--------------------------------------------------------------
	*/
	const double a = cos(gamma) * cos(beta);
	const auto b = cos(gamma) * sin(beta) * sin(alpha) + sin(gamma) * cos(alpha);
	const auto c = -cos(gamma) * sin(beta) * cos(alpha) + sin(gamma) * sin(alpha);
	const auto d = -sin(gamma) * cos(beta);
	const auto e = -sin(gamma) * sin(beta) * sin(alpha) + cos(gamma) * cos(alpha);
	const auto f = sin(gamma) * sin(beta) * cos(alpha) + cos(gamma) * sin(alpha);
	const auto g = sin(beta);
	const auto h = -cos(beta) * sin(alpha);
	const auto i = cos(beta) * cos(alpha);

	/*
	**--------------------------------------------------------------
	**    Calculate the elements of the vector input_point:
	**    point_2 = input_point - pivot_point
	**--------------------------------------------------------------
	*/
	const auto x = x_in - xa;
	const auto y = y_in - ya;
	const auto z = z_in - za;

	/*
	**--------------------------------------------------------------
	**    Calculate the elements of the output vector:
	**    output_point = scale * rotation_matrix * point_2 + translation_vector + pivot_point
	**--------------------------------------------------------------
	*/
	x_out = (1.0 + delta) * (a * x + b * y + c * z) + tx + xa;
	y_out = (1.0 + delta) * (d * x + e * y + f * z) + ty + ya;
	z_out = (1.0 + delta) * (g * x + h * y + i * z) + tz + za;
}

/*
**--------------------------------------------------------------
**    Function name: rd_projection
**    Description:   stereographic double projection
**
**    Parameter      Type        In/Out Req/Opt Default
**    phi            double      in     req     none
**    lambda         double      in     req     none
**    x_rd           double      out    -       none
**    y_rd           double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    phi         input Bessel latitude in degrees
**    lambda      input Bessel longitude in degrees
**    x_rd, rd_y  output RD coordinates
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void Helpers::rd_projection(double phi, double lambda,
	double& x_rd, double& y_rd)
{
	/*
	**--------------------------------------------------------------
	**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
	**            G. Strang van Hees, "Globale en lokale geodetische systemen". Delft: Nederlandse Commissie voor Geodesie (NCG), 1997.
	**--------------------------------------------------------------
	*/

	/*
	**--------------------------------------------------------------
	**    Explanation of the meaning of constants:
	**        f                         flattening of the ellipsoid
	**        ee                        first eccentricity squared (e squared in some notations)
	**        e                         first eccentricity
	**        eea                       second eccentricity squared (e' squared in some notations)
	**
	**        phi_amersfoort_sphere     latitude of projection base point Amersfoort on sphere in degrees
	**        lambda_amersfoort_sphere  longitude of projection base point Amersfoort on sphere in degrees
	**
	**        r1                        first (North South) principal radius of curvature in Amersfoort (M in some notations)
	**        r2                        second (East West) principal radius of curvature in Amersfoort (N in some notations)
	**        r_sphere                  radius of sphere
	**
	**        n                         constant of Gaussian projection n = 1.000475...
	**        q_amersfoort              isometric latitude of Amersfoort on ellipsiod
	**        w_amersfoort              isometric latitude of Amersfoort on sphere
	**        m                         constant of Gaussian projection m = 0.003773... (also named c in some notations)
	**--------------------------------------------------------------
	*/
	const auto f = 1 / Constants::INV_F_BESSEL;
	const auto ee = f * (2 - f);
	const auto e = sqrt(ee);
	const auto eea = ee / (1.0 - ee);

	const auto phi_amersfoort_sphere = deg_atan(deg_tan(Constants::PHI_AMERSFOORT_BESSEL) / sqrt(1 + eea * pow(deg_cos(Constants::PHI_AMERSFOORT_BESSEL), 2)));
	const auto lambda_amersfoort_sphere = Constants::LAMBDA_AMERSFOORT_BESSEL;

	const auto r1 = Constants::A_BESSEL * (1 - ee) / pow(sqrt(1 - ee * pow(deg_sin(Constants::PHI_AMERSFOORT_BESSEL), 2)), 3);
	const auto r2 = Constants::A_BESSEL / sqrt(1.0 - ee * pow(deg_sin(Constants::PHI_AMERSFOORT_BESSEL), 2));
	const auto r_sphere = sqrt(r1 * r2);

	const auto n = sqrt(1 + eea * pow(deg_cos(Constants::PHI_AMERSFOORT_BESSEL), 4));
	const auto q_amersfoort = atanh(deg_sin(Constants::PHI_AMERSFOORT_BESSEL)) - e * atanh(e * deg_sin(Constants::PHI_AMERSFOORT_BESSEL));
	const auto w_amersfoort = log(deg_tan(45 + 0.5 * phi_amersfoort_sphere));
	const auto m = w_amersfoort - n * q_amersfoort;

	/*
	**--------------------------------------------------------------
	**    Explanation of the meaning of variables:
	**        q                    isometric latitude on ellipsiod
	**        w                    isometric latitude on sphere
	**        phi_sphere           latitide on sphere in degrees
	**        delta_lambda_sphere  difference in longitude on sphere with Amersfoort in degrees
	**        psi                  distance angle from Amersfoort on sphere
	**        alpha                azimuth from Amersfoort
	**        r                    distance from Amersfoort in projection plane
	**--------------------------------------------------------------
	*/
	const auto q = atanh(deg_sin(phi)) - e * atanh(e * deg_sin(phi));
	const auto w = n * q + m;
	const auto phi_sphere = 2 * deg_atan(exp(w)) - 90;
	const auto delta_lambda_sphere = n * (lambda - lambda_amersfoort_sphere);
	const auto sin_half_psi_squared = pow(deg_sin(0.5 * (phi_sphere - phi_amersfoort_sphere)), 2) + pow(deg_sin(0.5 * delta_lambda_sphere), 2) * deg_cos(phi_sphere) * deg_cos(phi_amersfoort_sphere);
	const auto sin_half_psi = sqrt(sin_half_psi_squared);
	const auto cos_half_psi = sqrt(1 - sin_half_psi_squared);
	const auto tan_half_psi = sin_half_psi / cos_half_psi;
	const auto sin_psi = 2 * sin_half_psi * cos_half_psi;
	const auto cos_psi = 1 - 2 * sin_half_psi_squared;
	const auto sin_alpha = deg_sin(delta_lambda_sphere) * (deg_cos(phi_sphere) / sin_psi);
	const auto cos_alpha = (deg_sin(phi_sphere) - deg_sin(phi_amersfoort_sphere) * cos_psi) / (deg_cos(phi_amersfoort_sphere) * sin_psi);
	const auto r = 2 * Constants::SCALE_RD * r_sphere * tan_half_psi;

	x_rd = r * sin_alpha + Constants::X_AMERSFOORT_RD;
	y_rd = r * cos_alpha + Constants::Y_AMERSFOORT_RD;
}

/*
**--------------------------------------------------------------
**    Function name: inv_rd_projection
**    Description:   inverse stereographic double projection
**
**    Parameter      Type        In/Out Req/Opt Default
**    x_rd           double      in     req     none
**    y_rd           double      in     req     none
**    phi            double      out    -       none
**    lambda         double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_rd, rd_y  input RD coordinates
**    phi         output Bessel latitude in degrees
**    lambda      output Bessel longitude in degrees
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void Helpers::inv_rd_projection(double x_rd, double y_rd,
	double& phi, double& lambda)
{
	/*
	**--------------------------------------------------------------
	**    Source: G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning at Sea". Delft University of Technology, 1995.
	**            G. Strang van Hees, "Globale en lokale geodetische systemen". Delft: Nederlandse Commissie voor Geodesie (NCG), 1997.
	**--------------------------------------------------------------
	*/

	/*
	**--------------------------------------------------------------
	**    Explanation of the meaning of constants:
	**        f                         flattening of the ellipsoid
	**        ee                        first eccentricity squared (e squared in some notations)
	**        e                         first eccentricity
	**        eea                       second eccentricity squared (e' squared in some notations)
	**
	**        phi_amersfoort_sphere     latitude of projection base point Amersfoort on sphere in degrees
	**
	**        r1                        first (North South) principal radius of curvature in Amersfoort (M in some notations)
	**        r2                        second (East West) principal radius of curvature in Amersfoort (N in some notations)
	**        r_sphere                  radius of sphere
	**
	**        n                         constant of Gaussian projection n = 1.000475...
	**        q_amersfoort              isometric latitude of Amersfoort on ellipsiod
	**        w_amersfoort              isometric latitude of Amersfoort on sphere
	**        m                         constant of Gaussian projection m = 0.003773... (also named c in some notations)
	**--------------------------------------------------------------
	*/
	const auto f = 1 / Constants::INV_F_BESSEL;
	const auto ee = f * (2 - f);
	const auto e = sqrt(ee);
	const auto eea = ee / (1.0 - ee);

	const auto phi_amersfoort_sphere = deg_atan(deg_tan(Constants::PHI_AMERSFOORT_BESSEL) / sqrt(1 + eea * pow(deg_cos(Constants::PHI_AMERSFOORT_BESSEL), 2)));

	const auto r1 = Constants::A_BESSEL * (1 - ee) / pow(sqrt(1 - ee * pow(deg_sin(Constants::PHI_AMERSFOORT_BESSEL), 2)), 3);
	const auto r2 = Constants::A_BESSEL / sqrt(1.0 - ee * pow(deg_sin(Constants::PHI_AMERSFOORT_BESSEL), 2));
	const auto r_sphere = sqrt(r1 * r2);

	const auto n = sqrt(1 + eea * pow(deg_cos(Constants::PHI_AMERSFOORT_BESSEL), 4));
	const auto q_amersfoort = atanh(deg_sin(Constants::PHI_AMERSFOORT_BESSEL)) - e * atanh(e * deg_sin(Constants::PHI_AMERSFOORT_BESSEL));
	const auto w_amersfoort = log(deg_tan(45 + 0.5 * phi_amersfoort_sphere));
	const auto m = w_amersfoort - n * q_amersfoort;

	/*
	**--------------------------------------------------------------
	**    Explanation of the meaning of variables:
	**        r                    distance from Amersfoort in projection plane
	**        alpha                azimuth from Amersfoort
	**        psi                  distance angle from Amersfoort on sphere in degrees
	**        phi_sphere           latitide on sphere in degrees
	**        delta_lambda_sphere  difference in longitude on sphere with Amersfoort in degrees
	**        w                    isometric latitude on sphere
	**        q                    isometric latitude on ellipsoid
	**--------------------------------------------------------------
	*/
	const auto r = sqrt(pow(x_rd - Constants::X_AMERSFOORT_RD, 2) + pow(y_rd - Constants::Y_AMERSFOORT_RD, 2));
	auto sin_alpha = (x_rd - Constants::X_AMERSFOORT_RD) / r;
	if (r < Constants::PRECISION)
	{
		sin_alpha = 0;
	}
	auto cos_alpha = (y_rd - Constants::Y_AMERSFOORT_RD) / r;
	if (r < Constants::PRECISION)
	{
		cos_alpha = 1;
	}
	const auto psi = 2 * deg_atan(r / (2 * Constants::SCALE_RD * r_sphere));
	const double phi_sphere = deg_asin(cos_alpha * deg_cos(phi_amersfoort_sphere) * deg_sin(psi) + deg_sin(phi_amersfoort_sphere) * deg_cos(psi));
	const auto delta_lambda_sphere = deg_asin((sin_alpha * deg_sin(psi)) / deg_cos(phi_sphere));

	lambda = delta_lambda_sphere / n + Constants::LAMBDA_AMERSFOORT_BESSEL;

	const double w = atanh(deg_sin(phi_sphere));
	const double q = (w - m) / n;

	/*
	**--------------------------------------------------------------
	**    Iterative calculation of phi
	**--------------------------------------------------------------
	*/
	phi = 0;
	double diff = 90;
	while (diff > Constants::DEG_PRECISION)
	{
		const auto previous = phi;
		phi = 2 * deg_atan(exp(q + 0.5 * e * log((1 + e * deg_sin(phi)) / (1 - e * deg_sin(phi))))) - 90;
		diff = fabs(phi - previous);
	}
}

/*
**--------------------------------------------------------------
**    Function name: rd_correction
**    Description:   apply the modeled distortions in the RD coordinate system
**
**    Parameter      Type        In/Out Req/Opt Default
**    x_pseudo_rd    double      in     req     none
**    y_pseudo_rd    double      in     req     none
**    x_rd           double      out    -       none
**    y_rd           double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_pseudo_rd, y_pseudo_rd  input coordinates in undistorted pseudo RD
**    x_rd, y_rd                output coordinates in real RD
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
int Helpers::rd_correction(double x_pseudo_rd, double y_pseudo_rd,
	double& x_rd, double& y_rd)
{
	double dx, dy;

	auto error = GrdFile::grid_interpolation(x_pseudo_rd, y_pseudo_rd, GRID_FILE_DX, dx);
	if (error != 0)
	{
		return error;
	}
	error = GrdFile::grid_interpolation(x_pseudo_rd, y_pseudo_rd, GRID_FILE_DY, dy);
	x_rd = x_pseudo_rd - dx;
	y_rd = y_pseudo_rd - dy;
	return error;
}

/*
**--------------------------------------------------------------
**    Function name: inv_rd_correction
**    Description:   remove the modeled distortions in the RD coordinate system
**
**    Parameter      Type        In/Out Req/Opt Default
**    x_rd           double      in     req     none
**    y_rd           double      in     req     none
**    x_pseudo_rd    double      out    -       none
**    x_pseudo_rd    double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x_rd, y_rd                input coordinates in real RD
**    x_pseudo_rd, y_pseudo_rd  output coordinates in undistorted pseudo RD
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
int Helpers::inv_rd_correction(double x_rd, double y_rd,
	double& x_pseudo_rd, double& y_pseudo_rd)
{
	double dx, dy;

	/*
	**--------------------------------------------------------------
	**    The grid values are formally in pseudo RD. For the interpolation below the RD values are used. The introduced error is certainly smaller than 0.0001 m for the X2c.grd and Y2c.grd.
	**--------------------------------------------------------------
	*/
	x_pseudo_rd = x_rd;
	y_pseudo_rd = y_rd;
	auto error = GrdFile::grid_interpolation(x_pseudo_rd, y_pseudo_rd, GRID_FILE_DX, dx);
	if (error != 0)
	{
		return error;
	}
	error = GrdFile::grid_interpolation(x_pseudo_rd, y_pseudo_rd, GRID_FILE_DY, dy);
	x_pseudo_rd = x_rd + dx;
	y_pseudo_rd = y_rd + dy;
	return error;
}

