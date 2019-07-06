#pragma once
class Helpers
{
public:
	static double deg_sin(double alpha);
	static double deg_cos(double alpha);
	static double deg_tan(double alpha);
	static double deg_asin(double a);
	static double deg_atan(double a);
	static double atanh(double a);
	static void deg_min_sec2decimal(double deg, double min, double sec, double& dec_deg);
	static void decimal2deg_min_sec(double dec_deg, int& deg, int& min, double& sec);
	static void geographic2cartesian(double phi, double lambda, double h,
		double a, double inv_f,
		double& x, double& y, double& z);
	static void cartesian2geographic(double x, double y, double z,
		double a, double inv_f,
		double& phi, double& lambda, double& h);
	static void sim_trans(double x_in, double y_in, double z_in,
		double tx, double ty, double tz,
		double alpha, double beta, double gamma,
		double delta,
		double xa, double ya, double za,
		double& x_out, double& y_out, double& z_out);
	static void rd_projection(double phi, double lambda,
		double& x_rd, double& y_rd);
	static void inv_rd_projection(double x_rd, double y_rd,
		double& phi, double& lambda);
	static int rd_correction(double x_pseudo_rd, double y_pseudo_rd,
		double& x_rd, double& y_rd);
	static int inv_rd_correction(double x_rd, double y_rd,
		double& x_pseudo_rd, double& y_pseudo_rd);
};
