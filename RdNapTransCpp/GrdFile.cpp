// ***********************************************************************
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-06-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-07-2019
// ***********************************************************************
// C++ PORT from C version of RDNAPTRANS
// ***********************************************************************
#include "GrdFile.h"
#include <cmath>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include "Constants.h"
using namespace std;

/*
**--------------------------------------------------------------
**    Function name: read_grd_file_header
**    Description:   reads the header of a grd file
**
**    Parameter      Type        In/Out Req/Opt Default
**    filename       string      in     req     none
**    size_x         short int   out    -       none
**    size_y         short int   out    -       none
**    min_x          double      out    -       none
**    max_x          double      out    -       none
**    min_y          double      out    -       none
**    max_y          double      out    -       none
**    min_value      double      out    -       none
**    max_value      double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    filename   name of the to be read binary file
**    size_x     number of grid values in x direction (row)
**    size_y     number of grid values in y direction (col)
**    min_x      minimum of x
**    max_x      maximum of x
**    min_y      minimum of y
**    max_y      maximum of x
**    min_value  minimum value in grid (besides the error values)
**    max_value  maximum value in grid (besides the error values)
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/

	int GrdFile::read_grd_file_header(const string& file_name,
		short int& size_x, short int& size_y,
		double& min_x, double& max_x,
		double& min_y, double& max_y,
		double& min_value, double& max_value)
{
	/*
	**--------------------------------------------------------------
	**    Grd files are binary grid files in the format of the program Surfer(R)
	**--------------------------------------------------------------
	*/

	fstream file(file_name.c_str(), ios::in | ios::binary);

	/*
	**--------------------------------------------------------------
	**    Read file id
	**--------------------------------------------------------------
	*/
	char id[5];
	for (auto i = 0; i < 4; i = i + 1)
	{
		file.seekg(i, ios::beg);
		file.read(static_cast<char*>(&id[i]), 1);
	}
	id[4] = '\0';
	const auto id_string = string(id);

	/*
	**--------------------------------------------------------------
	**    Checks
	**--------------------------------------------------------------
	*/
	if (!file)
	{
		cerr << file_name << " does not exist" << endl;
		return -1;
	}

	if (id_string != "DSBB")
	{
		cerr << file_name << " is not a valid grd file" << endl;
		return -1;
	}

	/*
	**--------------------------------------------------------------
	**    Read output parameters
	**--------------------------------------------------------------
	*/
	file.seekg(4, ios::beg);
	file.read(reinterpret_cast<char*>(&size_x), 2);

	file.seekg(6, ios::beg);
	file.read(reinterpret_cast<char*>(&size_y), 2);


	file.seekg(8, ios::beg);
	file.read(reinterpret_cast<char*>(&min_x), 8);

	file.seekg(16, ios::beg);
	file.read(reinterpret_cast<char*>(&max_x), 8);

	file.seekg(24, ios::beg);
	file.read(reinterpret_cast<char*>(&min_y), 8);

	file.seekg(32, ios::beg);
	file.read(reinterpret_cast<char*>(&max_y), 8);

	file.seekg(40, ios::beg);
	file.read(reinterpret_cast<char*>(&min_value), 8);

	file.seekg(48, ios::beg);
	file.read(reinterpret_cast<char*>(&max_value), 8);

	return 0;
}


/*
**--------------------------------------------------------------
**    Function name: read_grd_file_body
**    Description:   reads a value from a grd file
**
**    Parameter      Type        In/Out Req/Opt Default
**    filename       string      in     req     none
**    number         long int    in     req     none
**    value          float       out    -       none
**
**    Additional explanation of the meaning of parameters
**    filename       name of the grd file to be read
**    record_number  number defining the position in the file
**    record_value   output of the read value
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/

int GrdFile::read_grd_file_body(const string& file_name, long int record_number, float& record_value)
{
	const auto record_length = 4;
	const auto header_length = 56;

	/*
	**--------------------------------------------------------------
	**    Read
	**    Grd files are binary grid files in the format of the program Surfer(R)
	**    The first "header_length" bytes are the header of the file
	**    The body of the file consists of records of "record_length" bytes
	**    The records have a "record_number", starting with 0,1,2,...
	**--------------------------------------------------------------
	*/
	fstream file(file_name.c_str(), ios::in | ios::binary);
	file.seekg(record_length * record_number + header_length, ios::beg);
	file.read(reinterpret_cast<char*>(&record_value), record_length);

	/*
	**--------------------------------------------------------------
	**    Checks
	**--------------------------------------------------------------
	*/
	if (!file)
	{
		cerr << file_name << " does not exist" << endl;
		return -1;
	}
	return 0;
}


/*
**--------------------------------------------------------------
**    Function name: grid_interpolation
**    Description:   grid interpolation using Overhauser splines
**
**    Parameter      Type        In/Out Req/Opt Default
**    x              double      in     req     none
**    y              double      in     req     none
**    grd_file       string      in     req     none
**    value          double      out    -       none
**
**    Additional explanation of the meaning of parameters
**    x, y           coordinates of the point for which a interpolated value is desired
**    grd_file       name of the grd file to be read
**    record_value   output of the interpolated value
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/

int GrdFile::grid_interpolation(double x, double y, const string& grd_file, double& value)
{
	short int size_x, size_y;
	double min_x, max_x, min_y, max_y, min_value, max_value;
	long int record_number[16];
	float record_value[16];
	double f[4], g[4];
	double gfac[16];

	/*
	**--------------------------------------------------------------
	**    Explanation of the meaning of variables:
	**    size_x     number of grid values in x direction (row)
	**    size_y     number of grid values in y direction (col)
	**    min_x      minimum of x
	**    max_x      maximum of x
	**    min_y      minimum of y
	**    max_y      maximum of x
	**    min_value  minimum value in grid (besides the error values)
	**    max_value  maximum value in grid (besides the error values)
	**--------------------------------------------------------------
	*/
	auto error = read_grd_file_header(grd_file, size_x, size_y, min_x, max_x, min_y, max_y, min_value, max_value);
	if (error != 0)
	{
		return -1;
	}

	const auto step_size_x = (max_x - min_x) / (size_x - 1);
	const auto step_size_y = (max_y - min_y) / (size_y - 1);

	/*
	**--------------------------------------------------------------
	**    Check for location safely inside the bounding box of grid
	**--------------------------------------------------------------
	*/
	if (x <= (min_x + step_size_x) || x >= (max_x - step_size_x) ||
		y <= (min_y + step_size_y) || y >= (max_y - step_size_y))
	{
		cerr << "Outside bounding box of " << grd_file << endl;
		if (grd_file == "x2c.grd") { error = 1; value = 0.0; }
		if (grd_file == "y2c.grd") { error = 2; value = 0.0; }
		if (grd_file == "nlgeo04.grd") error = 3;
		return error;
	}

	/*
	**--------------------------------------------------------------
	**    The selected grid points are situated around point X like this:
	**
	**        12  13  14  15
	**
	**         8   9  10  11
	**               X
	**         4   5   6   7
	**
	**         0   1   2   3
	**
	**    ddx and ddy (in parts of the grid interval) are defined relative to grid point 9, respectively to the right and down.
	**--------------------------------------------------------------
	*/
	const auto ddx = (x - min_x) / step_size_x - floor((x - min_x) / step_size_x);
	const auto ddy = 1 - ((y - min_y) / step_size_y - floor((y - min_y) / step_size_y));

	/*
	**--------------------------------------------------------------
	**    Calculate the record numbers of the selected grid points
	**    The records are numbered from lower left corner to the upper right corner starting with 0:
	**
	**    size_x*(size_y-1) . . size_x*size_y-1
	**                   .                    .
	**                   .                    .
	**                   0 . . . . . . size_x-1
	**--------------------------------------------------------------
	*/
	record_number[5] = int((x - min_x) / step_size_x + floor((y - min_y) / step_size_y) * size_x);
	record_number[0] = record_number[5] - size_x - 1;
	record_number[1] = record_number[5] - size_x;
	record_number[2] = record_number[5] - size_x + 1;
	record_number[3] = record_number[5] - size_x + 2;
	record_number[4] = record_number[5] - 1;
	record_number[6] = record_number[5] + 1;
	record_number[7] = record_number[5] + 2;
	record_number[8] = record_number[5] + size_x - 1;
	record_number[9] = record_number[5] + size_x;
	record_number[10] = record_number[5] + size_x + 1;
	record_number[11] = record_number[5] + size_x + 2;
	record_number[12] = record_number[5] + 2 * size_x - 1;
	record_number[13] = record_number[5] + 2 * size_x;
	record_number[14] = record_number[5] + 2 * size_x + 1;
	record_number[15] = record_number[5] + 2 * size_x + 2;

	/*
	**--------------------------------------------------------------
	**    Read the record values of the selected grid point
	**    Outside the validity area the records have a very large value (circa 1.7e38).
	**--------------------------------------------------------------
	*/
	for (auto i = 0; i < 16; i++)
	{
		error = read_grd_file_body(grd_file, record_number[i], record_value[i]);
		if (error != 0)
		{
			return -1;
		}
		if (record_value[i] > max_value + Constants::PRECISION || record_value[i] < min_value - Constants::PRECISION)
		{
			cerr << "Outside validity area of " << grd_file << endl;
			if (grd_file == "x2c.grd") { error = 1; value = 0.0; }
			if (grd_file == "y2c.grd") { error = 2; value = 0.0; }
			if (grd_file == "nlgeo04.grd") error = 3;
			return error;
		}
	}

	/*
	**--------------------------------------------------------------
	**    Calculation of the multiplication factors
	**--------------------------------------------------------------
	*/
	f[0] = -0.5 * ddx + ddx * ddx - 0.5 * ddx * ddx * ddx;
	f[1] = 1.0 - 2.5 * ddx * ddx + 1.5 * ddx * ddx * ddx;
	f[2] = 0.5 * ddx + 2.0 * ddx * ddx - 1.5 * ddx * ddx * ddx;
	f[3] = -0.5 * ddx * ddx + 0.5 * ddx * ddx * ddx;
	g[0] = -0.5 * ddy + ddy * ddy - 0.5 * ddy * ddy * ddy;
	g[1] = 1.0 - 2.5 * ddy * ddy + 1.5 * ddy * ddy * ddy;
	g[2] = 0.5 * ddy + 2.0 * ddy * ddy - 1.5 * ddy * ddy * ddy;
	g[3] = -0.5 * ddy * ddy + 0.5 * ddy * ddy * ddy;

	gfac[12] = f[0] * g[0];
	gfac[8] = f[0] * g[1];
	gfac[4] = f[0] * g[2];
	gfac[0] = f[0] * g[3];
	gfac[13] = f[1] * g[0];
	gfac[9] = f[1] * g[1];
	gfac[5] = f[1] * g[2];
	gfac[1] = f[1] * g[3];
	gfac[14] = f[2] * g[0];
	gfac[10] = f[2] * g[1];
	gfac[6] = f[2] * g[2];
	gfac[2] = f[2] * g[3];
	gfac[15] = f[3] * g[0];
	gfac[11] = f[3] * g[1];
	gfac[7] = f[3] * g[2];
	gfac[3] = f[3] * g[3];

	/*
	**--------------------------------------------------------------
	**    Calculation of the interpolated value
	**    Applying the multiplication factors on the selected grid values
	**--------------------------------------------------------------
	*/
	value = 0.0;
	for (auto i = 0; i < 16; i++)
	{
		value += gfac[i] * record_value[i];
	}

	return 0;
}