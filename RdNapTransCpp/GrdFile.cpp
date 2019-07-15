// ***********************************************************************
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-06-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-11-2019
// ***********************************************************************
// C++ PORT from C version of RDNAPTRANS
// ***********************************************************************
#include "GrdFile.h"
#include <cmath>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include "Constants.h"

#ifdef _MSC_VER  // If Visual C++ is used
#include "../out/build/x64-Debug/_cmrc/include/cmrc/cmrc.hpp"
#else // If GCC or CLang is used
#include "../cmake-build-debug/_cmrc/include/cmrc/cmrc.hpp"
#endif


CMRC_DECLARE(rdnaptransrc);
using namespace std;

// character string used for decoding the base64 encoded Grid strings
static const std::string base64_chars =
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz"
"0123456789+/";

/*
**--------------------------------------------------------------
**    Function name: get_decoded_string
**    Description:   extracts a decoded string from a resource file containing a base64 encoded string
**    (converted from a .grd file)
**
**    Parameter        Type                In/Out Req/Opt Default
**    file_name        std::string         in     req     none
**    result           std::string         out    req      -
**    Additional explanation of the meaning of parameters
**    file_name contains the name of the resource file
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
std::string GrdFile::get_decoded_string(const std::string& file_name)
{
	const auto fs = cmrc::rdnaptransrc::get_filesystem();
	const auto ff = fs.open(file_name);


	const std::stringstream string_stream(ff.cbegin());
	auto decoded = base64_decode(string_stream.str());
	return decoded;
}
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
**    decoded_string   the decoded string to be read
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
	int GrdFile::read_grd_file_header(const string& decoded_string,
		short int& size_x, short int& size_y,
		double& min_x, double& max_x,
		double& min_y, double& max_y,
		double& min_value, double& max_value)
{
	/*
	**--------------------------------------------------------------
	**    Grd files are binary grid files in the format of the program Surfer(R)
	**--------------------------------------------------------------
	**--------------------------------------------------------------
	**    Read file id
	**--------------------------------------------------------------
	*/
	char id[5];
	auto decoded_pos = 0;
	auto counter = 0;
	for (; decoded_pos < 4; decoded_pos++)
	{
		id[counter++] = decoded_string[decoded_pos];
	}
	id[4] = '\0';
	const auto id_string = string(id);

	/*
	**--------------------------------------------------------------
	**    Check
	**--------------------------------------------------------------
	*/
	if (id_string != "DSBB")
	{
		cerr << "not a valid grd resource" << endl;
		return -1;
	}

	/*
	**--------------------------------------------------------------
	**    Read output parameters
	**--------------------------------------------------------------
	*/
	extract_short(decoded_string, decoded_pos, size_x);
	extract_short(decoded_string, decoded_pos, size_y);
	extract_double(decoded_string, decoded_pos, min_x);
	extract_double(decoded_string, decoded_pos, max_x);
	extract_double(decoded_string, decoded_pos, min_y);
	extract_double(decoded_string, decoded_pos, max_y);
	extract_double(decoded_string, decoded_pos, min_value);
	extract_double(decoded_string, decoded_pos, max_value);

	return 0;
}

/*
**--------------------------------------------------------------
**    Function name: read_grd_file_body
**    Description:   reads a value from a grd file
**
**    Parameter      Type        In/Out Req/Opt Default
**    decoded_string std::string in     req     none
**    number         long int    in     req     none
**    value          float       out    -       none
**
**    Additional explanation of the meaning of parameters
**    decoded_string   the decoded string to be read
**    record_number  number defining the position in the file
**    record_value   output of the read value
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/

int GrdFile::read_grd_file_body(const string& decoded_string, long int record_number, float& record_value)
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

	const auto start = header_length + record_number * record_length;
	char record[record_length]{};
	auto counter = 0;
	for (auto i = start; i < start + record_length; i++)
	{
		record[counter++] = decoded_string[i];
	}
	memcpy(&record_value, &record, sizeof(record_value));

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
	const auto decoded = get_decoded_string(grd_file);
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
	auto error = read_grd_file_header(decoded, size_x, size_y, min_x, max_x, min_y, max_y, min_value, max_value);
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
		if (grd_file == "x2c.b64")
		{
			error = 1;
			value = 0.0;
		}
		if (grd_file == "y2c.b64")
		{
			error = 2;
			value = 0.0;
		}
		if (grd_file == "nlgeo04.b64")
		{
			error = 3;
		}
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
		error = read_grd_file_body(decoded, record_number[i], record_value[i]);
		if (error != 0)
		{
			return -1;
		}
		if (record_value[i] > max_value + Constants::PRECISION || record_value[i] < min_value - Constants::PRECISION)
		{
			cerr << "Outside validity area of " << grd_file << endl;
			if (grd_file == "x2c.b64")
			{
				error = 1; value = 0.0;
			}
			if (grd_file == "y2c.b64")
			{
				error = 2; value = 0.0;
			}
			if (grd_file == "nlgeo04.b64")
			{
				error = 3;
			}
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

	return error;
}

/*
**--------------------------------------------------------------
**    Function name: base64_decode
**    Description:   decodes a base64 encoded string
**
**    Parameter      Type              In/Out Req/Opt Default
**    encoded_string std::string       in     req     none
**    value          std::string       out    -       none
**
**    Additional explanation of the meaning of parameters
**    Returns the base64 decoded equivalent of the input string
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/

std::string GrdFile::base64_decode(std::string const& encoded_string) {
	auto in_len = encoded_string.size();
	auto i = 0;
	auto encoded_string_counter = 0;
	unsigned char char_array_4[4], char_array_3[3];
	std::string return_value;

	while (in_len-- && (encoded_string[encoded_string_counter] != '=') && is_base64(encoded_string[encoded_string_counter]))
	{
		char_array_4[i++] = encoded_string[encoded_string_counter++];
		if (i == 4) {
			for (i = 0; i < 4; i++)
			{
				char_array_4[i] = base64_chars.find(char_array_4[i]) & 0xff;
			}

			char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
			char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
			char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

			for (i = 0; (i < 3); i++)
			{
				return_value += char_array_3[i];
			}
			i = 0;
		}
	}

	if (i) {
		for (auto j = 0; j < i; j++)
		{
			char_array_4[j] = base64_chars.find(char_array_4[j]) & 0xff;
		}

		char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
		char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);

		for (auto j = 0; (j < i - 1); j++)
		{
			return_value += char_array_3[j];
		}
	}

	return return_value;
}

/*
**--------------------------------------------------------------
**    Function name: is_base64
**    Description:   checks if a character is base64 encoded
**
**    Parameter      Type                In/Out Req/Opt Default
**    c              unsigned char       in     req     none
**    value          bool                out    -       none
**
**    Additional explanation of the meaning of parameters
**
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
bool GrdFile::is_base64(unsigned char c) {
	return (isalnum(c) || (c == '+') || (c == '/'));
}

/*
**--------------------------------------------------------------
**    Function name: extract_short
**    Description:   extracts a short value from 2 bytes (chars) in a base64 decoded string
**
**    Parameter        Type                In/Out Req/Opt Default
**    decoded_string   std::string         in     req     none
**    decoded_pos      int&                in     req     none
**    result           short&              out    -       none
**    Additional explanation of the meaning of parameters
**    The position in the decoded file (decoded_pos) is a reference because it is updated in the method itself
**    The result variable is a reference because it contains the desired value after the method returns
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void GrdFile::extract_short(std::string decoded_string, int& decoded_pos, short& result)
{
	char int_array[2]{};
	auto counter = 0;
	const auto end = decoded_pos + 2;
	for (; decoded_pos < end; decoded_pos++)
	{
		int_array[counter++] = decoded_string[decoded_pos];
	}

	memcpy(&result, int_array, sizeof result);
}

/*
**--------------------------------------------------------------
**    Function name: extract_double
**    Description:   extracts a short value from 2 bytes (chars) in a base64 decoded string
**
**    Parameter        Type                In/Out Req/Opt Default
**    decoded_string   std::string         in     req     none
**    decoded_pos      int&                in     req     none
**    result           double&             out    -       none
**    Additional explanation of the meaning of parameters
**    The position in the decoded file (decoded_pos) is a reference because it is updated in the method itself
**    The result variable is a reference because it contains the desired value after the method returns
**    Return value: (besides the standard return values)
**    none
**--------------------------------------------------------------
*/
void GrdFile::extract_double(std::string decoded_string, int& decoded_pos, double& result)
{
	char double_array[8]{};
	auto counter = 0;
	const auto end = decoded_pos + 8;
	for (; decoded_pos < end; decoded_pos++)
	{
		double_array[counter++] = decoded_string[decoded_pos];
	}
	memcpy(&result, double_array, sizeof result);
}