// ***********************************************************************
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-06-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-11-2019
// ***********************************************************************
// C++ PORT from C version of RDNAPTRANS
// ***********************************************************************
#pragma once
#include <string>

class GrdFile
{
	static bool is_base64(unsigned char c);
	static std::string base64_decode(std::string const& encoded_string);
	static std::string get_decoded_string(const std::string& file_name);
	static void extract_short(std::string decoded_string, int& decoded_pos, short& result);
	static void extract_double(std::string decoded_string, int& decoded_pos, double& result);
public:
	static int read_grd_file_header(const std::string& decoded_string,
		short int& size_x, short int& size_y,
		double& min_x, double& max_x,
		double& min_y, double& max_y,
		double& min_value, double& max_value);
	static auto read_grd_file_body(const std::string& decoded_string, long int record_number, float& record_value) -> int;
	static auto grid_interpolation(double x, double y, const std::string& grd_file, double& value) -> int;

};
