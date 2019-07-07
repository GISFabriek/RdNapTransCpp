// ***********************************************************************
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-06-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-07-2019
// ***********************************************************************
// C++ PORT from C version of RDNAPTRANS
// ***********************************************************************
#include "Cartesian.h"

Cartesian::Cartesian(const double x, const double y, const double z) {
	m_x = x;
	m_y = y;
	m_z = z;
}

Cartesian::Cartesian(const double x, const double y) : Cartesian(x, y, 0){
}

// Copy constructor.
Cartesian::Cartesian(const Cartesian& other)
{
	m_x = other.m_x; m_y = other.m_y; m_z = other.m_z;
}

// Move constructor.
Cartesian::Cartesian(Cartesian&& other) noexcept
{
	m_x = other.m_x; m_y = other.m_y; m_z = other.m_z;
}

// Copy assignment operator.
Cartesian& Cartesian::operator=(const Cartesian& other)
= default;

// Move assignment operator.
Cartesian& Cartesian::operator=(Cartesian&& other) noexcept
{
	if (this != &other)
	{
		m_x = other.m_x;
		m_y = other.m_y;
		m_z = other.m_z;
		other.m_x = 0;
		other.m_y = 0;
		other.m_z = 0;
	}
	return *this;
}

void Cartesian::set_x(const double x)
{
	m_x = x;
}

void Cartesian::set_y(const double y)
{
	m_y = y;
}

void Cartesian::set_z(const double z)
{
	m_z = z;
}

double Cartesian::get_x() const
{
	return m_x;
}

double Cartesian::get_y() const
{
	return m_y;
}

double Cartesian::get_z() const
{
	return m_z;
}

Cartesian Cartesian::with_z(const double z)
{
	Cartesian cartesian(m_x, m_y, z);
	return cartesian;
}
