// ***********************************************************************
// Author           : Willem A. Ligtendag, De GISFabriek
// Created          : 07-06-2019
//
// Last Modified By : Willem A. Ligtendag, De GISFabriek
// Last Modified On : 07-07-2019
// ***********************************************************************
// C++ PORT from C version of RDNAPTRANS
// ***********************************************************************
#include "Geographic.h"

Geographic::~Geographic() = default;

Geographic::Geographic(const double phi, const double lambda, const double h) {
	m_phi = phi;
	m_lambda = lambda;
	m_h = h;
}

Geographic::Geographic(const double phi, const double lambda): Geographic(phi, lambda, 0)
{
}

// Copy constructor.
Geographic::Geographic(const Geographic& other)
{
	m_phi = other.m_phi; m_lambda = other.m_lambda; m_h = other.m_h;
}

// Move constructor.
Geographic::Geographic(Geographic&& other) noexcept
{
	m_phi = other.m_phi; m_lambda = other.m_lambda; m_h = other.m_h;
}

// Copy assignment operator.
Geographic& Geographic::operator=(const Geographic& other)
= default;

// Move assignment operator.
Geographic& Geographic::operator=(Geographic&& other) noexcept
{
	if (this != &other)
	{
		m_phi = other.m_phi;
		m_lambda = other.m_lambda;
		m_h = other.m_h;
		other.m_phi = 0;
		other.m_lambda = 0;
		other.m_h = 0;
	}
	return *this;
}

// Latitude.
void Geographic::set_phi(const double phi)
{
	m_phi = phi;
}

// Longitude.
void Geographic::set_lambda(const double lambda)
{
	m_lambda = lambda;
}

// Ellipsoidal height.
void Geographic::set_h(const double h)
{
	m_h = h;
}

// Latitude.
double Geographic::get_phi() const
{
	return m_phi;
}

// Longitude.
double Geographic::get_lambda() const
{
	return m_lambda;
}

// Ellipsoidal height.
double Geographic::get_h() const
{
	return m_h;
}

// Creates a copy containing the height h as well.
Geographic Geographic::with_h(const double h)
{
	Geographic geographic(m_phi, m_lambda, h);
	return geographic;
}