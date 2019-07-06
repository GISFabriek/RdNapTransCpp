#pragma once

class Cartesian {
	double m_x, m_y, m_z;
public:
	virtual ~Cartesian() = default;
	Cartesian(double x, double y, double z);
	Cartesian(double x, double y);
	Cartesian(const Cartesian& other);
	Cartesian(Cartesian&& other) noexcept;
	Cartesian& operator = (const Cartesian& other);
	auto operator =(Cartesian&& other) noexcept->Cartesian &;
	void set_x(double x);
	void set_y(double y);
	void set_z(double z);
	double get_x() const;
	double get_y() const;
	double get_z() const;
	virtual Cartesian with_z(double z);
};
