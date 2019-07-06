#pragma once
class Geographic {
	double m_phi, m_lambda, m_h;
public:
	virtual ~Geographic();
	Geographic(double phi, double lambda, double h);
	Geographic(double phi, double lambda);
	Geographic(const Geographic& other);
	Geographic(Geographic&& other) noexcept;	
	Geographic& operator = (const Geographic& other);
	auto operator =(Geographic&& other) noexcept -> Geographic&;
	void set_phi(double phi);
	void set_lambda(double lambda);
	void set_h(double h);
	double get_phi() const;
	double get_lambda() const;
	double get_h() const;
	virtual Geographic with_h(double h);
};