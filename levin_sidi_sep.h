#pragma once
#define DEF_UNDEFINED_SUM 0
#define ZERO 1e-10

#include "series_acceleration.h" // Include the series header
#include <vector> // Include the vector library

/**
  * @brief generalised Levi-Sidi class template.
  * @authors Naumov
  * @tparam T The type of the elements in the series, K The type of enumerating integer, series_templ is the type of series whose convergence we accelerate
  */

template <typename T, typename K>
constexpr const T minus_one_raised_to_power_n(K n)
{
	return n % 2 ? -1 : 1;
}
template <typename T, typename K>
constexpr const T binomial_coefficient(const T n, const K k)
{
	T b_c = 1;
	for (int i = 0; i < k; ++i)
		b_c = b_c * (n - static_cast<T>(i)) / (i + 1);
	return b_c;
}

template <typename T, typename K, typename series_templ>
class u_levi_sidi_algorithm : public series_acceleration<T, K, series_templ>
{
private:
	const T beta = 1;
public:
	/**
	* @brief Parameterized constructor to initialize the Levin Algorithm.
	* @authors Naumov
	* @param series The series class object to be accelerated
	*/
	u_levi_sidi_algorithm(const series_templ& series);

	/**
	* @brief Fast impimentation of Levin algorithm.
	* Computes the partial sum after the transformation using the Levin Algorithm.
	* For more information, 
	* @param k The number of terms in the partial sum.
	* @param n The order of transformation.
	* @param beta is arbitrary parametr, beta > 0
	* @param remainder is parametr to choose the form of a w_i
	* realized u,t and v remainders
	* 1 -> u
	* 2 -> t
	* 3 -> v
	* @return The partial sum after the transformation.
	*/
	T operator()(const K k, const int n) const;
};

template<typename T, typename K, typename series_templ>
u_levi_sidi_algorithm< T, K, series_templ> :: u_levi_sidi_algorithm(const series_templ& series) : series_acceleration<T, K, series_templ>(series) {}

template<typename T, typename K, typename series_templ>
T u_levi_sidi_algorithm<T, K, series_templ>::operator()(const K k, const int n) const
{
	if (n < 0)
		throw std::domain_error("negative integer in input");
	if (beta <= 0)
		throw std::domain_error("beta cannot be initiared by a negative number or a zero");


	T numerator = T(0), denominator = T(0);
	T w_n, rest;
	T up = T(1), down = T(1);

	for (int j = 0; j <= k; ++j) {

		rest = minus_one_raised_to_power_n<T,K>(j) * binomial_coefficient<T,K>(k, j);

		T up = 1, down = 1;
		for (int m = 0; m < k - 1; ++m) {
			up *= (beta + n + j + m);
			down *= (beta + n + k + m);
		}

		rest = rest * (up / down);

		w_n = 1 / ((beta + n)*this->series->operator()(n + j));

		numerator += rest * this->series->S_n(n + j) * w_n;
		denominator += rest * w_n;

	}

	if(std::abs(denominator)<ZERO)
		throw std::overflow_error("division by zero");

	if (!std::isfinite(numerator))
		throw std::overflow_error("numerator is infinite");

	numerator/= denominator;

	return numerator;
}

template <typename T, typename K, typename series_templ>
class t_levi_sidi_algorithm : public series_acceleration<T, K, series_templ>
{
private:
	const T beta = 1;
public:
	/**
	* @brief Parameterized constructor to initialize the Levin Algorithm.
	* @authors Naumov
	* @param series The series class object to be accelerated
	*/
	t_levi_sidi_algorithm(const series_templ& series);

	/**
	* @brief Fast impimentation of Levin algorithm.
	* Computes the partial sum after the transformation using the Levin Algorithm.
	* For more information,
	* @param k The number of terms in the partial sum.
	* @param n The order of transformation.
	* @param beta is arbitrary parametr, beta > 0
	* @param remainder is parametr to choose the form of a w_i
	* realized u,t and v remainders
	* 1 -> u
	* 2 -> t
	* 3 -> v
	* @return The partial sum after the transformation.
	*/
	T operator()(const K k, const int n) const;
};

template<typename T, typename K, typename series_templ>
t_levi_sidi_algorithm < T, K, series_templ> ::t_levi_sidi_algorithm(const series_templ& series) : series_acceleration<T, K, series_templ>(series) {}

template<typename T, typename K, typename series_templ>
T t_levi_sidi_algorithm<T, K, series_templ>::operator()(const K k, const int n) const
{
	if (n < 0)
		throw std::domain_error("negative integer in input");
	if (beta <= 0)
		throw std::domain_error("beta cannot be initiared by a negative number or a zero");

	T numerator = T(0), denominator = T(0);
	T w_n, rest;
	T up = T(1), down = T(1);

	for (int j = 0; j <= k; ++j) {

		rest = minus_one_raised_to_power_n<T, K>(j) * binomial_coefficient<T, K>(k, j);

		T up = 1, down = 1;
		for (int m = 0; m < k - 1; ++m) {
			up *= (beta + n + j + m);
			down *= (beta + n + k + m);
		}

		rest = rest * (up / down);

		w_n = 1 / this->series->operator()(n + j);

		numerator += rest * this->series->S_n(n + j) * w_n;
		denominator += rest * w_n;

	}

	if (std::abs(denominator) < ZERO)
	throw std::overflow_error("division by zero");

	if (!std::isfinite(numerator))
		throw std::overflow_error("numerator is infinite o");

	numerator /= denominator;

	return numerator;
}

template <typename T, typename K, typename series_templ>
class v_levi_sidi_algorithm : public series_acceleration<T, K, series_templ>
{
private:
	/**
	* @brief
	* @param beta 
	*/
	const T beta = 1;
public:
	/**
	* @brief Parameterized constructor to initialize the Levin Algorithm.
	* @authors Naumov
	* @param series The series class object to be accelerated
	*/
	v_levi_sidi_algorithm(const series_templ& series);

	/**
	* @brief Fast impimentation of Levin algorithm.
	* Computes the partial sum after the transformation using the Levin Algorithm.
	* For more information,
	* @param k The number of terms in the partial sum.
	* @param n The order of transformation.
	* @param beta is arbitrary parametr, beta > 0
	* @param remainder is parametr to choose the form of a w_i
	* realized u,t and v remainders
	* 1 -> u
	* 2 -> t
	* 3 -> v
	* @return The partial sum after the transformation.
	*/
	T operator()(const K k, const int n) const;
};


template<typename T, typename K, typename series_templ>
v_levi_sidi_algorithm < T, K, series_templ> ::v_levi_sidi_algorithm(const series_templ& series) : series_acceleration<T, K, series_templ>(series) {}

template<typename T, typename K, typename series_templ>
T v_levi_sidi_algorithm<T, K, series_templ>::operator()(const K k, const int n) const
{
	if (n < 0)
		throw std::domain_error("negative integer in input");
	if (beta <= 0)
		throw std::domain_error("beta cannot be initiared by a negative number or a zero");

	T numerator = T(0), denominator = T(0);
	T w_n, rest;
	T up = T(1), down = T(1);

	for (int j = 0; j <= k; ++j) {

		rest = minus_one_raised_to_power_n<T, K>(j) * binomial_coefficient<T, K>(k, j);

		T up = 1, down = 1;
		for (int m = 0; m < k - 1; ++m) {
			up *= (beta + n + j + m);
			down *= (beta + n + k + m);
		}

		rest = rest * (up / down);

		T a0 = this->series->operator()(n + j), a1 = this->series->operator()(n + j + 1);
		w_n = (a1 - a0) / (a1 * a0);

		numerator += rest * this->series->S_n(n + j) * w_n;
		denominator += rest * w_n;

	}

	if (std::abs(denominator) < ZERO)
	throw std::overflow_error("division by zero");

	if (!std::isfinite(numerator))
		throw std::overflow_error("numerator is infinite");

	numerator /= denominator;

	return numerator;
}