#pragma once
#define DEF_UNDEFINED_SUM 0
#define ZERO 1e-10
#define BETA 1

#include "series_acceleration.h" // Include the series header
#include <vector> // Include the vector library

/**
  * @brief generalised Levi-Sidi class template.
  * @authors Venajalainen
  * @tparam T The type of the elements in the series, K The type of enumerating integer, series_templ is the type of series whose convergence we accelerate
  */

template<typename T, typename K, typename series_templ> class levi_sidi_algorithm;

template<typename T, typename K, typename series_templ> class u_levi_sidi_algorithm;

template<typename T, typename K, typename series_templ> class t_levi_sidi_algorithm;

template<typename T, typename K, typename series_templ> class v_levi_sidi_algorithm;

class u_transform {
public:
	u_transform() {}

	template<typename T, typename K>
	T operator()(const int& n, const int& j, const series_base<T, K>* series) const;
};

template<typename T, typename K>
T u_transform::operator()(const int& n, const int& j, const series_base<T, K>* series) const {
	return 1 / ((BETA + n) * series->operator()(n + j));
}


class t_transform {
public:
	t_transform() {}

	template<typename T, typename K>
	T operator()(const int& n, const int& j, const series_base<T, K>* series) const;
};

template<typename T, typename K>
T t_transform::operator()(const int& n, const int& j, const series_base<T, K>* series) const {
	return 1 / series->operator()(n + j);
}


class v_transform {
public:
	v_transform() {}

	template<typename T, typename K>
	T operator()(const int& n, const int& j, const series_base<T, K>* series) const;
};

template<typename T, typename K>
T v_transform::operator()(const int& n, const int& j, const series_base<T, K>* series) const {
	return (series->operator()(n + j + 1) - series->operator()(n + j)) / (series->operator()(n + j + 1) * series->operator()(n + j));
}


template<typename T, typename K, typename series_templ>
class levi_sidi_algorithm : public series_acceleration<T, K, series_templ>
{
protected:
	T minus_one_raised_to_power_n(K n) const
	{
		return n % 2 ? -1 : 1;
	}

	T binomial_coefficient(const T n, const K k) const
	{
		T b_c = 1;
		for (int i = 0; i < k; ++i)
			b_c = b_c * (n - static_cast<T>(i)) / (i + 1);
		return b_c;
	}

	/**
	* @brief General function to calculate S-tranformation. Implemented u,t and v transformations
	* @authors Venajalainen
	* @param k The number of terms in the partial sum.
	* @param n the order of transformation
	* @param id the type of remainder to use in transformation
	*/

	template<class remainderType>
	T calculate(const K& k, const int& n, remainderType remainder_func) const {

		if (n < 0)
			throw std::domain_error("negative integer in input");
		if (BETA <= 0)
			throw std::domain_error("beta cannot be initiared by a negative number or a zero");


		T numerator = T(0), denominator = T(0);
		T w_n, rest;
		T up = T(1), down = T(1);

		for (int j = 0; j <= k; ++j) {

			rest = minus_one_raised_to_power_n(j) * binomial_coefficient(k, j);

			T up = 1, down = 1;
			for (int m = 0; m < k - 1; ++m) {
				up *= (BETA + n + j + m);
				down *= (BETA + n + k + m);
			}

			rest = rest * (up / down);

			w_n = remainder_func(n, j, this->series);

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

public:

	friend class u_levi_sidi_algorithm<T, K, series_templ>;
	friend class t_levi_sidi_algorithm<T, K, series_templ>;
	friend class v_levi_sidi_algorithm<T, K, series_templ>;


	/**
	* @brief Parameterized constructor to initialize the Levin Algorithm.
	* @authors Venajalainen
	* @param series The series class object to be accelerated
	*/

	levi_sidi_algorithm(const series_templ& series);

	/**
	* @brief Abstract method for the inherited classes below
	* Computes the partial sum after the transformation using the Levin Algorithm.
	* For more information,
	* @param k The number of terms in the partial sum.
	* @param n The order of transformation.
	* @return The partial sum after the transformation.
	*/

	T operator()(const K k, const int n) const = 0;

};

template<typename T, typename K, typename series_templ>
levi_sidi_algorithm< T, K, series_templ>::levi_sidi_algorithm(const series_templ& series) : series_acceleration<T, K, series_templ>(series) {}

template <typename T, typename K, typename series_templ>
class u_levi_sidi_algorithm : public levi_sidi_algorithm<T, K, series_templ>
{
public:
	/**
	* @brief Parameterized constructor to initialize the Levin Algorithm.
	* @authors Venajalainen
	* @param series The series class object to be accelerated
	*/
	u_levi_sidi_algorithm(const series_templ& series);

	/**
	* @brief Fast impimentation of Levin algorithm.
	* Computes the partial sum after the transformation using the Levin Algorithm.
	* For more information,
	* @param k The number of terms in the partial sum.
	* @param n The order of transformation.
	* @return The partial sum after the transformation.
	*/
	T operator()(const K k, const int n) const;

};

template<typename T, typename K, typename series_templ>
u_levi_sidi_algorithm< T, K, series_templ> ::u_levi_sidi_algorithm(const series_templ& series) : levi_sidi_algorithm<T, K, series_templ>(series) {}

template<typename T, typename K, typename series_templ>
T u_levi_sidi_algorithm<T, K, series_templ>::operator()(const K k, const int n) const { return levi_sidi_algorithm<T, K, series_templ>::calculate(k, n, u_transform{}); }

template <typename T, typename K, typename series_templ>
class t_levi_sidi_algorithm : public levi_sidi_algorithm<T, K, series_templ>
{
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
	* @return The partial sum after the transformation.
	*/
	T operator()(const K k, const int n) const;
};

template<typename T, typename K, typename series_templ>
t_levi_sidi_algorithm < T, K, series_templ> ::t_levi_sidi_algorithm(const series_templ& series) : levi_sidi_algorithm<T, K, series_templ>(series) {}

template<typename T, typename K, typename series_templ>
T t_levi_sidi_algorithm<T, K, series_templ>::operator()(const K k, const int n) const { return levi_sidi_algorithm<T, K, series_templ>::calculate(k, n, t_transform{}); }

template <typename T, typename K, typename series_templ>
class v_levi_sidi_algorithm : public levi_sidi_algorithm<T, K, series_templ>
{
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
	* @param remainder is parametr to choose the form of a w_i
	* @return The partial sum after the transformation.
	*/
	T operator()(const K k, const int n) const;
};


template<typename T, typename K, typename series_templ>
v_levi_sidi_algorithm < T, K, series_templ> ::v_levi_sidi_algorithm(const series_templ& series) : levi_sidi_algorithm<T, K, series_templ>(series) {}

template<typename T, typename K, typename series_templ>
T v_levi_sidi_algorithm<T, K, series_templ>::operator()(const K k, const int n) const { return levi_sidi_algorithm<T, K, series_templ>::calculate(k, n, v_transform{}); }