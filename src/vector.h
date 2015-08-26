#include <vector>

//! Type definition for a standard vector of doubles
typedef std::vector<double> double_vect;

/*! Type definition for a 2D vector of doubles.
 *
 * Implemented as a standard vector of standard vectors of doubles.
 */
typedef std::vector<double_vect> double_2d;

//! @brief Centre a vector by subtracting the mean.
double_vect centre_vector(double_vect vect);

//! @brief Square a vector.
double_vect square_vector(double_vect vect);

//! @brief Sum a vector.
double sum_vector(double_vect vect);

//! @brief Elementwise multiplication of two vectors.
double_vect mult_vectors(double_vect vect1, double_vect vect2);

//! @brief Elementwise division of two vectors.
double_vect div_vectors(double_vect vect1, double_vect vect2);

//! @brief Calculate correlation between vectors.
double_vect correl_vectors(double_vect vect1, double_vect vect2,
                           double_vect vect3);

//! @brief Calculate rm values from two vectors.
double_vect rm_vectors(double_vect vect1, double_vect vect2);

//! @brief Calculate f values from two vectors.
double_vect f_vectors(double_vect correl_vect, double_vect rm_vect);

//! @brief Calculate h values from two vectors.
double_vect h_vectors(double_vect f_vect, double_vect rm_vect);

//! @brief Calculate z values.
double_vect z_vectors(double_vect cor1, double_vect cor2, double_vect sqrtn,
                      double_vect cross_cor, double_vect h_vect);

//! @brief Apply function to each element of a vector.
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect, F func);

//! @brief Apply function that combines elements of two vectors.
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect1, std::vector<T> vect2,
                                                                    F func);
//! @brief Apply function that reduces a 2D vector to a 1D vector.
template <typename T, typename F>
std::vector<T> reduce_2D_vect (std::vector<std::vector<T>> vect2D, F func);
