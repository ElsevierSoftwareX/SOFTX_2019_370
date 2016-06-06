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

/*! Template function that takes a 2D vector of any type and applies a function
 * to each element vector to give a 1D vector of the same type.
 *
 * @param vect2D 2D vector to apply the function too.
 * @param func Function to apply to each element of the 2D vector. Must take a
 * vector and return a single value of same type as the vector.
 *
 * @returns Vector containing results from applying _func_ to _vect2D_.
 */
template <typename T, typename F>
std::vector<T> reduce_2D_vect (std::vector<std::vector<T>> vect2D, F func)
{
    std::vector<T> reduced;

    for (auto vect : vect2D) {
        reduced.push_back(func(vect));
    }

    return reduced;
}


/*! Template function that takes two vectors of any type and applies a function
 * to each pair of elements to give another vector of the same type.
 *
 * @param vect1 First vector to apply the function too.
 * @param vect2 Second vector to apply the function too
 * @param func Function to apply to each pair of elements. Must return a single
 * value of the same type as the vectors
 *
 * @returns Vector containing results from applying _func_ to _vect1_ and
 * _vect2_.
 */
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect1, std::vector<T> vect2,
                                                                     F func)
{
    std::vector<T> applied;

    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        applied.push_back(func(vect1[idx], vect2[idx]));
    }

    return applied;
}

/*! Template function that takes a vector of any type and applies a function
 * to each element to give another vector of the same type.
 *
 * @param vect Vector to apply the function too.
 * @param func Function to apply to each element of the vector. Must return a
 * value of same type as the vector.
 *
 * @returns Vector containing results from applying _func_ to _vect_.
 */
template <typename T, typename F>
std::vector<T> apply_vect_func(std::vector<T> vect, F func)
{
    std::vector<T> applied;

    for (auto v : vect) {
        applied.push_back(func(v));
    }

    return applied;
}

