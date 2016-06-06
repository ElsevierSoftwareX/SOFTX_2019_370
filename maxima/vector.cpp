#include <numeric>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "vector.h"

/*! Calculate the mean of all values in the vector and subtract from each
 * individual value.
 *
 * @param vect The vector to centre.
 *
 * @return The centred vector.
 */
double_vect centre_vector(double_vect vect)
{
    double sum  = std::accumulate(vect.begin(), vect.end(), 0.0);
    double mean = sum / vect.size();
    double_vect centered;

    for (auto v : vect) {
        centered.push_back(v - mean);
    }

    return centered;
}

/*! Multiply each value in a vector by itself.
 *
 * @param vect The vector to square.
 *
 * @return The squared vector.
 */
double_vect square_vector(double_vect vect)
{
    double_vect squared;

    for (auto v : vect) {
        squared.push_back(v * v);
    }

    return squared;
}

/*! Add up all the values in a vector.
 *
 * @param vect The vector to sum.
 *
 * @return The summed value.
 */
double sum_vector(double_vect vect)
{
    double sum = std::accumulate(vect.begin(), vect.end(), 0.0);

    return sum;
}

/*! Multiply pairs of values from two vectors of equal length.
 *
 * @param vect1 First vector to be multiplied.
 * @param vect2 Second vector to be multiplied.
 *
 * @return Vector of pairwise multiples.
 */
double_vect mult_vectors(double_vect vect1, double_vect vect2)
{
    double_vect mult;

    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Pairwise multiplication
    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        mult.push_back(vect1[idx] * vect2[idx]);
    }

    return mult;
}

/*! Divide values from one vector by values from a second vector of equal
 * length.
 *
 * @param vect1 The vector to be divided (numerator).
 * @param vect2 The vector to divide by (denominator).
 *
 * @return Vector of vect1[i] / vect2[i] values.
 */
double_vect div_vectors(double_vect vect1, double_vect vect2)
{
    double_vect divided;

    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Pairwise division
    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        divided.push_back(vect1[idx] / vect2[idx]);
    }

    return divided;

}

/*!@todo Explain input vectors and return correl vector
 */
double_vect correl_vectors(double_vect vect1, double_vect vect2,
                           double_vect vect3)
{
    double_vect correlated;
    double_vect mult;

    mult = mult_vectors(vect2, vect3);
    std::transform(mult.begin(), mult.end(), mult.begin(),
                                            (double(*)(double)) std::sqrt);
    correlated = div_vectors(vect1, mult);

    for (auto& c : correlated) {
        if(std::isnan(c)) {
            c = 0;
        }
        if(c < 0) {
            c = 0;
        }
    }

    return correlated;
}

/*! Calculate values equal to 0.5 * (vect1[i]^2 + vect2[i]^2)
 *
 * @param vect1 The first values in the equation.
 * @param vect2 The second values in the equation.
 *
 * @return Vector of calculated RM values.
 */
double_vect rm_vectors(double_vect vect1, double_vect vect2)
{
    double_vect rm;

    // Throw exception if vectors of different lengths
    if (vect1.size() != vect2.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Calculate RM values.
    for (size_t idx = 0; idx < vect1.size(); ++idx) {
        double tmp = (vect1[idx] * vect1[idx]) + (vect2[idx] * vect2[idx]);
        rm.push_back(0.5 * tmp);
    }

    return rm;

}

/*!F values are bounded at an upper value of 1.0.
 *
 * @param correl_vect Vector returned by correl_vectors.
 * @param rm_vect Vector return by rm_vectors.
 *
 * @return Vector of calculated f values.
 *
 * @todo Explain what f value is.
 */
double_vect f_vectors(double_vect correl_vect, double_vect rm_vect)
{
    double_vect f_vect;

    // Throw exception if vectors of different lengths
    if (correl_vect.size() != rm_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Calculate f values
    for (size_t idx = 0; idx < correl_vect.size(); ++idx) {
        double correl = correl_vect[idx];
        double rm     = rm_vect[idx];

        f_vect.push_back((1.0 - correl) / (2.0 * (1.0 - rm)));
    }

    // Set upper bound to 1.0
    for (auto& f : f_vect) {
        if(f > 1.0) {
            f = 1.0;
        }
    }

    return f_vect;
}

/*! @param f_vect Vector of values returned by f_vectors.
 * @param rm_vect Vector of values returned by rm_vectors.
 *
 * @return Vector of calculated h values.
 *
 * @todo Explain what h value is.
 */
double_vect h_vectors(double_vect f_vect, double_vect rm_vect)
{
    double_vect h_vect;

    // Throw exception if vectors of different lengths
    if (f_vect.size() != rm_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Calculate h value
    for (size_t idx = 0; idx < f_vect.size(); ++idx) {
        double f  = f_vect[idx];
        double rm = rm_vect[idx];

        h_vect.push_back((1.0 - f * rm) / (1.0 - rm));
    }

    return h_vect;
}

/*! @todo Explain z score calculation and params
 */
double_vect z_vectors(double_vect cor1, double_vect cor2, double_vect sqrtn,
                      double_vect cross_cor, double_vect h_vect)
{
    double_vect z_vect;

    // Throw exception if vectors of different lengths
    if (cor1.size() != cor2.size()      ||
        cor1.size() != sqrtn.size()     ||
        cor1.size() != cross_cor.size() ||
        cor1.size() != h_vect.size()) {
        throw std::invalid_argument("Vectors have different lengths");
    }

    // Calculate z score
    for (size_t idx = 0; idx < cor1.size(); ++idx) {

        double z1  = std::atanh(cor1[idx]);
        double z2  = std::atanh(cor2[idx]);

        double num   = (z1 - z2) * sqrtn[idx];
        double denom = 2.0 * (1.0 - cross_cor[idx]) * h_vect[idx];

        z_vect.push_back(num / std::sqrt(denom));
    }

    return z_vect;
}
