#ifndef POLY_H
#define POLY_H

#include <vector>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <complex>
#include <thread>
#include <mutex>
#include <cmath>

using power = size_t;
using coeff = int;

class polynomial
{

public:
    /**
     * @brief Construct a new polynomial object that is the number 0 (ie. 0x^0)
     *
     */
    polynomial();

    /**
     * @brief Construct a new polynomial object from an iterator to pairs of <power,coeff>
     *
     * @tparam Iter
     *  An iterator that points to a std::pair<power, coeff>
     * @param begin
     *  The start of the container to copy elements from
     * @param end
     *  The end of the container to copy elements from
     */
    template <typename Iter>
    polynomial(Iter begin, Iter end) {
        if (begin == end) {
            polyData[0] = 0; // Handle empty case
            return;
        }
        size_t total_elements = std::distance(begin, end);
        if (total_elements > 1000) { // Threshold for multithreading
            size_t num_threads = std::thread::hardware_concurrency();
            size_t chunk_size = total_elements / num_threads;
            std::vector<std::thread> threads;
            Iter it = begin;
            for (size_t i = 0; i < num_threads; ++i) {
                Iter chunk_start = it;
                Iter chunk_end = chunk_start;
                size_t steps = (i == num_threads - 1) ? total_elements - chunk_size * i : chunk_size;
                for (size_t j = 0; j < steps && chunk_end != end; ++j) {
                    ++chunk_end;
                }
                threads.emplace_back([this, chunk_start, chunk_end]() {
                    for (auto it = chunk_start; it != chunk_end; ++it) {
                        std::lock_guard<std::mutex> lock(this->mutex_);
                        this->polyData[it->first] = it->second;
                    }
                });
                it = chunk_end;
            }
            for (auto& thread : threads) {
                thread.join();
            }
        } else {
        while (begin != end) {
            polyData[begin->first] = begin->second;
            ++begin;
            }
        }
    }

    /**
     * @brief Construct a new polynomial object from an existing polynomial object
     *
     * @param other
     *  The polynomial to copy
     */
    polynomial(const polynomial &other);

    /**
     * @brief Prints the polynomial.
     *
     * Only used for debugging, isn't graded.
     *
     */
    void print() const;

    /**
     * @brief Turn the current polynomial instance into a deep copy of another
     * polynomial
     *
     * @param other
     * The polynomial to copy
     * @return
     * A reference to the copied polynomial
     */
    polynomial &operator=(const polynomial &other);


    /**
     * Overload the +, * and % operators. The function prototypes are not
     * provided.
     * 
     * Addition (+) should support
     * 1. polynomial + polynomial
     * 2. polynomial + int
     * 3. int + polynomial
     * 
     * Multiplication (*) should support
     * 1. polynomial * polynomial
     * 2. polynomial * int
     * 3. int * polynomial
     * 
     * Modulo (%) should support
     * 1. polynomial % polynomial
     */
    polynomial operator+(const polynomial &other) const;
    polynomial operator+(int val) const;
    friend polynomial operator+(int val, const polynomial &other);

    polynomial operator*(const polynomial &other) const;
    polynomial operator*(int val) const;
    friend polynomial operator*(int val, const polynomial &other);

    polynomial operator-(const polynomial &other) const;

    polynomial operator%(const polynomial &other) const;
    

    /**
     * @brief Returns the degree of the polynomial
     *
     * @return size_t
     *  The degree of the polynomial
     */
    size_t find_degree_of() const;

    /**
     * @brief Returns a vector that contains the polynomial is canonical form. This
     *        means that the power at index 0 is the largest power in the polynomial,
     *        the power at index 1 is the second largest power, etc.
     *
     *        ie. x^2 + 7x^4 + 1 would be returned as [(4,7),(2,1),(0,1)]
     *
     *        Note: any terms that have a coefficient of zero aren't returned in
     *        in the canonical form, with one exception.
     *        See the above example (there's no x^3 term, so
     *        there's no entry in the vector for x^3)
     *
     *        The only exception to this is the polynomial 0. If your polynomial is
     *        just the constant 0 then you would return a single entry of [(0,0)]
     *
     *        ie. y = 0 would be returned as [(0,0)]
     *
     * @return std::vector<std::pair<power, coeff>>
     *  A vector of pairs representing the canonical form of the polynomial
     */
    std::vector<std::pair<power, coeff>> canonical_form() const;

private:
    std::unordered_map<power, coeff> polyData;
    mutable std::mutex mutex_;
    
    // FFT helper functions
    static void fft(std::vector<std::complex<double>> &a, bool inverse = false);
    polynomial multiply_fft(const polynomial &other) const;
    static size_t next_power_of_two(size_t n);
};

#endif
