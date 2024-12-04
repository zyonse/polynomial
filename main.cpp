#include <iostream>
#include <chrono>
#include <optional>
#include <vector>
#include <fstream>
#include <sstream>

#include "poly.h"

std::optional<double> poly_test(polynomial& p1,
                                polynomial& p2,
                                std::vector<std::pair<power, coeff>> solution)

{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    polynomial p3 = p1 * p2;

    auto p3_can_form = p3.canonical_form();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    p3.print();

    if (p3_can_form != solution)
    {
        return std::nullopt;
    }

    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
}

// Convert polynomial string to polynomial object
polynomial parsePolynomial (std::string poly_str) {
    std::vector<std::pair<power, coeff>> poly_input;
    std::stringstream ss(poly_str);
    std::string term;
    while (std::getline(ss, term, '\n')) {
        size_t x_pos = term.find('x');
        size_t power_pos = term.find('^');
        if (x_pos != std::string::npos && power_pos != std::string::npos) {
            std::string coeff_str = term.substr(0, x_pos);
            std::string power_str = term.substr(power_pos + 1);
            coeff c = std::stod(coeff_str);
            power p = std::stoi(power_str);
            poly_input.push_back({p, c});
        }
    }
    return polynomial(poly_input.begin(), poly_input.end());
}

bool test_file() {
        // Read polynomials from simple_poly.txt
    std::ifstream infile("simple_poly.txt");
    if (!infile) {
        std::cerr << "Failed to open simple_poly.txt" << std::endl;
        return 1;
    }
    std::string poly1_str, poly2_str;
    getline(infile, poly1_str, ';');
    getline(infile, poly2_str, ';');
    infile.close();

    // Parse the polynomials
    polynomial poly1 = parsePolynomial(poly1_str);
    polynomial poly2 = parsePolynomial(poly2_str);

    // Addition with timing
    auto start_add = std::chrono::high_resolution_clock::now();
    polynomial sum = poly1 + poly2;
    auto end_add = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_add = end_add - start_add;
    std::cout << "Addition time: " << duration_add.count() << " seconds" << std::endl;

    // Multiplication with timing
    auto start_mul = std::chrono::high_resolution_clock::now();
    polynomial product = poly1 * poly2;
    auto end_mul = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_mul = end_mul - start_mul;
    std::cout << "Multiplication time: " << duration_mul.count() << " seconds" << std::endl;

    // Modulus with timing
    auto start_mod = std::chrono::high_resolution_clock::now();
    //polynomial remainder = poly1 % poly2;
    auto end_mod = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_mod = end_mod - start_mod;
    std::cout << "Modulus time: " << duration_mod.count() << " seconds" << std::endl;

    // Verify multiplication result with expected result from result.txt
    std::ifstream result_file("result.txt");
    if (result_file) {
        std::string expected_str;
        getline(result_file, expected_str, ';');
        result_file.close();
        polynomial expected = parsePolynomial(expected_str);
        if (product.canonical_form() == expected.canonical_form()) {
            std::cout << "Multiplication result matches expected result." << std::endl;
        } else {
            std::cout << "Multiplication result does not match expected result." << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Failed to open result.txt" << std::endl;
        return 1;
    }
    return 0;
}

// Function to test multiplication of sparse polynomials
void test_sparse_polynomials() {
    std::vector<std::pair<power, coeff>> sparse_poly1 = {
        {1000, 5}, {5000, 10}, {100000, 3}
    };
    std::vector<std::pair<power, coeff>> sparse_poly2 = {
        {2000, 2}, {5000, 4}, {150000, 6}
    };

    polynomial sp1(sparse_poly1.begin(), sparse_poly1.end());
    polynomial sp2(sparse_poly2.begin(), sparse_poly2.end());

    auto start_sparse_mul = std::chrono::high_resolution_clock::now();
    polynomial sparse_product = sp1 * sp2;
    auto end_sparse_mul = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_sparse_mul = end_sparse_mul - start_sparse_mul;

    std::cout << "Sparse Multiplication time: " << duration_sparse_mul.count() << " seconds" << std::endl;
}

std::optional<double> test_polynomial_modulo(polynomial& dividend, 
                                           polynomial& divisor,
                                           std::vector<std::pair<power, coeff>> expected_result) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    polynomial remainder = dividend % divisor;
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    
    std::cout << "Remainder: ";
    remainder.print();
    
    if (remainder.canonical_form() != expected_result) {
        return std::nullopt;
    }
    
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
}

int main()
{
    /** We're doing (x+1)^2, so solution is x^2 + 2x + 1*/
    std::vector<std::pair<power, coeff>> solution = {{2,1}, {1,2}, {0,1}};

    /** This holds (x+1), which we'll pass to each polynomial */
    std::vector<std::pair<power, coeff>> poly_input = {{1,1}, {0,1}};

    polynomial p1(poly_input.begin(), poly_input.end());
    polynomial p2(poly_input.begin(), poly_input.end());

    std::optional<double> result = poly_test(p1, p2, solution);

    // Test polynomial in simple_poly.txt and compare to expected result in result.txt
    if (result.has_value())
    {
        std::cout << "Passed test, took " << result.value()/1000 << " seconds" << std::endl;
    } 
    else 
    {
        std::cout << "Failed test" << std::endl;
    }

    if (!test_file()) {
        std::cout << "Passed simple_poly.txt test" << std::endl;
    } else {
        std::cout << "Failed simple_poly.txt test" << std::endl;
        return 1;
    }

    // Test polynomial modulo
    std::vector<std::pair<power, coeff>> modulo_expected = {{0,0}}; // Expected result for (x+1) % (x+1)
    std::optional<double> modulo_result = test_polynomial_modulo(p1, p2, modulo_expected);
    
    if (modulo_result.has_value()) {
        std::cout << "Passed modulo test, took " << modulo_result.value() << " ms" << std::endl;
    } else {
        std::cout << "Failed modulo test" << std::endl;
    }

    // Test polynomial modulo with complex polynomials
    // Testing (2x^3 + 3x^2 + 1) % (x^2 + 1)
    std::vector<std::pair<power, coeff>> dividend_input = {
        {9, 1},    // \( x^9 \)
        {8, -2},   // \( -2x^8 \)
        {7, 3},    // \( 3x^7 \)
        {6, -4},   // \( -4x^6 \)
        {5, 5},    // \( 5x^5 \)
        {4, -6},   // \( -6x^4 \)
        {3, 7},    // \( 7x^3 \)
        {2, -8},   // \( -8x^2 \)
        {1, 9},    // \( 9x \)
        {0, -10}   // \( -10 \)
    };
    std::vector<std::pair<power, coeff>> divisor_input = {
        {3, 1},    // \( x^3 \)
        {2, -1},   // \( -x^2 \)
        {0, 1}     // \( 1 \)
    };
    std::vector<std::pair<power, coeff>> modulo_expected_complex = {
        {2, -6},   // \( -6x^2 \)
        {1, 14},   // \( 14x \)
        {0, -15}   // \( -15 \)
    };
    polynomial complex_dividend(dividend_input.begin(), dividend_input.end());
    polynomial complex_divisor(divisor_input.begin(), divisor_input.end());
    
    std::optional<double> modulo_result_complex = test_polynomial_modulo(complex_dividend, complex_divisor, modulo_expected_complex);
    
    if (modulo_result_complex.has_value()) {
        std::cout << "Passed complex modulo test, took " << modulo_result_complex.value() << " ms" << std::endl;
    } else {
        std::cout << "Failed complex modulo test" << std::endl;
    }

    // Call the sparse polynomial test function
    test_sparse_polynomials();

    return 0;
}