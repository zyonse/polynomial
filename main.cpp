#include <iostream>
#include <chrono>
#include <optional>
#include <vector>
#include <fstream>

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

    return 0;
}