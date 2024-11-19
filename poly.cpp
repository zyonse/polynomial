#include "poly.h"

polynomial::polynomial() : polyData(1, 0) {}

polynomial::polynomial(const polynomial &other) {
    polyData.resize(other.polyData.size(), 0);
    for (size_t i = 0; i < other.polyData.size(); i++) {
        polyData[i] = other.polyData[i];
    }
}

void polynomial::print() const {
    for (int i = 0; i < polyData.size(); i++) {
        std::cout << polyData[i] << "x^" << i << " ";
    }
    std::cout << std::endl;
}

polynomial &polynomial::operator=(const polynomial &other) {
    size_t size = other.polyData.size();
    polyData.resize(size, 0);
    for (int i = 0; i < size; i++) {
        polyData[i] = other.polyData[i];
    }
    return *this;
}

polynomial polynomial::operator+(const polynomial &other) const {
    polynomial result;
    // Loop through the larger polynomial
    for (int i = 0; i < std::max(polyData.size(), other.polyData.size()); i++) {
        result.polyData[i] = polyData[i] + other.polyData[i];
    }
    return result;
}

polynomial polynomial::operator+(int val) const {
    // Use copy constructor to create a new polynomial
    polynomial result = *this;
    result.polyData[0] = polyData[0] + val;
    return result;
}

polynomial operator+(int val, const polynomial &other) {
    // Use copy constructor to create a new polynomial
    polynomial result = other;
    result.polyData[0] = other.polyData[0] + val;
    return result;
}

polynomial polynomial::operator*(const polynomial &other) const {
    polynomial result;
    const size_t len1 = polyData.size();
    const size_t len2 = other.polyData.size();
    result.polyData.resize(len1 + len2 - 1, 0);
    // Multiply each term in the first polynomial by each term in the second polynomial
    for (int i = 0; i < len1; i++) {
        for (int j = 0; j < len2; j++) {
            result.polyData[i + j] += polyData[i] * other.polyData[j];
        }
    }
    return result;
}

polynomial polynomial::operator*(int val) const {
    // Use copy constructor to create a new polynomial
    polynomial result = *this;
    for (int i = 0; i < polyData.size(); i++) {
        result.polyData[i] = polyData[i] * val;
    }
    return result;
}

polynomial operator*(int val, const polynomial &other) {
    // Use copy constructor to create a new polynomial
    polynomial result = other;
    for (int i = 0; i < other.polyData.size(); i++) {
        result.polyData[i] = other.polyData[i] * val;
    }
    return result;
}

polynomial polynomial::operator-(const polynomial &other) const {
    polynomial result;
    const size_t max_size = std::max(polyData.size(), other.polyData.size());
    result.polyData.resize(max_size, 0);

    // Subtract terms
    for (size_t i = 0; i < max_size; ++i) {
        const coeff term1 = i < polyData.size() ? polyData[i] : 0;
        const coeff term2 = i < other.polyData.size() ? other.polyData[i] : 0;
        result.polyData[i] = term1 - term2;
    }
    return result;
}

polynomial polynomial::operator%(const polynomial &divisor) const {
    polynomial dividend = *this;
    polynomial remainder = *this;
    int divisor_degree = divisor.find_degree_of();
    int dividend_degree = dividend.find_degree_of();

    while (dividend_degree >= divisor_degree) {
        int degree_diff = dividend_degree - divisor_degree;
        int coeff = dividend.polyData[dividend_degree] / divisor.polyData[divisor_degree];
        polynomial temp;
        temp.polyData.resize(degree_diff + 1, 0);
        temp.polyData[degree_diff] = coeff;

        remainder = remainder - (temp * divisor);
        dividend = remainder;
        dividend_degree = dividend.find_degree_of();
    }

    return remainder;
}

size_t polynomial::find_degree_of() const {
    for (int i = polyData.size() - 1; i >= 0; i--) {
        if (polyData[i] != 0) {
            return i;
        }
    }
    return 0;
}

std::vector<std::pair<power, coeff>> polynomial::canonical_form() const {
    std::vector<std::pair<power, coeff>> canonical;
    
    // Find the actual degree (highest non-zero coefficient)
    int actual_degree = -1;
    for (int i = polyData.size() - 1; i >= 0; i--) {
        if (polyData[i] != 0) {
            actual_degree = i;
            break;
        }
    }

    // If polynomial is zero
    if (actual_degree == -1) {
        return {std::make_pair(0, 0)};
    }

    // Build canonical form from highest to lowest degree
    canonical.reserve(actual_degree + 1);
    for (int i = actual_degree; i >= 0; i--) {
        if (polyData[i] != 0) {
            canonical.emplace_back(i, polyData[i]);
        }
    }
    
    return canonical;
}
