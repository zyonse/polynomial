#include "poly.h"

polynomial::polynomial() : polyData(1, 0) {}

polynomial::polynomial(const polynomial &other) {
    for (int i = 0; i < other.polyData.size(); i++) {
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
    for (int i = 0; i < other.polyData.size(); i++) {
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
    result.polyData.resize(polyData.size() + other.polyData.size() - 1, 0);
    // Multiply each term in the first polynomial by each term in the second polynomial
    for (int i = 0; i < polyData.size(); i++) {
        for (int j = 0; j < other.polyData.size(); j++) {
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

size_t polynomial::find_degree_of() {
    for (int i = polyData.size() - 1; i >= 0; i--) {
        if (polyData[i] != 0) {
            return i;
        }
    }
    return 0;
}

std::vector<std::pair<power, coeff>> polynomial::canonical_form() const {
    std::vector<std::pair<power, coeff>> canonical;
    for (int i = polyData.size() - 1; i >= 0; i--) {
        if (polyData[i] != 0) {
            canonical.push_back(std::make_pair(i, polyData[i]));
        }
    }
    return canonical;
}
