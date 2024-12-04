#include "poly.h"

polynomial::polynomial() {
    polyData[0] = 0;
}

polynomial::polynomial(const polynomial &other) : polyData(other.polyData) {}

void polynomial::print() const {
    auto terms = canonical_form();
    for (const auto& term : terms) {
        std::cout << term.second << "x^" << term.first << " ";
    }
    std::cout << std::endl;
}

polynomial &polynomial::operator=(const polynomial &other) {
    polyData = other.polyData;
    return *this;
}

polynomial polynomial::operator+(const polynomial &other) const {
    polynomial result = *this;
    for (const auto& term : other.polyData) {
        result.polyData[term.first] += term.second;
        if (result.polyData[term.first] == 0) {
            result.polyData.erase(term.first);
        }
    }
    return result;
}

polynomial polynomial::operator+(int val) const {
    polynomial result = *this;
    result.polyData[0] += val;
    if (result.polyData[0] == 0) {
        result.polyData.erase(0);
    }
    return result;
}

polynomial operator+(int val, const polynomial &other) {
    polynomial result = other;
    // Use .find() to check if power 0 exists, then add val
    auto it = other.polyData.find(0);
    result.polyData[0] = (it != other.polyData.end() ? it->second : 0) + val;
    return result;
}

polynomial polynomial::operator*(const polynomial &other) const {
    polynomial result;
    for (const auto& term1 : polyData) {
        for (const auto& term2 : other.polyData) {
            power new_power = term1.first + term2.first;
            coeff new_coeff = term1.second * term2.second;
            result.polyData[new_power] += new_coeff;
            if (result.polyData[new_power] == 0) {
                result.polyData.erase(new_power);
            }
        }
    }
    if (result.polyData.empty()) {
        result.polyData[0] = 0;
    }
    return result;
}

polynomial polynomial::operator*(int val) const {
    if (val == 0) {
        polynomial result;
        result.polyData[0] = 0;
        return result;
    }
    polynomial result;
    for (const auto& term : polyData) {
        result.polyData[term.first] = term.second * val;
    }
    return result;
}

polynomial operator*(int val, const polynomial &other) {
    if (val == 0) {
        polynomial result;
        return result;
    }
    polynomial result;
    for (const auto& term : other.polyData) {
        result.polyData[term.first] = term.second * val;
    }
    return result;
}

polynomial polynomial::operator-(const polynomial &other) const {
    polynomial result = *this;
    for (const auto& term : other.polyData) {
        result.polyData[term.first] -= term.second;
        if (result.polyData[term.first] == 0) {
            result.polyData.erase(term.first);
        }
    }
    if (result.polyData.empty()) {
        result.polyData[0] = 0;
    }
    return result;
}

polynomial polynomial::operator%(const polynomial &divisor) const {
    if (divisor.polyData.empty() || (divisor.polyData.size() == 1 && divisor.polyData.at(0) == 0)) {
        throw std::runtime_error("Division by zero polynomial");
    }

    polynomial remainder = *this;
    size_t divisor_degree = divisor.find_degree_of();

    while (!remainder.polyData.empty()) {
        size_t remainder_degree = remainder.find_degree_of();
        if (remainder_degree < divisor_degree) {
            break;
        }

        // Find the leading terms
        auto div_lead_term = divisor.polyData.find(divisor_degree);
        auto rem_lead_term = remainder.polyData.find(remainder_degree);
        if (div_lead_term == divisor.polyData.end() || rem_lead_term == remainder.polyData.end()) {
            break;
        }

        // Calculate the quotient term
        size_t new_power = remainder_degree - divisor_degree;
        coeff new_coeff = rem_lead_term->second / div_lead_term->second;

        // Create temporary polynomial for subtraction
        polynomial temp;
        temp.polyData[new_power] = new_coeff;

        // Subtract
        remainder = remainder - (temp * divisor);

        // Remove zero terms
        for (auto it = remainder.polyData.begin(); it != remainder.polyData.end();) {
            if (it->second == 0) {
                it = remainder.polyData.erase(it);
            } else {
                ++it;
            }
        }
    }

    if (remainder.polyData.empty()) {
        remainder.polyData[0] = 0;
    }
    return remainder;
}

size_t polynomial::find_degree_of() const {
    if (polyData.empty() || (polyData.size() == 1 && polyData.count(0) && polyData.at(0) == 0)) {
        return 0;
    }
    size_t max_degree = 0;
    for (const auto& term : polyData) {
        max_degree = std::max(max_degree, term.first);
    }
    return max_degree;
}

std::vector<std::pair<power, coeff>> polynomial::canonical_form() const {
    // Special case for zero polynomial
    if (polyData.empty() || 
        (polyData.size() == 1 && polyData.find(0) != polyData.end() && polyData.at(0) == 0)) {
        return {{0, 0}};  // Return {{0,0}} directly
    }
    
    std::vector<std::pair<power, coeff>> canonical;
    for (const auto& term : polyData) {
        if (term.second != 0) {
            canonical.emplace_back(term.first, term.second);
        }
    }
    
    std::sort(canonical.begin(), canonical.end(), 
        [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // If after filtering we have no terms, it's a zero polynomial
    if (canonical.empty()) {
        return {{0, 0}};
    }
    
    return canonical;
}
