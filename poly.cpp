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

void polynomial::fft(std::vector<std::complex<double>> &a, bool invert, int depth) const {
    size_t n = a.size();
    if (n <= 1)
        return;

    std::vector<std::complex<double>> a0(n / 2), a1(n / 2);
    for (size_t i = 0; 2 * i < n; ++i) {
        a0[i] = a[2 * i];
        a1[i] = a[2 * i + 1];
    }

    if (active_threads.load() < max_threads) {
        active_threads += 2;
        std::thread t1(&polynomial::fft, this, std::ref(a0), invert, depth + 1);
        std::thread t2(&polynomial::fft, this, std::ref(a1), invert, depth + 1);
        t1.join();
        t2.join();
        active_threads -= 2;
    } else {
        fft(a0, invert, depth + 1);
        fft(a1, invert, depth + 1);
    }

    double ang = 2 * M_PI / n * (invert ? -1 : 1);
    std::complex<double> w(1), wn(cos(ang), sin(ang));
    for (size_t i = 0; 2 * i < n; ++i) {
        std::complex<double> u = a0[i];
        std::complex<double> v = w * a1[i];
        a[i] = u + v;
        a[i + n / 2] = u - v;
        if (invert) {
            a[i] /= 2;
            a[i + n / 2] /= 2;
        }
        w *= wn;
    }
}

polynomial polynomial::operator*(const polynomial &other) const {
    // For sparse polynomials with few terms, use direct multiplication
    if (polyData.size() * other.polyData.size() < 1000) {
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

    // For dense polynomials, use FFT
    size_t n = 1;
    size_t deg1 = find_degree_of();
    size_t deg2 = other.find_degree_of();
    while (n < deg1 + deg2 + 1) n <<= 1;

    auto fa = to_fft_format(polyData, n);
    auto fb = to_fft_format(other.polyData, n);

    fft(fa, false, 0);
    fft(fb, false, 0);

    for (size_t i = 0; i < n; i++) {
        fa[i] *= fb[i];
    }

    fft(fa, true, 0);

    polynomial result;
    for (size_t i = 0; i < n; i++) {
        int rounded = static_cast<int>(std::round(fa[i].real()));
        if (rounded != 0) {
            result.polyData[i] = rounded;
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
    std::vector<std::pair<power, coeff>> canonical;
    if (polyData.empty() || (polyData.size() == 1 && polyData.count(0) && polyData.at(0) == 0)) {
        canonical.emplace_back(0, 0);
        return canonical;
    }
    
    for (const auto& term : polyData) {
        if (term.second != 0) {
            canonical.emplace_back(term.first, term.second);
        }
    }
    
    std::sort(canonical.begin(), canonical.end(), 
        [](const auto& a, const auto& b) { return a.first > b.first; });
    
    return canonical;
}

// Convert polynomial to vector for FFT
std::vector<std::complex<double>> polynomial::to_fft_format(const std::unordered_map<power, coeff>& data, size_t size) {
    std::vector<std::complex<double>> result(size);
    for (const auto& term : data) {
        if (term.first < size) {
            result[term.first] = std::complex<double>(term.second, 0);
        }
    }
    return result;
}
