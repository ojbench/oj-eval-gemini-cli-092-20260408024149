#ifndef SRC_HPP
#define SRC_HPP

#include "fraction.hpp"
#include <utility>

#define IGNORE_MATRIX 1

#if IGNORE_MATRIX

class matrix {
private:
    int m, n;
    fraction **data;

public:
    matrix() {
        m = n = 0;
        data = nullptr;
    }

    matrix(int m_, int n_) {
        m = m_;
        n = n_;
        if (m == 0 || n == 0) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; ++i) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; ++j) {
                    data[i][j] = fraction(0);
                }
            }
        }
    }

    matrix(const matrix &obj) {
        m = obj.m;
        n = obj.n;
        if (m == 0 || n == 0) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; ++i) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; ++j) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
    }

    matrix(matrix &&obj) noexcept {
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
    }

    ~matrix() {
        if (data != nullptr) {
            for (int i = 0; i < m; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }
    }

    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        if (data != nullptr) {
            for (int i = 0; i < m; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }
        m = obj.m;
        n = obj.n;
        if (m == 0 || n == 0) {
            data = nullptr;
        } else {
            data = new fraction*[m];
            for (int i = 0; i < m; ++i) {
                data[i] = new fraction[n];
                for (int j = 0; j < n; ++j) {
                    data[i][j] = obj.data[i][j];
                }
            }
        }
        return *this;
    }

    matrix &operator=(matrix &&obj) noexcept {
        if (this == &obj) return *this;
        if (data != nullptr) {
            for (int i = 0; i < m; ++i) {
                delete[] data[i];
            }
            delete[] data;
        }
        m = obj.m;
        n = obj.n;
        data = obj.data;
        obj.m = 0;
        obj.n = 0;
        obj.data = nullptr;
        return *this;
    }

    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n) {
            throw matrix_error();
        }
        return data[i - 1][j];
    }

    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n != rhs.m) {
            throw matrix_error();
        }
        matrix res(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; ++i) {
            for (int j = 0; j < rhs.n; ++j) {
                fraction sum(0);
                for (int k = 0; k < lhs.n; ++k) {
                    sum = sum + lhs.data[i][k] * rhs.data[k][j];
                }
                res.data[i][j] = sum;
            }
        }
        return res;
    }

    matrix transposition() {
        if (m == 0 || n == 0) {
            throw matrix_error();
        }
        matrix res(n, m);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                res.data[j][i] = data[i][j];
            }
        }
        return res;
    }

    fraction determination() {
        if (m != n || m == 0) {
            throw matrix_error();
        }
        matrix tmp(*this);
        fraction det(1);
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i; j < n; ++j) {
                if (!(tmp.data[j][i] == fraction(0))) {
                    pivot = j;
                    break;
                }
            }
            if (tmp.data[pivot][i] == fraction(0)) {
                return fraction(0);
            }
            if (pivot != i) {
                std::swap(tmp.data[i], tmp.data[pivot]);
                det = det * fraction(-1);
            }
            det = det * tmp.data[i][i];
            fraction inv = fraction(1) / tmp.data[i][i];
            for (int j = i + 1; j < n; ++j) {
                fraction factor = tmp.data[j][i] * inv;
                for (int k = i; k < n; ++k) {
                    tmp.data[j][k] = tmp.data[j][k] - factor * tmp.data[i][k];
                }
            }
        }
        return det;
    }
};

#endif

class resistive_network {
private:
    int interface_size, connection_size;
    matrix adjacency, conduction, laplace;

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        interface_size = interface_size_;
        connection_size = connection_size_;
        adjacency = matrix(connection_size, interface_size);
        conduction = matrix(connection_size, connection_size);
        for (int i = 0; i < connection_size; ++i) {
            adjacency(i + 1, from[i] - 1) = fraction(1);
            adjacency(i + 1, to[i] - 1) = fraction(-1);
            conduction(i + 1, i) = fraction(1) / resistance[i];
        }
        laplace = adjacency.transposition() * conduction * adjacency;
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) return fraction(0);
        int n = interface_size;
        matrix Mi(n - 1, n - 1);
        int r = 1;
        for (int i = 1; i <= n; ++i) {
            if (i == interface_id1) continue;
            int c = 0;
            for (int j = 1; j <= n; ++j) {
                if (j == interface_id1) continue;
                Mi(r, c) = laplace(i, j - 1);
                c++;
            }
            r++;
        }
        fraction detMi = Mi.determination();

        fraction detMij;
        if (n == 2) {
            detMij = fraction(1);
        } else {
            matrix Mij(n - 2, n - 2);
            r = 1;
            for (int i = 1; i <= n; ++i) {
                if (i == interface_id1 || i == interface_id2) continue;
                int c = 0;
                for (int j = 1; j <= n; ++j) {
                    if (j == interface_id1 || j == interface_id2) continue;
                    Mij(r, c) = laplace(i, j - 1);
                    c++;
                }
                r++;
            }
            detMij = Mij.determination();
        }

        return detMij / detMi;
    }

    fraction get_voltage(int id, fraction current[]) {
        if (id == interface_size) return fraction(0);
        int n = interface_size;
        matrix Mn(n - 1, n - 1);
        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - 1; ++j) {
                Mn(i, j) = laplace(i, j);
            }
        }
        fraction detMn = Mn.determination();

        matrix Mni(n - 1, n - 1);
        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < n - 1; ++j) {
                if (j == id - 1) {
                    Mni(i, j) = current[i - 1];
                } else {
                    Mni(i, j) = laplace(i, j);
                }
            }
        }
        fraction detMni = Mni.determination();

        return detMni / detMn;
    }

    fraction get_power(fraction voltage[]) {
        fraction power(0);
        for (int i = 0; i < connection_size; ++i) {
            fraction u = fraction(0);
            for (int j = 0; j < interface_size; ++j) {
                u = u + adjacency(i + 1, j) * voltage[j];
            }
            power = power + u * u * conduction(i + 1, i);
        }
        return power;
    }
};

#endif //SRC_HPP