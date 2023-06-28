
#include "spline.h"
#include <cassert>


namespace Spline {
    struct dmatrix {
        std::vector<double> data;
        int rows, cols;

        dmatrix() : rows(0), cols(0) {}
        dmatrix(int rows, int cols ) : rows(rows), cols(cols), data(rows*cols, 0) {}

        inline double& operator()(int row, int col) {
            return data[row + col*rows];
        }
        inline double operator()(int row, int col) const {
            return data[row + col*rows];
        }
    };
    struct dvector {
        std::vector<double> data;

        dvector() {}
        dvector(const std::vector<double>& d) : data(d) {}
        dvector(int n) : data(n, 0) {}

        inline std::size_t length() const {
            return data.size();
        }
        inline double& operator()(int i) {
            return data[i];
        }
        inline const double& operator()(int i) const {
            return data[i];
        }
    };
    struct cvector {
        std::vector<std::complex<double>> data;

        cvector() {}
        cvector(const std::vector<std::complex<double>>& d) : data(d) {}
        cvector(int n) : data(n, 0) {}

        inline std::size_t length() const {
            return data.size();
        }
        inline std::complex<double>& operator()(int i) {
            return data[i];
        }
        inline const std::complex<double>& operator()(int i) const {
            return data[i];
        }
    };

    static dvector operator*(const dmatrix& m, const dvector& x) {
        assert(x.length() == m.cols);

        dvector b(m.rows);
        for (int i = 0; i < m.rows; i++) {
            b(i) = 0;
            for (int j = 0; j < m.cols; j++) {
                b(i) += m(i, j)*x(j);
            }
        }
        return b;
    }
    static cvector operator*(const dmatrix& m, const cvector& x) {
        assert(x.length() == m.cols);

        cvector b(m.rows);
        for (int i = 0; i < m.rows; i++) {
            b(i) = 0;
            for (int j = 0; j < m.cols; j++) {
                b(i) += m(i, j)*x(j);
            }
        }
        return b;
    }
    static dmatrix operator*(const dmatrix& m, const dmatrix& n) {
        assert(m.cols == n.rows);

        dmatrix b(m.rows, n.cols);
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < n.cols; j++) {
                b(i, j) = 0;
                for (int k = 0; k < m.cols; k++) {
                    b(i,j) += m(i,k)*n(k,j);
                }
            }
        }
        return b;
    }
    static void diff(const cvector& a, cvector& b) {
        for (int i = 0; i < a.length() - 1; i++) {
            b(i) = a(i+1) - a(i);
        }
    }
    static void diff(const dvector& a, dvector& b) {
        for (int i = 0; i < a.length() - 1; i++) {
            b(i) = a(i+1) - a(i);
        }
    }


    static void Thomas(const dmatrix& A, const cvector& r, cvector& x) {
        size_t N = A.rows;
        std::vector<double> a(N, 0);
        std::vector<double> b(N, 0);
        std::vector<double> c(N, 0);
        x = cvector(N);

        for (int i = 0; i < N; i++) {
            if (i > 0)
                a[i] = A(i,i-1); 
            
            b[i] = A(i,i);

            if (i < N-1)
                c[i] = A(i,i+1); 
        }

        thomas_algorithm(a, b, c, r.data, x.data);
    }










    void Cubic::Initialize(const std::vector<double>& x, const std::vector<std::complex<double>>& y) {
        assert(x.size() == y.size());
        
        std::size_t n = x.size();

        dvector Dx(n-1);				// derivative of x
        cvector yp(n-1);				// derivative of y (dy/dx)
        cvector r(n-2);     			// r-side

        dmatrix Tsp(n - 2, n - 2);			// sparse matrix
        //Tsp.reserve(Eigen::VectorXi::Constant(n - 2, 3));		//three nonzero entries in each column


        diff(x, Dx);
        diff(y, yp);
        for (int i = 0; i < yp.length(); ++i)
            yp(i) /= Dx(i);

        for(auto i = 1;i<n-3;++i) {
            Tsp(i, i - 1)   = 1./Dx(i + 1);
            Tsp(i, i)       = 2. * (1./Dx(i) + 1./Dx(i + 1));
            Tsp(i, i + 1)   = 1./Dx(i);

            r(i) = 3. * (yp(i) / Dx(i + 1) + yp(i + 1) / Dx(i));
        }

        // not-a-knot computation of slopes
        Tsp(0, 0) =  1./Dx(0) + Dx(0)/Dx(1)/Dx(1) + 2./Dx(1);
        Tsp(0, 1) =  1./Dx(1) + Dx(0)/Dx(1)/Dx(1);

        Tsp(n - 3, n - 3) =  2./Dx(n-3) + Dx(n-2)/Dx(n-3)/Dx(n-3) + 1./Dx(n-2);
        Tsp(n - 3, n - 4) = Dx(n-2) / Dx(n-3)/Dx(n-3) + 1./Dx(n-3);

        r(0) = yp(0) / Dx(0) + 3. * yp(1) / Dx(1) + 2.*yp(1)*Dx(0)/Dx(1)/Dx(1);
        r(n - 3) = yp(n-2) / Dx(n-2) + 3. * yp(n-3) / Dx(n-3) + 2.*yp(n-3)*Dx(n-2) / Dx(n-3) / Dx(n-3);


        cvector stilde(n-2);

        Thomas(Tsp, r, stilde);

        // first and last slopes
        std::complex<double> s0 = -stilde(0) + 2.*yp(0);
        s0 = 2.*yp(0) - stilde(0) + (Dx(0)*Dx(0)/Dx(1)/Dx(1)) * (stilde(0) + stilde(1)  - 2.*yp(1)); 

        std::complex<double> sn = -stilde(n - 3) + 2. * yp(n - 2);
        sn = 2. * yp(n - 2) - stilde(n - 3) + (Dx(n - 2) * Dx(n - 2) / Dx(n - 3) / Dx(n - 3)) * (stilde(n - 4) + stilde(n - 3) - 2. * yp(n - 3));

        cvector s(n);
        for (auto i = 1; i < n - 1; ++i)
            s(i) = stilde(i - 1);
        s(0) = s0;
        s(n - 1) = sn;


        // vectors of coefficients a,b,c,d
        a = std::vector<std::complex<double>>(&y[0], &y[n - 1]);
        b = std::vector<std::complex<double>>(&s(0), &s(n - 1));
        c = std::vector<std::complex<double>>(n - 1);
        d = std::vector<std::complex<double>>(n - 1);
        for (auto i = 0; i < c.size(); ++i) {
            c[i] = (yp(i) - s(i)) / Dx(i); // I can't remember why this is right
            d[i] = (s(i + 1) + s(i) - 2. * yp(i)) / (Dx(i) * Dx(i));
        }
    }

    std::complex<double> Cubic::operator ()(double x_s, const std::vector<double>& x) const {
        std::vector<double>::const_iterator it;
        it = std::lower_bound(x.begin(), x.end(), x_s);
        int idx = (std::max)(int(it - x.begin()) - 1, 0);
        idx = (std::min)(idx,int(x.size()-2));

        std::complex<double> y_s = d[idx] * (x_s - x[idx + 1]) + c[idx];
        y_s = y_s * (x_s - x[idx]) + b[idx];
        y_s = y_s * (x_s - x[idx]) + a[idx];

        return y_s;
    }





    
    void thomas_algorithm(  const std::vector<double>& a,
                            const std::vector<double>& b,
                            const std::vector<double>& c,
                            const std::vector<std::complex<double>>& d,
                            std::vector<std::complex<double>>& f) {
        size_t N = d.size();
        std::vector<double> c_star(N, 0.0);
        std::vector<std::complex<double>> d_star(N, 0.0);
        c_star[0] = c[0] / b[0];
        d_star[0] = d[0] / b[0];
        for (int i=1; i<N; i++) {
            double m = 1.0 / (b[i] - a[i] * c_star[i-1]);
            c_star[i] = c[i] * m;
            d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
        }
        f[N-1] = d_star[N-1];
        for (int i=N-2; i >= 0; i--) {
            f[i] = d_star[i] - c_star[i] * f[i+1];
        }
    }
}