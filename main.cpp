#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <pybind11/complex.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>

using namespace std;
const double PI = std::acos(-1);

std::vector<std::complex<double>> dft(const std::vector<double>& input) {
	size_t N = input.size();
	std::vector<std::complex<double>> output(N);

	for (size_t k = 0; k < N; ++k) {
		std::complex<double> sum(0.0, 0.0);
		for (size_t n = 0; n < N; ++n) {
			double angle = -2.0 * PI * k * n / N;
			sum += input[n] * std::exp(std::complex<double>(0.0, angle));
		}
		output[k] = sum;
	}

	return output;
}

std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& input) {
	size_t N = input.size();
	std::vector<std::complex<double>> output(N);

	for (size_t n = 0; n < N; ++n) {
		std::complex<double> sum(0.0, 0.0);
		for (size_t k = 0; k < N; ++k) {
			double angle = 2.0 * PI * k * n / N;
			sum += input[k] * std::exp(std::complex<double>(0.0, angle));
		}
		output[n] = sum / static_cast<double>(N);
	}

	return output;
}

vector<double> sine(double A, double w, double x0, double xk, double sample_s) {
	vector<double> sygnal;
	vector<double> os_x;
	double krok = (xk - x0) / sample_s;
	for (double i = x0; i < xk; i += krok) {
		sygnal.push_back(A * sin(w * i));
		os_x.push_back(i);
	}
	matplot::plot(os_x, sygnal, "o");
	matplot::show();
	return sygnal;
}

vector<double> cosine(double A, double w, double x0, double xk, double sample_s) {
	vector<double> sygnal;
	vector<double> os_x;
	double krok = (xk - x0) / sample_s;
	for (double i = x0; i < xk; i += krok) {
		sygnal.push_back(A * cos(w * i));
		os_x.push_back(i);
	}
	matplot::plot(os_x, sygnal, "o");
	matplot::show();
	return sygnal;
}

vector<double> saw(double A, double freq, double x0, double xk, double sample_s) {
	vector<double> sygnal;
	vector<double> os_x;
	double krok = (xk - x0) / sample_s;
	double okres = 1 / freq;
	for (double i = x0;i < xk;i += krok) {
		double hel = fmod((i + okres / 2), okres) - okres / 2;
		sygnal.push_back(hel);
		os_x.push_back(i);

	}
	matplot::plot(os_x, sygnal, "o");
	matplot::show();
	return sygnal;


}

vector<double> rect(double A, double T, double fill, double x0, double xk, double sample_s) {
	double step = (xk - x0) / sample_s;
	vector<double> sygnal;
	vector<double> os_x;
	for (x0;x0 < xk;x0 += step) {
		double helper = fmod(x0, T);
		if (helper > T * fill)
			sygnal.push_back(A);
		if (helper < T * fill)
			sygnal.push_back(0);
		if (helper == T * fill || helper == 0)
			sygnal.push_back(A / 2);
		os_x.push_back(x0);

	}
	matplot::plot(os_x, sygnal, "o");
	matplot::show();
	return sygnal;
}

vector < double> OD(vector<double> sygnal, int deg, bool type) {
	vector<double> fsyg;
	vector<double> os_x;
	int mian = 0;
	for (int i = 0;i < sygnal.size();i++) {
		double val = 0;
		for (int j = 1;j <= deg;j++) {
			if (i + j < sygnal.size()) {
				val += sygnal[i + j];
				mian += 1;
			}
		}
		if (type) {
			for (int k = 1; k <= deg; k++)
				if (i - k >= 0) {
					val += sygnal[i - k];
					mian += 1;
				}
		}
		val /= mian;

		fsyg.push_back(val);
		os_x.push_back(i);

	}
	matplot::plot(os_x, fsyg, "o");
	matplot::show();
	return fsyg;
}

vector<double> DD(vector<double>sygnal, vector<double> filter) {
	vector<double> fsyg;
	vector<double> os_x;
	int deg = filter.size();
	double val = 0;
	if (deg % 2 == 0) {
		cout << "invalid filter size" << endl;
	}
	for (int i = 0; i < sygnal.size();i++) {
		val = 0;
		int h = i - deg / 2;
		for (int j = 0; j < filter.size();j++) {
			if (h >= 0 && h < sygnal.size()) {
				val += sygnal[h] * filter[j];
			}
			h++;
		}
		fsyg.push_back(val);
		os_x.push_back(i);

	}
	matplot::plot(os_x, fsyg, "o");
	matplot::show();
	return fsyg;
}

std::vector<double> corr(const std::vector<double>& x, const std::vector<double>& y) {
	int N = x.size();
	int M = y.size();
	int len = N + M - 1;

	std::vector<double> wynik(len, 0.0);
	vector<double> os_x;
	for (int lag = -(M - 1); lag < N; ++lag) {
		double suma = 0.0;
		for (int i = 0; i < M; ++i) {
			int j = i + lag;
			if (j >= 0 && j < N) {
				suma += x[j] * y[i];
			}
		}
		wynik[lag + M - 1] = suma;
		os_x.push_back(lag + M - 1);
	}
	matplot::plot(os_x, wynik, "o");
	matplot::show();
	return wynik;
}


PYBIND11_MODULE(signals, m) {
	m.def("sine", &sine, "A function that generates sine");
	m.def("cosine", &cosine, "A function that generates sine");
	m.def("saw", &saw, "A function that generates sawsine");
	m.def("rect", &rect, "A function that generates rectasine");
	m.def("OD", &OD, "Apply a 1D filter to a vector of doubles");
	m.def("dft", &dft, "Apply Fouriers transform");
	m.def("idft", &idft, "Apply reverse Fouriers transform");
	m.def("DD", &DD, "Apply a 2D filter with to a vector of doubles");
	m.def("corr", &corr, "Calculate correlation of two signals");
}


