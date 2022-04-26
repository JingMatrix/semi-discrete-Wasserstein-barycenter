#include "integration.hpp"
#include "power-diagram.hpp"

int main() {
  std::list<K::Point_2> bbox = {K::Point_2(0, 0), K::Point_2(1, 0),
                                K::Point_2(1, 1), K::Point_2(0, 1)};

  CGAL::IO::set_pretty_mode(std::cout);
  CGAL::IO::set_pretty_mode(std::cerr);
  auto rt = triangulation_from_data("data/weight_points");
  auto pd = generate_power_diagram(rt);
  Integral_power_diagram data_area = area(pd, bbox);
  double check_area = 0;
  for (auto cit = data_area.begin(); cit != data_area.end(); ++cit) {
    /* std::cout << cit->first << " with cell area: " << cit->second <<
     * std::endl; */
    check_area += cit->second;
  }
  std::cout << "In total, they have area " << check_area << "." << std::endl;
  /* rt.draw_dual(std::cout); */
  return 0;
}

double exact = 1.3932039296856768591842462603255;

double g(double *k, size_t dim, void *params) {
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);
  double A = 1.0 / (M_PI * M_PI * M_PI);
  return A / (1.0 - cos(k[0]) * cos(k[1]) * cos(k[2]));
}

void display_results(const char *title, double result, double error) {
  printf("%s ==================\n", title);
  printf("result = % .6f\n", result);
  printf("sigma  = % .6f\n", error);
  printf("exact  = % .6f\n", exact);
  printf("error  = % .6f = %.2g sigma\n", result - exact,
         fabs(result - exact) / error);
}

void test_integration(void) {
  double res, err;

  double xl[3] = {0, 0, 0};
  double xu[3] = {M_PI, M_PI, M_PI};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G = {&g, 3, 0};

  size_t calls = 500000;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc(3);
    gsl_monte_miser_integrate(&G, xl, xu, 3, calls, r, s, &res, &err);
    gsl_monte_miser_free(s);

    display_results("miser", res, err);
  }

  gsl_rng_free(r);
}

/* int main() { test_integration(); } */
