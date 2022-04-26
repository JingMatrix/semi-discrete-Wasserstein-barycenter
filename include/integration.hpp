#pragma once

#include "power-diagram.hpp"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

typedef std::map<Regular_triangulation::Weighted_point, double>
    Integral_power_diagram;

Integral_power_diagram area(Power_diagram &pd,
                            std::list<K::Point_2> &support);
Integral_power_diagram integral(Power_diagram &pd, gsl_monte_function &f);
