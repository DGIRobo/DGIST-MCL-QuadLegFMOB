#include "filter.h"

double Ts = 0.001;
double pi = 3.141592;

double tustin_derivative(double input, double input_old, double output_old, double cutoff_freq)
{
    double time_const = 1 / (2 * pi * cutoff_freq);
    double output = 0;

    output = (2 * (input - input_old) - (Ts - 2 * time_const) * output_old) / (Ts + 2 * time_const);

    return output;
}

double tustin_integrate(double input, double input_old, double output_old)
{
  double output = 0;

  output = (Ts * (input + input_old) + 2 * output_old) / 2;

  return output;
}

double lowpassfilter(double input, double input_old, double output_old, double cutoff_freq)
{
    double time_const = 1 / (2 * pi * cutoff_freq);
    double output = 0;

    output = (Ts * (input + input_old) - (Ts - 2 * time_const) * output_old) / (Ts + 2 * time_const);

    return output;
}
