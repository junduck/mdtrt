#pragma once

#include "RcppArmadillo.h"
#include "boost/circular_buffer.hpp"

#include <cmath>
#include <algorithm>
#include <array>
//#include <boost/circular_buffer.hpp>
//#include "include/armadillo"

namespace sloth
{
  namespace updater
  {
    // Kahan compensated sum
    inline void csum(double &s, double x, double &c)
    {
      double y = x - c;
      volatile double tmp = s + y;
      c = (tmp - s) - y;
      s = tmp;
    }

    template <bool kahan>
    struct accumulator
    {
    };

    template<>
    struct accumulator<false>
    {
      void sum(double &s, double x)
      {
        s += x;
      }
    };

    template<>
    struct accumulator<true>
    {
      double c;

      accumulator()
        : c(0.0)
      {
      }

      void sum(double &s, double x)
      {
        double y = x - c;
        volatile double tmp = s + x;
        c = (tmp - s) - y;
        s = tmp;
      }
    };

  } // namespace updater

} // namespace sloth
