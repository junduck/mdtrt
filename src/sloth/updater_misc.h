#pragma once

#include "updater_common.h"

namespace sloth
{
  namespace updater
  {

    struct lindecay_updater
    {
      double a;
      double decay;

      explicit lindecay_updater(int period)
          : a(period),
            decay(0.0)
      {
      }

      void update(double x)
      {
        decay = std::max(std::max(x, decay - a), 0.0);
      }
    };

    struct expdecay_updater
    {
      double a;
      double decay;

      explicit expdecay_updater(int period)
          : a(1.0 - 1.0 / period),
            decay(0.0)
      {
      }

      void update(double x)
      {
        decay = std::max(x, decay * a);
      }
    };

    template <class T>
    struct lag_updater
    {
      T last;

      explicit lag_updater(const T &init)
          : T(init)
      {
      }

      T update(const T &x)
      {
        T tmp{last};
        last = x;
        return tmp;
      }

      void update_byref(T &x)
      {
        std::swap(x, last);
      }
    };

    template <class T>
    struct diff_updater
    {
      T last;

      explicit diff_updater(const T &init)
          : T(init)
      {
      }

      T update(const T &x)
      {
        T d{x - last};
        last = x;
        return d;
      }

      void update_byref(T &x)
      {
        T tmp{x};
        x -= last;
        std::swap(tmp, last);
      }
    };

  } // namespace updater

} // namespace sloth
