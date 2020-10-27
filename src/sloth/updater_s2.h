#pragma once

#include "updater_common.h"

namespace sloth
{
  namespace updater
  {
    template <bool kahan>
    struct s2_updater
    {
    };

    template <>
    struct s2_updater<false>
    {
      size_t n;
      double m, s2;

      s2_updater()
          : n(0),
            m(0.0),
            s2(0.0)
      {
      }

      void insert(double x)
      {
        double d = x - m;
        n += 1;
        m += d / n;
        s2 += (x - m) * d;
      }

      void remove(double x)
      {
        double d = x - m;
        n -= 1;
        m -= d / n;
        s2 -= (x - m) * d;
      }
    };

    template <>
    struct s2_updater<true>
    {
      size_t n;
      double m, s2;
      double cm, cs2;

      s2_updater()
          : n(0),
            m(0.0),
            s2(0.0),
            cm(0.0),
            cs2(0.0)
      {
      }

      void insert(double x)
      {
        double d = x - m;
        n += 1;
        csum(m, d / n, cm);
        csum(s2, (x - m) * d, cs2);
      }

      void remove(double x)
      {
        double d = m - x;
        n -= 1;
        csum(m, d / n, cm);
        csum(s2, (x - m) * d, cs2);
      }
    };

    // West, D. H. D. (1979). Updating mean and variance estimates: An improved method. Communications of the ACM, 22(9), 532-535.
    template <bool kahan>
    struct s2weighted_updater
    {
    };

    template <>
    struct s2weighted_updater<false>
    {
      size_t n;
      double wsum, w2sum, m, s2;

      s2weighted_updater()
          : n(0),
            wsum(0.0),
            w2sum(0.0),
            m(0.0),
            s2(0.0)
      {
      }

      void insert(double x, double w)
      {
        double d = x - m;
        wsum += w;
        w2sum += w * w;
        m += (w / wsum) * d;
        s2 += w * (x - m) * d;
      }

      void remove(double x, double w)
      {
        double d = x - m;
        wsum -= w;
        w2sum -= w * w;
        m -= (w / wsum) * d;
        s2 -= w * (x - m) * d;
      }
    };

    template <>
    struct s2weighted_updater<true>
    {
      double wsum, w2sum, m, s2;
      double cwsum, cw2sum, cm, cs2;

      s2weighted_updater()
          : wsum(0.0),
            w2sum(0.0),
            m(0.0),
            s2(0.0),
            cwsum(0.0),
            cw2sum(0.0),
            cm(0.0),
            cs2(0.0)
      {
      }

      void insert(double x, double w)
      {
        double d = x - m;
        csum(wsum, w, cwsum);
        csum(w2sum, w * w, cw2sum);
        csum(m, (w / wsum) * d, cm);
        csum(s2, w * (x - m) * d, cs2);
      }

      void remove(double x, double w)
      {
        double d = m - x;
        csum(wsum, -w, cwsum);
        csum(w2sum, -w * w, cw2sum);
        csum(m, (w / wsum) * d, cm);
        csum(s2, w * (x - m) * d, cs2);
      }
    };

    template <bool kahan>
    struct s2ew_updater
    {
    };

    template <>
    struct s2ew_updater<false>
    {
      bool init;
      double a, b, m, s2;

      explicit s2ew_updater(int period)
          : s2ew_updater(2.0 / (period + 1))
      {
      }

      explicit s2ew_updater(double alpha)
          : init(false),
            a(alpha),
            b(1.0 - a),
            m(0.0),
            s2(0.0)
      {
      }

      void insert(double x)
      {
        if (init)
        {
          double d = x - m;
          double i = a * d;
          m += i;
          s2 = b * (s2 + i * d);
        }
        else
        {
          m = x;
          init = true;
        }
      }

      void remove(double x)
      {
      }
    };

    template <>
    struct s2ew_updater<true>
    {
      bool init;
      double a, b, m, s2, mc; // , s2c;

      explicit s2ew_updater(int period)
          : s2ew_updater(2.0 / (period + 1))
      {
      }

      explicit s2ew_updater(double alpha)
          : init(false),
            a(alpha),
            b(1.0 - a),
            m(0.0),
            s2(0.0),
            mc(0.0) // ,
                    // s2c(0.0)
      {
      }

      void insert(double x)
      {
        if (init)
        {
          double d = x - m;
          double i = a * d;
          csum(m, i, mc);
          // TODO: Does compensation even scale with b?
          // csum(s2, i * d, s2c);
          // s2 *= b;
          // s2c *= b;
          s2 = b * (s2 + i * d);
        }
        else
        {
          m = x;
          init = true;
        }
      }

      void remove(double x)
      {
      }
    };

    // Naiive impl
    struct cov_updater
    {
      size_t n;
      arma::rowvec d, m;
      arma::mat s2;

      explicit cov_updater(int nvar)
          : n(0),
            d(nvar, arma::fill::zeros),
            m(nvar, arma::fill::zeros),
            s2(nvar, nvar, arma::fill::zeros)
      {
      }

      void insert(const arma::rowvec &x)
      {
        d = x - m;
        n += 1;
        m += d / n;
        s2 += (x - m).t() * d;
      }

      void remove(const arma::rowvec &x)
      {
        d = x - m;
        n -= 1;
        m -= d / n;
        s2 -= (x - m).t() * d;
      }
    };

    struct ewcov_updater
    {
      bool init;
      double a, b;
      arma::rowvec m;
      arma::mat s2;

      ewcov_updater(int nvar, double alpha)
          : init(false),
            a(alpha),
            b(1 - a),
            m(nvar, arma::fill::zeros),
            s2(nvar, nvar, arma::fill::zeros)
      {
      }

      ewcov_updater(int nvar, int period)
          : ewcov_updater(nvar, 2.0 / (period + 1))
      {
      }

      void insert(const arma::rowvec &x)
      {
        if (init)
        {
          arma::rowvec d{x - m};
          arma::rowvec i{a * d};
          m += i;
          s2 += i.t() * d;
          s2 *= b;
        }
        else
        {
          m = x;
          init = true;
        }
      }

      void remove(const arma::rowvec &x)
      {
      }
    };

  } // namespace updater

} // namespace sloth
