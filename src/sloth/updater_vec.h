#pragma once

#include "updater_common.h"

namespace sloth
{
  namespace updater
  {

    // adapted updater algo to support vector (multiple variables) input for convenience

    // s2

    struct s2_vec_updater
    {
      size_t n;
      arma::colvec m, s2;

      explicit s2_vec_updater(int nvar)
          : n(0),
            m(nvar, arma::fill::zeros),
            s2(nvar, arma::fill::zeros)
      {
      }

      void insert(const arma::colvec &x)
      {
        arma::colvec d = x - m;
        n += 1;
        m += d / n;
        s2 += (x - m) % d;
      }

      void remove(const arma::colvec &x)
      {
        arma::colvec d = x - m;
        n -= 1;
        m -= d / n;
        s2 -= (x - m) % d;
      }
    };

    struct s2weighted_vec_updater
    {
      arma::colvec wsum, w2sum, m, s2;

      void insert(const arma::colvec &x, const arma::colvec &w)
      {
        arma::colvec d = x - m;
        wsum += w;
        w2sum += w % w;
        m += (w / wsum) % d;
        s2 += w % (x - m) % d;
      }

      void remove(const arma::colvec &x, const arma::colvec &w)
      {
        arma::colvec d = x - m;
        wsum -= w;
        w2sum -= w % w;
        m -= (w / wsum) % d;
        s2 -= w % (x - m) % d;
      }
    };

    struct s2ew_vec_updater
    {
      bool init;
      double a, b;
      arma::colvec m, s2;

      s2ew_vec_updater(int nvar, int period)
          : s2ew_vec_updater(nvar, 2.0 / (period + 1))
      {
      }

      s2ew_vec_updater(int nvar, double alpha)
          : init(false),
            a(alpha),
            b(1.0 - a),
            m(nvar, arma::fill::zeros),
            s2(nvar, arma::fill::zeros)
      {
      }

      void insert(const arma::colvec &x)
      {
        if (init)
        {
          arma::colvec d = x - m;
          arma::colvec i = a * d;
          m += i;
          s2 += i % d;
          s2 *= b;
        }
        else
        {
          m = x;
          init = true;
        }
      }

      void remove(const arma::colvec &x)
      {
      }
    };

    // ma

    struct sma_vec_updater
    {
      size_t n;
      arma::colvec m;

      explicit sma_vec_updater(int nvar)
          : n(0),
            m(nvar, arma::fill::zeros)
      {
      }

      void insert(const arma::colvec &x)
      {
        n += 1;
        m += (x - m) / n;
      }

      void remove(const arma::colvec &x)
      {
        n -= 1;
        m -= (x - m) / n;
      }

      void roll(const arma::colvec &x, const arma::colvec &y)
      {
        m += (x - y) / n;
      }
    };

    struct wma_vec_updater
    {
      size_t n, w;
      arma::colvec xsum, wsum;

      explicit wma_vec_updater(int nvar)
          : n(0),
            w(0),
            xsum(nvar, arma::fill::zeros),
            wsum(nvar, arma::fill::zeros)
      {
      }

      void insert(const arma::colvec &x)
      {
        n += 1;
        w += n;
        wsum += x * n;
        xsum += x;
      }

      void remove(const arma::colvec &x)
      {
        n -= 1;
        w -= n;
        wsum -= xsum;
        xsum -= x;
      }

      void roll(const arma::colvec &x, const arma::colvec &y)
      {
        wsum += x * n - xsum;
        xsum += x - y;
      }
    };

    struct ema_vec_updater
    {
      bool init;
      double a;
      arma::colvec m;

      ema_vec_updater(int nvar, int period)
          : ema_vec_updater(nvar, 2.0 / (period + 1))
      {
      }

      ema_vec_updater(int nvar, double alpha)
          : init(false),
            a(alpha),
            m(nvar, arma::fill::zeros)
      {
      }

      void insert(const arma::colvec &x)
      {
        if (init)
        {
          m += a * (x - m);
        }
        else
        {
          m = x;
          init = true;
        }
      }
    };

    // misc



  } // namespace updater

} // namespace sloth
