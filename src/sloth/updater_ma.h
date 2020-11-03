#pragma once

#include "updater_common.h"

namespace sloth
{
  namespace updater
  {

    template <bool kahan>
    struct msum_updater
    {
    };

    template <>
    struct msum_updater<false>
    {
      double s;

      msum_updater()
          : s(0.0)
      {
      }

      void insert(double x)
      {
        s += x;
      }

      void remove(double x)
      {
        s -= x;
      }

      void roll(double x, double y)
      {
        s += (x - y);
      }
    };

    template <>
    struct msum_updater<true>
    {
      double s, sc;

      msum_updater()
          : s(0.0),
            sc(0.0)
      {
      }

      void insert(double x)
      {
        csum(s, x, sc);
      }

      void remove(double x)
      {
        csum(s, -x, sc);
      }

      void roll(double x, double y)
      {
        csum(s, x - y, sc);
      }
    };

    template <bool kahan>
    struct sma_updater
    {
    };

    template <>
    struct sma_updater<false>
    {
      size_t n;
      double m;

      sma_updater()
          : n(0),
            m(0.0)
      {
      }

      void insert(double x)
      {
        double d = x - m;
        n += 1;
        m += d / n;
      }

      void remove(double x)
      {
        double d = x - m;
        n -= 1;
        m -= d / n;
      }

      void roll(double x, double y)
      {
        double d = x - y;
        m += d / n;
      }
    };

    template <>
    struct sma_updater<true>
    {
      size_t n;
      double m, cm;

      sma_updater()
          : n(0),
            m(0.0),
            cm(0.0)
      {
      }

      void insert(double x)
      {
        double d = x - m;
        n += 1;
        csum(m, d / n, cm);
      }

      void remove(double x)
      {
        double d = m - x;
        n -= 1;
        csum(m, d / n, cm);
      }

      void roll(double x, double y)
      {
        double d = x - y;
        csum(m, d / n, cm);
      }
    };

    template <bool kahan>
    struct wma_updater
    {
    };

    template <>
    struct wma_updater<false>
    {
      double xsum, wsum;
      size_t n, w;

      wma_updater()
          : xsum(0.0),
            wsum(0.0),
            n(0),
            w(0)
      {
      }

      void insert(double x)
      {
        n += 1;
        w += n;
        wsum += x * n;
        xsum += x;
      }

      void remove(double x)
      {
        n -= 1;
        w -= n;
        wsum -= xsum;
        xsum -= x;
      }

      void roll(double x, double y)
      {
        wsum += x * n - xsum;
        xsum += x - y;
      }
    };

    template <>
    struct wma_updater<true>
    {
      double xsum, wsum, xsumc, wsumc;
      size_t n, w;

      wma_updater()
          : xsum(0.0),
            wsum(0.0),
            xsumc(0.0),
            wsumc(0.0),
            n(0),
            w(0)
      {
      }

      void insert(double x)
      {
        n += 1;
        w += n;
        csum(wsum, x * n, wsumc);
        csum(xsum, x, xsumc);
      }

      void remove(double x)
      {
        n -= 1;
        w -= n;
        csum(wsum, -xsum, wsumc);
        csum(xsum, -x, xsumc);
      }

      void roll(double x, double y)
      {
        csum(wsum, x * n - xsum, wsumc);
        csum(xsum, x - y, xsumc);
      }
    };

    template <bool kahan>
    struct wsma_updater
    {
    };

    template <>
    struct wsma_updater<false>
    {
      double m, wsum;

      wsma_updater()
          : m(0.0),
            wsum(0.0)
      {
      }

      void insert(double x, double w)
      {
        double d = x - m;
        wsum += w;
        m += (w / wsum) * d;
      }

      void remove(double x, double w)
      {
        double d = x - m;
        wsum -= w;
        m -= (w / wsum) * d;
      }

      void roll(double x, double wx, double y, double wy)
      {
        insert(x, wx);
        remove(y, wy);
      }
    };

    template <>
    struct wsma_updater<true>
    {
      double m, wsum, mc, wsumc;

      wsma_updater()
          : m(0.0),
            wsum(0.0),
            mc(0.0),
            wsumc(0.0)
      {
      }

      void insert(double x, double w)
      {
        double d = x - m;
        csum(wsum, w, wsumc);
        csum(m, (w / wsum) * d, mc);
      }

      void remove(double x, double w)
      {
        double d = m - x;
        csum(wsum, -w, wsumc);
        csum(m, (w / wsum) * d, mc);
      }

      void roll(double x, double wx, double y, double wy)
      {
        insert(x, wx);
        remove(y, wy);
      }
    };

    template <bool kahan>
    struct ema_updater
    {
    };

    template <>
    struct ema_updater<false>
    {
      bool init;
      double a, m;

      explicit ema_updater(double alpha)
          : init(false),
            a(alpha),
            m(0.0)
      {
      }

      void insert(double x)
      {
        if (init)
        {
          m += (x - m) * a;
        }
        else
        {
          m = x;
          init = true;
        }
      }
    };

    template <>
    struct ema_updater<true>
    {
      bool init;
      double a, m, mc;

      explicit ema_updater(double alpha)
          : init(false),
            a(alpha),
            m(0.0),
            mc(0.0)
      {
      }

      void insert(double x)
      {
        if (init)
        {
          csum(m, (x - m) * a, mc);
        }
        else
        {
          m = x;
          init = true;
        }
      }
    };

    template <size_t order, bool kahan>
    struct xema_updater
    {
    };

    template <size_t order>
    struct xema_updater<order, false>
    {
      bool init;
      double a;
      std::array<double, order> m;

      explicit xema_updater(double alpha)
          : init(false),
            a(alpha)
      {
        m.fill(0.0);
      }

      void insert(double x)
      {
        if (init)
        {
          m[0] += (x - m[0]) * a;
          for (size_t i = 1; i < order; ++i)
          {
            m[i] += (m[i - 1] - m[i]) * a;
          }
        }
        else
        {
          m.fill(x);
          init = true;
        }
      }
    };

    template <size_t order>
    struct xema_updater<order, true>
    {
      bool init;
      double a;
      std::array<double, order> m, mc;

      explicit xema_updater(double alpha)
          : init(false),
            a(alpha)
      {
        m.fill(0.0);
        mc.fill(0.0);
      }

      void insert(double x)
      {
        if (init)
        {
          csum(m[0], (x - m[0]) * a, mc[0]);
          for (size_t i = 1; i < order; ++i)
          {
            csum(m[i], (m[i - 1] - m[i]) * a, mc[i]);
          }
        }
        else
        {
          m.fill(x);
          init = true;
        }
      }
    };

  } // namespace updater

} // namespace sloth
