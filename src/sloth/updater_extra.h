#pragma once

#include "updater_common.h"

namespace sloth
{
  namespace updater
  {
    struct s4_updater
    {
      size_t n, n2, n3;
      double m, s2, s3, s4;

      s4_updater()
        : n(0), n2(0), n3(0),
          m(0.0), s2(0.0), s3(0.0), s4(0.0)
      {
      }

      void insert(double x)
      {
        n += 1;
        double d = x - m;
        double dn = d / n;
        double dn2 = dn * dn;
        double d_dn_n1 = d * dn * (n - 1);
        s4 += d_dn_n1 * dn2 * (n * n - 3 * n + 3) + 6 * dn2 * s2 - 4 * dn * s3;
        s3 += d_dn_n1 * dn * (n - 2) - 3 * dn * s2;
        s2 += d_dn_n1;
        m += dn;
      }

      void roll(double x, double y)
      {
        n2 = n * n;
        n3 = n2 * n;
        double d = x - m;
        double d2 = d * d;
        double d3 = d2 * d;
        double d0 = y - m;
        double d02 = d0 * d0;
        double d03 = d02 * d0;
        double dd0 = d - d0;
        double dd02 = dd0 * dd0;
        s4 += d2 * d2 - d02 * d02 
          - 4.0 * dd0 * (s3 + d3 - d03) / n 
          + 6.0 * (s2 + d2 - d02) * dd02 / n2 
          - 3.0 * dd02 * dd02 / n3;
        s3 += d3 - d03 
          - 3.0 * dd0 * (s2 + d2 - d02) / n 
          + 2.0 * dd02 * dd0 / n2;
        s2 += d2 - d02 - dd02 / n;
        m += (x - y) / n;
      }

      std::array<double, 4> moments() const
      {
        return {m, s2 / n, s3 / n, s4 / n};
      }

      std::array<double, 4> statistics() const
      {
        auto y = moments();
        y[3] = y[3] / std::pow(y[1], 2) - 3;
        y[2] = std::sqrt(n * (n - 1)) / (n - 2) * y[2] / std::pow(y[1], 1.5); // G1
        y[1] = std::sqrt(s2 / (n - 1));
        return y;
      }
    };
    
    template <class T>
    struct minmax_updater
    {
      boost::circular_buffer<T> mindeq, maxdeq;

      explicit minmax_updater(int deque_size)
          : mindeq(deque_size),
            maxdeq(deque_size)
      {
      }

      void insert(T x)
      {
        while (!mindeq.empty() && mindeq.back() > x)
          mindeq.pop_back();
        while (!maxdeq.empty() && maxdeq.back() < x)
          maxdeq.pop_back();
        mindeq.push_back(x);
        maxdeq.push_back(x);
      }

      void remove(T x)
      {
        if (mindeq.front() == x)
          mindeq.pop_front();
        if (maxdeq.front() == x)
          maxdeq.pop_front();
      }
    };

    template <class T>
    struct argminmax_updater
    {
      boost::circular_buffer<T> mindeq, maxdeq;
      boost::circular_buffer<int> minidx, maxidx;

      explicit argminmax_updater(int deque_size)
          : mindeq(deque_size),
            maxdeq(deque_size),
            minidx(deque_size),
            maxidx(deque_size)
      {
      }

      void insert(T x)
      {
        while (!mindeq.empty() && mindeq.back() > x)
        {
          mindeq.pop_back();
          minidx.pop_back();
        }
        while (!maxdeq.empty() && maxdeq.back() < x)
        {
          maxdeq.pop_back();
          minidx.pop_back();
        }
        for (auto &it : minidx)
          it += 1;
        for (auto &it : maxidx)
          it += 1;
        mindeq.push_back(x);
        maxdeq.push_back(x);
        minidx.push_back(0);
        maxidx.push_back(0);
      }

      void remove(T x)
      {
        if (mindeq.front() == x)
        {
          mindeq.pop_front();
          minidx.pop_front();
        }
        if (maxdeq.front() == x)
        {
          maxdeq.pop_front();
          maxidx.pop_front();
        }
      }
    };

    struct psquare_updater
    {
      bool init;
      char n;
      std::array<double, 5> h, mpos, dpos, incr;

      explicit psquare_updater(double proba)
          : init(false),
            n(0),
            mpos{1.0, 2.0, 3.0, 4.0, 5.0},
            dpos{1.0, 1.0 + 2.0 * proba, 1.0 + 4.0 * proba, 3.0 + 2.0 * proba, 5.0},
            incr{0.0, proba / 2.0, proba, (1.0 + proba) / 2.0, 1.0}
      {
      }

      void insert(double x)
      {
        if (!init)
        {
          h[n] = x;
          n += 1;
          if (n == 5)
          {
            init = true;
            std::sort(h.begin(), h.end());
          }
        }
        else
        {
          int k;
          if (x < h[0])
          {
            k = 0;
            h[0] = x;
          }
          else if (x < h[1])
          {
            k = 0;
          }
          else if (x < h[2])
          {
            k = 1;
          }
          else if (x < h[3])
          {
            k = 2;
          }
          else if (x <= h[4])
          {
            k = 3;
          }
          else
          {
            k = 3;
            h[4] = x;
          }
          for (auto i = k + 1; i < 5; ++i)
          {
            mpos[i] += 1.0;
          }
          for (auto i = 0; i < 5; ++i)
          {
            dpos[i] += incr[i];
          }
          for (auto i = 1; i < 4; ++i)
          {
            const double delta_dpos = dpos[i] - mpos[i];
            if ((delta_dpos >= 1.0 && mpos[i + 1] - mpos[i] > 1.0) ||
                (delta_dpos <= -1.0 && mpos[i - 1] - mpos[i] < -1.0))
            {
              const int isgn = delta_dpos < 0 ? -1 : 1;
              const double sgn = static_cast<double>(isgn);
              auto t1 = sgn / (mpos[i + 1] - mpos[i - 1]);
              auto t2 = (mpos[i] - mpos[i - 1] + sgn) * (h[i + 1] - h[i]) / (mpos[i + 1] - mpos[i]);
              auto t3 = (mpos[i + 1] - mpos[i] - sgn) * (h[i] - h[i - 1]) / (mpos[i] - mpos[i - 1]);
              auto adj = h[i] + t1 * (t2 + t3);
              if (h[i - 1] < adj && h[i + 1] > adj)
              {
                h[i] = adj;
              }
              else
              {
                h[i] += sgn * (h[i + isgn] - h[i]) / (mpos[i + isgn] - mpos[i]);
              }
              mpos[i] += sgn;
            }
          }
        }
      }
    };

    // Naiive impl
    struct rls_updater
    {
      arma::colvec w;
      arma::mat P;
      double lamb, ilamb;

      rls_updater(int size, double lambda = 0.98, double sigma = 1.0)
          : w(size, arma::fill::zeros),
            P(size, size, arma::fill::zeros),
            lamb(lambda),
            ilamb(1.0 / lambda)
      {
        P.diag().fill(sigma);
      }

      void fit(const arma::rowvec &x, double d)
      {
        double a = d - arma::dot(x, w);
        arma::mat Px{P * x.t()};
        arma::colvec g{Px / (lamb + arma::dot(x, Px))};
        P *= ilamb * (1 - arma::dot(g, x));
        w += a * g;
      }

      double predict(const arma::rowvec &x)
      {
        return arma::dot(x, w);
      }
    };

    // Naiive impl
    struct nlms_updater
    {
      double _mu, _eps;
      arma::colvec w;

      nlms_updater(int size, double mu, double eps = 1e-11)
          : _mu(mu),
            _eps(eps),
            w(size, arma::fill::zeros)
      {
      }

      void fit(const arma::rowvec &x, double d)
      {
        double a = d - arma::dot(x, w);
        w += _mu * a * x.t() / (arma::dot(x, x) + _eps);
      }

      double predict(const arma::rowvec &x)
      {
        return arma::dot(x, w);
      }
    };

    // Naiive impl
    struct kmeans_updater
    {
      size_t k, n;
      double hl;
      arma::mat c;
      arma::colvec d;

      kmeans_updater(int ncluster, int ndim, double halflife)
          : k(ncluster),
            n(0),
            hl(halflife),
            c(k, ndim, arma::fill::zeros),
            d(k, arma::fill::zeros)
      {
      }

      kmeans_updater(int ncluster, int ndim, int halflife_period)
          : kmeans_updater(ncluster, ndim, 2.0 / (halflife_period + 1))
      {
      }

      size_t update(const arma::rowvec &x)
      {
        size_t idx = n;
        if (n < k)
        {
          arma::rowvec tmp(x.n_cols, arma::fill::randu);
          c.row(n) = x % tmp;
          n += 1;
        }
        else
        {
          idx = predict(x);
          c.row(idx) += hl * (x - c.row(idx));
        }
        return idx;
      }

      size_t predict(const arma::rowvec &x)
      {
        for (size_t i = 0; i < n; ++i)
        {
          d(i) = arma::norm(c.row(i) - x);
        }
        return d.index_min();
      }
    };

  } // namespace updater

} // namespace sloth
