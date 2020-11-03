#include "RcppArmadillo.h"
#include "sloth/updater.h"

using namespace Rcpp;

class tumbling_window_clock
{
private:
  sloth::updater::twclk_updater<double, double> tclk;
  double vlag;

  void _update(double t, double p, double v, NumericMatrix m, R_xlen_t& cnt)
  {
    tclk.insert(t, p, v);
      if (tclk.bar.nbin)
      {
        m(cnt, 0) = tclk.bar.nbin;
        m(cnt, 1) = tclk.bar.ibin;
        m(cnt, 2) = tclk.bar.open;
        m(cnt, 3) = tclk.bar.high;
        m(cnt, 4) = tclk.bar.low;
        m(cnt, 5) = tclk.bar.close;
        m(cnt, 6) = tclk.bar.vol;
        m(cnt, 7) = tclk.bar.vwap;
        m(cnt, 8) = tclk.bar.time;
        ++cnt;
      }
  }

public:
  explicit tumbling_window_clock(double window_size)
      : tclk(window_size, true),
        vlag(0.0)
  {
  }

  NumericMatrix update(NumericVector t, NumericVector p, NumericVector v)
  {
    R_xlen_t npt = t.length();
    NumericMatrix m(static_cast<R_xlen_t>((t[npt - 1] - t[0]) / tclk.bsize) + 1, 9);

    R_xlen_t count = 0;
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(t[i], p[i], v[i], m, count);
    }
    if (count)
    {
      return m(Range(0, count - 1), _);
    }
    else
    {
      return m(Range(0, 0), _);
    }
  }

  NumericMatrix update_tvol(NumericVector t, NumericVector p, NumericVector v)
  {
    R_xlen_t npt = t.length();
    NumericMatrix m(static_cast<R_xlen_t>((t[npt - 1] - t[0]) / tclk.bsize) + 1, 9);

    R_xlen_t count = 0;
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(t[i], p[i], v[i] - vlag, m, count);
      vlag = v[i];
    }
    if (count)
    {
      return m(Range(0, count - 1), _);
    }
    else
    {
      return m(Range(0, 0), _);
    }
  }
};

class sliding_window_clock
{
private:
  double w;
  sloth::updater::swclk_updater<double> ohlcv;
  double lagv;

  void _update(double t, double p, double v)
  {
    ohlcv.insert(t, p, v);
    while ((t - ohlcv.tbuf.back()) > w)
    {
      ohlcv.remove();
    }
  }

  void _assign(NumericMatrix m, R_xlen_t idx)
  {
      m(idx, 0) = ohlcv.open;
      m(idx, 1) = ohlcv.high;
      m(idx, 2) = ohlcv.low;
      m(idx, 3) = ohlcv.close;
      m(idx, 4) = ohlcv.tnvr;
      m(idx, 5) = ohlcv.vol;
      m(idx, 6) = ohlcv.vwap;
      m(idx, 7) = ohlcv.time;
  }

public:
  sliding_window_clock(double window_size, int max_buffer_size)
      : w(window_size),
        ohlcv(max_buffer_size),
        lagv(0.0)
  {
  }

  // v: snapshot volume
  NumericMatrix update(NumericVector t, NumericVector p, NumericVector v)
  {
    R_xlen_t npt = t.length();
    NumericMatrix m(npt, 8);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(t[i], p[i], v[i]);
      _assign(m, i);
    }
    return m;
  }

  // v: total volume
  NumericMatrix update_tvol(NumericMatrix t, NumericMatrix p, NumericVector v)
  {
    R_xlen_t npt = t.length();
    NumericMatrix m(npt, 8);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(t[i], p[i], v[i] - lagv);
      lagv = v[i];
      _assign(m, i);
    }
    return m;
  }
};

class sliding_window_clock_market
{
private:
  size_t n;
  double w;
  std::vector<sloth::updater::swclk_updater<double>> ohlcv;
  std::vector<double> lagv;

  void _update(double t, double p, double v, R_xlen_t idx)
  {
    ohlcv[idx].insert(t, p, v);
    while ((t - ohlcv[idx].tbuf.back()) > w)
    {
      ohlcv[idx].remove();
    }
  }

  void _assign(NumericMatrix m, R_xlen_t idx)
  {
        m(idx, 0) = ohlcv[idx].open;
        m(idx, 1) = ohlcv[idx].high;
        m(idx, 2) = ohlcv[idx].low;
        m(idx, 3) = ohlcv[idx].close;
        m(idx, 4) = ohlcv[idx].tnvr;
        m(idx, 5) = ohlcv[idx].vol;
        m(idx, 6) = ohlcv[idx].vwap;
        m(idx, 7) = ohlcv[idx].time;
  }

public:
  sliding_window_clock_market(int num_obs, double window_size, int max_buffer_size)
      : n(num_obs),
        w(window_size),
        ohlcv(),
        lagv(n, 0.0)
  {
    ohlcv.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
      ohlcv.emplace_back(max_buffer_size);
    }
  }

  // v: snapshot volume
  bool update(NumericVector t, NumericVector p, NumericVector v, NumericMatrix m)
  {
    bool cnt = false;
    for (R_xlen_t i = 0; i < n; ++i)
    {
      if (t[i] > ohlcv[i].time)
      {
        _update(t[i], p[i], v[i], i);
        _assign(m, i);
        cnt = true;
      }
    }
    return cnt;
  }

  // v: total volume
  bool update_tvol(NumericVector t, NumericVector p, NumericVector v, NumericMatrix m)
  {
    bool cnt = false;
    for (R_xlen_t i = 0; i < n; ++i)
    {
      const double dv = v[i] - lagv[i];
      lagv[i] = v[i];
      if (dv > 0)
      {
        _update(t[i], p[i], dv, i);
        _assign(m, i);
        cnt = true;
      }
    }
    return cnt;
  }
};

class volume_clock
{
private:
  sloth::updater::volclk_updater vclk;
  double vlag;

  void _update(double p, double v, NumericMatrix m, R_xlen_t& cnt)
  {
    vclk.insert(p, v);
    if (vclk.bar.nbin)
    {
        m(cnt, 0) = vclk.bar.nbin;
        m(cnt, 1) = vclk.bar.ibin;
        m(cnt, 2) = vclk.bar.open;
        m(cnt, 3) = vclk.bar.high;
        m(cnt, 4) = vclk.bar.low;
        m(cnt, 5) = vclk.bar.close;
        m(cnt, 6) = vclk.bar.vwap;
        ++cnt;
    }
  }

public:
  explicit volume_clock(double bin_size)
      : vclk(bin_size),
        vlag(0.0)
  {
  }

  NumericMatrix update(NumericVector p, NumericVector v)
  {
    R_xlen_t count = 0;
    double tot = std::accumulate(v.begin(), v.end(), vclk.res.vol);
    NumericMatrix m(static_cast<R_xlen_t>(tot / vclk.bsize) + 1, 7);
    for (R_xlen_t i = 0; i < p.length(); ++i)
    {
      _update(p[i], v[i], m, count);
    }
    if (count)
    {
      return m(Range(0, count - 1), _);
    }
    else
    {
      return m(Range(0, 0), _);
    }
  }

  NumericMatrix update_tvol(NumericVector p, NumericVector v)
  {
    R_xlen_t count = 0;
    NumericMatrix m(static_cast<R_xlen_t>((v[v.length() - 1] + vclk.res.vol) / vclk.bsize) + 1, 7);
    for (R_xlen_t i = 0; i < p.length(); ++i)
    {
      _update(p[i], v[i] - vlag, m, count);
      vlag = v[i];
    }
    if (count)
    {
      return m(Range(0, count - 1), _);
    }
    else
    {
      return m(Range(0, 0), _);
    }
  }
};

class turnover_clock
{
private:
  sloth::updater::tnvrclk_updater vclk;
  double vlag;

  void _update(double p, double v, NumericMatrix m, R_xlen_t& cnt)
  {
    vclk.insert(p, v);
    if (vclk.bar.nbin)
    {
        m(cnt, 0) = vclk.bar.nbin;
        m(cnt, 1) = vclk.bar.ibin;
        m(cnt, 2) = vclk.bar.open;
        m(cnt, 3) = vclk.bar.high;
        m(cnt, 4) = vclk.bar.low;
        m(cnt, 5) = vclk.bar.close;
        m(cnt, 6) = vclk.bar.vwap;
        ++cnt;
    }
  }

public:
  explicit turnover_clock(double bin_size)
      : vclk(bin_size),
        vlag(0.0)
  {
  }

  NumericMatrix update(NumericVector p, NumericVector v)
  {
    R_xlen_t count = 0;
    auto pv = p * v;
    auto tot = std::accumulate(pv.begin(), pv.end(), vclk.res.tnvr);
    NumericMatrix m(static_cast<R_xlen_t>(tot / vclk.bsize) + 1, 7);
    for (R_xlen_t i = 0; i < p.length(); ++i)
    {
      _update(p[i], v[i], m, count);
    }
    if (count)
    {
      return m(Range(0, count - 1), _);
    }
    else
    {
      return m(Range(0, 0), _);
    }
  }

  NumericMatrix update_tvol(NumericVector p, NumericVector v)
  {
    R_xlen_t count = 0;
    R_xlen_t npt = p.length();
    NumericMatrix m(npt, 7);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(p[i], v[i] - vlag, m, count);
      vlag = v[i];
    }
    if (count)
    {
      return m(Range(0, count - 1), _);
    }
    else
    {
      return m(Range(0, 0), _);
    }
  }
};

RCPP_MODULE(basic_clock)
{
  using namespace Rcpp;

  class_<tumbling_window_clock>("tumbling_window_clock")
      .constructor<double>()
      .method("update", &tumbling_window_clock::update, "Update state")
      .method("update_tvol", &tumbling_window_clock::update_tvol, "Update state, cumulative volume");

  class_<sliding_window_clock>("sliding_window_clock")
      .constructor<double, int>()
      .method("update", &sliding_window_clock::update, "Update state")
      .method("update_tvol", &sliding_window_clock::update_tvol, "Update state, cumulative volume");

  class_<sliding_window_clock_market>("sliding_window_clock_market")
      .constructor<int, double, int>()
      .method("update", &sliding_window_clock_market::update, "Update state")
      .method("update_tvol", &sliding_window_clock_market::update_tvol, "Update state, cumulative volume");

  class_<volume_clock>("volume_clock")
      .constructor<double>()
      .method("update", &volume_clock::update, "Update state")
      .method("update_tvol", &volume_clock::update_tvol, "Update state, cumulative volume");

  class_<turnover_clock>("turnover_clock")
      .constructor<double>()
      .method("update", &turnover_clock::update, "Update state")
      .method("update_tvol", &turnover_clock::update_tvol, "Update state, cumulative volume");
}
