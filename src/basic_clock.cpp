#include "RcppArmadillo.h"
#include "sloth/updater.h"

using namespace Rcpp;

class tumbling_window_clock
{
private:
  sloth::updater::twclk_updater<double, double> tclk;

public:
  explicit tumbling_window_clock(double window_size)
      : tclk(window_size, true)
  {
  }

  NumericMatrix update(NumericVector t, NumericVector p, NumericVector v)
  {
    R_xlen_t count = 0;
    R_xlen_t npt = t.length();
    NumericMatrix m(npt, 9);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      tclk.insert(t[i], p[i], v[i]);
      if (tclk.bar.nbin)
      {
        m(count, 0) = tclk.bar.nbin;
        m(count, 1) = tclk.bar.ibin;
        m(count, 2) = tclk.bar.open;
        m(count, 3) = tclk.bar.high;
        m(count, 4) = tclk.bar.low;
        m(count, 5) = tclk.bar.close;
        m(count, 6) = tclk.bar.vol;
        m(count, 7) = tclk.bar.vwap;
        m(count, 8) = tclk.bar.time;
        ++count;
      }
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

public:
  sliding_window_clock(double window_size, int max_buffer_size)
      : w(window_size),
        ohlcv(max_buffer_size)
  {
  }

  // v: snapshot volume
  NumericMatrix update(NumericVector t, NumericVector p, NumericVector v)
  {
    R_xlen_t npt = t.length();
    NumericMatrix m(npt, 8);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      ohlcv.insert(t[i], p[i], v[i]);
      if ((t[i] - ohlcv.tbuf.back()) > w)
      {
        ohlcv.remove();
      }
      m(i, 0) = ohlcv.open;
      m(i, 1) = ohlcv.high;
      m(i, 2) = ohlcv.low;
      m(i, 3) = ohlcv.close;
      m(i, 4) = ohlcv.tnvr;
      m(i, 5) = ohlcv.vol;
      m(i, 6) = ohlcv.vwap;
      m(i, 7) = ohlcv.time;
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
      ohlcv.insert(t[i], p[i], v[i] - ohlcv.vol);
      if ((t[i] - ohlcv.tbuf.back()) > w)
      {
        ohlcv.remove();
      }
      m(i, 0) = ohlcv.open;
      m(i, 1) = ohlcv.high;
      m(i, 2) = ohlcv.low;
      m(i, 3) = ohlcv.close;
      m(i, 4) = ohlcv.tnvr;
      m(i, 5) = ohlcv.vol;
      m(i, 6) = ohlcv.vwap;
      m(i, 7) = ohlcv.time;
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

public:
  sliding_window_clock_market(int num_obs, double window_size, int max_buffer_size)
      : n(num_obs),
        w(window_size),
        ohlcv()
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
        ohlcv[i].insert(t[i], p[i], v[i]);
        if ((t[i] - ohlcv[i].tbuf.back()) > w)
        {
          ohlcv[i].remove();
        }
        m(i, 0) = ohlcv[i].open;
        m(i, 1) = ohlcv[i].high;
        m(i, 2) = ohlcv[i].low;
        m(i, 3) = ohlcv[i].close;
        m(i, 4) = ohlcv[i].tnvr;
        m(i, 5) = ohlcv[i].vol;
        m(i, 6) = ohlcv[i].vwap;
        m(i, 7) = ohlcv[i].time;
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
      const double dv = v[i] - ohlcv[i].vol;
      if (dv > 0)
      {
        ohlcv[i].insert(t[i], p[i], dv);
        if ((t[i] - ohlcv[i].tbuf.back()) > w)
        {
          ohlcv[i].remove();
        }
        m(i, 0) = ohlcv[i].open;
        m(i, 1) = ohlcv[i].high;
        m(i, 2) = ohlcv[i].low;
        m(i, 3) = ohlcv[i].close;
        m(i, 4) = ohlcv[i].tnvr;
        m(i, 5) = ohlcv[i].vol;
        m(i, 6) = ohlcv[i].vwap;
        m(i, 7) = ohlcv[i].time;
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

public:
  explicit volume_clock(double bin_size)
      : vclk(bin_size)
  {
  }

  NumericMatrix update(NumericVector p, NumericVector v)
  {
    R_xlen_t count = 0;
    R_xlen_t npt = p.length();
    NumericMatrix m(npt, 7);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      vclk.insert(p[i], v[i]);
      if (vclk.bar.nbin)
      {
        m(count, 0) = vclk.bar.nbin;
        m(count, 1) = vclk.bar.ibin;
        m(count, 2) = vclk.bar.open;
        m(count, 3) = vclk.bar.high;
        m(count, 4) = vclk.bar.low;
        m(count, 5) = vclk.bar.close;
        m(count, 6) = vclk.bar.vwap;
        ++count;
      }
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

public:
  explicit turnover_clock(double bin_size)
      : vclk(bin_size)
  {
  }

  NumericMatrix update(NumericVector p, NumericVector v)
  {
    R_xlen_t count = 0;
    R_xlen_t npt = p.length();
    NumericMatrix m(npt, 7);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      vclk.insert(p[i], v[i]);
      if (vclk.bar.nbin)
      {
        m(count, 0) = vclk.bar.nbin;
        m(count, 1) = vclk.bar.ibin;
        m(count, 2) = vclk.bar.open;
        m(count, 3) = vclk.bar.high;
        m(count, 4) = vclk.bar.low;
        m(count, 5) = vclk.bar.close;
        m(count, 6) = vclk.bar.vwap;
        ++count;
      }
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
      .method("update", &tumbling_window_clock::update, "Update state");

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
      .method("update", &volume_clock::update, "Update state");

  class_<turnover_clock>("turnover_clock")
      .constructor<double>()
      .method("update", &turnover_clock::update, "Update state");
}
