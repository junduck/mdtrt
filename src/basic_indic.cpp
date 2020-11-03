#include "RcppArmadillo.h"
#include "boost/circular_buffer.hpp"
#include "sloth/updater.h"

using namespace Rcpp;

class zlema
{
private:
  size_t lag;
  boost::circular_buffer<double> xbuf;
  sloth::updater::ema_updater<false> _ema;

  double _update(double x)
  {
    double datum;
    if (xbuf.size() == lag)
    {
      datum = x + (x - xbuf.back());
    }
    else
    {
      datum = x;
    }
    xbuf.push_front(x);
    _ema.insert(datum);
    return _ema.m;
  }

public:
  zlema(int period)
      : lag((period - 1) / 2),
        xbuf(lag),
        _ema(2.0 / (period - 1))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class wma
{
private:
  size_t n;
  boost::circular_buffer<double> xbuf;
  sloth::updater::wma_updater<false> _wma;

  double _update(double x)
  {
    if (xbuf.size() < n)
    {
      _wma.insert(x);
    }
    else
    {
      _wma.roll(x, xbuf.back());
    }
    xbuf.push_front(x);
    return _wma.wsum / _wma.w;
  }

public:
  friend class hma;

  wma(int period)
      : n(period),
        xbuf(n),
        _wma()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class willr
{
private:
  size_t n;
  boost::circular_buffer<double> hbuf, lbuf;
  sloth::updater::minmax_updater<double> minmax;

  double _update(double high, double low, double close)
  {
    minmax.insert(high);
    minmax.insert(low);
    if (hbuf.size() == n)
    {
      minmax.remove(hbuf.back());
      minmax.remove(lbuf.back());
    }
    hbuf.push_front(high);
    lbuf.push_front(low);
    double min = minmax.mindeq.front();
    double max = minmax.maxdeq.front();
    return -100.0 * ((min - close) / (max - min));
  }

public:
  willr(int period)
      : n(period),
        hbuf(period),
        lbuf(period),
        minmax(period)
  {
  }

  NumericVector run(NumericVector high, NumericVector low, NumericVector close)
  {
    R_xlen_t npt = high.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(high[i], low[i], close[i]);
    }
    return val;
  }

  void run_inplace(NumericVector high, NumericVector low, NumericVector close, NumericVector val)
  {
    for (R_xlen_t i = 0; i < high.length(); ++i)
    {
      val[i] = _update(high[i], low[i], close[i]);
    }
  }
};

class wilders
{
private:
  sloth::updater::ema_updater<false> _ema;

  double _update(double x)
  {
    _ema.insert(x);
    return _ema.m;
  }

public:
  wilders(int period)
      : _ema(1.0 / period)
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class wad
{
private:
  double preclose, _wad;

  double _update(double high, double low, double close)
  {
    double ad = 0.0;
    if (close > preclose)
    {
      ad = close - std::min(preclose, low);
    }
    else if (close < preclose)
    {
      ad = close - std::max(preclose, high);
    }
    _wad += ad;
    return _wad;
  }

public:
  wad()
      : preclose(0.0),
        _wad(0.0)
  {
  }

  NumericVector run(NumericVector high, NumericVector low, NumericVector close)
  {
    R_xlen_t npt = high.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(high[i], low[i], close[i]);
    }
    return val;
  }

  void run_inplace(NumericVector high, NumericVector low, NumericVector close, NumericVector val)
  {
    for (R_xlen_t i = 0; i < high.length(); ++i)
    {
      val[i] = _update(high[i], low[i], close[i]);
    }
  }
};

class vwma
{
private:
  size_t n;
  boost::circular_buffer<double> cbuf, vbuf;
  sloth::updater::wsma_updater<false> _vwma;

  double _update(double close, double vol)
  {
    _vwma.insert(close, vol);
    if (cbuf.size() == n)
    {
      _vwma.remove(cbuf.back(), vbuf.back());
    }
    cbuf.push_front(close);
    vbuf.push_front(vol);
    return _vwma.m;
  }

public:
  vwma(int period)
      : n(period),
        cbuf(n),
        vbuf(n),
        _vwma()
  {
  }

  NumericVector run(NumericVector close, NumericVector vol)
  {
    R_xlen_t npt = close.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(close[i], vol[i]);
    }
    return val;
  }

  void run_inplace(NumericVector close, NumericVector vol, NumericVector val)
  {
    for (R_xlen_t i = 0; i < close.length(); ++i)
    {
      val[i] = _update(close[i], vol[i]);
    }
  }
};

class smosc
{
private:
  size_t fast_n, slow_n;
  boost::circular_buffer<double> fast_buf, slow_buf;
  sloth::updater::sma_updater<false> fast, slow;

  double _update(double x)
  {
    fast.insert(x);
    slow.insert(x);
    if (slow_buf.size() == slow_n)
    {
      fast.roll(x, fast_buf.back());
      slow.roll(x, slow_buf.back());
    }
    else if (fast_buf.size() == fast_n)
    {
      fast.roll(x, fast_buf.back());
      slow.insert(x);
    }
    else
    {
      fast.insert(x);
      slow.insert(x);
    }
    fast_buf.push_back(x);
    slow_buf.push_back(x);
    return (fast.m - slow.m) / slow.m;
  }

public:
  smosc(int fast_period, int slow_period)
      : fast_n(std::min(fast_period, slow_period)),
        slow_n(std::max(fast_period, slow_period)),
        fast_buf(fast_n),
        slow_buf(slow_n),
        fast(),
        slow()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class emosc
{
private:
  int fast_n, slow_n;
  sloth::updater::ema_updater<false> fast, slow;

  double _update(double x)
  {
    fast.insert(x);
    slow.insert(x);
    return (fast.m - slow.m) / slow.m;
  }

public:
  emosc(int fast_period, int slow_period)
      : fast_n(std::min(fast_period, slow_period)),
        slow_n(std::max(fast_period, slow_period)),
        fast(2.0 / (fast_n + 1)),
        slow(2.0 / (slow_n + 1))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class volatility_hisotry
{
private:
  size_t n;
  double ann;
  boost::circular_buffer<double> xbuf;
  sloth::updater::s2_updater<false> _s2;

  double _update(double x)
  {
    auto _n = xbuf.size();
    if (_n)
    {
      _s2.insert(x / xbuf.front() - 1.0);
    }
    if (_n == n)
    {
      _s2.remove(xbuf.back());
    }
    xbuf.push_front(x);
    return ann * std::sqrt(_s2.s2 / _s2.n);
  }

public:
  volatility_hisotry(int period, int annualise = 250)
      : n(period),
        ann(std::sqrt(annualise)),
        xbuf(n),
        _s2()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class vidya
{
private:
  int fperiod, speriod;
  double factor, _vidya;
  boost::circular_buffer<double> fbuf, sbuf;
  sloth::updater::s2_updater<false> fast, slow;

  double _update(double x)
  {
    fast.insert(x);
    slow.insert(x);
    if (sbuf.full())
    {
      fast.remove(fbuf.back());
      slow.remove(sbuf.back());
    }
    else if (fbuf.full())
    {
      fast.remove(fbuf.back());
    }
    fbuf.push_front(x);
    sbuf.push_front(x);
    if (fbuf.empty())
    {
      _vidya = x;
    }
    else
    {
      const double a = factor * std::sqrt(fast.s2 / slow.s2);
      _vidya += a * (x - _vidya);
    }
    return _vidya;
  }

public:
  vidya(int fast_period, int slow_period)
      : fperiod(std::min(fast_period, slow_period)),
        speriod(std::max(fast_period, slow_period)),
        factor(std::sqrt(static_cast<double>(speriod) / fperiod)),
        _vidya(0.0),
        fbuf(fperiod),
        sbuf(speriod),
        fast(),
        slow()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class vhf
{
private:
  boost::circular_buffer<double> xbuf, dbuf;
  sloth::updater::minmax_updater<double> minmax;
  sloth::updater::sma_updater<false> sma;

  double _update(double x)
  {
    double d = 0.0;
    if (!xbuf.empty())
    {
      d = std::abs(x - xbuf.front());
    }
    minmax.insert(x);
    if (xbuf.full())
    {
      minmax.remove(xbuf.back());
      sma.roll(d, dbuf.back());
    }
    else
    {
      sma.insert(d);
    }
    xbuf.push_front(x);
    dbuf.push_front(d);
    return (minmax.maxdeq.front() - minmax.mindeq.front()) / sma.m;
  }

public:
  vhf(int period)
      : xbuf(period),
        dbuf(period),
        minmax(period),
        sma()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class variance
{
  bool roll;
  boost::circular_buffer<double> xbuf;
  sloth::updater::s2_updater<false> _s2;

  double _update(double x)
  {
    _s2.insert(x);
    if (roll)
    {
      if (xbuf.full())
      {
        _s2.remove(xbuf.back());
      }
      xbuf.push_front(x);
    }
    return _s2.s2 / _s2.n;
  }

public:
  variance(int period)
      : roll(period),
        xbuf(period),
        _s2()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class sma
{
private:
  bool roll;
  boost::circular_buffer<double> xbuf;
  sloth::updater::sma_updater<false> _sma;

  double _update(double x)
  {
    if (roll)
    {
      if (xbuf.full())
      {
        _sma.roll(x, xbuf.back());
      }
      else
      {
        _sma.insert(x);
      }
      xbuf.push_front(x);
    }
    else
    {
      _sma.insert(x);
    }
    return _sma.m;
  }

public:
  sma(int period)
      : roll(period),
        xbuf(period),
        _sma()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class rsi
{
private:
  bool init;
  double alpha, last, _ema_up, _ema_down;

  double _update(double x)
  {
    double up = 0.0;
    double down = 0.0;
    if (init)
    {
      if (x > last)
      {
        up = x - last;
      }
      else
      {
        down = last - x;
      }
      _ema_up += alpha * (up - _ema_up);
      _ema_down += alpha * (down - _ema_down);
    }
    else
    {
      init = true;
    }
    last = x;
    return 100.0 - (100.0 / (1.0 + _ema_up / _ema_down));
  }

public:
  explicit rsi(int period)
      : init(false),
        alpha(1.0 / period),
        last(0.0),
        _ema_up(1e-10),
        _ema_down(1e-10)
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

// lagged rate

// lagged diff

class rollsum
{
private:
  boost::circular_buffer<double> xbuf;
  sloth::updater::msum_updater<false> xsum;

  double _update(double x)
  {
    if (xbuf.full())
    {
      xsum.roll(x, xbuf.back());
    }
    else
    {
      xsum.insert(x);
    }
    xbuf.push_front(x);
    return xsum.s;
  }

public:
  friend class cmo;

  rollsum(int period)
      : xbuf(period)
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class volidx
{
private:
  bool init;
  double _pvi, _nvi;
  double lv, lp;

  void _update(double p, double v)
  {
    if (init)
    {
      if (lv > v)
      {
        _pvi += _pvi * ((p - lp) / lp);
      }
      else if (lv < v)
      {
        _nvi += _nvi * ((p - lp) / lp);
      }
    }
    else
    {
      init = true;
    }
    lp = p;
    lv = v;
  }

public:
  volidx(int period, double initialise = 1000.0)
      : init(false),
        _pvi(initialise),
        _nvi(initialise),
        lv(0.0),
        lp(0.0)
  {
  }

  List run(NumericVector p, NumericVector v)
  {
    R_xlen_t npt = p.length();
    NumericVector pvi(npt);
    NumericVector nvi(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(p[i], v[i]);
      pvi[i] = _pvi;
      nvi[i] = _nvi;
    }
    return List::create(
        Named("pvi") = pvi,
        Named("nvi") = nvi);
  }

  void run_inplace(NumericVector p, NumericVector v, NumericVector pvi, NumericVector nvi)
  {
    R_xlen_t npt = p.length();
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(p[i], v[i]);
      pvi[i] = _pvi;
      nvi[i] = _nvi;
    }
  }
};

// PSAR: one pass algo!

class obv
{
private:
  bool init;
  double _obv, lp;

  double _update(double p, double v)
  {
    if (init)
    {
      if (p = lp)
      {
        _obv = 0.0;
      }
      else if (p < lp)
      {
        _obv -= v;
      }
      else
      {
        _obv += v;
      }
    }
    else
    {
      init = true;
    }
    lp = p;
    return _obv;
  }

public:
  obv()
      : init(false),
        _obv(0.0),
        lp(0.0)
  {
  }

  NumericVector run(NumericVector p, NumericVector v)
  {
    R_xlen_t npt = p.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(p[i], v[i]);
    }
    return val;
  }

  void run_inplace(NumericVector p, NumericVector v, NumericVector val)
  {
    for (R_xlen_t i = 0; i < p.length(); ++i)
    {
      val[i] = _update(p[i], v[i]);
    }
  }
};

class minmax
{
private:
  boost::circular_buffer<double> xbuf;
  sloth::updater::minmax_updater<double> mm;

  void _update(double x)
  {
    mm.insert(x);
    if (xbuf.full())
    {
      mm.remove(xbuf.back());
    }
    xbuf.push_front(x);
  }

public:
  minmax(int period)
      : xbuf(period),
        mm(period)
  {
  }

  List run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector min(npt);
    NumericVector max(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(x[i]);
      min[i] = mm.mindeq.front();
      max[i] = mm.maxdeq.front();
    }
    return List::create(
        Named("min") = min,
        Named("max") = max);
  }

  void run_inplace(NumericVector x, NumericVector min, NumericVector max)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      _update(x[i]);
      min[i] = mm.mindeq.front();
      max[i] = mm.maxdeq.front();
    }
  }
};

class mfi
{
private:
  bool init;
  boost::circular_buffer<double> up, down;
  double uptnvr, downtnvr, lp;

  double _update(double p, double v)
  {
    double tnvr = p * v;
    if (!init)
    {
      lp = p;
      init = true;
    }
    if (up.full())
    {
      if (p > lp)
      {
        uptnvr += tnvr - up.back();
        up.push_front(tnvr);
        down.push_front(0.0);
      }
      else if (p < lp)
      {
        downtnvr += tnvr - down.back();
        up.push_front(0.0);
        down.push_front(tnvr);
      }
      else
      {
        up.push_front(0.0);
        down.push_front(0.0);
      }
    }
    else
    {
      if (p > lp)
      {
        uptnvr += tnvr;
        up.push_front(tnvr);
        down.push_front(0.0);
      }
      else if (p < lp)
      {
        downtnvr += tnvr;
        up.push_front(0.0);
        down.push_front(tnvr);
      }
      else
      {
        up.push_front(0.0);
        down.push_front(0.0);
      }
    }
    lp = p;
    return 100.0 - (100.0 / (1.0 + uptnvr / downtnvr));
  }

public:
  mfi(int period)
      : init(false),
        up(period),
        down(period),
        uptnvr(0.0),
        downtnvr(0.0),
        lp(0.0)
  {
  }

  NumericVector run(NumericVector p, NumericVector v)
  {
    R_xlen_t npt = p.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(p[i], v[i]);
    }
    return val;
  }

  void run_inplace(NumericVector p, NumericVector v, NumericVector val)
  {
    for (R_xlen_t i = 0; i < p.length(); ++i)
    {
      val[i] = _update(p[i], v[i]);
    }
  }
};

// mean absolute error

class mass
{
private:
  boost::circular_buffer<double> ebuf;
  sloth::updater::xema_updater<2, false> _ema;
  double _mass;

  double _update(double x)
  {
    _ema.insert(x);
    double e = _ema.m[0] / _ema.m[1];
    if (ebuf.full())
    {
      _mass += e - ebuf.back();
    }
    else
    {
      _mass += e;
    }
    ebuf.push_front(e);
    return _mass;
  }

public:
  mass(int period)
      : ebuf(period),
        _ema(2.0 / (period + 1))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class macd
{
private:
  size_t n, m, p, t;
  sloth::updater::ema_updater<false> ema_n, ema_m, ema_p;
  double _macd, _sig, _hist;

  void _update(double x)
  {
    ema_n.insert(x);
    ema_m.insert(x);
    if (t == m)
    {
      _macd = ema_n.m - ema_m.m;
      ema_p.insert(_macd);
      _sig = ema_p.m;
      _hist = _macd - _sig;
    }
    else
    {
      t += 1;
    }
  }

public:
  macd(int short_period = 12, int long_period = 26, int sig_period = 9)
      : n(std::min(short_period, long_period)),
        m(std::max(short_period, long_period)),
        p(sig_period),
        ema_n(2.0 / (n + 1)),
        ema_m(2.0 / (m + 1)),
        ema_p(2.0 / (p + 1)),
        _macd(0.0),
        _sig(0.0),
        _hist(0.0)
  {
  }

  List run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector vmacd(npt);
    NumericVector vsig(npt);
    NumericVector vhist(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(x[i]);
      vmacd[i] = _macd;
      vsig[i] = _sig;
      vhist[i] = _hist;
    }
    return List::create(
        Named("macd") = vmacd,
        Named("signal") = vsig,
        Named("histo") = vhist);
  }

  void run_inplace(NumericVector x, NumericVector vmacd, NumericVector vsig, NumericVector vhist)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      _update(x[i]);
      vmacd[i] = _macd;
      vsig[i] = _sig;
      vhist[i] = _hist;
    }
  }
};

// linear regression: one pass

class kama
{
private:
  bool init;
  boost::circular_buffer<double> xbuf, dbuf;
  double _kama, dsum, f, s, fs;

  double _update(double x)
  {
    double d, dd, e, a;
    if (init)
    {
      d = std::abs(x - xbuf.front());
      if (dbuf.full())
      {
        dsum += d - dbuf.back();
        dd = x - xbuf.back();
        e = dd / dsum;
        a = std::pow(e, 2);
        _kama += a * (x - _kama);
      }
      else
      {
        dsum += d;
        _kama = x;
      }
      dbuf.push_front(d);
    }
    else
    {
      init = true;
      _kama = x;
    }
    xbuf.push_front(x);
    return _kama;
  }

public:
  kama(int period, int fast_period = 2, int slow_period = 30)
      : init(false),
        xbuf(period),
        dbuf(period),
        _kama(0.0),
        dsum(0.0),
        f(2.0 / (fast_period + 1)),
        s(2.0 / (slow_period - +1)),
        fs(std::abs(f - s))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class hma
{
private:
  wma wma_m, wma_n, wma_s;

  double _update(double x)
  {
    double mm = wma_m._update(x);
    double mn = wma_n._update(x);
    double dd = mm * 2.0 - mn;
    return wma_s._update(dd);
  }

public:
  hma(int period)
      : wma_m(period / 2),
        wma_n(static_cast<int>(std::sqrt(period))),
        wma_s(period)
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

// fisher transform: one pass norm of p?

class ema
{
private:
  sloth::updater::ema_updater<false> _ema;

  double _update(double x)
  {
    _ema.insert(x);
    return _ema.m;
  }

public:
  ema(int period)
      : _ema(2.0 / (period + 1))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

// DMI?

class ppo
{
private:
  size_t fast, slow;
  sloth::updater::ema_updater<false> _ema_fast, _ema_slow;

  double _update(double x)
  {
    _ema_fast.insert(x);
    _ema_slow.insert(x);
    return _ema_fast.m / _ema_slow.m - 1.0;
  }

public:
  ppo(int fast_period = 12, int slow_period = 26)
      : fast(std::min(fast_period, slow_period)),
        slow(std::max(fast_period, slow_period)),
        _ema_fast(2.0 / (fast + 1)),
        _ema_slow(2.0 / (slow + 1))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class dpo
{
private:
  size_t m;
  boost::circular_buffer<double> xbuf;
  sloth::updater::sma_updater<false> sma;

  double _update(double x)
  {
    double val;
    if (xbuf.full())
    {
      sma.roll(x, xbuf.back());
      val = xbuf[m] - sma.m;
    }
    else
    {
      sma.insert(x);
      val = 0.0;
    }
    xbuf.push_front(x);
    return val;
  }

public:
  dpo(int period)
      : m(period / 2),
        xbuf(period),
        sma()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class dema
{
private:
  sloth::updater::xema_updater<2, false> _ema;

  double _update(double x)
  {
    _ema.insert(x);
    return 2.0 * _ema.m[0] - _ema.m[1];
  }

public:
  dema(int period)
      : _ema(2.0 / (period + 1))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

class tema
{
private:
  sloth::updater::xema_updater<3, false> _ema;

  double _update(double x)
  {
    _ema.insert(x);
    return 3.0 * _ema.m[0] - 3 * _ema.m[1] + _ema.m[2];
  }

public:
  tema(int period)
      : _ema(2.0 / (period + 1))
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

// crossover, crossany

class cmo
{
private:
  bool init;
  rollsum up, down;
  double last;

  double _update(double x)
  {
    double sup, sdown;
    if (init)
    {
      if (x > last)
      {
        sup = up._update(x - last);
        sdown = down._update(0.0);
      }
      else if (x < last)
      {
        sup = up._update(0.0);
        sdown = down._update(last - x);
      }
      else
      {
        sup = up._update(0.0);
        sdown = down._update(0.0);
      }
    }
    else
    {
      init = true;
    }
    last = x;
    return 100.0 * ((sup - sdown) / (sup + sdown));
  }

public:
  cmo(int period)
      : init(false),
        up(period),
        down(period),
        last(0.0)
  {
    up._update(1e-10);
    down._update(1e-10);
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

// WARNING: O(N) update
class cmi
{
private:
  boost::circular_buffer<double> xbuf;
  sloth::updater::sma_updater<false> sma;

  double _update(double x)
  {
    if (xbuf.full())
    {
      sma.roll(x, xbuf.back());
    }
    else
    {
      sma.insert(x);
    }
    xbuf.push_front(x);
    // MAE: O(N)
    double m = sma.m;
    double md = 0.0;
    for (auto xx : xbuf)
    {
      md += std::abs(xx - m);
    }
    md /= xbuf.size();
    return (x - m) / (0.015 * md);
  }

public:
  cmi(int period)
      : xbuf(period),
        sma()
  {
  }

  NumericVector run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(x[i]);
    }
    return val;
  }

  void run_inplace(NumericVector x)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      x[i] = _update(x[i]);
    }
  }
};

// Bollinger alternative

class bband
{
private:
  boost::circular_buffer<double> xbuf;
  sloth::updater::s2_updater<false> sigma;
  double z, m, band;

  void _update(double x)
  {
    sigma.insert(x);
    if (xbuf.full())
    {
      sigma.remove(xbuf.back());
    }
    xbuf.push_front(x);
    m = sigma.m;
    band = z * std::sqrt(sigma.s2 / (sigma.n - 1));
  }

public:
  bband(int period, double z)
      : xbuf(period),
        sigma(),
        z(z)
  {
  }

  List run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector upper(npt);
    NumericVector lower(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(x[i]);
      upper[i] = m + band;
      lower[i] = m - band;
    }
    return List::create(
        Named("upper") = upper,
        Named("lower") = lower);
  }

  void run_inplace(NumericVector x, NumericVector upper, NumericVector lower)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      _update(x[i]);
      upper[i] = m + band;
      lower[i] = m - band;
    }
  }
};

class ebband
{
private:
  sloth::updater::s2ew_updater<false> sigma;
  double z, m, band;

  void _update(double x)
  {
    sigma.insert(x);
    m = sigma.m;
    band = z * std::sqrt(sigma.s2);
  }

public:
  ebband(int period, double z)
      : sigma(2.0 / (period + 1)),
        z(z)
  {
  }

  List run(NumericVector x)
  {
    R_xlen_t npt = x.length();
    NumericVector upper(npt);
    NumericVector lower(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      _update(x[i]);
      upper[i] = m + band;
      lower[i] = m - band;
    }
    return List::create(
        Named("upper") = upper,
        Named("lower") = lower);
  }

  void run_inplace(NumericVector x, NumericVector upper, NumericVector lower)
  {
    for (R_xlen_t i = 0; i < x.length(); ++i)
    {
      _update(x[i]);
      upper[i] = m + band;
      lower[i] = m - band;
    }
  }
};

class atr
{
private:
  bool init;
  sloth::updater::ema_updater<false> ema;
  double last;

  double _update(double high, double low, double close)
  {
    double tr;
    if (init)
    {
      tr = std::max(std::max(high - low, std::abs(high - last)), std::abs(low - last));
    }
    else
    {
      tr = high - low;
      init = true;
    }
    last = close;
    ema.insert(tr);
    return ema.m;
  }

public:
  atr(int period)
      : init(false),
        ema(1.0 / period)
  {
  }

  NumericVector run(NumericVector high, NumericVector low, NumericVector close)
  {
    R_xlen_t npt = high.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(high[i], low[i], close[i]);
    }
    return val;
  }

  void run_inplace(NumericVector high, NumericVector low, NumericVector close, NumericVector val)
  {
    for (R_xlen_t i = 0; i < high.length(); ++i)
    {
      val[i] = _update(high[i], low[i], close[i]);
    }
  }
};

class aroonosc
{
private:
  boost::circular_buffer<double> hbuf, lbuf;
  sloth::updater::argminmax_updater<double> minmax;

  double _update(double high, double low)
  {
    if (hbuf.full())
    {
      minmax.remove(hbuf.back());
      minmax.remove(lbuf.back());
    }
    minmax.insert(high);
    minmax.insert(low);
    hbuf.push_front(high);
    lbuf.push_front(low);
    double up = 1.0 - static_cast<double>(minmax.maxidx.front()) / hbuf.size();
    double down = 1.0 - static_cast<double>(minmax.minidx.front()) / hbuf.size();
    return 100.0 * (up - down);
  }

public:
  aroonosc(int period)
      : hbuf(period),
        lbuf(period),
        minmax(period)
  {
  }

  NumericVector run(NumericVector high, NumericVector low)
  {
    R_xlen_t npt = high.length();
    NumericVector val(npt);
    for (R_xlen_t i = 0; i < npt; ++i)
    {
      val[i] = _update(high[i], low[i]);
    }
    return val;
  }

  void run_inplace(NumericVector high, NumericVector low, NumericVector val)
  {
    for (R_xlen_t i = 0; i < high.length(); ++i)
    {
      val[i] = _update(high[i], low[i]);
    }
  }
};

// ADX? 
