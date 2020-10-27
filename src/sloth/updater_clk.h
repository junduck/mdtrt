#pragma once

#include "updater_common.h"

namespace sloth
{
  namespace updater
  {

    // Volume/Turnover clock bar
    struct vclk_bar
    {
      int nbin;
      long ibin;
      double open, high, low, close, vol, tnvr, vwap;
    };

    // Time clock bar
    template <class TT>
    struct tclk_bar
    {
      TT time;
      int nbin;
      long ibin;
      double open, high, low, close, vol, tnvr, vwap;
    };

    struct volclk_updater
    {
      bool init;
      double bsize;
      vclk_bar bar, res;

      explicit volclk_updater(double bucket_size)
          : init(false),
            bsize(bucket_size),
            bar(), res()
      {
      }

      void insert(double p, double v)
      {
        if (!init)
        {
          res.open = res.high = res.low = p;
          init = true;
        }
        res.close = p;
        double tot_vol = res.vol + v;
        if (tot_vol < bsize)
        {
          res.nbin = 0;
          res.high = std::max(p, res.high);
          res.low = std::min(p, res.low);
          res.tnvr += v * p;
          res.vol = tot_vol;
        }
        else
        {
          bar.open = res.open;
          bar.high = res.high;
          bar.low = res.low;
          bar.close = res.close;
          bar.vwap = (res.tnvr + (bsize - res.vol) * p) / bsize;
          // reservoir
          res.vol = std::fmod(tot_vol, bsize);
          res.nbin = static_cast<int>(tot_vol / bsize);
          res.tnvr = res.vol * p;
          res.open = res.high = res.low = p;
        }
        res.ibin += res.nbin;
        bar.nbin = res.nbin;
        bar.ibin = res.ibin;
      }
    };

    struct tnvrclk_updater
    {
      bool init;
      double bsize;
      vclk_bar bar, res;

      explicit tnvrclk_updater(double bucket_size)
          : init(false),
            bsize(bucket_size),
            bar(), res()
      {
      }

      void insert(double p, double v)
      {
        if (!init)
        {
          res.open = res.high = res.low = p;
          init = true;
        }
        res.close = p;
        double tot_tnvr = res.tnvr + v * p;
        if (tot_tnvr < bsize)
        {
          bar.nbin = 0;
          res.high = std::max(p, res.high);
          res.low = std::min(p, res.low);
          res.vol += v;
          res.tnvr = tot_tnvr;
        }
        else
        {
          bar.open = res.open;
          bar.high = res.high;
          bar.low = res.low;
          bar.close = res.close;
          bar.vwap = bsize / (res.vol + (bsize - res.tnvr) / p);
          // reservoir
          res.tnvr = std::fmod(tot_tnvr, bsize);
          res.nbin = static_cast<int>(tot_tnvr / bsize);
          res.vol = res.tnvr / p;
          res.open = res.high = res.low = p;
        }
        res.ibin += res.nbin;
        bar.nbin = res.nbin;
        bar.ibin = res.ibin;
      }
    };

    // Naiive tumbling window impl, left open: (0, dt]
    template <class TT, class TDT = TT>
    struct twclk_updater
    {
      bool init, align;
      TDT bsize;
      tclk_bar<TT> bar, res;

      twclk_updater(TDT window_size, bool align_to_window)
          : init(false),
            align(align_to_window),
            bsize(window_size),
            bar(), res()
      {
      }

      void insert(TT t, double p, double v)
      {
        if (!init)
        {
          res.open = res.high = res.low = p;
          res.time = t;
          if (align)
          {
            TT remain = res.time;
            while (remain >= bsize)
            {
              remain -= bsize;
            }
            res.time -= remain;
          }
          init = true;
        }
        TDT dt = t - res.time;
        res.close = p;
        res.nbin = 0;
        if (dt < bsize)
        {
          res.high = std::max(p, res.high);
          res.low = std::min(p, res.low);
          res.vol += v;
          res.tnvr += p * v;
        }
        else
        {
          while (dt >= bsize)
          {
            dt -= bsize;
            res.time += bsize;
            ++res.nbin;
          }
          if (dt)
          {
            // p, v belong to next bar, do not merge
            bar.open = res.open;
            bar.high = res.high;
            bar.low = res.low;
            bar.close = res.close;
            bar.vol = res.vol;
            bar.vwap = res.tnvr / res.vol;
            res.vol = v;
            res.tnvr = p * v;
          }
          else
          {
            // merge p, v into bar
            bar.open = res.open;
            bar.high = std::max(p, res.high);
            bar.low = std::min(p, res.low);
            bar.close = p;
            bar.vol = res.vol + v;
            bar.vwap = (res.tnvr + p * v) / bar.vol;
            res.vol = 0;
            res.tnvr = 0;
          }
          res.open = res.high = res.low = p;
        }
        res.ibin += res.nbin;
        bar.ibin = res.ibin;
        bar.nbin = res.nbin;
      }
    };

    // Naiive dynamic sliding window impl
    template <class TT>
    struct swclk_updater
    {
      TT time;
      double open, high, low, close, tnvr, vol, vwap;
      boost::circular_buffer<double> pbuf, vbuf, hbuf, lbuf;
      boost::circular_buffer<TT> tbuf;

      swclk_updater(int max_buffer_size)
          : open(0.0), high(0.0), low(0.0), close(0.0), vol(0.0),
            pbuf(max_buffer_size),
            vbuf(max_buffer_size),
            hbuf(max_buffer_size),
            lbuf(max_buffer_size),
            tbuf(max_buffer_size)
      {
      }

      void insert(TT t, double p, double v)
      {
        // close, vol : O(1)
        close = p;
        vol += v;
        tnvr += p * v;
        vwap = tnvr / vol;
        pbuf.push_front(p);
        vbuf.push_front(v);
        // high, low : O(log N)
        while (!hbuf.empty() && hbuf.back() < p)
          hbuf.pop_back();
        while (!lbuf.empty() && lbuf.back() > p)
          lbuf.pop_back();
        hbuf.push_back(p);
        lbuf.push_back(p);
        high = hbuf.front();
        low = lbuf.front();
        // time : O(1)
        time = t;
        tbuf.push_front(t);
      }

      void remove()
      {
        // open, vol : O(1)
        double p = pbuf.back();
        double v = vbuf.back();
        pbuf.pop_back();
        vbuf.pop_back();
        open = pbuf.back();
        vol -= v;
        tnvr -= p * v;
        vwap = tnvr / vol;
        // high, low : O(1)
        if (hbuf.front() == p)
          hbuf.pop_front();
        if (lbuf.front() == p)
          lbuf.pop_front();
        high = hbuf.front();
        low = lbuf.front();
        // time : O(1)
        tbuf.pop_back();
      }
    };

  } // namespace updater

} // namespace sloth
