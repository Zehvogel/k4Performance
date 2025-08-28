#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_all_hists.py

Plots validation histograms written by PFOtoMCviaClusterLink.

Features:
  - Sanity check: fails if required histograms are missing (unless --no-strict).
  - Looks only in a single directory (--hist-dir, default PFOCluster).
  - Handles single-particle pgun as well as generic samples.
  - Produces spectra, efficiency, response/resolution, PID confusion plots.
  - Outputs images and optionally a combined multi-page PDF.
"""

import os, sys, argparse, math
from collections import namedtuple
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# -------------------- Categories --------------------
CAT_TO_PNAME = {0: "mu", 1: "gamma", 2: "electron", 3: "nh", 4: "ch", 5: "other"}
DEFAULT_HIST_DIR = "PFOCluster"

REQUIRED_HISTS = [
    "mcE", "pfoE", "mcCos", "pfoCos", "respE",
    "effE_den", "effE_num", "effCos_den", "effCos_num", "pidConf"
]

EPOINTS = [5, 10, 20, 50, 100]
Region = namedtuple("Region", "name tmin tmax")

# -------------------- Utils --------------------
def mkdir_p(p):
    os.makedirs(p, exist_ok=True)

def open_root(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise IOError("Cannot open %s" % path)
    return f

def geth(f, name, hist_dir, debug=False):
    """Get histogram by name inside hist_dir only"""
    base = hist_dir.strip("/")
    obj = f.Get(base + "/" + name)
    if obj and debug:
        print(f"  [geth] Found {base}/{name}")
    return obj

# -------------------- Physics helpers --------------------
def infer_theta_range_deg_from_mcCos(f, hist_dir, debug=False, thr=1):
    h = geth(f, "mcCos", hist_dir, debug)
    if not h or h.GetEntries() == 0:
        return (None, None)
    ax = h.GetXaxis()
    i_lo, i_hi = None, None
    for i in range(1, ax.GetNbins() + 1):
        if h.GetBinContent(i) > thr:
            i_lo = i
            break
    for i in range(ax.GetNbins(), 0, -1):
        if h.GetBinContent(i) > thr:
            i_hi = i
            break
    if i_lo is None or i_hi is None:
        return (None, None)
    lo_c = ax.GetBinLowEdge(i_lo)
    hi_c = ax.GetBinUpEdge(i_hi)
    def c2d(c):
        c = max(-1, min(1, c))
        return math.degrees(math.acos(c))
    return min(c2d(lo_c), c2d(hi_c)), max(c2d(lo_c), c2d(hi_c))

def pid_purity_from_pidConf(f, hist_dir, debug=False):
    h = geth(f, "pidConf", hist_dir, debug)
    if not h:
        return 5, None
    x = h.GetXaxis(); y = h.GetYaxis()
    sums = [sum(h.GetBinContent(ix, iy) for iy in range(1, y.GetNbins() + 1))
            for ix in range(1, x.GetNbins() + 1)]
    tot = sum(sums)
    if tot <= 0:
        return 5, 0.0
    dom = int(max(range(len(sums)), key=lambda i: sums[i]))
    return dom, sums[dom]/tot

# -------------------- Plot helpers --------------------
def save_canvas(c, png, pdf=None):
    c.SaveAs(png)
    if pdf:
        if not os.path.exists(pdf):
            c.SaveAs(pdf + "(")
        else:
            c.SaveAs(pdf)

def teff_graph(num, den):
    if not num or not den or den.Integral() <= 0:
        return None
    eff = ROOT.TEfficiency(num, den)
    return eff.CreateGraph()

def project_resp_y_at_E(h2, e_center, de=0.5):
    ax = h2.GetXaxis()
    i1 = max(1, min(ax.FindBin(e_center - de), ax.GetNbins()))
    i2 = max(1, min(ax.FindBin(e_center + de), ax.GetNbins()))
    return h2.ProjectionY(f"{h2.GetName()}_E{e_center}", i1, i2)

def iterative_gauss_sigma(h1, max_iter=8):
    if not h1 or h1.GetEntries() < 20:
        return h1.GetRMS() if h1 else 0.0
    mean = h1.GetMean(); sigma = h1.GetRMS()
    for _ in range(max_iter):
        lo, hi = mean - 3*sigma, mean + 3*sigma
        f = ROOT.TF1("g", "gaus", lo, hi)
        r = h1.Fit(f, "QRN")
        if int(r) != 0:
            break
        mean = f.GetParameter(1)
        sigma = abs(f.GetParameter(2))
    return sigma

def resp_from_Erec_hist(hErec, Etrue):
    if not hErec or Etrue <= 0:
        return None
    nb = hErec.GetNbinsX()
    out = ROOT.TH1F(hErec.GetName() + "_resp", ";(Erec/Etrue-1);Entries", nb, -1.0, 1.0)
    for i in range(1, nb + 1):
        Ecen = hErec.GetXaxis().GetBinCenter(i)
        out.Fill((Ecen/Etrue)-1.0, hErec.GetBinContent(i))
    return out

# -------------------- Sanity check --------------------
def sanity_check_required(f, hist_dir, debug, strict):
    missing = []
    for key in REQUIRED_HISTS:
        if not geth(f, key, hist_dir, debug):
            missing.append(key)
    if missing and strict:
        print(f"ERROR: Missing histograms in {hist_dir}: {', '.join(missing)}")
        sys.exit(2)
    return missing

# -------------------- CLI and driver --------------------
def main():
    ap = argparse.ArgumentParser(description="Plot ROOT histograms from validation.")
    ap.add_argument("files", nargs="+", help="Input ROOT files.")
    ap.add_argument("--out", default="plots", help="Output directory.")
    ap.add_argument("--pdf", default=None, help="Optional multi-page PDF output.")
    ap.add_argument("--hist-dir", default=DEFAULT_HIST_DIR,
                    help="Directory in ROOT file containing histograms.")
    ap.add_argument("--pgun", action="store_true", help="Force pgun behavior.")
    ap.add_argument("--species", choices=list(CAT_TO_PNAME.values()), help="Force particle species.")
    ap.add_argument("--purity-thr", type=float, default=0.85,
                    help="Purity threshold to assume pgun.")
    ap.add_argument("--min-proj", type=int, default=50, help="Min entries per E slice for resolution.")
    ap.add_argument("--no-plots", action="store_true", help="Do not produce plots.")
    ap.add_argument("--debug", action="store_true")
    strict_grp = ap.add_mutually_exclusive_group()
    strict_grp.add_argument("--strict", dest="strict", action="store_true")
    strict_grp.add_argument("--no-strict", dest="strict", action="store_false")
    ap.set_defaults(strict=True)
    args = ap.parse_args()

    mkdir_p(args.out)

    # Process files
    summary = []
    for path in args.files:
        print(f"Processing {path}")
        f = open_root(path)
        sanity_check_required(f, args.hist_dir, args.debug, args.strict)

        # Infer particle type and purity
        dom_cat, purity = pid_purity_from_pidConf(f, args.hist_dir, args.debug)
        pname = args.species if args.species else CAT_TO_PNAME.get(dom_cat, "other")
        tmin, tmax = infer_theta_range_deg_from_mcCos(f, args.hist_dir, args.debug)

        summary.append((os.path.basename(path), pname, purity, tmin, tmax))

        if args.no_plots:
            continue

        # Example: Efficiency vs E
        h_num = geth(f, "effE_num", args.hist_dir)
        h_den = geth(f, "effE_den", args.hist_dir)
        g = teff_graph(h_num, h_den)
        if g:
            c = ROOT.TCanvas("c_effE", "", 600, 500)
            g.Draw("AP")
            g.GetXaxis().SetTitle("MC E [GeV]")
            g.GetYaxis().SetTitle("Efficiency")
            save_canvas(c, os.path.join(args.out, f"{pname}_effE.png"), args.pdf)

        # PID confusion
        h_pid = geth(f, "pidConf", args.hist_dir)
        if h_pid:
            c = ROOT.TCanvas("c_pid", "", 600, 500)
            h_pid.Draw("COLZ TEXT")
            save_canvas(c, os.path.join(args.out, f"{pname}_pidConf.png"), args.pdf)

    # Print summary table
    headers = ["File", "Particle", "Purity", "theta range [deg]"]
    rows = []
    for fname, pname, pur, tmin, tmax in summary:
        purtxt = "n/a" if pur is None else f"{pur:.2f}"
        thtxt = "n/a" if tmin is None else f"{tmin:.1f}-{tmax:.1f}"
        rows.append([fname, pname, purtxt, thtxt])

    # compute column widths
    cols = list(zip(headers, *rows))
    widths = [max(len(str(c)) for c in col) for col in cols]

    def fmt_row(r):
        return " | ".join(str(val).ljust(w) for val, w in zip(r, widths))

    print("\nInput summary:")
    print(fmt_row(headers))
    print("-+-".join("-" * w for w in widths))
    for r in rows:
        print(fmt_row(r))

    if args.pdf:
        c = ROOT.TCanvas()
        c.SaveAs(args.pdf + "]")

if __name__ == "__main__":
    main()
