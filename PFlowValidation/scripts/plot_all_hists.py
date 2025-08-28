#!/usr/bin/env python3

import os, sys, math, argparse, collections
from array import array
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ---------------------------
# Small helpers
# ---------------------------

def parse_region(fname: str) -> str:
    n = os.path.basename(fname).lower()
    if 'barrel' in n:     return 'barrel'
    if 'transition' in n: return 'transition'
    if 'endcap' in n:     return 'endcap'
    return 'all'

def _try_get(f, full):
    obj = f.Get(full)
    return obj if obj else None

def get_hist(f, path, name):
    """
    Robust histogram getter:
      - tries <path>/<name>
      - tries /<path>/<name>
      - tries /MC/<path>/<name>
    Returns a detached clone or None.
    """
    trials = [
        f'{path}/{name}',
        f'/{path}/{name}',
        f'/MC/{path}/{name}',
    ]
    for full in trials:
        full = full.replace('//','/')
        obj = _try_get(f, full)
        if obj:
            h = obj.Clone(f'{name}__clone')
            h.SetDirectory(0)
            return h
    return None

def energy_from_mcE_sel(h):
    if not h or h.GetEntries() <= 0: return None
    ib = h.GetMaximumBin()
    return h.GetXaxis().GetBinCenter(ib)

def integral_and_err(h):
    """Use array('d') for PyROOT 6.36 compatibility."""
    if not h: return 0.0, 0.0
    err = array('d', [0.0])
    val = h.IntegralAndError(1, h.GetNbinsX(), err, "")
    return float(val), float(err[0])

def efficiency_overall(hnum, hden):
    nnum, _ = integral_and_err(hnum)
    nden, _ = integral_and_err(hden)
    if nden <= 0:
        return None, None
    eff = nnum/nden
    err = math.sqrt(max(eff*(1.0-eff)/nden, 0.0))
    return eff, err

def efficiency_vs_costheta(hnum, hden):
    out = []
    if (not hnum) or (not hden): return out
    xax = hden.GetXaxis()
    for i in range(1, hden.GetNbinsX()+1):
        den = hden.GetBinContent(i)
        num = hnum.GetBinContent(i) if i <= hnum.GetNbinsX() else 0.0
        x   = xax.GetBinCenter(i)
        ex  = xax.GetBinWidth(i)/2.0
        if den <= 0:
            out.append((x, 0.0, ex, 0.0))
            continue
        eff = num/den
        err = math.sqrt(max(eff*(1.0-eff)/den, 0.0))
        out.append((x, eff, ex, err))
    return out

def fit_gauss(h, rangemin=-0.3, rangemax=0.3):
    """Fit respR_sel with a Gaussian; return (sigma%, err%, TF1 or None, ok)."""
    if (not h) or h.GetEntries() < 10:
        return None, None, None, False
    h.GetXaxis().SetRangeUser(rangemin, rangemax)
    rms = h.GetRMS() or 0.05
    fn  = ROOT.TF1(f'gaus_{id(h)}', 'gaus', rangemin, rangemax)
    fn.SetParameters(max(h.GetMaximum(), 1.0), 0.0, rms)
    res = h.Fit(fn, 'RQSN')
    ok  = (int(res) == 0)
    if not ok:
        fn.SetRange(-0.5, 0.5)
        res = h.Fit(fn, 'RQSN')
        ok  = (int(res) == 0)
    if not ok:
        return None, None, None, False
    sig  = fn.GetParameter(2)
    esig = fn.GetParError(2)
    return 100.0*sig, 100.0*esig, fn, True

def style_graph(gr, color, marker=20):
    gr.SetLineColor(color); gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker); gr.SetLineWidth(2)

COL = {
    'barrel':     ROOT.kBlue+1,
    'transition': ROOT.kOrange+7,
    'endcap':     ROOT.kRed+1,
    'all':        ROOT.kBlack,
}

# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Minimal plotting: resolution (respR_sel) & efficiencies (selected histos only)."
    )
    ap.add_argument('files', nargs='+', help='Input ROOT files (one particle type per run).')
    ap.add_argument('-o','--outdir', default='figs', help='Output directory (default: figs)')
    ap.add_argument('--histpath', default='PFOCluster', help='Histogram directory (default: PFOCluster)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fitdir = os.path.join(args.outdir, 'fit'); os.makedirs(fitdir, exist_ok=True)

    print(f'[info] Using histpath="{args.histpath}" (will also try /{args.histpath} and /MC/{args.histpath})')

    File = collections.namedtuple('File','path region E h')
    files = []

    # ---------------- load files ----------------
    for p in args.files:
        f = ROOT.TFile.Open(p)
        if not f or f.IsZombie():
            print(f'[warn] cannot open: {p} — skipping.')
            continue

        region = parse_region(p)
        # Only "selected" particles are used (= fullfilling criteria: match in type, angles, pT)
        h_mcE_sel    = get_hist(f, args.histpath, 'mcE_sel')
        h_respR_sel  = get_hist(f, args.histpath, 'respR_sel')
        h_effE_num   = get_hist(f, args.histpath, 'effE_sel_num')
        h_effE_den   = get_hist(f, args.histpath, 'effE_sel_den')
        h_effCos_num = get_hist(f, args.histpath, 'effCos_sel_num')
        h_effCos_den = get_hist(f, args.histpath, 'effCos_sel_den')

        bn = os.path.basename(p)
        if not h_mcE_sel:
            print(f'[warn] {bn}: missing mcE_sel under {args.histpath} — skipping file.')
            f.Close(); continue

        E = energy_from_mcE_sel(h_mcE_sel)
        if E is None:
            print(f'[warn] {bn}: mcE_sel is empty — cannot infer energy; skipping file.')
            f.Close(); continue
        print(f'[info] {bn}: region={region}, E≈{E:g} GeV')

        if not h_respR_sel or h_respR_sel.GetEntries() == 0:
            print(f'[warn] {bn}: respR_sel empty — no resolution point will be plotted.')
        if (not h_effE_den) or h_effE_den.Integral() == 0:
            print(f'[warn] {bn}: effE_sel_den empty/missing — efficiency vs E point will be skipped.')
        if (not h_effCos_den) or h_effCos_den.Integral() == 0:
            print(f'[warn] {bn}: effCos_sel_den empty/missing — efficiency vs cosθ will be skipped.')

        files.append(File(
            path=p, region=region, E=E,
            h={'mcE_sel':h_mcE_sel, 'respR_sel':h_respR_sel,
               'effE_num':h_effE_num, 'effE_den':h_effE_den,
               'effCos_num':h_effCos_num, 'effCos_den':h_effCos_den}
        ))

        # Fit preview (if available)
        if h_respR_sel and h_respR_sel.GetEntries() > 0:
            c = ROOT.TCanvas('c_fit', 'fit', 700, 500); c.SetGrid()
            hdraw = h_respR_sel.Clone(f'respR_sel__{bn}'); hdraw.SetDirectory(0)
            hdraw.SetTitle(f'{region} | E={E:g} GeV; (E_{{rec}}/E_{{true}}-1); Entries')
            sigp, esigp, fn, ok = fit_gauss(hdraw)
            hdraw.Draw('HIST')
            if ok and fn:
                fn.SetLineColor(ROOT.kRed+1); fn.Draw('SAME')
                box = ROOT.TPaveText(0.55, 0.70, 0.88, 0.88, 'NDC')
                box.SetFillStyle(0); box.SetBorderSize(0)
                box.AddText(f'#sigma = {sigp:.2f} #pm {esigp:.2f} %'); box.Draw()
                print(f'[info] {bn}: fit σ={sigp:.2f}±{esigp:.2f} %')
            else:
                print(f'[warn] {bn}: Gaussian fit failed.')
            out = os.path.join(fitdir, f'{region}_E{E:g}GeV_fit.pdf')
            c.SaveAs(out); print(f'[write] {out}')
            del c

        f.Close()

    if not files:
        print('[error] no usable files — nothing to plot.')
        sys.exit(1)

    # group by region
    by_region = collections.defaultdict(list)
    for fi in files: by_region[fi.region].append(fi)

    # ---------------- Resolution vs E ----------------
    c_res = ROOT.TCanvas('c_res', 'Resolution vs E', 900, 700); c_res.SetGrid()
    leg   = ROOT.TLegend(0.58, 0.18, 0.88, 0.38); leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetHeader('Regions','C')
    graphs = {}; xmax, ymax = 0.0, 0.0

    for region, fis in by_region.items():
        pts = []
        for fi in sorted(fis, key=lambda x: x.E):
            h = fi.h['respR_sel']
            if not h or h.GetEntries()==0:
                print(f'[info] {os.path.basename(fi.path)}: no respR_sel; skipping in resolution.')
                continue
            sigp, esigp, _, ok = fit_gauss(h)
            if not ok:
                print(f'[warn] {os.path.basename(fi.path)}: fit failed; skipping point.')
                continue
            pts.append((fi.E, sigp, 0.0, esigp))
            xmax = max(xmax, fi.E); ymax = max(ymax, sigp + esigp)
        if not pts:
            print(f'[info] region {region}: no resolution points.')
            continue
        gr = ROOT.TGraphErrors(len(pts))
        for i,(x,y,ex,ey) in enumerate(pts): gr.SetPoint(i,x,y); gr.SetPointError(i,ex,ey)
        style_graph(gr, COL.get(region, ROOT.kGray+2), marker=20)
        graphs[region]=gr; leg.AddEntry(gr, region, 'lp')

    frame = c_res.DrawFrame(0.0, 0.0, max(1.0, 1.05*xmax), max(5.0, 1.2*ymax if ymax>0 else 10.0))
    frame.SetTitle('Energy resolution (selected);E [GeV];#sigma(E)/E [%]')
    for gr in graphs.values(): gr.Draw('P SAME')
    leg.Draw()
    out_res = os.path.join(args.outdir, 'resolution_vs_energy_by_region.pdf')
    c_res.SaveAs(out_res); print(f'[write] {out_res}')

    # ---------------- Efficiency vs E ----------------
    c_effE = ROOT.TCanvas('c_effE', 'Efficiency vs E', 900, 700); c_effE.SetGrid()
    legE   = ROOT.TLegend(0.58, 0.18, 0.88, 0.38); legE.SetFillStyle(0); legE.SetBorderSize(0); legE.SetHeader('Regions','C')
    graphsE = {}; xmax, ymax = 0.0, 0.0

    for region, fis in by_region.items():
        pts = []
        for fi in sorted(fis, key=lambda x: x.E):
            eff, err = efficiency_overall(fi.h['effE_num'], fi.h['effE_den'])
            if eff is None:
                print(f'[info] {os.path.basename(fi.path)}: missing/empty effE_sel_den — skipping point.')
                continue
            pts.append((fi.E, 100.0*eff, 0.0, 100.0*err))
            xmax = max(xmax, fi.E); ymax = max(ymax, 100.0*(eff+err))
        if not pts:
            print(f'[info] region {region}: no efficiency(E) points.')
            continue
        gr = ROOT.TGraphErrors(len(pts))
        for i,(x,y,ex,ey) in enumerate(pts): gr.SetPoint(i,x,y); gr.SetPointError(i,ex,ey)
        style_graph(gr, COL.get(region, ROOT.kGray+2), marker=21)
        graphsE[region]=gr; legE.AddEntry(gr, region, 'lp')

    frame = c_effE.DrawFrame(0.0, 0.0, max(1.0, 1.05*xmax), min(105.0, max(10.0, 1.2*ymax)))
    frame.SetTitle('Reconstruction efficiency (selected);E [GeV];Efficiency [%]')
    for gr in graphsE.values(): gr.Draw('P SAME')
    legE.Draw()
    out_effE = os.path.join(args.outdir, 'efficiency_vs_energy_by_region_sel.pdf')
    c_effE.SaveAs(out_effE); print(f'[write] {out_effE}')

    # ---------------- Efficiency vs cosTheta (one canvas per region) ----------------
    for region, fis in by_region.items():
        energies = sorted({fi.E for fi in fis})
        if not energies:
            print(f'[info] region {region}: no energies for eff vs cosTheta.'); continue

        c = ROOT.TCanvas(f'c_effCos_{region}', f'eff vs cos#theta [{region}]', 900, 700); c.SetGrid()
        legC = ROOT.TLegend(0.58, 0.18, 0.88, 0.38); legC.SetFillStyle(0); legC.SetBorderSize(0); legC.SetHeader(region,'C')
        graphsC = []; ymax, ymin = 0.0, 1.0

        # x-range from any denominator
        any_den = next((fi.h['effCos_den'] for fi in fis if fi.h['effCos_den']), None)
        if not any_den:
            print(f'[info] region {region}: no effCos_sel_den — skipping canvas.')
            continue
        xmin = any_den.GetXaxis().GetXmin(); xmaxx = any_den.GetXaxis().GetXmax()

        emin, emax = min(energies), max(energies)
        for E in energies:
            fi = next((x for x in fis if abs(x.E - E) < 1e-9), None)
            if not fi: continue
            pts = efficiency_vs_costheta(fi.h['effCos_num'], fi.h['effCos_den'])
            if not pts:
                print(f'[info] {os.path.basename(fi.path)}: no per-bin denom — skip E={E:g} GeV.'); continue
            gr = ROOT.TGraphErrors(len(pts))
            for i,(x,y,ex,ey) in enumerate(pts):
                gr.SetPoint(i, x, 100.0*y); gr.SetPointError(i, ex, 100.0*ey)
                ymax = max(ymax, 100.0*(y+ey)); ymin = min(ymin, 100.0*max(0.0, y-ey))
            # palette by energy
            color = ROOT.kBlack if len(energies)==1 else ROOT.TColor.GetColorPalette(
                int(50 + 200 * (0.0 if emax==emin else (E-emin)/(emax-emin)))
            )
            style_graph(gr, color, marker=20)
            graphsC.append(gr); legC.AddEntry(gr, f'E={E:g} GeV', 'lp')

        if not graphsC:
            print(f'[info] region {region}: nothing to draw for eff vs cosθ.'); continue

        frame = c.DrawFrame(xmin, max(0.0, 0.9*ymin), xmaxx, min(105.0, max(60.0, 1.1*ymax)))
        frame.SetTitle('Reconstruction efficiency vs cos#theta (selected);cos#theta;Efficiency [%]')
        for gr in graphsC: gr.Draw('P SAME')
        legC.Draw()
        out_effC = os.path.join(args.outdir, f'efficiency_vs_costheta_{region}_sel.pdf')
        c.SaveAs(out_effC); print(f'[write] {out_effC}')

    print('[done] all plots written.')

if __name__ == '__main__':
    main()
