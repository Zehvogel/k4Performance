#!/usr/bin/env python3
import os
import sys
import re
import math
import argparse
import collections
import csv
from array import array
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# ---------------------------
# Parsing helpers
# ---------------------------

def parse_region(fname: str) -> str:
    n = os.path.basename(fname).lower()
    if 'barrel' in n:
        return 'barrel'
    if 'transition' in n:
        return 'transition'
    if 'endcap' in n:
        return 'endcap'
    return 'all'  # filenames without a region tag

def parse_particle(fname: str) -> str:
    n = os.path.basename(fname).lower()
    if any(k in n for k in ['gamma', 'photon', 'gam']):
        return 'gamma'
    if any(k in n for k in ['kaon0l', 'k0l', 'kaon0long', 'k0long', 'k_0^l', 'kl']):
        return 'kaon0l'
    if any(k in n for k in ['pi-', 'pim', 'pion-', 'pi+', 'pip', 'pion+']):
        return 'pion'
    if any(k in n for k in ['mu-', 'mum', 'muon-', 'mu+', 'mup', 'muon+']):
        return 'muon'
    return 'particle'

def parse_energy_from_name(fname: str):
    n = os.path.basename(fname)
    m = re.search(r'(\d+(?:\.\d+)?)\s*GeV', n, flags=re.I)
    return float(m.group(1)) if m else None

# ---------------------------
# ROOT helpers
# ---------------------------

def _try_get(f, full):
    return f.Get(full) if f else None

def get_hist(f, path, name):
    trials = [
        f'{path}/{name}',
        f'/{path}/{name}',
        f'/MC/{path}/{name}',
        f'PFOCluster/{name}',
        f'/PFOCluster/{name}',
        f'/MC/PFOCluster/{name}',
    ]
    for full in trials:
        full = full.replace('//', '/')
        obj = _try_get(f, full)
        if obj:
            h = obj.Clone(f'{name}__clone')
            h.SetDirectory(0)
            return h
    return None

def energy_from_mcE_sel(h):
    if not h or h.GetEntries() <= 0:
        return None
    ib = h.GetMaximumBin()
    return h.GetXaxis().GetBinCenter(ib)

def integral_and_err(h):
    if not h:
        return 0.0, 0.0
    err = array('d', [0.0])
    val = h.IntegralAndError(1, h.GetNbinsX(), err, "")
    return float(val), float(err[0])

def efficiency_overall(hnum, hden):
    nnum, _ = integral_and_err(hnum)
    nden, _ = integral_and_err(hden)
    if nden <= 0:
        return None, None, 0.0, 0.0
    eff = nnum / nden
    err = math.sqrt(max(eff * (1.0 - eff) / nden, 0.0))
    return eff, err, nnum, nden

def efficiency_vs_costheta(hnum, hden):
    out = []
    if (not hnum) or (not hden):
        return out
    xax = hden.GetXaxis()
    for i in range(1, hden.GetNbinsX() + 1):
        den = hden.GetBinContent(i)
        num = hnum.GetBinContent(i) if i <= hnum.GetNbinsX() else 0.0
        x = xax.GetBinCenter(i)
        ex = xax.GetBinWidth(i) / 2.0
        if den <= 0:
            out.append((x, ex, 0.0, 0.0, num, den))
            continue
        eff = num / den
        err = math.sqrt(max(eff * (1.0 - eff) / den, 0.0))
        out.append((x, ex, eff, err, num, den))
    return out

# ---------------------------
# Two-step Gaussian fit (no rebin)
# ---------------------------

def fit_gauss_twostep(h, default_range=(-0.3, 0.3)):
    """
    Two-step Gaussian fit on 'h' without any rebinning:
      1) Fit in a general window (default_range clipped to axis).
      2) Refit in +/- 3*sigma around the first fit's mean.
    Returns dict with:
      ok, h (clone), fn (final TF1), xmin, xmax, sigp, esigp, mu, sigma, stage
    """
    if (not h) or h.GetEntries() < 10:
        return {'ok': False}

    # Work on a detached clone (keep original binning, no Rebin)
    hfit = h.Clone(h.GetName() + '_fitclone')
    hfit.SetDirectory(0)

    # Axis limits and initial window
    ax = hfit.GetXaxis()
    xlo_tot, xhi_tot = ax.GetXmin(), ax.GetXmax()
    x1 = max(default_range[0], xlo_tot)
    x2 = min(default_range[1], xhi_tot)
    if x2 <= x1:
        x1, x2 = xlo_tot, xhi_tot

    # 1) General fit
    f1 = ROOT.TF1('g1_{}'.format(id(hfit)), 'gaus', x1, x2)
    amp0 = max(hfit.GetBinContent(hfit.GetMaximumBin()), 1.0)
    mean0 = hfit.GetMean()
    rms0  = hfit.GetRMS() or (x2 - x1) / 6.0
    f1.SetParameters(amp0, mean0, rms0)
    r1 = hfit.Fit(f1, 'RQSN')
    ok1 = (int(r1) == 0)
    if not ok1:
        # Retry on full axis once
        f1.SetRange(xlo_tot, xhi_tot)
        f1.SetParameters(amp0, 0.0 if abs(mean0) > 0.5 else mean0, rms0)
        if int(hfit.Fit(f1, 'RQSN')) != 0:
            return {'ok': False}

    mu1 = f1.GetParameter(1)
    sg1 = abs(f1.GetParameter(2)) or (x2 - x1) / 6.0

    # 2) Second fit in +/- 3*sigma from the first fit
    w2_lo = max(xlo_tot, mu1 - 3.0 * sg1)
    w2_hi = min(xhi_tot, mu1 + 3.0 * sg1)
    if w2_hi - w2_lo < max(ax.GetBinWidth(1) * 3.0, 1e-3):
        # Window collapsed; use the first fit as the result
        sigp = 100.0 * abs(f1.GetParameter(2))
        esigp = 100.0 * f1.GetParError(2)
        return {'ok': True, 'h': hfit, 'fn': f1, 'xmin': x1, 'xmax': x2,
                'sigp': sigp, 'esigp': esigp, 'mu': mu1, 'sigma': sg1, 'stage': 'first'}

    f2 = ROOT.TF1('g2_{}'.format(id(hfit)), 'gaus', w2_lo, w2_hi)
    f2.SetParameters(f1.GetParameter(0), mu1, sg1)
    r2 = hfit.Fit(f2, 'RQSN')
    ok2 = (int(r2) == 0)

    if ok2:
        sigp  = 100.0 * abs(f2.GetParameter(2))
        esigp = 100.0 * f2.GetParError(2)
        return {'ok': True, 'h': hfit, 'fn': f2, 'xmin': w2_lo, 'xmax': w2_hi,
                'sigp': sigp, 'esigp': esigp, 'mu': f2.GetParameter(1),
                'sigma': abs(f2.GetParameter(2)), 'stage': 'second'}
    else:
        # Fall back to first fit
        sigp  = 100.0 * abs(f1.GetParameter(2))
        esigp = 100.0 * f1.GetParError(2)
        return {'ok': True, 'h': hfit, 'fn': f1, 'xmin': x1, 'xmax': x2,
                'sigp': sigp, 'esigp': esigp, 'mu': mu1, 'sigma': sg1, 'stage': 'first'}

# ---------------------------
# Styling
# ---------------------------

def style_graph(gr, color, marker=20):
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetLineWidth(2)

COL = {
    'barrel': ROOT.kBlue + 1,
    'transition': ROOT.kOrange + 7,
    'endcap': ROOT.kRed + 1,
    'all': ROOT.kViolet + 1,
}

# high-contrast, colorblind-friendly palette for energies
ENERGY_PALETTE_HEX = [
    "#4E79A7", "#E15759", "#59A14F", "#B07AA1",
    "#F28E2B", "#76B7B2", "#EDC948", "#FF9DA7",
    "#9C755F", "#8CD17D", "#499894", "#D37295",
    "#86BCB6", "#FABFD2", "#C49C94", "#1F77B4",
]

def energy_color(idx: int) -> int:
    hexcol = ENERGY_PALETTE_HEX[idx % len(ENERGY_PALETTE_HEX)]
    return ROOT.TColor.GetColor(hexcol)

def draw_note_ignored_all(canvas):
    t = ROOT.TLatex()
    t.SetNDC(True)
    t.SetTextSize(0.032)
    t.SetTextColor(ROOT.kRed + 1)
    t.DrawLatex(0.13, 0.92, "NOTE: unlabeled files present; 'all' curves ignored")
    canvas.Modified()

# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Minimal plotting + CSV: resolution (respR_sel) and selected efficiencies."
    )
    ap.add_argument('files', nargs='+', help='Input ROOT files (one particle type per run).')
    ap.add_argument('-o', '--outdir', default='figs', help='Output directory (default: figs)')
    ap.add_argument('--histpath', default='PFOCluster', help='Histogram directory (default: PFOCluster)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fitdir = os.path.join(args.outdir, 'fit')
    os.makedirs(fitdir, exist_ok=True)

    print('[info] Using histpath hint "{}" (also trying /{} and /MC/{})'.format(
        args.histpath, args.histpath, args.histpath))

    File = collections.namedtuple('File', 'path region E particle h')
    files = []

    # ---------------- load files ----------------
    for p in args.files:
        f = ROOT.TFile.Open(p)
        if not f or f.IsZombie():
            print('[warn] cannot open: {} -- skipping.'.format(p))
            continue

        region = parse_region(p)
        particle = parse_particle(p)

        # Selected histos
        h_mcE_sel = get_hist(f, args.histpath, 'mcE_sel')
        h_respR_sel = get_hist(f, args.histpath, 'respR_sel')
        h_effE_num = get_hist(f, args.histpath, 'effE_sel_num')
        h_effE_den = get_hist(f, args.histpath, 'effE_sel_den')
        h_effCos_num = get_hist(f, args.histpath, 'effCos_sel_num')
        h_effCos_den = get_hist(f, args.histpath, 'effCos_sel_den')

        bn = os.path.basename(p)
        if not h_mcE_sel:
            print('[warn] {}: missing mcE_sel under {} -- skipping file.'.format(bn, args.histpath))
            f.Close()
            continue

        # prefer energy from filename; fall back to mcE_sel
        E_nam = parse_energy_from_name(p)
        E_mc = energy_from_mcE_sel(h_mcE_sel)
        if E_nam is not None:
            E = E_nam
        elif E_mc is not None:
            E = E_mc
        else:
            print('[warn] {}: cannot determine energy (no name match, mcE_sel empty) -- skipping.'.format(bn))
            f.Close()
            continue

        if (E_nam is not None) and (E_mc is not None) and (abs(E_mc - E_nam) / max(E_nam, 1e-9) > 0.05):
            print('[warn] {}: energy mismatch mcE_sel~{} vs name~{} GeV (using name).'.format(bn, E_mc, E_nam))

        print('[info] {}: particle={}, region={}, E~{} GeV (from {})'.format(
            bn, particle, region, E, 'name' if E_nam is not None else 'mcE_sel'))

        if not h_respR_sel or h_respR_sel.GetEntries() == 0:
            print('[warn] {}: respR_sel empty -- no resolution point will be plotted.'.format(bn))
        if (not h_effE_den) or h_effE_den.Integral() == 0:
            print('[warn] {}: effE_sel_den empty/missing -- efficiency vs E point will be skipped.'.format(bn))
        if (not h_effCos_den) or h_effCos_den.Integral() == 0:
            print('[warn] {}: effCos_sel_den empty/missing -- efficiency vs cosTheta will be skipped.'.format(bn))

        files.append(File(
            path=p, region=region, E=E, particle=particle,
            h={
                'mcE_sel': h_mcE_sel, 'respR_sel': h_respR_sel,
                'effE_num': h_effE_num, 'effE_den': h_effE_den,
                'effCos_num': h_effCos_num, 'effCos_den': h_effCos_den
            }
        ))

        # --- Fit preview with two-step Gaussian (no rebin) ---
        if h_respR_sel and h_respR_sel.GetEntries() > 0:
            fit = fit_gauss_twostep(h_respR_sel)
            c = ROOT.TCanvas('c_fit', 'fit', 700, 500)
            c.SetGrid()
            if fit['ok']:
                hdraw = fit['h'].Clone('respR_sel__{}'.format(bn))
                hdraw.SetDirectory(0)

                ax = hdraw.GetXaxis()
                ax_lo, ax_hi = ax.GetXmin(), ax.GetXmax()

                # Candidate plot window: mu +/- 6 sigma (twice wider than the +/-3 sigma fit)
                mu = float(fit.get('mu', 0.0))
                sg = abs(float(fit.get('sigma', 0.0)))
                use_full_axis = False
                xmin_plot, xmax_plot = fit['xmin'], fit['xmax']  # fallback to fit window

                inwin = 0.0
                total = hdraw.Integral(1, hdraw.GetNbinsX())

                if sg > 0.0:
                    cand_lo = max(ax_lo, mu - 6.0 * sg)
                    cand_hi = min(ax_hi, mu + 6.0 * sg)
                    # Ensure we at least include the actual fit window
                    cand_lo = min(cand_lo, fit['xmin'])
                    cand_hi = max(cand_hi, fit['xmax'])

                    # Check fraction inside candidate window
                    thr = 0.80
                    nb = hdraw.GetNbinsX()
                    i1 = max(1, ax.FindBin(cand_lo))
                    i2 = min(nb, ax.FindBin(cand_hi))
                    inwin = hdraw.Integral(i1, i2)
                    frac = (inwin / total) if total > 0 else 0.0

                    if frac >= thr and (cand_hi - cand_lo) > max(ax.GetBinWidth(1) * 3.0, 1e-3):
                        xmin_plot, xmax_plot = cand_lo, cand_hi
                    else:
                        use_full_axis = True
                        xmin_plot, xmax_plot = ax_lo, ax_hi
                else:
                    # No sensible sigma -> full axis
                    use_full_axis = True
                    xmin_plot, xmax_plot = ax_lo, ax_hi

                # y max over the chosen x-range
                i1 = max(1, ax.FindBin(xmin_plot))
                i2 = min(hdraw.GetNbinsX(), ax.FindBin(xmax_plot))
                ymax = 0.0
                for ib in range(i1, i2 + 1):
                    ymax = max(ymax, hdraw.GetBinContent(ib))
                if ymax <= 0.0:
                    ymax = max(1.0, hdraw.GetMaximum())

                frame = c.DrawFrame(xmin_plot, 0.0, xmax_plot, 1.15 * ymax)
                title_suffix = 'full axis' if use_full_axis else 'mu +/- 6 sigma'
                frame.SetTitle('{}, {}, E={} GeV; (E_rec/E_true-1) [{}]; Entries'.format(
                    particle, region, E, title_suffix))

                hdraw.Draw('HIST SAME')
                fit['fn'].SetLineColor(ROOT.kRed + 1)
                fit['fn'].Draw('SAME')

                box = ROOT.TPaveText(0.55, 0.70, 0.88, 0.88, 'NDC')
                box.SetFillStyle(0)
                box.SetBorderSize(0)
                box.AddText('#sigma = {:.2f} #pm {:.2f} %'.format(fit['sigp'], fit['esigp']))
                box.Draw()

                frac_print = 100.0 * (inwin / total) if total > 0 else 0.0
                print('[info] {}: first fit mu={:.4f}, sigma={:.4f}; used stage={}'.format(
                    bn, fit.get('mu', 0.0), fit.get('sigma', 0.0), fit.get('stage', '?')))
                print('[info] {}: preview window={}, entries fraction inside={:.1f}%'.format(
                    bn, ('full axis' if use_full_axis else '+/- 6 sigma'), frac_print))
                print('[info] {}: fit sigma={:.2f} +/- {:.2f} %'.format(bn, fit['sigp'], fit['esigp']))
            else:
                frame = c.DrawFrame(-0.3, 0.0, 0.3, 1.0)
                frame.SetTitle('{}, {}, E={} GeV; (E_rec/E_true-1); Entries'.format(particle, region, E))
                print('[warn] {}: Gaussian fit failed.'.format(bn))
            out = os.path.join(fitdir, '{}_{}_E{}GeV_fit.pdf'.format(particle, region, E))
            c.SaveAs(out)
            print('[write] {}'.format(out))
            del c

        f.Close()

    if not files:
        print('[error] no usable files -- nothing to plot.')
        sys.exit(1)

    # Determine particle label (assume homogeneous set)
    particle_label = files[0].particle if files else 'particle'

    # Presence of labeled regions and unlabeled
    has_region_tags = any(fi.region in ('barrel', 'transition', 'endcap') for fi in files)
    unlabeled_files = [fi for fi in files if fi.region == 'all']
    if unlabeled_files:
        for fi in unlabeled_files:
            print('[warn] unlabeled file (no barrel/transition/endcap): {}'.format(os.path.basename(fi.path)))
    if has_region_tags and unlabeled_files:
        print('[warn] Found labeled regions; will ignore "all" files in plots.')

    # group by region
    by_region = collections.defaultdict(list)
    for fi in files:
        by_region[fi.region].append(fi)

    # CSV collectors
    res_csv_rows = []
    effE_csv_rows = []
    effCt_csv_rows = []

    # ---------------- Resolution vs E ----------------
    c_res = ROOT.TCanvas('c_res', 'Resolution vs E', 900, 700)
    c_res.SetGrid()
    leg = ROOT.TLegend(0.58, 0.18, 0.88, 0.38)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    legend_has_entries = False
    graphs = {}
    xmax = 0.0

    # Target bands by particle type
    if particle_label in ('gamma', 'photon'):
        target_lo, target_hi = 0.0, 8.0
    elif particle_label in ('kaon0l', 'k0l', 'kaon', 'pion', 'pi', 'neutron'):
        target_lo, target_hi = 5.0, 25.0  # hadrons
    else:
        target_lo, target_hi = 0.0, 25.0

    # Collect all y values to decide whether to use the target band
    all_y = []

    for region, fis in by_region.items():
        if has_region_tags and region == 'all':
            continue
        pts = []
        for fi in sorted(fis, key=lambda x: x.E):
            h = fi.h['respR_sel']
            if not h or h.GetEntries() == 0:
                print('[info] {}: no respR_sel; skipping in resolution.'.format(os.path.basename(fi.path)))
                continue
            fit = fit_gauss_twostep(h)
            if not fit['ok']:
                print('[warn] {}: fit failed; skipping point.'.format(os.path.basename(fi.path)))
                continue
            sigp = fit['sigp']
            esigp = fit['esigp']
            pts.append((fi.E, sigp, 0.0, esigp))
            all_y.append(sigp)
            res_csv_rows.append([fi.particle, region, '{:g}'.format(fi.E),
                                 '{:.4f}'.format(sigp), '{:.4f}'.format(esigp),
                                 os.path.basename(fi.path)])
            xmax = max(xmax, fi.E)
        if not pts:
            print('[info] region {}: no resolution points.'.format(region))
            continue
        gr = ROOT.TGraphErrors(len(pts))
        for i, (x, y, ex, ey) in enumerate(pts):
            gr.SetPoint(i, x, y)
            gr.SetPointError(i, ex, ey)
        style_graph(gr, COL.get(region, ROOT.kGray + 2), marker=20)
        graphs[region] = gr
        if region != 'all':
            leg.AddEntry(gr, region, 'lp')
            legend_has_entries = True

    # Decide y-axis based on fraction of points inside the target band
    ylo_plot, yhi_plot = target_lo, target_hi
    thr_frac = 0.80  # "most of it" threshold
    if all_y:
        inside = sum(1 for y in all_y if (target_lo <= y <= target_hi))
        frac = inside / float(len(all_y))
        if frac >= thr_frac:
            print('[info] resolution y-axis: using target band [{:.1f}, {:.1f}]% ({} of {} points inside, {:.0f}%)'
                  .format(target_lo, target_hi, inside, len(all_y), 100*frac))
            ylo_plot, yhi_plot = target_lo, target_hi
        else:
            ymin = min(all_y)
            ymax = max(all_y)
            span = max(1e-6, ymax - ymin)
            pad = 0.15 * span
            ylo_plot = max(0.0, math.floor((ymin - pad) / 0.5) * 0.5)
            yhi_plot = math.ceil((ymax + pad) / 0.5) * 0.5
            print('[info] resolution y-axis: auto range [{:.2f}, {:.2f}]% (only {:.0f}% inside target band)'
                  .format(ylo_plot, yhi_plot, 100*frac))
    else:
        print('[info] resolution y-axis: no points; keeping default band [{:.1f}, {:.1f}]%'.format(target_lo, target_hi))

    frame = c_res.DrawFrame(0.0, ylo_plot, max(1.0, 1.05 * xmax), yhi_plot)
    frame.SetTitle('{}: energy resolution;E [GeV];#sigma(E)/E [%]'.format(particle_label))
    for region, gr in graphs.items():
        if has_region_tags and region == 'all':
            continue
        gr.Draw('P SAME')
    if legend_has_entries:
        leg.SetHeader('Regions', 'C')
        leg.Draw()
    if has_region_tags and unlabeled_files:
        draw_note_ignored_all(c_res)
    out_res = os.path.join(args.outdir, '{}_resolution_vs_energy_by_region.pdf'.format(particle_label))
    c_res.SaveAs(out_res)
    print('[write] {}'.format(out_res))

    # ---------------- Efficiency vs E ----------------
    c_effE = ROOT.TCanvas('c_effE', 'Efficiency vs E', 900, 700)
    c_effE.SetGrid()
    legE = ROOT.TLegend(0.58, 0.18, 0.88, 0.38)
    legE.SetFillStyle(0)
    legE.SetBorderSize(0)
    legendE_has_entries = False
    graphsE = {}
    xmax, ymin, ymax = 0.0, 100.0, 0.0

    for region, fis in by_region.items():
        if has_region_tags and region == 'all':
            continue
        pts = []
        for fi in sorted(fis, key=lambda x: x.E):
            eff, err, nnum, nden = efficiency_overall(fi.h['effE_num'], fi.h['effE_den'])
            if eff is None:
                print('[info] {}: missing/empty effE_sel_den -- skipping eff(E) point.'.format(os.path.basename(fi.path)))
                continue
            y = 100.0 * eff
            ye = 100.0 * err
            pts.append((fi.E, y, 0.0, ye))
            effE_csv_rows.append([fi.particle, region, '{:g}'.format(fi.E),
                                  '{:.4f}'.format(y), '{:.4f}'.format(ye),
                                  '{:.0f}'.format(nnum), '{:.0f}'.format(nden),
                                  os.path.basename(fi.path)])
            xmax = max(xmax, fi.E)
            ymin = min(ymin, y - ye)
            ymax = max(ymax, y + ye)
            print('[info] eff(E) {} @ {:g} GeV = {:.1f} +/- {:.1f} % (num={:.0f}, den={:.0f})'.format(
                region, fi.E, y, ye, nnum, nden))
        if not pts:
            print('[info] region {}: no efficiency(E) points.'.format(region))
            continue
        gr = ROOT.TGraphErrors(len(pts))
        for i, (x, y, ex, ey) in enumerate(pts):
            gr.SetPoint(i, x, y)
            gr.SetPointError(i, ex, ey)
        style_graph(gr, COL.get(region, ROOT.kGray + 2), marker=21)
        graphsE[region] = gr
        if region != 'all':
            legE.AddEntry(gr, region, 'lp')
            legendE_has_entries = True

    ylo, yhi = 85.0, 100.5
    if ymin < ylo:
        ylo = max(0.0, math.floor((ymin - 2.5) / 5.0) * 5.0)
    if ymax > yhi:
        yhi = min(120.0, math.ceil((ymax + 2.5) / 5.0) * 5.0)

    frame = c_effE.DrawFrame(0.0, ylo, max(1.0, 1.05 * xmax), yhi)
    frame.SetTitle('{}: reconstruction efficiency;E [GeV];Efficiency [%]'.format(particle_label))
    for region, gr in graphsE.items():
        if has_region_tags and region == 'all':
            continue
        gr.Draw('P SAME')
    if legendE_has_entries:
        legE.SetHeader('Regions', 'C')
        legE.Draw()
    if has_region_tags and unlabeled_files:
        draw_note_ignored_all(c_effE)
    out_effE = os.path.join(args.outdir, '{}_efficiency_vs_energy_by_region_sel.pdf'.format(particle_label))
    c_effE.SaveAs(out_effE)
    print('[write] {}'.format(out_effE))

    # ---------------- Efficiency vs cosTheta (one canvas per region) ----------------
    for region, fis in by_region.items():
        if has_region_tags and region == 'all':
            continue

        energies = sorted({fi.E for fi in fis})
        if not energies:
            print('[info] region {}: no energies for efficiency vs cosTheta.'.format(region))
            continue

        c = ROOT.TCanvas('c_effCos_{}'.format(region), 'efficiency vs cosTheta [{}]'.format(region), 900, 700)
        c.SetGrid()
        legC = ROOT.TLegend(0.58, 0.18, 0.88, 0.38)
        legC.SetFillStyle(0)
        legC.SetBorderSize(0)
        if region != 'all':
            legC.SetHeader(region, 'C')
        graphsC = []
        ymin, ymax = 100.0, 0.0

        any_den = next((fi.h['effCos_den'] for fi in fis if fi.h['effCos_den']), None)
        if not any_den:
            print('[info] region {}: no effCos_sel_den -- skipping canvas.'.format(region))
            continue
        xmin = any_den.GetXaxis().GetXmin()
        xmaxx = any_den.GetXaxis().GetXmax()

        energy_to_color = {e: energy_color(i) for i, e in enumerate(energies)}

        for E in energies:
            fi = next((x for x in fis if abs(x.E - E) < 1e-9), None)
            if not fi:
                continue
            pts = efficiency_vs_costheta(fi.h['effCos_num'], fi.h['effCos_den'])
            if not pts:
                print('[info] {}: no per-bin denom -- skip E={} GeV.'.format(os.path.basename(fi.path), E))
                continue
            gr = ROOT.TGraphErrors(len(pts))
            for i, (xv, ex, ev, ee, num, den) in enumerate(pts):
                yv = 100.0 * ev
                ey = 100.0 * ee
                gr.SetPoint(i, xv, yv)
                gr.SetPointError(i, ex, ey)
                ymin = min(ymin, yv - ey)
                ymax = max(ymax, yv + ey)
                effCt_csv_rows.append([fi.particle, region, '{:g}'.format(fi.E),
                                       '{:.6f}'.format(xv), '{:.6f}'.format(ex),
                                       '{:.4f}'.format(yv), '{:.4f}'.format(ey),
                                       '{:.0f}'.format(num), '{:.0f}'.format(den),
                                       os.path.basename(fi.path)])
            color = energy_to_color[E]
            style_graph(gr, color, marker=20)
            graphsC.append(gr)
            legC.AddEntry(gr, 'E={} GeV'.format(E), 'lp')

        if not graphsC:
            print('[info] region {}: nothing to draw for efficiency vs cosTheta.'.format(region))
            continue

        ylo, yhi = 85.0, 100.5
        if ymin < ylo:
            ylo = max(0.0, math.floor((ymin - 2.5) / 5.0) * 5.0)
        if ymax > yhi:
            yhi = min(120.0, math.ceil((ymax + 2.5) / 5.0) * 5.0)

        frame = c.DrawFrame(xmin, ylo, xmaxx, yhi)
        frame.SetTitle('{}: efficiency vs cos#theta;cos#theta;Efficiency [%]'.format(particle_label))
        for gr in graphsC:
            gr.Draw('P SAME')
        legC.Draw()
        if has_region_tags and unlabeled_files:
            draw_note_ignored_all(c)
        out_effC = os.path.join(args.outdir, '{}_efficiency_vs_costheta_{}_sel.pdf'.format(particle_label, region))
        c.SaveAs(out_effC)
        print('[write] {}'.format(out_effC))

    # ---------------- Write CSVs ----------------
    csv_res = os.path.join(args.outdir, '{}_resolution_vs_energy.csv'.format(particle_label))
    csv_effE = os.path.join(args.outdir, '{}_efficiency_vs_energy_sel.csv'.format(particle_label))
    csv_effCt = os.path.join(args.outdir, '{}_efficiency_vs_costheta_sel.csv'.format(particle_label))

    with open(csv_res, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['particle', 'region', 'energy_GeV', 'sigma_percent', 'sigma_err_percent', 'file'])
        w.writerows(res_csv_rows)
    print('[write] {}  ({} rows)'.format(csv_res, len(res_csv_rows)))

    with open(csv_effE, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['particle', 'region', 'energy_GeV', 'eff_percent', 'eff_err_percent',
                    'num_integral', 'den_integral', 'file'])
        w.writerows(effE_csv_rows)
    print('[write] {}  ({} rows)'.format(csv_effE, len(effE_csv_rows)))

    with open(csv_effCt, 'w', newline='') as fp:
        w = csv.writer(fp)
        w.writerow(['particle', 'region', 'energy_GeV', 'cosTheta', 'cosTheta_halfwidth',
                    'eff_percent', 'eff_err_percent', 'num_bin', 'den_bin', 'file'])
        w.writerows(effCt_csv_rows)
    print('[write] {}  ({} rows)'.format(csv_effCt, len(effCt_csv_rows)))

    print('[done] all plots and CSVs written.')

if __name__ == '__main__':
    main()
