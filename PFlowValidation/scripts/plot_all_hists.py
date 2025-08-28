#!/usr/bin/env python3
import sys, os, argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(1110)

def mkdir_p(p): os.makedirs(p, exist_ok=True)

def is_hist(obj):
    return isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.TH2) or isinstance(obj, ROOT.TH3)

def list_keys_recursive(dirobj, prefix=""):
    """Yield (path, obj) for all objects under a TDirectory (recursively)."""
    out = []
    for key in dirobj.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        path = f"{prefix}/{name}" if prefix else name
        if isinstance(obj, ROOT.TDirectory):
            out.extend(list_keys_recursive(obj, path))
        else:
            out.append((path, obj))
    return out

def draw_th1(h, png, pdf=None):
    c = ROOT.TCanvas("c1","",900,700)
    h.SetLineWidth(2)
    h.Draw("hist")
    c.SaveAs(png)
    if pdf: c.SaveAs(pdf+"(" if not os.path.exists(pdf) else pdf)
    c.Close()

def draw_th2(h, png, pdf=None):
    c = ROOT.TCanvas("c2","",900,800)
    c.SetRightMargin(0.15)
    h.Draw("COLZ")
    c.SaveAs(png)
    if pdf: c.SaveAs(pdf+"(" if not os.path.exists(pdf) else pdf)
    c.Close()

def draw_teff(num, den, title, png, pdf=None):
    # Ensure num<=den bin-by-bin for TEfficiency
    numc = num.Clone(num.GetName()+"_clip")
    denc = den.Clone(den.GetName()+"_copy")
    for ib in range(1, denc.GetNbinsX()+1):
        if numc.GetBinContent(ib) > denc.GetBinContent(ib):
            numc.SetBinContent(ib, denc.GetBinContent(ib))
    eff = ROOT.TEfficiency(numc, denc)
    c = ROOT.TCanvas("ceff","",900,700)
    eff.SetTitle(title)
    eff.Draw()
    c.SaveAs(png)
    if pdf: c.SaveAs(pdf+"(" if not os.path.exists(pdf) else pdf)
    c.Close()
    numc.Delete(); denc.Delete()

def draw_th3(h3, base_name, outdir, pdf=None, max_slices=30):
    # 1) XY projection (sum over Z)
    hxy = h3.Project3D("xy")
    out_xy = os.path.join(outdir, base_name + "_xy.png")
    draw_th2(hxy, out_xy, pdf)

    # 2) per-Z-bin XY slices (cap to max_slices for huge histos)
    nz = h3.GetNbinsZ()
    step = max(1, nz // max_slices)
    for iz in range(1, nz+1, step):
        h3.GetZaxis().SetRange(iz, iz)
        hslice = h3.Project3D("xy")
        zlow = h3.GetZaxis().GetBinLowEdge(iz)
        zup  = h3.GetZaxis().GetBinUpEdge(iz)
        out = os.path.join(outdir, f"{base_name}_z{iz}_[{zlow:.3g},{zup:.3g}].png")
        draw_th2(hslice, out, pdf)
    h3.GetZaxis().SetRange(1, nz)  # reset

def main():
    ap = argparse.ArgumentParser(description="Plot TH1/TH3 from k4Performance ROOT files and build efficiencies from *_num/_den.")
    ap.add_argument("rootfile", help="Input ROOT file (from THistSvc)")
    ap.add_argument("outdir",   help="Output directory for PNGs")
    ap.add_argument("--pdf",    help="Optional multi-page PDF to collect all plots", default=None)
    ap.add_argument("--max-slices", type=int, default=30, help="Max Z-slices to draw per TH3 (default: 30)")
    args = ap.parse_args()

    if not os.path.isfile(args.rootfile):
        print(f"ERROR: file not found: {args.rootfile}")
        sys.exit(2)
    mkdir_p(args.outdir)

    f = ROOT.TFile.Open(args.rootfile)
    if not f or f.IsZombie():
        print(f"ERROR: cannot open {args.rootfile}")
        sys.exit(3)

    # Gather histograms
    entries = [(p, o) for (p, o) in list_keys_recursive(f) if is_hist(o)]

    # Map path->hist for easy lookup
    name2hist = {p: o for (p, o) in entries}

    # First pass: draw all TH1 and TH3 (and TH2 if present), keep *_num/_den too
    for path, h in entries:
        safe_dir = os.path.join(args.outdir, os.path.dirname(path))
        mkdir_p(safe_dir)
        base = os.path.basename(path)
        stem = base  # no extension
        outpng = os.path.join(safe_dir, stem + ".png")

        if isinstance(h, ROOT.TH3):
            draw_th3(h, stem, safe_dir, args.pdf, args.max_slices)
        elif isinstance(h, ROOT.TH2):
            draw_th2(h, outpng, args.pdf)
        else:  # TH1
            # Always save a PNG for reference (even for *_num/_den)
            draw_th1(h, outpng, args.pdf)

    # Second pass: build efficiencies from *_num / *_den pairs in the SAME directory
    # Pairing rules:
    #   <dir>/<stem>_num  with  <dir>/<stem>_den
    pairs = []
    for path in list(name2hist.keys()):
        if not path.endswith("_num"): continue
        den_path = path[:-4] + "den"
        if den_path in name2hist:
            # Only TH1 pairs make sense here
            if isinstance(name2hist[path], ROOT.TH1) and not isinstance(name2hist[path], ROOT.TH2) and not isinstance(name2hist[path], ROOT.TH3):
                if isinstance(name2hist[den_path], ROOT.TH1) and not isinstance(name2hist[den_path], ROOT.TH2) and not isinstance(name2hist[den_path], ROOT.TH3):
                    pairs.append((path, den_path))

    for num_path, den_path in pairs:
        num = name2hist[num_path]
        den = name2hist[den_path]
        outdir = os.path.join(args.outdir, os.path.dirname(num_path))
        mkdir_p(outdir)
        stem = os.path.basename(num_path)[:-4]  # remove "_num"
        title = den.GetTitle() if den.GetTitle() else stem
        outpng = os.path.join(outdir, stem + "_eff.png")
        draw_teff(num, den, title, outpng, args.pdf)

    # Close multi-page PDF if requested (ROOT uses ")" on the last SaveAs)
    if args.pdf and os.path.exists(args.pdf):
        c = ROOT.TCanvas()
        c.SaveAs(args.pdf + ")")
        c.Close()

    print(f"Done. PNGs in: {args.outdir}" + (f" and PDF: {args.pdf}" if args.pdf else ""))

if __name__ == "__main__":
    main()
