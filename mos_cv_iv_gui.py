#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MOS C–V (HF/LF) & I–V GUI — Shichman–Hodges nMOS + simple MOS-C model
Author: ChatGPT (for educational use)
Updates:
- 圖上自動標註：Accumulation/Depletion/Inversion區、Vth 與三極/飽和分界點（每條 Vgs 的拐點以 ● 標記）。
"""
import sys, math
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from dataclasses import dataclass
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import ticker
# 字型與負號
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = ['Microsoft JhengHei','Noto Sans CJK TC','PingFang TC','Arial Unicode MS','DejaVu Sans']
mpl.rcParams['axes.unicode_minus'] = False
from matplotlib.figure import Figure

q = 1.602176634e-19
k_B = 1.380649e-23
T = 300.0
Vt_th = k_B*T/q

eps0 = 8.854187817e-12     # F/m
eps_ox = 3.9*eps0
eps_si = 11.7*eps0
ni_Si = 1.0e16  # 1/m^3  (≈1e10 cm^-3)


def _beautify_linear(ax, xlim=None, ylim=None, nx=6, ny=6, sci_y=True, title=None, xlabel=None, ylabel=None):
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nx))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(ny))
    if sci_y:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    ax.grid(True, which='major', alpha=0.35)
    ax.grid(True, which='minor', alpha=0.12)
    if title: ax.set_title(title)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.margins(x=0.02, y=0.05)

def _beautify_logy(ax, xlim=None, ylim=None, nx=6, base=10, title=None, xlabel=None, ylabel=None):
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)
    ax.yaxis.set_major_locator(ticker.LogLocator(base=base))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=base, subs=np.arange(2, base), numdecs=base-2))
    ax.yaxis.set_major_formatter(ticker.LogFormatterExponent(base=base))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nx))
    ax.grid(True, which='major', alpha=0.35)
    ax.grid(True, which='minor', alpha=0.12)
    if title: ax.set_title(title)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.margins(x=0.02, y=0.05)

@dataclass
class Params:
    NA_cm3: float = 1e16      # cm^-3
    tox_nm: float = 5.0       # nm
    mu_n_cm: float = 450.0    # cm^2/Vs
    W_um: float = 10.0        # micrometer
    L_um: float = 1.0         # micrometer
    VFB: float = -0.9         # V
    Vg_min: float = -2.0
    Vg_max: float = 2.0
    Vds_fixed: float = 1.0    # for transfer curve
    vds_sweep_max: float = 2.0
    n_vgs_curves: int = 5

def um_to_m(x_um): return x_um*1e-6
def nm_to_m(x_nm): return x_nm*1e-9

def build_device(p: Params):
    NA = p.NA_cm3*1e6
    tox = nm_to_m(p.tox_nm)
    mu = p.mu_n_cm*1e-4
    W = um_to_m(p.W_um)
    L = um_to_m(p.L_um)

    Cox = eps_ox/tox
    phiF = Vt_th*math.log(NA/ni_Si)
    gamma = math.sqrt(2*q*eps_si*NA)/Cox
    Vth = p.VFB + 2*phiF + gamma*math.sqrt(2*phiF)
    beta = mu*Cox*(W/L)
    return dict(NA=NA, tox=tox, mu=mu, W=W, L=L, Cox=Cox, phiF=phiF, gamma=gamma, Vth=Vth, beta=beta)

def cv_curve(p: Params, npts=401):
    d = build_device(p)
    Vg = np.linspace(p.Vg_min, p.Vg_max, npts)
    Cox, phiF, gamma = d["Cox"], d["phiF"], d["gamma"]
    psi_grid = np.linspace(-3*abs(phiF), 2*abs(phiF), npts)
    Vg_map = lambda ps: p.VFB + ps + gamma*np.sqrt(np.maximum(ps,0.0))
    Vg_ps = Vg_map(psi_grid)
    psi_of_Vg = np.interp(Vg, Vg_ps, psi_grid, left=psi_grid[0], right=psi_grid[-1])
    Wd = np.sqrt(np.maximum(2*eps_si*np.maximum(psi_of_Vg,0.0)/(q*d["NA"]), 1e-30))
    Cs = eps_si/np.maximum(Wd,1e-30)
    C_HF = (Cox*Cs)/(Cox+Cs)
    C_LF = C_HF.copy()
    mask_inv = psi_of_Vg >= 2*abs(phiF)
    C_LF[mask_inv] = Cox
    return Vg, C_HF, C_LF, d

def iv_output(p: Params, nVds=161):
    d = build_device(p)
    beta, Vth = d["beta"], d["Vth"]
    Vds = np.linspace(0, p.vds_sweep_max, nVds)
    curves = []
    Vgs_list = np.linspace(max(Vth-0.5, p.Vg_min), p.Vg_max, p.n_vgs_curves)
    for Vgs in Vgs_list:
        Id = np.zeros_like(Vds)
        Vov = Vgs - Vth
        for i, vd in enumerate(Vds):
            if Vov <= 0:
                Id[i] = 0.0
            elif vd < Vov:
                Id[i] = beta*(Vov*vd - 0.5*vd*vd)
            else:
                Id[i] = 0.5*beta*Vov*Vov
        curves.append((Vgs, Id, Vov))
    return Vds, curves, d

def iv_transfer(p: Params, nVgs=201):
    d = build_device(p)
    beta, Vth = d["beta"], d["Vth"]
    Vgs = np.linspace(p.Vg_min, p.Vg_max, nVgs)
    Id = np.zeros_like(Vgs)
    for i, vg in enumerate(Vgs):
        Vov = vg - Vth
        if Vov <= 0:
            Id[i] = 1e-18
        elif p.Vds_fixed < Vov:
            Id[i] = beta*(Vov*p.Vds_fixed - 0.5*p.Vds_fixed*p.Vds_fixed)
        else:
            Id[i] = 0.5*beta*Vov*Vov
    return Vgs, Id, d

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MOS C–V & I–V 模擬器 (HF/LF + SH nMOS) — 含標註")
        self.p = Params()
        self._build_ui()
        self.update_plots()

    def _build_ui(self):
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        frm = ttk.Frame(self, padding=6)
        frm.grid(row=0, column=0, sticky="ns")

        def add_slider(row, label, frm_parent, from_, to, get, set, fmt="{:.2f}", length=180):
            ttk.Label(frm_parent, text=label, width=22, anchor="w").grid(row=row, column=0, sticky="w")
            var = tk.DoubleVar(value=get())
            s = ttk.Scale(frm_parent, from_=from_, to=to, value=get(), orient="horizontal", length=length,
                          command=lambda v:(set(float(v)), self.update_plots()))
            s.grid(row=row, column=1, sticky="ew", padx=6)
            val_lbl = ttk.Label(frm_parent, text=fmt.format(get()))
            val_lbl.grid(row=row, column=2, sticky="e")
            def on_move(v): val_lbl.config(text=fmt.format(float(v)))
            s.configure(command=lambda v:(set(float(v)), on_move(v), self.update_plots()))
            return s

        r = 0
        add_slider(r:=r+1, "Doping NA (1e15–1e18 cm⁻³)", frm, 1e15, 1e18, lambda:self.p.NA_cm3, lambda v:setattr(self.p,"NA_cm3",v), fmt="{:.3g}")
        add_slider(r:=r+1, "Oxide thickness tox (nm)", frm, 1.0, 20.0, lambda:self.p.tox_nm, lambda v:setattr(self.p,"tox_nm",v))
        add_slider(r:=r+1, "Electron mobility μn (cm²/V·s)", frm, 100.0, 800.0, lambda:self.p.mu_n_cm, lambda v:setattr(self.p,"mu_n_cm",v), fmt="{:.1f}")
        add_slider(r:=r+1, "W (µm)", frm, 1.0, 100.0, lambda:self.p.W_um, lambda v:setattr(self.p,"W_um",v), fmt="{:.1f}")
        add_slider(r:=r+1, "L (µm)", frm, 0.1, 10.0, lambda:self.p.L_um, lambda v:setattr(self.p,"L_um",v), fmt="{:.2f}")
        add_slider(r:=r+1, "VFB (V)", frm, -1.5, 0.5, lambda:self.p.VFB, lambda v:setattr(self.p,"VFB",v))
        add_slider(r:=r+1, "Vg min (V)", frm, -5.0, 0.0, lambda:self.p.Vg_min, lambda v:setattr(self.p,"Vg_min",v))
        add_slider(r:=r+1, "Vg max (V)", frm, 0.0, 5.0, lambda:self.p.Vg_max, lambda v:setattr(self.p,"Vg_max",v))
        add_slider(r:=r+1, "Vds for transfer (V)", frm, 0.05, 3.0, lambda:self.p.Vds_fixed, lambda v:setattr(self.p,"Vds_fixed",v))
        add_slider(r:=r+1, "Vds sweep max (output) (V)", frm, 0.2, 5.0, lambda:self.p.vds_sweep_max, lambda v:setattr(self.p,"vds_sweep_max",v))

        self.btn_save = ttk.Button(frm, text="Export three figures (PNG)", command=self.export_pngs)
        self.btn_save.grid(row=r+2, column=0, columnspan=3, sticky="ew", pady=(8,4))

        self.fig = Figure(figsize=(8.5, 6.0), dpi=100)
        self.ax_cv = self.fig.add_subplot(2,2,1)
        self.ax_out = self.fig.add_subplot(2,2,2)
        self.ax_tr = self.fig.add_subplot(2,1,2)
        self.fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.15, wspace=0.35, hspace=0.35)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.16, wspace=0.35, hspace=0.38)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=1, sticky="nsew", padx=6, pady=6)

        self.rowconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

    def update_plots(self):
        # C–V with region shading
        Vg, C_HF, C_LF, d = cv_curve(self.p)
        Cox, Vth, VFB = d["Cox"], d["Vth"], self.p.VFB
        self.ax_cv.clear()
        scale = 1e5  # F/m^2 -> nF/cm^2
        self.ax_cv.plot(Vg, C_HF*scale, label="C-V (HF)")
        self.ax_cv.plot(Vg, C_LF*scale, label="C-V (LF)", linestyle="--")
        self.ax_cv.axhline(Cox*scale, linestyle=":", label="Cox")
        # 區域：Accumulation (Vg<<VFB)、Depletion (VFB~Vth)、Inversion (Vg>Vth)
        ymin, ymax = self.ax_cv.get_ylim()

        self.ax_cv.axvline(VFB, linestyle=":", linewidth=1)
        self.ax_cv.axvline(Vth, linestyle="--", linewidth=1)
        # Arrow annotations for VFB and Vth
        self.ax_cv.annotate("VFB", xy=(VFB, ymax*0.98), xytext=(VFB-0.1*(self.p.Vg_max-self.p.Vg_min), ymax*1.05),
                            textcoords="data", ha="right", va="bottom",
                            arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=5, lw=0.8))
        self.ax_cv.annotate("Vth (= VFB + 2φF + γ√(2φF))", 
                            xy=(Vth, ymax*0.98), 
                            xytext=(Vth+0.1*(self.p.Vg_max-self.p.Vg_min), ymax*1.05),
                            textcoords="data", ha="left", va="bottom",
                            arrowprops=dict(arrowstyle="->", shrinkA=0, shrinkB=5, lw=0.8))
        # 2φF label near Vth top
        self.ax_cv.text(Vth, ymax*0.90, "2φF", ha="center", va="top", fontsize=9)
        # Region captions
        self.ax_cv.text((min(Vg)+VFB)/2, ymax*0.80, "Accumulation", ha="center", fontsize=9)
        self.ax_cv.text((VFB+Vth)/2,      ymax*0.80, "Depletion",    ha="center", fontsize=9)
        self.ax_cv.text((Vth+max(Vg))/2,  ymax*0.80, "Inversion",    ha="center", fontsize=9)

        Cv_title = 'MOS C–V (per area)'
        _beautify_linear(self.ax_cv, xlim=(self.p.Vg_min, self.p.Vg_max), ylim=(0, Cox*1e5*1.15), title=Cv_title, xlabel='Vg (V)', ylabel='C (nF/cm²)', sci_y=False)
        self.ax_cv.legend(loc='lower left', framealpha=0.85, fontsize=8)

        # I–V output with knee markers
        Vds, curves, d2 = iv_output(self.p)
        self.ax_out.clear()
        self.ax_out.set_facecolor('white')
        for Vgs, Id, Vov in curves:
            Id_mA = Id*1e3
            self.ax_out.plot(Vds, Id_mA, label=f"Vgs={Vgs:.2f} V")
            if Vov > 0:
                knee = min(Vov, Vds[-1])
                idx = np.argmin(np.abs(Vds - knee))
                self.ax_out.plot([Vds[idx]], [Id_mA[idx]], marker="o")
        max_id = max([np.max(Id*1e3) for _, Id, _ in curves] + [1e-9])
        _beautify_linear(self.ax_out, xlim=(0, self.p.vds_sweep_max), ylim=(0, max_id*1.15), title='nMOS 輸出特性 Id–Vds(●: saturation knee)', xlabel='Vds (V)', ylabel='Id (mA)', sci_y=False)
        self.ax_out.legend(fontsize=8, ncols=2, loc='upper left', framealpha=0.85)

        # I–V transfer with Vth line + region label
        Vgs, Idtr, d3 = iv_transfer(self.p)
        self.ax_tr.clear()
        self.ax_tr.set_facecolor('white')
        self.ax_tr.semilogy(Vgs, Idtr*1e3)
        self.ax_tr.axvline(d3["Vth"], linestyle="--", linewidth=1)
        y_min, y_max = self.ax_tr.get_ylim()
        self.ax_tr.text(self.p.Vg_min + 0.05*(self.p.Vg_max-self.p.Vg_min), y_max/3, "Subthreshold", fontsize=9)
        self.ax_tr.text(d3["Vth"]+0.05*(self.p.Vg_max-self.p.Vg_min), y_max/3, "強Inversion", fontsize=9)
        _beautify_logy(self.ax_tr, xlim=(self.p.Vg_min, self.p.Vg_max), title=f'nMOS 轉移特性 Id–Vgs @ Vds={self.p.Vds_fixed:.2f} V（dashed: Vth）', xlabel='Vgs (V)', ylabel='Id (mA)')

        footer = f"Vth = {d['Vth']:.3f} V, β = {d['beta']:.3e} A/V², Cox = {d['Cox']:.3e} F/m²"
        self.fig.suptitle(footer, fontsize=9)
        self.canvas.draw()

    def export_pngs(self):
        try:
            folder = filedialog.askdirectory(title="Choose output folder")
            if not folder: return
            self.update_plots()
            import matplotlib.pyplot as plt
            # CV
            fig1 = Figure(figsize=(6,4), dpi=150); ax1 = fig1.add_subplot(1,1,1)
            Vg, C_HF, C_LF, d = cv_curve(self.p)
            scale = 1e5
            ax1.plot(Vg, C_HF*scale, label="C-V (HF)")
            ax1.plot(Vg, C_LF*scale, label="C-V (LF)", linestyle="--")
            Cox, Vth, VFB = d["Cox"], d["Vth"], self.p.VFB
            ax1.axhline(Cox*scale, linestyle=":", label="Cox")
            ymin, ymax = ax1.get_ylim()
            ax1.axvline(VFB, linestyle=":", linewidth=1); ax1.axvline(Vth, linestyle="--", linewidth=1)
            _beautify_linear(ax1, xlim=(self.p.Vg_min, self.p.Vg_max), ylim=(0, Cox*1e5*1.15), xlabel='Vg (V)', ylabel='C (nF/cm²)', sci_y=False); ax1.legend(loc='upper left', fontsize=8, framealpha=0.85)
            fig1.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.16); fig1.savefig(folder + "/cv_curve.png")
            # IV output
            fig2 = Figure(figsize=(6,4), dpi=150); ax2 = fig2.add_subplot(1,1,1)
            ax1.set_facecolor('white')
            ax2.set_facecolor('white')
            Vds, curves, d2 = iv_output(self.p)
            for Vgs, Id, Vov in curves:
                Id_mA = Id*1e3
                ax2.plot(Vds, Id_mA, label=f"Vgs={Vgs:.2f} V")
                if Vov>0:
                    knee = min(Vov, Vds[-1]); idx = np.argmin(np.abs(Vds-knee)); ax2.plot([Vds[idx]],[Id[idx]],"o")
            max_id = max([np.max(Id) for _, Id, _ in curves] + [1e-12])
            _beautify_linear(ax2, xlim=(0, self.p.vds_sweep_max), ylim=(0, max_id*1.15), xlabel='Vds (V)', ylabel='Id (mA)'); ax2.legend(fontsize=8, ncols=2, loc='upper left', framealpha=0.85)
            fig2.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.16); fig2.savefig(folder + "/iv_output.png")
            # IV transfer
            fig3 = Figure(figsize=(6,4), dpi=150); ax3 = fig3.add_subplot(1,1,1)
            ax3.set_facecolor('white')
            Vgs, Idtr, d3 = iv_transfer(self.p)
            ax3.semilogy(Vgs, Idtr*1e3); ax3.axvline(d3['Vth'], linestyle='--', linewidth=1)
            ax3.text(0.05, 0.18, 'Subthreshold', transform=ax3.transAxes, fontsize=9)
            ax3.text(0.55, 0.18, 'Strong inversion', transform=ax3.transAxes, fontsize=9)
            _beautify_logy(ax3, xlim=(self.p.Vg_min, self.p.Vg_max), xlabel='Vgs (V)', ylabel='Id (mA)')
            fig3.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.16); fig3.savefig(folder + "/iv_transfer.png")
            messagebox.showinfo("Done", "Exported cv_curve.png, iv_output.png, iv_transfer.png")
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    App().mainloop()
