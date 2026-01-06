#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Notebook Launcher（不改原始 code）
- 6 個分頁（1~6）
- monkey patch：把 tk.Tk() 替換成可嵌入的 Frame
- Lazy load：點到分頁才載入
- 2_MOS：載入後修正 ParamRow 的 grid 權重/Entry sticky，減少灰色空白區（不改 2_MOS.py 原檔）
- ✅ PyInstaller 相容：可在 onefile(_MEIPASS) / onedir / 原始 .py 直接跑
"""

import os
import sys
import runpy
import tkinter as tk
from tkinter import ttk, messagebox


# ============================================================
#  PyInstaller 資源路徑工具
# ============================================================
def _runtime_base_dir() -> str:
    """
    - 原始跑 .py：回傳 final.py 所在資料夾
    - PyInstaller onedir/onefile：回傳 sys._MEIPASS（onefile 解壓資料夾）或 exe 所在資料夾（保險）
    """
    if getattr(sys, "frozen", False):
        # onefile 會有 _MEIPASS；onedir 也可能有
        base = getattr(sys, "_MEIPASS", None)
        if base and os.path.isdir(base):
            return base
        # 保險：回 exe 所在資料夾
        return os.path.dirname(os.path.abspath(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def _exe_dir() -> str:
    if getattr(sys, "frozen", False):
        return os.path.dirname(os.path.abspath(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def find_experiment_script(filename: str) -> str:
    """
    依序嘗試常見位置：
    1) base/experiments/filename          (建議 add-data 放這)
    2) base/filename
    3) exe_dir/experiments/filename      (你若把檔案放 exe 同層)
    4) exe_dir/filename
    5) cwd/experiments/filename
    6) cwd/filename
    """
    base = _runtime_base_dir()
    exedir = _exe_dir()
    cwd = os.getcwd()

    candidates = [
        os.path.join(base, "experiments", filename),
        os.path.join(base, filename),
        os.path.join(exedir, "experiments", filename),
        os.path.join(exedir, filename),
        os.path.join(cwd, "experiments", filename),
        os.path.join(cwd, filename),
    ]

    for p in candidates:
        if os.path.exists(p):
            return p

    # 找不到 → 回傳詳細錯誤
    msg = "找不到實驗檔案：{}\n\n已嘗試路徑：\n- {}".format(
        filename, "\n- ".join(candidates)
    )
    raise FileNotFoundError(msg)


# ============================================================
#  Tk -> Embedded Frame monkey patch（核心）
# ============================================================
_EMBED_PARENT = None
_ORIG_TK = None


class EmbeddedTk(ttk.Frame):
    """取代 tk.Tk()：實際上是嵌入到 Notebook tab 的 Frame。"""

    def __init__(self, *args, **kwargs):
        global _EMBED_PARENT
        if _EMBED_PARENT is None:
            raise RuntimeError("EmbeddedTk: _EMBED_PARENT is None (no tab parent).")

        super().__init__(_EMBED_PARENT)

        try:
            _EMBED_PARENT.rowconfigure(0, weight=1)
            _EMBED_PARENT.columnconfigure(0, weight=1)
        except Exception:
            pass

        self.grid(row=0, column=0, sticky="nsew")
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

    # tk.Tk 常用視窗方法：no-op
    def title(self, *a, **k): return None
    def geometry(self, *a, **k): return None
    def resizable(self, *a, **k): return None
    def iconbitmap(self, *a, **k): return None
    def iconphoto(self, *a, **k): return None
    def attributes(self, *a, **k): return None
    def wm_attributes(self, *a, **k): return None
    def protocol(self, *a, **k): return None
    def minsize(self, *a, **k): return None
    def maxsize(self, *a, **k): return None
    def state(self, *a, **k): return None
    def wm_state(self, *a, **k): return None
    def lift(self, *a, **k): return None
    def focus_force(self, *a, **k): return None
    def withdraw(self, *a, **k): return None
    def deiconify(self, *a, **k): return None

    def mainloop(self, *a, **k):
        return None

    def quit(self, *a, **k):
        return None

    def destroy(self):
        try:
            super().destroy()
        except Exception:
            pass


class TkPatch:
    """context manager：暫時把 tk.Tk 換成 EmbeddedTk"""

    def __init__(self, parent):
        self.parent = parent

    def __enter__(self):
        global _EMBED_PARENT, _ORIG_TK
        _EMBED_PARENT = self.parent
        _ORIG_TK = tk.Tk
        tk.Tk = EmbeddedTk
        return self

    def __exit__(self, exc_type, exc, tb):
        global _EMBED_PARENT, _ORIG_TK
        try:
            tk.Tk = _ORIG_TK
        except Exception:
            pass
        _EMBED_PARENT = None
        _ORIG_TK = None
        return False


def embed_script_into_tab(tab: ttk.Frame, script_filename: str):
    """把原始 .py 執行並嵌入到 tab（不修改原始檔）"""
    script_path = find_experiment_script(script_filename)
    with TkPatch(tab):
        runpy.run_path(script_path, run_name="__main__")


# ============================================================
#  2_MOS：載入後修正灰色空白（不改 2_MOS.py 原檔）
# ============================================================
def _walk_widgets(root_widget):
    stack = [root_widget]
    while stack:
        w = stack.pop()
        yield w
        try:
            stack.extend(list(w.winfo_children()))
        except Exception:
            pass


def post_patch_mos(tab: ttk.Frame):
    for w in _walk_widgets(tab):
        try:
            if w.__class__.__name__ == "ParamRow":
                try:
                    w.columnconfigure(0, weight=0)
                    w.columnconfigure(1, weight=1)
                except Exception:
                    pass

                for ch in w.winfo_children():
                    if isinstance(ch, ttk.Entry):
                        try:
                            ch.grid_configure(sticky="ew")
                        except Exception:
                            pass
                        try:
                            ch.configure(width=10)
                        except Exception:
                            pass
                    if isinstance(ch, ttk.Scale):
                        try:
                            ch.grid_configure(sticky="ew")
                        except Exception:
                            pass
        except Exception:
            pass


# ============================================================
#  Lazy tab
# ============================================================
class LazyLoadTab(ttk.Frame):
    def __init__(self, master, script_filename, title, post_hook=None):
        super().__init__(master)
        self.script_filename = script_filename
        self.title = title
        self.post_hook = post_hook
        self.loaded = False

        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        self.hint = ttk.Label(self, text=f"{title}\n\n(點選此分頁後載入…)", justify="center")
        self.hint.grid(row=0, column=0, sticky="nsew", padx=20, pady=20)

    def load_once(self):
        if self.loaded:
            return
        self.loaded = True
        try:
            self.hint.destroy()
        except Exception:
            pass

        try:
            embed_script_into_tab(self, self.script_filename)
            if callable(self.post_hook):
                self.after(10, lambda: self.post_hook(self))
        except Exception as e:
            err = ttk.Label(
                self,
                text=f"載入失敗：{self.title}\n\n{type(e).__name__}: {e}",
                justify="left"
            )
            err.grid(row=0, column=0, sticky="nsew", padx=16, pady=16)


# ============================================================
#  主程式：6 tabs（1→6）
# ============================================================
def main():
    root = tk.Tk()
    root.title("光電子學作業整合 GUI（6 Experiments / Notebook）")
    root.geometry("1350x860")
    root.minsize(1200, 720)

    try:
        style = ttk.Style()
        if "vista" in style.theme_names():
            style.theme_use("vista")
    except Exception:
        pass

    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1)

    nb = ttk.Notebook(root)
    nb.grid(row=0, column=0, sticky="nsew")

    scripts = [
        ("1 Fiber",     "1_fiber.py", None),
        ("2 MOS",       "2_MOS.py", post_patch_mos),
        ("3 LED/MQW",   "3_MQW.py", None),
        ("4 Laser/MQW", "4_Laser.py", None),
        ("5 PIN/APD",   "5_PIN_APD.py", None),
        ("6 SolarCell", "6_SolarCell.py", None),
    ]

    tabs = []
    for title, fname, hook in scripts:
        t = LazyLoadTab(nb, fname, title, post_hook=hook)
        nb.add(t, text=title)
        tabs.append(t)

    def on_tab_changed(_evt=None):
        try:
            sel = nb.index(nb.select())
            tabs[sel].load_once()
        except Exception:
            pass

    nb.bind("<<NotebookTabChanged>>", on_tab_changed)
    root.after(50, on_tab_changed)
    root.mainloop()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        messagebox.showerror("Launcher Error", f"{type(e).__name__}: {e}")
