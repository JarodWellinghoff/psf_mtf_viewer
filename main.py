import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from scipy.interpolate import interp1d


def psf_gaussian_core(i, j, p_i, p_j, t_i, t_j):
    core = np.exp(-p_i * (i - t_i) ** 2) * np.exp(-p_j * (j - t_j) ** 2)
    return core / (p_i * p_j)


def psf_gaussian_ring(x_mm, y_mm, r0_mm, sigma_mm):
    R = np.sqrt(x_mm**2 + y_mm**2)
    return np.exp(-0.5 * ((R - r0_mm) / max(1e-12, sigma_mm)) ** 2)


def psf_gaussian_isotropic(x_mm, y_mm, sigma_mm):
    R2 = x_mm**2 + y_mm**2
    return np.exp(-0.5 * (R2 / max(1e-12, sigma_mm**2)))


def _gauss_iso(x_mm, y_mm, sigma_mm):
    R2 = x_mm**2 + y_mm**2
    s2 = max(1e-12, sigma_mm**2)
    return np.exp(-0.5 * R2 / s2) / (2 * np.pi * s2)


def psf_dog_unsharp(x_mm, y_mm, alpha, sigma0_mm, sigma1_mm):
    """(1+alpha)G(sigma0) - alpha*G(sigma1), sums to 1 by construction."""
    return (1.0 + alpha) * _gauss_iso(x_mm, y_mm, sigma0_mm) - alpha * _gauss_iso(
        x_mm, y_mm, sigma1_mm
    )


def normalize_psf(psf, background=0.0, mode="unit_energy"):
    if mode in ("unit_energy", "unit", "energy"):
        sig = psf - background
        s = sig.sum()
        return (sig / s + background) if s > 0 else psf
    if mode in ("peak_one", "peak", "max1", "max"):
        m = psf.max()
        return psf / m if m > 0 else psf
    return psf


def simulate_psf_gaussian(
    size_xy, p_i, p_j, t_i, t_j, background=0.0, normalize="unit_energy"
):
    """
    Gaussian PSF with optional background and normalization.

    h(i,j) = (1/(p_i*p_j)) * exp(-p_i*(i-t_i)^2) * exp(-p_j*(j-t_j)^2) + B
      - i = row index (0..H-1), j = col index (0..W-1)
      - p_i, p_j in pixel^-1
      - t_i, t_j can be subpixel (float)
      - normalize in {"unit_energy","peak_one","none"}
    """
    H, W = size_xy
    i = np.arange(H)[:, None]
    j = np.arange(W)[None, :]

    core = np.exp(-p_i * (i - t_i) ** 2) * np.exp(-p_j * (j - t_j) ** 2)
    psf = core / (p_i * p_j) + background

    if normalize.lower() in ("unit_energy", "unit", "energy"):
        s = (psf - background).sum()
        if s > 0:
            psf = (psf - background) / s + background
    elif normalize.lower() in ("peak_one", "peak", "max1", "max"):
        m = psf.max()
        if m > 0:
            psf = psf / m
    # else: "none" -> leave scale as-is

    return psf


class PSFMTFViewer:
    def __init__(self):
        # Grid / display params
        self.size = 513  # force odd
        self.wire_fov_mm = 50.0  # only used for axis extents

        # PSF (paper) params & defaults
        self.p_i = 0.4522  # pixel^-1 (row direction)
        self.p_j = 0.4599  # pixel^-1 (col direction)
        self.t_i = None  # set to center on first update
        self.t_j = None
        self.background = 0.0
        self.normalize = "unit_energy"  # "unit_energy" | "peak_one" | "none"
        self.r0_mm = 2.5  # ring radius (mm)
        self.sigma_mm = 0.6  # ring width (mm)
        self.core_weight = 0.12  # 0..1 (e.g., 0.05–0.2)
        self.core_sigma_mm = 0.25  # narrow bump; try 0.15–0.4 mm
        self.mode = "dog"  # "gaussian", "ring", or "dog"
        self.alpha = 0.5  # sharpening amount (0.05–0.2 typical)
        self.sigma0_mm = 0.6  # narrow core width
        self.sigma1_mm = 1.2  # broad halo width (>= 3*sigma0)

        # Create figure and axes
        self.fig = plt.figure(figsize=(14, 6))
        self.ax_psf = plt.subplot(1, 3, 1)
        self.ax_mtf2d = plt.subplot(1, 3, 2)
        self.ax_mtf = plt.subplot(1, 3, 3)

        # Make room for text boxes
        plt.subplots_adjust(bottom=0.42)

        # Create text boxes
        self.create_textboxes()

        # Initial plot
        self.update_plots()
        plt.show()

    def create_textboxes(self):

        col_w = 0.18  # width of each textbox
        row_h = 0.05  # height of each textbox
        col_x = [0.08, 0.30, 0.52, 0.74]  # left positions for columns 1-4
        start_y = 0.32  # top row y-position
        row_spacing = 0.06  # vertical spacing (distance between rows)

        # Helper to compute textbox positions
        def tb_axis(col, row):
            return plt.axes([col_x[col], start_y - row * row_spacing, col_w, row_h])

        # === Column 1: Basic Setup ===
        ax_size = tb_axis(0, 0)
        ax_fov = tb_axis(0, 1)
        ax_mode = tb_axis(0, 2)  # “gaussian / ring / dog”
        ax_bg = tb_axis(0, 3)
        ax_norm = tb_axis(0, 4)
        self.tb_size = TextBox(ax_size, "size:", initial=str(self.size))
        self.tb_fov = TextBox(ax_fov, "wire_fov_mm:", initial=str(self.wire_fov_mm))
        self.tb_bg = TextBox(ax_bg, "background B:", initial=f"{self.background:.6f}")
        self.tb_mode = TextBox(
            ax_mode, "mode:", initial=self.mode
        )  # "ring" or "gaussian"
        self.tb_norm = TextBox(ax_norm, "normalize:", initial=self.normalize)

        # === Column 2: Mode Specific Settings ===
        # Gaussian
        # if self.mode == "gaussian":
        #     ax_pi = tb_axis(1, 0)
        #     ax_pj = tb_axis(1, 1)
        # self.tb_pi = TextBox(ax_pi, "p_i (pix^-1):", initial=f"{self.p_i:.6f}")
        # self.tb_pj = TextBox(ax_pj, "p_j (pix^-1):", initial=f"{self.p_j:.6f}")

        # DoG
        # elif self.mode == "dog":
        # ax_alpha = tb_axis(1, 0)
        # ax_s0 = tb_axis(1, 1)
        # ax_s1 = tb_axis(1, 2)
        # self.tb_alpha = TextBox(ax_alpha, "alpha:", initial=f"{self.alpha:.3f}")
        # self.tb_s0 = TextBox(ax_s0, "sigma0_mm:", initial=f"{self.sigma0_mm:.3f}")
        # self.tb_s1 = TextBox(ax_s1, "sigma1_mm:", initial=f"{self.sigma1_mm:.3f}")

        # Ring ===
        # elif self.mode == "ring":
        ax_r0 = tb_axis(1, 0)
        ax_sigma = tb_axis(1, 1)
        ax_cw = tb_axis(1, 2)
        ax_cs = tb_axis(1, 3)
        self.tb_r0 = TextBox(ax_r0, "r0_mm:", initial=f"{self.r0_mm:.3f}")
        self.tb_sigma = TextBox(ax_sigma, "sigma_mm:", initial=f"{self.sigma_mm:.3f}")
        self.tb_cw = TextBox(ax_cw, "core_weight:", initial=f"{self.core_weight:.3f}")
        self.tb_cs = TextBox(
            ax_cs, "core_sigma_mm:", initial=f"{self.core_sigma_mm:.3f}"
        )

        # t_i/t_j default to center; show "center" initially
        # self.tb_ti = TextBox(ax_ti, "t_i (row, px):", initial="center")
        # self.tb_tj = TextBox(ax_tj, "t_j (col, px):", initial="center")

        # Bind
        self.tb_size.on_submit(self.on_size_change)
        self.tb_fov.on_submit(self.on_fov_change)
        # self.tb_pi.on_submit(self.on_pi_change)
        # self.tb_pj.on_submit(self.on_pj_change)
        # self.tb_ti.on_submit(self.on_ti_change)
        # self.tb_tj.on_submit(self.on_tj_change)
        self.tb_bg.on_submit(self.on_bg_change)
        self.tb_norm.on_submit(self.on_norm_change)
        self.tb_mode.on_submit(self.on_mode_change)
        self.tb_r0.on_submit(self.on_r0_change)
        self.tb_sigma.on_submit(self.on_sigma_change)
        self.tb_cw.on_submit(self.on_cw_change)
        self.tb_cs.on_submit(self.on_cs_change)
        self.tb_mode.on_submit(
            lambda t: (setattr(self, "mode", t.strip().lower()), self.update_plots())
        )
        # self.tb_alpha.on_submit(
        #     lambda t: (setattr(self, "alpha", float(t)), self.update_plots())
        # )
        # self.tb_s0.on_submit(
        #     lambda t: (
        #         setattr(self, "sigma0_mm", max(1e-6, float(t))),
        #         self.update_plots(),
        #     )
        # )
        # self.tb_s1.on_submit(
        #     lambda t: (
        #         setattr(self, "sigma1_mm", max(1e-6, float(t))),
        #         self.update_plots(),
        #     )
        # )

    # Handlers
    def on_size_change(self, text):
        try:
            self.size = int(text)
            self.update_plots()
        except ValueError:
            pass

    def on_fov_change(self, text):
        try:
            self.wire_fov_mm = float(text)
            self.update_plots()
        except ValueError:
            pass

    def on_pi_change(self, text):
        try:
            self.p_i = float(text)
            self.update_plots()
        except ValueError:
            pass

    def on_pj_change(self, text):
        try:
            self.p_j = float(text)
            self.update_plots()
        except ValueError:
            pass

    def _parse_center_or_float(self, text, fallback):
        if isinstance(text, str) and text.strip().lower() in (
            "center",
            "centre",
            "mid",
            "middle",
            "",
        ):
            return None  # means "center"
        try:
            return float(text)
        except Exception:
            return fallback

    def on_ti_change(self, text):
        self.t_i = self._parse_center_or_float(text, self.t_i)
        self.update_plots()

    def on_tj_change(self, text):
        self.t_j = self._parse_center_or_float(text, self.t_j)
        self.update_plots()

    def on_bg_change(self, text):
        try:
            self.background = float(text)
            self.update_plots()
        except ValueError:
            pass

    def on_norm_change(self, text):
        val = (text or "").strip().lower()
        if val in (
            "unit_energy",
            "peak_one",
            "none",
            "unit",
            "energy",
            "peak",
            "max1",
            "max",
        ):
            self.normalize = val
            self.update_plots()

    def on_mode_change(self, text):
        val = (text or "").strip().lower()
        if val in ("ring", "gaussian", "dog"):
            self.mode = val
            self.update_plots()

    def on_r0_change(self, text):
        try:
            self.r0_mm = float(text)
            self.update_plots()
        except ValueError:
            pass

    def on_sigma_change(self, text):
        try:
            self.sigma_mm = float(text)
            self.update_plots()
        except ValueError:
            pass

    def on_cw_change(self, text):
        try:
            v = float(text)
            self.core_weight = min(1.0, max(0.0, v))
            self.update_plots()
        except ValueError:
            pass

    def on_cs_change(self, text):
        try:
            self.core_sigma_mm = max(1e-6, float(text))
            self.update_plots()
        except ValueError:
            pass

    def compute_psf(self):
        size = int(self.size)
        if size % 2 == 0:
            size += 1
        H = W = size

        t_i = (H - 1) / 2.0 if self.t_i is None else float(self.t_i)
        t_j = (W - 1) / 2.0 if self.t_j is None else float(self.t_j)
        dx_mm = self.wire_fov_mm / size

        i = np.arange(H)[:, None]
        j = np.arange(W)[None, :]
        x_mm = (j - t_j) * dx_mm
        y_mm = (i - t_i) * dx_mm

        if self.mode == "dog":
            psf = psf_dog_unsharp(
                x_mm, y_mm, self.alpha, self.sigma0_mm, self.sigma1_mm
            )
        elif self.mode == "ring":
            ring = psf_gaussian_ring(x_mm, y_mm, self.r0_mm, self.sigma_mm)
            bump = psf_gaussian_isotropic(x_mm, y_mm, self.core_sigma_mm)
            psf = (1.0 - float(self.core_weight)) * ring + float(
                self.core_weight
            ) * bump
        else:  # "gaussian"
            core = np.exp(-self.p_i * (i - t_i) ** 2) * np.exp(
                -self.p_j * (j - t_j) ** 2
            )
            psf = core / (self.p_i * self.p_j)

        # add background ONCE, then normalize according to mode
        psf = psf + float(self.background)
        psf = normalize_psf(psf, background=float(self.background), mode=self.normalize)
        return psf, dx_mm, size

    def compute_mtf(self, psf, dx, size):
        half = size // 2

        # frequency grids in cycles/mm
        fx = np.fft.fftshift(np.fft.fftfreq(size, d=dx))
        fy = np.fft.fftshift(np.fft.fftfreq(size, d=dx))
        Fx, Fy = np.meshgrid(fx, fy)
        F = np.sqrt(Fx**2 + Fy**2)

        # 1) unit-energy PSF (no clipping, no min-subtraction)
        psf_n = psf / psf.sum()

        # 2) OTF/MTF, normalized so DC == 1
        OT = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(psf_n)))
        mtf2d = np.abs(OT)
        dc = mtf2d[half, half]
        if dc > 0:
            mtf2d = mtf2d / dc

        # radial slice from center to Nyquist (x-axis)
        f_rad = F[half, half:]
        mtf_rad = mtf2d[half, half:]
        return fx, fy, f_rad, mtf_rad, mtf2d

    def update_plots(self):
        psf, dx, size = self.compute_psf()
        fx, fy, f_rad, mtf_rad, mtf2d = self.compute_mtf(psf, dx, size)

        # PSF image (mm)
        self.ax_psf.clear()
        self.ax_mtf2d.clear()
        self.ax_mtf.clear()
        extent_mm = [
            -self.wire_fov_mm / 2,
            self.wire_fov_mm / 2,
            -self.wire_fov_mm / 2,
            self.wire_fov_mm / 2,
        ]
        self.ax_psf.imshow(psf, extent=extent_mm, cmap="hot", origin="lower")
        self.ax_psf.set_xlabel("x (mm)")
        self.ax_psf.set_ylabel("y (mm)")
        self.ax_psf.set_title("PSF")

        # 2D MTF (cycles/cm) — use frequency extent
        extent_freq = [fx.min() * 10, fx.max() * 10, fy.min() * 10, fy.max() * 10]
        self.ax_mtf2d.imshow(
            mtf2d, extent=extent_freq, cmap="viridis", origin="lower", aspect="auto"
        )
        self.ax_mtf2d.set_xlabel("Spatial Frequency x (cycles/cm)")
        self.ax_mtf2d.set_ylabel("Spatial Frequency y (cycles/cm)")
        self.ax_mtf2d.set_title("MTF (2D)")

        # 1D MTF
        self.ax_mtf.plot(f_rad * 10, mtf_rad, linewidth=2)
        self.ax_mtf.axhline(y=0.5, linestyle="--", alpha=0.5, label="MTF50")
        self.ax_mtf.set_xlabel("Spatial Frequency (cycles/cm)")
        self.ax_mtf.set_ylabel("MTF")
        self.ax_mtf.set_title("Modulation Transfer Function")
        self.ax_mtf.grid(True, alpha=0.3)
        self.ax_mtf.set_ylim([0, 1.05])
        nyq_cm = (1.0 / (2.0 * dx)) * 10
        self.ax_mtf.set_xlim([0, 10])
        self.ax_mtf.legend()
        self.fig.canvas.draw_idle()


if __name__ == "__main__":
    viewer = PSFMTFViewer()
