#!/usr/bin/env python3
"""
Modern Interactive PSF/MTF Viewer using Plotly Dash
A sleek, fast, and interactive replacement for the matplotlib version
"""

import numpy as np
from scipy.interpolate import interp1d
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc


# ===== PSF/MTF Calculation Functions (same as original) =====

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


def compute_psf(size, wire_fov_mm, mode, background, normalize_mode,
                p_i, p_j, r0_mm, sigma_mm, core_weight, core_sigma_mm,
                alpha, sigma0_mm, sigma1_mm, t_i=None, t_j=None):
    """Compute PSF based on selected mode and parameters"""
    if size % 2 == 0:
        size += 1
    H = W = size

    t_i = (H - 1) / 2.0 if t_i is None else float(t_i)
    t_j = (W - 1) / 2.0 if t_j is None else float(t_j)
    dx_mm = wire_fov_mm / size

    i = np.arange(H)[:, None]
    j = np.arange(W)[None, :]
    x_mm = (j - t_j) * dx_mm
    y_mm = (i - t_i) * dx_mm

    if mode == "dog":
        psf = psf_dog_unsharp(x_mm, y_mm, alpha, sigma0_mm, sigma1_mm)
    elif mode == "ring":
        ring = psf_gaussian_ring(x_mm, y_mm, r0_mm, sigma_mm)
        bump = psf_gaussian_isotropic(x_mm, y_mm, core_sigma_mm)
        psf = (1.0 - float(core_weight)) * ring + float(core_weight) * bump
    else:  # "gaussian"
        core = np.exp(-p_i * (i - t_i) ** 2) * np.exp(-p_j * (j - t_j) ** 2)
        psf = core / (p_i * p_j)

    psf = psf + float(background)
    psf = normalize_psf(psf, background=float(background), mode=normalize_mode)
    return psf, dx_mm, size


def compute_mtf(psf, dx, size):
    """Compute MTF from PSF"""
    half = size // 2

    # frequency grids in cycles/mm
    fx = np.fft.fftshift(np.fft.fftfreq(size, d=dx))
    fy = np.fft.fftshift(np.fft.fftfreq(size, d=dx))
    Fx, Fy = np.meshgrid(fx, fy)
    F = np.sqrt(Fx**2 + Fy**2)

    # unit-energy PSF
    psf_n = psf / psf.sum()

    # OTF/MTF, normalized so DC == 1
    OT = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(psf_n)))
    mtf2d = np.abs(OT)
    dc = mtf2d[half, half]
    if dc > 0:
        mtf2d = mtf2d / dc

    # radial slice from center to Nyquist
    f_rad = F[half, half:]
    mtf_rad = mtf2d[half, half:]
    return fx, fy, f_rad, mtf_rad, mtf2d


# ===== Dash App Setup =====

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "PSF/MTF Interactive Viewer"

# ===== Layout =====

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("PSF/MTF Interactive Viewer", className="text-center mb-4 mt-3"),
            html.P("Modern, interactive visualization of Point Spread Function and Modulation Transfer Function",
                   className="text-center text-muted mb-4")
        ])
    ]),

    # Control Panel
    dbc.Card([
        dbc.CardHeader(html.H4("Parameters")),
        dbc.CardBody([
            dbc.Row([
                # Column 1: Basic Settings
                dbc.Col([
                    html.H5("Basic Settings", className="mb-3"),
                    dbc.Label("Grid Size:"),
                    dcc.Slider(id='size-slider', min=129, max=1025, step=128, value=513,
                               marks={129: '129', 257: '257', 513: '513', 769: '769', 1025: '1025'},
                               tooltip={"placement": "bottom", "always_visible": True}),
                    html.Br(),

                    dbc.Label("Field of View (mm):"),
                    dcc.Slider(id='fov-slider', min=10, max=100, step=5, value=50,
                               marks={10: '10', 30: '30', 50: '50', 70: '70', 100: '100'},
                               tooltip={"placement": "bottom", "always_visible": True}),
                    html.Br(),

                    dbc.Label("PSF Mode:"),
                    dcc.Dropdown(
                        id='mode-dropdown',
                        options=[
                            {'label': 'Gaussian', 'value': 'gaussian'},
                            {'label': 'Ring (Airy-like)', 'value': 'ring'},
                            {'label': 'DoG (Unsharp Mask)', 'value': 'dog'}
                        ],
                        value='dog',
                        clearable=False
                    ),
                    html.Br(),

                    dbc.Label("Normalization:"),
                    dcc.Dropdown(
                        id='normalize-dropdown',
                        options=[
                            {'label': 'Unit Energy', 'value': 'unit_energy'},
                            {'label': 'Peak = 1', 'value': 'peak_one'},
                            {'label': 'None', 'value': 'none'}
                        ],
                        value='unit_energy',
                        clearable=False
                    ),
                    html.Br(),

                    dbc.Label("Background:"),
                    dcc.Input(id='background-input', type='number', value=0.0,
                             step=0.000001, className="form-control")
                ], md=3),

                # Column 2: Gaussian Parameters
                dbc.Col([
                    html.Div(id='gaussian-params', children=[
                        html.H5("Gaussian Parameters", className="mb-3"),
                        dbc.Label("p_i (pixel⁻¹):"),
                        dcc.Input(id='pi-input', type='number', value=0.4522,
                                 step=0.0001, className="form-control"),
                        html.Br(),
                        dbc.Label("p_j (pixel⁻¹):"),
                        dcc.Input(id='pj-input', type='number', value=0.4599,
                                 step=0.0001, className="form-control"),
                    ])
                ], md=3),

                # Column 3: Ring Parameters
                dbc.Col([
                    html.Div(id='ring-params', children=[
                        html.H5("Ring Parameters", className="mb-3"),
                        dbc.Label("Ring Radius (mm):"),
                        dcc.Slider(id='r0-slider', min=0.5, max=10, step=0.1, value=2.5,
                                   marks={0.5: '0.5', 3: '3', 5: '5', 8: '8', 10: '10'},
                                   tooltip={"placement": "bottom", "always_visible": True}),
                        html.Br(),

                        dbc.Label("Ring Width (mm):"),
                        dcc.Slider(id='sigma-slider', min=0.1, max=3, step=0.1, value=0.6,
                                   marks={0.1: '0.1', 1: '1', 2: '2', 3: '3'},
                                   tooltip={"placement": "bottom", "always_visible": True}),
                        html.Br(),

                        dbc.Label("Core Weight:"),
                        dcc.Slider(id='core-weight-slider', min=0, max=1, step=0.01, value=0.12,
                                   marks={0: '0', 0.25: '0.25', 0.5: '0.5', 0.75: '0.75', 1: '1'},
                                   tooltip={"placement": "bottom", "always_visible": True}),
                        html.Br(),

                        dbc.Label("Core Sigma (mm):"),
                        dcc.Slider(id='core-sigma-slider', min=0.1, max=1, step=0.05, value=0.25,
                                   marks={0.1: '0.1', 0.3: '0.3', 0.5: '0.5', 0.8: '0.8', 1: '1'},
                                   tooltip={"placement": "bottom", "always_visible": True}),
                    ])
                ], md=3),

                # Column 4: DoG Parameters
                dbc.Col([
                    html.Div(id='dog-params', children=[
                        html.H5("DoG Parameters", className="mb-3"),
                        dbc.Label("Alpha (sharpening):"),
                        dcc.Slider(id='alpha-slider', min=0, max=1, step=0.01, value=0.5,
                                   marks={0: '0', 0.25: '0.25', 0.5: '0.5', 0.75: '0.75', 1: '1'},
                                   tooltip={"placement": "bottom", "always_visible": True}),
                        html.Br(),

                        dbc.Label("Sigma0 (mm):"),
                        dcc.Slider(id='sigma0-slider', min=0.1, max=2, step=0.1, value=0.6,
                                   marks={0.1: '0.1', 0.6: '0.6', 1: '1', 1.5: '1.5', 2: '2'},
                                   tooltip={"placement": "bottom", "always_visible": True}),
                        html.Br(),

                        dbc.Label("Sigma1 (mm):"),
                        dcc.Slider(id='sigma1-slider', min=0.2, max=4, step=0.1, value=1.2,
                                   marks={0.2: '0.2', 1: '1', 2: '2', 3: '3', 4: '4'},
                                   tooltip={"placement": "bottom", "always_visible": True}),
                    ])
                ], md=3),
            ])
        ])
    ], className="mb-4"),

    # Visualization Area
    dbc.Card([
        dbc.CardHeader(html.H4("Visualizations")),
        dbc.CardBody([
            dcc.Loading(
                id="loading",
                type="default",
                children=[
                    dcc.Graph(id='psf-mtf-plot', style={'height': '600px'})
                ]
            )
        ])
    ]),

    html.Footer([
        html.Hr(),
        html.P("Interactive PSF/MTF Viewer | Powered by Plotly Dash",
               className="text-center text-muted small mb-3")
    ])
], fluid=True)


# ===== Callbacks =====

@app.callback(
    [Output('gaussian-params', 'style'),
     Output('ring-params', 'style'),
     Output('dog-params', 'style')],
    [Input('mode-dropdown', 'value')]
)
def update_param_visibility(mode):
    """Show/hide parameter panels based on selected mode"""
    gaussian_style = {'display': 'block' if mode == 'gaussian' else 'none'}
    ring_style = {'display': 'block' if mode == 'ring' else 'none'}
    dog_style = {'display': 'block' if mode == 'dog' else 'none'}
    return gaussian_style, ring_style, dog_style


@app.callback(
    Output('psf-mtf-plot', 'figure'),
    [Input('size-slider', 'value'),
     Input('fov-slider', 'value'),
     Input('mode-dropdown', 'value'),
     Input('normalize-dropdown', 'value'),
     Input('background-input', 'value'),
     Input('pi-input', 'value'),
     Input('pj-input', 'value'),
     Input('r0-slider', 'value'),
     Input('sigma-slider', 'value'),
     Input('core-weight-slider', 'value'),
     Input('core-sigma-slider', 'value'),
     Input('alpha-slider', 'value'),
     Input('sigma0-slider', 'value'),
     Input('sigma1-slider', 'value')]
)
def update_plots(size, fov, mode, normalize_mode, background,
                p_i, p_j, r0_mm, sigma_mm, core_weight, core_sigma_mm,
                alpha, sigma0_mm, sigma1_mm):
    """Update all three plots when parameters change"""

    # Compute PSF and MTF
    psf, dx, actual_size = compute_psf(
        size, fov, mode, background, normalize_mode,
        p_i, p_j, r0_mm, sigma_mm, core_weight, core_sigma_mm,
        alpha, sigma0_mm, sigma1_mm
    )
    fx, fy, f_rad, mtf_rad, mtf2d = compute_mtf(psf, dx, actual_size)

    # Create subplot figure
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('PSF (Point Spread Function)',
                       'MTF 2D (Frequency Domain)',
                       'MTF 1D (Radial Profile)'),
        specs=[[{'type': 'heatmap'}, {'type': 'heatmap'}, {'type': 'scatter'}]],
        horizontal_spacing=0.1
    )

    # PSF Heatmap
    extent_mm = [-fov/2, fov/2]
    fig.add_trace(
        go.Heatmap(
            z=psf,
            x=np.linspace(extent_mm[0], extent_mm[1], actual_size),
            y=np.linspace(extent_mm[0], extent_mm[1], actual_size),
            colorscale='Hot',
            name='PSF',
            hovertemplate='x: %{x:.2f} mm<br>y: %{y:.2f} mm<br>Intensity: %{z:.4f}<extra></extra>',
            colorbar=dict(x=0.29, len=0.8)
        ),
        row=1, col=1
    )

    # 2D MTF Heatmap
    extent_freq = [fx.min() * 10, fx.max() * 10, fy.min() * 10, fy.max() * 10]
    fig.add_trace(
        go.Heatmap(
            z=mtf2d,
            x=np.linspace(extent_freq[0], extent_freq[1], actual_size),
            y=np.linspace(extent_freq[2], extent_freq[3], actual_size),
            colorscale='Viridis',
            name='MTF 2D',
            hovertemplate='fx: %{x:.2f} cyc/cm<br>fy: %{y:.2f} cyc/cm<br>MTF: %{z:.4f}<extra></extra>',
            colorbar=dict(x=0.63, len=0.8)
        ),
        row=1, col=2
    )

    # 1D MTF Line Plot
    f_rad_cm = f_rad * 10  # Convert to cycles/cm
    fig.add_trace(
        go.Scatter(
            x=f_rad_cm,
            y=mtf_rad,
            mode='lines',
            name='MTF',
            line=dict(color='royalblue', width=3),
            hovertemplate='Frequency: %{x:.2f} cyc/cm<br>MTF: %{y:.4f}<extra></extra>'
        ),
        row=1, col=3
    )

    # Add MTF50 reference line
    fig.add_trace(
        go.Scatter(
            x=[0, f_rad_cm[-1]],
            y=[0.5, 0.5],
            mode='lines',
            name='MTF50',
            line=dict(color='red', width=2, dash='dash'),
            hovertemplate='MTF50 Reference<extra></extra>'
        ),
        row=1, col=3
    )

    # Update axes
    fig.update_xaxes(title_text="x (mm)", row=1, col=1)
    fig.update_yaxes(title_text="y (mm)", row=1, col=1)

    fig.update_xaxes(title_text="Spatial Frequency x (cycles/cm)", row=1, col=2)
    fig.update_yaxes(title_text="Spatial Frequency y (cycles/cm)", row=1, col=2)

    fig.update_xaxes(title_text="Spatial Frequency (cycles/cm)", range=[0, 10], row=1, col=3)
    fig.update_yaxes(title_text="MTF", range=[0, 1.05], row=1, col=3)

    # Add grid to 1D MTF plot
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray', row=1, col=3)
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray', row=1, col=3)

    # Update layout
    fig.update_layout(
        height=600,
        showlegend=True,
        legend=dict(x=0.88, y=0.5),
        hovermode='closest',
        template='plotly_white'
    )

    return fig


# ===== Run App =====

if __name__ == '__main__':
    print("\n" + "="*60)
    print("PSF/MTF Interactive Viewer - Starting...")
    print("="*60)
    print("\nOpen your browser and navigate to: http://127.0.0.1:8050")
    print("\nPress Ctrl+C to stop the server\n")
    app.run_server(debug=True, host='127.0.0.1', port=8050)
