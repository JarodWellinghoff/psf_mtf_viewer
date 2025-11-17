# PSF/MTF Viewer

Interactive visualization tool for Point Spread Function (PSF) and Modulation Transfer Function (MTF) analysis of optical systems.

## 🚀 Quick Start

### Modern Interactive Version (Recommended)

The new **Plotly Dash** version provides a sleek, modern, and highly interactive experience:

```bash
# Install dependencies
pip install -r requirements.txt

# Run the interactive dashboard
python3 main_dash.py
```

Then open your browser to: **http://127.0.0.1:8050**

### Classic Matplotlib Version

For the original matplotlib implementation:

```bash
python3 main.py
```

## ✨ Features

### Modern Dash Version (`main_dash.py`) - ⭐ Recommended

- **Highly Interactive**: Zoom, pan, and hover for detailed information
- **Real-time Updates**: Parameters update instantly as you adjust sliders
- **Modern UI**: Clean, professional dashboard layout with Bootstrap styling
- **Fast Rendering**: Much faster than matplotlib for complex visualizations
- **Responsive Design**: Works great on different screen sizes
- **Better Performance**: Optimized for smooth interactions even with large grids

### Classic Version (`main.py`)

- Traditional matplotlib interface
- TextBox-based parameter input
- Static plot updates

## 📊 Visualization Components

The tool provides three synchronized views:

1. **PSF Heatmap** (Left)
   - Spatial intensity distribution
   - Shows how a point source spreads through the optical system
   - Units: millimeters (mm)

2. **2D MTF Heatmap** (Center)
   - Frequency domain representation
   - Shows modulation at different spatial frequencies
   - Units: cycles/cm

3. **1D MTF Profile** (Right)
   - Radial slice from center to edge
   - Includes MTF50 reference line (50% modulation point)
   - Critical metric for optical system performance

## 🎛️ PSF Models

### 1. Gaussian Mode
Simple 2D Gaussian blur model:
- **p_i**: Precision in row direction (pixel⁻¹)
- **p_j**: Precision in column direction (pixel⁻¹)
- Default: p_i=0.4522, p_j=0.4599

### 2. Ring Mode (Airy-like)
Diffraction ring pattern with optional central bump:
- **r0_mm**: Ring radius (mm)
- **sigma_mm**: Ring width (mm)
- **core_weight**: Central bump weight (0-1)
- **core_sigma_mm**: Central bump width (mm)

Models Airy disk-like patterns common in optical systems with circular apertures.

### 3. DoG (Difference of Gaussians)
Unsharp masking filter:
- **alpha**: Sharpening strength (0-1)
- **sigma0_mm**: Narrow core width (mm)
- **sigma1_mm**: Broad halo width (mm)
- Formula: (1+α)G(σ₀) - αG(σ₁)

Used to model edge enhancement and sharpening effects.

## 🔧 Parameters

### Basic Settings
- **Grid Size**: Resolution of the PSF computation (129-1025 pixels)
- **Field of View**: Physical size of the visualization area (mm)
- **Normalization**:
  - `unit_energy`: Total PSF energy = 1
  - `peak_one`: Peak PSF value = 1
  - `none`: No normalization
- **Background**: Constant background intensity

## 📐 Understanding MTF

**MTF (Modulation Transfer Function)** describes how well an optical system reproduces contrast at different spatial frequencies:

- **MTF = 1.0**: Perfect reproduction of that frequency
- **MTF = 0.5**: 50% contrast reduction (MTF50 - key metric)
- **MTF = 0.0**: Complete loss of that frequency

**Key Metric:**
- **MTF50**: The spatial frequency at which modulation drops to 50%
- Higher MTF50 = sharper optical system
- Displayed as dashed red line in the 1D plot

## 🎨 Why Use the Dash Version?

### Performance Comparison

| Feature | Matplotlib (`main.py`) | Dash (`main_dash.py`) |
|---------|----------------------|---------------------|
| **Interactive** | Limited | Highly interactive |
| **Zoom/Pan** | Basic | Smooth, native |
| **Hover Info** | No | Yes, detailed |
| **Update Speed** | Slow (~1-2s) | Fast (<0.5s) |
| **UI Design** | Basic textboxes | Modern dashboard |
| **Parameter Control** | Text input only | Sliders + tooltips |
| **Responsiveness** | Static | Responsive |
| **Web Export** | No | Yes |

### User Experience

**Matplotlib Version:**
- Enter values manually in text boxes
- Press Enter to update
- Slow redraws
- Clunky parameter adjustment

**Dash Version:**
- Drag sliders for instant feedback
- Hover over plots for precise values
- Zoom into regions of interest
- Pan around the visualization
- Professional, modern interface

## 📦 Dependencies

```
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.4.0          # For classic version
dash>=2.14.0               # For modern version
dash-bootstrap-components>=1.5.0
plotly>=5.17.0
```

## 🛠️ Installation

```bash
# Clone the repository
git clone <repository-url>
cd psf_mtf_viewer

# Install dependencies
pip install -r requirements.txt

# Run the modern version (recommended)
python3 main_dash.py

# Or run the classic version
python3 main.py
```

## 📖 Usage Examples

### Quick PSF Analysis

1. **Launch the Dash app**: `python3 main_dash.py`
2. **Select a PSF mode**: Choose from Gaussian, Ring, or DoG
3. **Adjust parameters**: Use sliders to see real-time updates
4. **Analyze MTF**: Look at the MTF50 value to assess sharpness
5. **Export**: Use Plotly's built-in tools to save plots

### Comparing Different Systems

To compare two optical configurations:
1. Note the MTF50 value for the first configuration
2. Adjust parameters to the second configuration
3. Compare MTF50 values - higher is sharper

### Finding Optimal Parameters

Use the interactive sliders to explore the parameter space:
- Adjust one parameter at a time
- Watch how the MTF changes in real-time
- Look for configurations that maximize MTF50

## 🌐 Web Version

There's also a pure JavaScript web application version (`index.html`) that runs entirely in the browser without Python. See `README_WEBAPP.md` for details.

## 📚 Technical Background

### PSF (Point Spread Function)
The PSF describes how a point source of light is imaged by an optical system. It represents the impulse response of the system in the spatial domain.

### MTF (Modulation Transfer Function)
The MTF is the magnitude of the Fourier transform of the PSF. It quantifies the contrast transfer at different spatial frequencies and is a standard metric for characterizing optical system performance.

### Computation Method
1. Generate PSF based on selected model
2. Normalize PSF to unit energy
3. Apply 2D FFT to obtain OTF (Optical Transfer Function)
4. Take magnitude to get MTF
5. Normalize so DC component = 1
6. Extract radial profile for 1D MTF curve

## 🤝 Contributing

Contributions are welcome! Feel free to:
- Report bugs
- Suggest new PSF models
- Improve the UI
- Add new features

## 📄 License

[Your License Here]

## 🔗 References

- Born & Wolf, "Principles of Optics"
- Goodman, "Introduction to Fourier Optics"
- ISO 9334: Optics and optical instruments - Optical transfer function

## 🙋 Support

For questions or issues:
- Check the code comments for implementation details
- Review the web app version for JavaScript implementation
- Open an issue on GitHub

---

**Built with:** Python, NumPy, SciPy, Plotly, Dash, and Bootstrap

**Made modern with:** Interactive Plotly Dash replacing slow matplotlib ⚡
