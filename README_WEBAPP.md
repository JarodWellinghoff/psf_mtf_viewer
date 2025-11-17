# PSF/MTF Viewer - Web Application

An interactive web-based tool for analyzing Point Spread Functions (PSF) and Modulation Transfer Functions (MTF) of optical systems. This application runs entirely in the browser with no server required.

## Features

- **Three PSF Models**:
  - **Gaussian**: Simple 2D Gaussian distribution for basic optical blur analysis
  - **Ring**: Simulates diffraction patterns with ring structures (Airy disk patterns)
  - **DoG (Difference of Gaussians)**: Edge enhancement and unsharp masking filters

- **Real-time Visualization**:
  - Point Spread Function (PSF) heatmap
  - 2D Modulation Transfer Function (MTF)
  - 1D Radial MTF profile with MTF50 reference

- **Interactive Parameters**: Adjust all parameters in real-time and see immediate results

- **Flexible Normalization**: Unit energy, peak normalization, or no normalization

## Usage

### Running Locally

1. Clone this repository
2. Open `index.html` in a modern web browser (Chrome, Firefox, Safari, Edge)
3. No installation or build process required!

### GitHub Pages Deployment

This application is designed to be hosted on GitHub Pages:

1. Push the repository to GitHub
2. Go to repository Settings → Pages
3. Select the branch to deploy (e.g., `main` or `gh-pages`)
4. Set the source folder to `/` (root)
5. Save and wait a few minutes for deployment
6. Access your app at `https://<username>.github.io/<repository-name>/`

## Parameters Guide

### Basic Setup
- **Size**: Image size in pixels (must be odd for symmetry)
- **Wire FOV**: Field of view in millimeters
- **Background**: Background intensity level
- **PSF Mode**: Choose between Gaussian, Ring, or DoG
- **Normalize**: Normalization method for the PSF

### Gaussian Mode
- **p_i**: Precision parameter in row direction (pixel⁻¹)
- **p_j**: Precision parameter in column direction (pixel⁻¹)

### Ring Mode
- **r0**: Ring radius in millimeters
- **sigma**: Ring width in millimeters
- **Core Weight**: Blend factor for central bump (0-1)
- **Core Sigma**: Width of the central peak in millimeters

### DoG Mode
- **Alpha**: Sharpening strength (typical: 0.05-0.2)
- **Sigma 0**: Narrow core width in millimeters
- **Sigma 1**: Broad halo width in millimeters

## Technical Details

- **Frontend**: Pure HTML, CSS, and JavaScript
- **Visualization**: Plotly.js for interactive plots
- **Computation**: All calculations performed client-side in JavaScript
- **FFT**: Custom FFT implementation for MTF calculation

## Browser Compatibility

Tested and working on:
- Chrome/Edge (latest)
- Firefox (latest)
- Safari (latest)

## Original Python Version

This web application is a port of the original Python-based PSF/MTF viewer. The mathematical calculations are identical, ensuring consistent results between both versions.

## License

See LICENSE file for details.
