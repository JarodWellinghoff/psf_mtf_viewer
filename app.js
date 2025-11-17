// Main Application Logic

// Current parameters
let currentParams = {
    size: 513,
    wire_fov_mm: 50.0,
    background: 0.0,
    mode: 'ring',
    normalize: 'unit_energy',
    // Gaussian parameters
    p_i: 0.2,
    p_j: 0.2,
    // Ring parameters
    r0_mm: 5.0,
    sigma_mm: 1.0,
    core_weight: 0.5,
    core_sigma_mm: 0.5,
    // DoG parameters
    alpha: 0.1,
    sigma0_mm: 1.0,
    sigma1_mm: 2.0
};

// Initialize the application
function init() {
    setupEventListeners();
    updateParameterVisibility();
    updateVisualization();
}

// Setup event listeners for all input elements
function setupEventListeners() {
    // Basic parameters
    document.getElementById('size').addEventListener('change', updateFromInputs);
    document.getElementById('wire_fov_mm').addEventListener('change', updateFromInputs);
    document.getElementById('background').addEventListener('change', updateFromInputs);
    document.getElementById('mode').addEventListener('change', function() {
        updateFromInputs();
        updateParameterVisibility();
    });
    document.getElementById('normalize').addEventListener('change', updateFromInputs);

    // Gaussian parameters
    document.getElementById('p_i').addEventListener('change', updateFromInputs);
    document.getElementById('p_j').addEventListener('change', updateFromInputs);

    // Ring parameters
    document.getElementById('r0_mm').addEventListener('change', updateFromInputs);
    document.getElementById('sigma_mm').addEventListener('change', updateFromInputs);
    document.getElementById('core_weight').addEventListener('change', updateFromInputs);
    document.getElementById('core_sigma_mm').addEventListener('change', updateFromInputs);

    // DoG parameters
    document.getElementById('alpha').addEventListener('change', updateFromInputs);
    document.getElementById('sigma0_mm').addEventListener('change', updateFromInputs);
    document.getElementById('sigma1_mm').addEventListener('change', updateFromInputs);

    // Update button
    document.getElementById('update-btn').addEventListener('click', updateVisualization);

    // Allow Enter key to trigger update
    document.querySelectorAll('input, select').forEach(element => {
        element.addEventListener('keypress', function(e) {
            if (e.key === 'Enter') {
                updateVisualization();
            }
        });
    });
}

// Update parameters from input elements
function updateFromInputs() {
    currentParams.size = parseInt(document.getElementById('size').value);
    currentParams.wire_fov_mm = parseFloat(document.getElementById('wire_fov_mm').value);
    currentParams.background = parseFloat(document.getElementById('background').value);
    currentParams.mode = document.getElementById('mode').value;
    currentParams.normalize = document.getElementById('normalize').value;

    // Gaussian parameters
    currentParams.p_i = parseFloat(document.getElementById('p_i').value);
    currentParams.p_j = parseFloat(document.getElementById('p_j').value);

    // Ring parameters
    currentParams.r0_mm = parseFloat(document.getElementById('r0_mm').value);
    currentParams.sigma_mm = parseFloat(document.getElementById('sigma_mm').value);
    currentParams.core_weight = parseFloat(document.getElementById('core_weight').value);
    currentParams.core_sigma_mm = parseFloat(document.getElementById('core_sigma_mm').value);

    // DoG parameters
    currentParams.alpha = parseFloat(document.getElementById('alpha').value);
    currentParams.sigma0_mm = parseFloat(document.getElementById('sigma0_mm').value);
    currentParams.sigma1_mm = parseFloat(document.getElementById('sigma1_mm').value);
}

// Show/hide parameters based on selected mode
function updateParameterVisibility() {
    const mode = document.getElementById('mode').value;

    // Hide all mode-specific parameters
    document.querySelectorAll('.mode-params').forEach(el => {
        el.classList.remove('active');
    });

    // Show parameters for selected mode
    const modeParamsId = mode + '-params';
    const modeParams = document.getElementById(modeParamsId);
    if (modeParams) {
        modeParams.classList.add('active');
    }
}

// Main visualization update function
function updateVisualization() {
    try {
        updateFromInputs();

        // Show loading indicator (optional)
        console.log('Computing PSF/MTF...');

        // Compute PSF
        const { psf, dx_mm, size } = computePsf(currentParams);

        // Compute MTF
        const { fx, fy, f_rad, mtf_rad, mtf2d } = computeMtf(psf, dx_mm, size);

        // Plot results
        plotPsf(psf, currentParams.wire_fov_mm);
        plotMtf2d(mtf2d, fx, fy);
        plotMtf1d(f_rad, mtf_rad, dx_mm);

        console.log('Visualization updated successfully');
    } catch (error) {
        console.error('Error updating visualization:', error);
        alert('Error computing PSF/MTF. Please check your parameters.');
    }
}

// Plot PSF heatmap
function plotPsf(psf, wire_fov_mm) {
    const size = psf.length;
    const extent = wire_fov_mm / 2;

    // Create coordinate arrays
    const x = [];
    const y = [];
    for (let i = 0; i < size; i++) {
        x.push(-extent + (i / (size - 1)) * wire_fov_mm);
        y.push(-extent + (i / (size - 1)) * wire_fov_mm);
    }

    const data = [{
        z: psf,
        x: x,
        y: y,
        type: 'heatmap',
        colorscale: 'Hot',
        colorbar: {
            title: 'Intensity',
            titleside: 'right'
        }
    }];

    const layout = {
        title: '',
        xaxis: {
            title: 'x (mm)',
            zeroline: false
        },
        yaxis: {
            title: 'y (mm)',
            zeroline: false,
            scaleanchor: 'x'
        },
        margin: { t: 30, b: 50, l: 60, r: 80 }
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        displaylogo: false
    };

    Plotly.newPlot('psf-plot', data, layout, config);
}

// Plot 2D MTF
function plotMtf2d(mtf2d, fx, fy) {
    // Convert frequencies to cycles/cm
    const fx_cm = fx.map(f => f * 10);
    const fy_cm = fy.map(f => f * 10);

    const data = [{
        z: mtf2d,
        x: fx_cm,
        y: fy_cm,
        type: 'heatmap',
        colorscale: 'Viridis',
        colorbar: {
            title: 'MTF',
            titleside: 'right'
        }
    }];

    const layout = {
        title: '',
        xaxis: {
            title: 'Spatial Frequency x (cycles/cm)',
            zeroline: false
        },
        yaxis: {
            title: 'Spatial Frequency y (cycles/cm)',
            zeroline: false,
            scaleanchor: 'x'
        },
        margin: { t: 30, b: 50, l: 60, r: 80 }
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        displaylogo: false
    };

    Plotly.newPlot('mtf-2d-plot', data, layout, config);
}

// Plot 1D MTF
function plotMtf1d(f_rad, mtf_rad, dx_mm) {
    // Convert to cycles/cm
    const f_rad_cm = f_rad.map(f => f * 10);

    // MTF trace
    const trace1 = {
        x: f_rad_cm,
        y: mtf_rad,
        mode: 'lines',
        name: 'MTF',
        line: {
            color: 'rgb(31, 119, 180)',
            width: 2
        }
    };

    // MTF50 reference line
    const trace2 = {
        x: [0, Math.max(...f_rad_cm)],
        y: [0.5, 0.5],
        mode: 'lines',
        name: 'MTF50',
        line: {
            color: 'rgba(255, 0, 0, 0.5)',
            width: 1,
            dash: 'dash'
        }
    };

    const data = [trace1, trace2];

    const layout = {
        title: '',
        xaxis: {
            title: 'Spatial Frequency (cycles/cm)',
            range: [0, 10],
            zeroline: false
        },
        yaxis: {
            title: 'MTF',
            range: [0, 1.05],
            zeroline: false
        },
        showlegend: true,
        legend: {
            x: 0.7,
            y: 0.95
        },
        margin: { t: 30, b: 50, l: 60, r: 30 },
        plot_bgcolor: 'rgba(240, 240, 240, 0.5)',
        paper_bgcolor: 'white'
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        displaylogo: false
    };

    Plotly.newPlot('mtf-1d-plot', data, layout, config);
}

// Initialize when DOM is ready
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
} else {
    init();
}
