// PSF/MTF Calculation Functions
// Ported from Python implementation

// Helper function: Create a 2D array filled with a value
function create2DArray(rows, cols, initialValue = 0) {
    return Array.from({ length: rows }, () =>
        Array.from({ length: cols }, () => initialValue)
    );
}

// Helper function: Create a range array
function range(start, end) {
    return Array.from({ length: end - start }, (_, i) => start + i);
}

// Gaussian core function
function psfGaussianCore(i, j, p_i, p_j, t_i, t_j) {
    const core = Math.exp(-p_i * Math.pow(i - t_i, 2)) *
                 Math.exp(-p_j * Math.pow(j - t_j, 2));
    return core / (p_i * p_j);
}

// Ring-shaped Gaussian PSF
function psfGaussianRing(x_mm, y_mm, r0_mm, sigma_mm) {
    const R = Math.sqrt(x_mm * x_mm + y_mm * y_mm);
    const sigma = Math.max(1e-12, sigma_mm);
    return Math.exp(-0.5 * Math.pow((R - r0_mm) / sigma, 2));
}

// Isotropic Gaussian
function psfGaussianIsotropic(x_mm, y_mm, sigma_mm) {
    const R2 = x_mm * x_mm + y_mm * y_mm;
    const sigma2 = Math.max(1e-12, sigma_mm * sigma_mm);
    return Math.exp(-0.5 * (R2 / sigma2));
}

// Normalized isotropic Gaussian (used in DoG)
function gaussIso(x_mm, y_mm, sigma_mm) {
    const R2 = x_mm * x_mm + y_mm * y_mm;
    const s2 = Math.max(1e-12, sigma_mm * sigma_mm);
    return Math.exp(-0.5 * R2 / s2) / (2 * Math.PI * s2);
}

// Difference of Gaussians (DoG) unsharp masking
function psfDogUnsharp(x_mm, y_mm, alpha, sigma0_mm, sigma1_mm) {
    return (1.0 + alpha) * gaussIso(x_mm, y_mm, sigma0_mm) -
           alpha * gaussIso(x_mm, y_mm, sigma1_mm);
}

// Normalize PSF
function normalizePsf(psf, background = 0.0, mode = "unit_energy") {
    const H = psf.length;
    const W = psf[0].length;

    if (mode === "unit_energy" || mode === "unit" || mode === "energy") {
        // Calculate sum of signal (psf - background)
        let sum = 0;
        for (let i = 0; i < H; i++) {
            for (let j = 0; j < W; j++) {
                sum += psf[i][j] - background;
            }
        }

        if (sum > 0) {
            const result = create2DArray(H, W);
            for (let i = 0; i < H; i++) {
                for (let j = 0; j < W; j++) {
                    result[i][j] = (psf[i][j] - background) / sum + background;
                }
            }
            return result;
        }
    } else if (mode === "peak_one" || mode === "peak" || mode === "max1" || mode === "max") {
        // Find maximum value
        let max = -Infinity;
        for (let i = 0; i < H; i++) {
            for (let j = 0; j < W; j++) {
                if (psf[i][j] > max) max = psf[i][j];
            }
        }

        if (max > 0) {
            const result = create2DArray(H, W);
            for (let i = 0; i < H; i++) {
                for (let j = 0; j < W; j++) {
                    result[i][j] = psf[i][j] / max;
                }
            }
            return result;
        }
    }

    return psf;
}

// Main PSF computation function
function computePsf(params) {
    let size = parseInt(params.size);
    if (size % 2 === 0) size += 1; // Ensure odd size

    const H = size;
    const W = size;
    const t_i = (H - 1) / 2.0;
    const t_j = (W - 1) / 2.0;
    const dx_mm = params.wire_fov_mm / size;

    let psf = create2DArray(H, W);

    // Compute PSF based on mode
    if (params.mode === "dog") {
        // Difference of Gaussians
        for (let i = 0; i < H; i++) {
            for (let j = 0; j < W; j++) {
                const x_mm = (j - t_j) * dx_mm;
                const y_mm = (i - t_i) * dx_mm;
                psf[i][j] = psfDogUnsharp(x_mm, y_mm, params.alpha, params.sigma0_mm, params.sigma1_mm);
            }
        }
    } else if (params.mode === "ring") {
        // Ring mode
        for (let i = 0; i < H; i++) {
            for (let j = 0; j < W; j++) {
                const x_mm = (j - t_j) * dx_mm;
                const y_mm = (i - t_i) * dx_mm;
                const ring = psfGaussianRing(x_mm, y_mm, params.r0_mm, params.sigma_mm);
                const bump = psfGaussianIsotropic(x_mm, y_mm, params.core_sigma_mm);
                psf[i][j] = (1.0 - params.core_weight) * ring + params.core_weight * bump;
            }
        }
    } else {
        // Gaussian mode
        for (let i = 0; i < H; i++) {
            for (let j = 0; j < W; j++) {
                const core = Math.exp(-params.p_i * Math.pow(i - t_i, 2)) *
                            Math.exp(-params.p_j * Math.pow(j - t_j, 2));
                psf[i][j] = core / (params.p_i * params.p_j);
            }
        }
    }

    // Add background
    for (let i = 0; i < H; i++) {
        for (let j = 0; j < W; j++) {
            psf[i][j] += params.background;
        }
    }

    // Normalize
    psf = normalizePsf(psf, params.background, params.normalize);

    return { psf, dx_mm, size };
}

// FFT implementation using a simple 2D FFT
// For simplicity, we'll use a straightforward approach
function fft2(data) {
    const N = data.length;
    const M = data[0].length;

    // Allocate output array (complex numbers represented as [real, imag])
    const result = Array.from({ length: N }, () =>
        Array.from({ length: M }, () => [0, 0])
    );

    // Apply 1D FFT to each row
    for (let i = 0; i < N; i++) {
        const row = data[i].map(val => [val, 0]);
        const fftRow = fft1d(row);
        for (let j = 0; j < M; j++) {
            result[i][j] = fftRow[j];
        }
    }

    // Apply 1D FFT to each column
    for (let j = 0; j < M; j++) {
        const col = result.map(row => row[j]);
        const fftCol = fft1d(col);
        for (let i = 0; i < N; i++) {
            result[i][j] = fftCol[i];
        }
    }

    return result;
}

// 1D FFT (Cooley-Tukey algorithm)
function fft1d(x) {
    const N = x.length;

    if (N <= 1) return x;

    // Check if N is power of 2, if not, pad with zeros
    const nextPow2 = Math.pow(2, Math.ceil(Math.log2(N)));
    if (N !== nextPow2) {
        const padded = [...x];
        while (padded.length < nextPow2) {
            padded.push([0, 0]);
        }
        return fft1d(padded).slice(0, N);
    }

    // Divide
    const even = fft1d(x.filter((_, idx) => idx % 2 === 0));
    const odd = fft1d(x.filter((_, idx) => idx % 2 === 1));

    // Conquer
    const result = Array(N);
    for (let k = 0; k < N / 2; k++) {
        const angle = -2 * Math.PI * k / N;
        const tReal = Math.cos(angle) * odd[k][0] - Math.sin(angle) * odd[k][1];
        const tImag = Math.cos(angle) * odd[k][1] + Math.sin(angle) * odd[k][0];

        result[k] = [even[k][0] + tReal, even[k][1] + tImag];
        result[k + N / 2] = [even[k][0] - tReal, even[k][1] - tImag];
    }

    return result;
}

// FFT shift (move zero-frequency component to center)
function fftshift(data) {
    const N = data.length;
    const M = data[0].length;
    const halfN = Math.floor(N / 2);
    const halfM = Math.floor(M / 2);

    const result = create2DArray(N, M);

    for (let i = 0; i < N; i++) {
        for (let j = 0; j < M; j++) {
            const newI = (i + halfN) % N;
            const newJ = (j + halfM) % M;
            result[newI][newJ] = data[i][j];
        }
    }

    return result;
}

// Inverse FFT shift
function ifftshift(data) {
    const N = data.length;
    const M = data[0].length;
    const halfN = Math.ceil(N / 2);
    const halfM = Math.ceil(M / 2);

    const result = Array.from({ length: N }, (_, i) =>
        Array.from({ length: M }, (_, j) => {
            const newI = (i + halfN) % N;
            const newJ = (j + halfM) % M;
            return Array.isArray(data[newI][newJ]) ? [...data[newI][newJ]] : data[newI][newJ];
        })
    );

    return result;
}

// Compute MTF from PSF
function computeMtf(psf, dx, size) {
    const half = Math.floor(size / 2);

    // Normalize PSF to unit energy
    let sum = 0;
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
            sum += psf[i][j];
        }
    }

    const psf_n = create2DArray(size, size);
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
            psf_n[i][j] = psf[i][j] / sum;
        }
    }

    // Apply ifftshift
    const psf_shifted = ifftshift(psf_n);

    // Compute 2D FFT
    const fftResult = fft2(psf_shifted);

    // Apply fftshift to center zero frequency
    const fftShifted = fftshift(fftResult);

    // Compute magnitude (MTF)
    const mtf2d = create2DArray(size, size);
    for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
            const real = fftShifted[i][j][0];
            const imag = fftShifted[i][j][1];
            mtf2d[i][j] = Math.sqrt(real * real + imag * imag);
        }
    }

    // Normalize so DC = 1
    const dc = mtf2d[half][half];
    if (dc > 0) {
        for (let i = 0; i < size; i++) {
            for (let j = 0; j < size; j++) {
                mtf2d[i][j] /= dc;
            }
        }
    }

    // Create frequency arrays (cycles/mm)
    const fx = Array(size);
    const fy = Array(size);
    for (let i = 0; i < size; i++) {
        const freq = (i - half) / (size * dx);
        fx[i] = freq;
        fy[i] = freq;
    }

    // Extract radial MTF (from center to right edge)
    const f_rad = [];
    const mtf_rad = [];
    for (let j = half; j < size; j++) {
        f_rad.push(fx[j]);
        mtf_rad.push(mtf2d[half][j]);
    }

    return { fx, fy, f_rad, mtf_rad, mtf2d };
}

// Calculate MTF50 (frequency at which MTF = 0.5)
function calculateMtf50(f_rad, mtf_rad) {
    for (let i = 0; i < mtf_rad.length - 1; i++) {
        if (mtf_rad[i] >= 0.5 && mtf_rad[i + 1] < 0.5) {
            // Linear interpolation
            const t = (0.5 - mtf_rad[i]) / (mtf_rad[i + 1] - mtf_rad[i]);
            return f_rad[i] + t * (f_rad[i + 1] - f_rad[i]);
        }
    }
    return null;
}
