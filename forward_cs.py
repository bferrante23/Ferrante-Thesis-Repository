import sys
sys.path.append('/raid/class239/bferrante/pyffit-main')

import os
import pyffit
import numpy as np
import matplotlib.pyplot as plt


# -------------------------------- Model parameters -----------------------------
#                    x        y       z        l        w     strike     dip      ss      ds
m_forward = [3.5681, -7.7641, 4.2947, 43.8623, 17.4530, 232.0857, 60.6675, 1.4604, 1.9999]
# -------------------------------------------------------------------------------

# ----------------------------- Plot selection flags ----------------------------
vector            = 0  # Full-resolution vector field
mis               = 1  # Full-resolution misfit
full_forward_data = 1  # Save forward model data to disk
forward           = 0  # Quadtree-subsampled forward model
forward_full      = 1  # Full-resolution forward model
full_data         = 1  # Full-resolution data plot
subsampled_data   = 0  # Subsampled data
# -------------------------------------------------------------------------------


def main():
    try:
        inversion()
    except ValueError:
        print('Value Error occurred. Try it again')
        exit()


def inversion():

    # --- Input files ---
    insar_files = [
        '/raid/class239/bferrante/A2/53/2800/merge/2021028_2024081/unwrap_mask_ll_deramp_cut_D53.grd',
        '/raid/class239/bferrante/A2/54/2800/merge/2023269_2024268/unwrap_mask_ll_deramp_cut_D54.grd',
        '/raid/class239/bferrante/A2/159/800/merge/2023276_2024219/unwrap_mask_ll_deramp_cut_A159.grd',
        '/raid/class239/bferrante/A2/160/800/merge/2023225_2024280/unwrap_mask_ll_deramp_cut_A160.grd',
    ]

    smooth_files = [
        '/raid/class239/bferrante/A2/53/2800/merge/2021028_2024081/unwrap_mask_ll_deramp_filt_cut_crop.grd',
        '/raid/class239/bferrante/A2/54/2800/merge/2023269_2024268/unwrap_mask_ll_deramp_filt_cut_crop.grd',
        '/raid/class239/bferrante/A2/159/800/merge/2023276_2024219/unwrap_mask_ll_deramp_filt_cut_crop.grd',
        '/raid/class239/bferrante/A2/160/800/merge/2023225_2024280/unwrap_mask_ll_deramp_filt_cut_crop.grd',
    ]

    look_dirs = [
        '/raid/class239/bferrante/A2/53/2800/merge',
        '/raid/class239/bferrante/A2/54/2800/merge',
        '/raid/class239/bferrante/A2/159/800/merge',
        '/raid/class239/bferrante/A2/160/800/merge',
    ]

    weights = [1, 1, 1, 1]
    model_dirs = []  # Paths to models, used when sampling_mode == 2

    sampling_mode = 1  # 0: no smoothing; 1: smooth input; 2: sample from forward model

    # --- Geographic parameters ---
    # These must match the setup.py file used to generate the MCMC parameters.
    ref_point = [78.63, 41.26]
    EPSG      = '32644'

    # --- Elastic parameters ---
    poisson_ratio = 0.25
    shear_modulus = 30e9
    lmda  = 2 * shear_modulus * poisson_ratio / (1 - 2 * poisson_ratio)
    alpha = (lmda + shear_modulus) / (lmda + 2 * shear_modulus)

    # --- Quadtree parameters ---
    # These are matched to the setup.py values used for the inversion.
    rms_min      = 0.01
    rms_min2     = 0.01
    nan_frac_max = 0.7
    width_min    = 0.5
    width_max    = 10
    mean_low     = 0.3

    # --- Output ---
    out_dir = 'forward2026cs_test'
    os.makedirs(out_dir, exist_ok=True)

    # --- Prior bounds, for reference only ---
    priors = {
        'x':           [0,  8],
        'y':           [-10,  -5],
        'z':           [2,  7],
        'l':           [35,  52],
        'w':           [14,  24],
        'strike':      [200, 250],
        'dip':         [50, 70],
        'strike_slip': [0,   2],
        'dip_slip':    [0,   4],
    }

    labels = ['x',   'y',  'z',  'l',  'w', 'strike', 'dip', 'strike_slip', 'dip_slip']
    units  = ['km', 'km', 'km', 'km', 'km',    'deg', 'deg',           'm',        'm']
    scales = [1,      1,    1,    1,    1,         1,     1,             1,          1]

    # --- Plotting parameters ---
    vlim_disp = [-0.4, 0.4]

    # -------------------------------------------------------------------------
    # Forward model functions
    # -------------------------------------------------------------------------

    def elliptical_tapering(row, col, N, u0):
        """Scale slip by an elliptical taper across an N by N sub-patch grid."""
        central = (N + 1) / 2 - 1
        axis    = N / 2
        scale1  = np.sqrt(1 - (row - central) ** 2 / axis ** 2)
        scale2  = np.sqrt(1 - (col - central) ** 2 / axis ** 2)
        return np.array(u0) * scale1 * scale2

    def patch_slip_ellip(m, coords, look):
        """
        Compute LOS displacement with elliptical slip taper.
        """
        X, Y, Z, L, W, strike, dip, strike_slip, dip_slip = m
        x, y = coords
        u0 = [strike_slip, dip_slip, 0]

        N       = 7
        central = (N + 1) / 2 - 1
        disp_summed = np.zeros((np.size(x), 3))

        for row in range(N):
            for col in range(N):
                x_sub = (
                    X
                    + (central - col) * L / N * np.sin(np.deg2rad(strike))
                    + row * W / N * np.cos(np.deg2rad(dip)) * np.cos(np.deg2rad(strike))
                )

                y_sub = (
                    Y
                    + (central - col) * L / N * np.cos(np.deg2rad(strike))
                    - row * W / N * np.cos(np.deg2rad(dip)) * np.sin(np.deg2rad(strike))
                )

                z_sub = Z + row * W / N * np.sin(np.deg2rad(dip))

                slip_sub = elliptical_tapering(row, col, N, u0)

                patch_sub = pyffit.finite_fault.Patch()
                patch_sub.add_self_geometry(
                    (x_sub, y_sub, z_sub),
                    strike,
                    dip,
                    L / N,
                    W / N,
                    slip=slip_sub
                )

                disp_summed += patch_sub.disp(x, y, 0, alpha, slip=slip_sub)

        disp_summed = disp_summed.reshape(x.size, 3, 1, 1)
        disp_LOS = pyffit.finite_fault.proj_greens_functions(disp_summed, look)[:, :, 0].reshape(-1)

        if -np.inf in disp_summed:
            return np.ones_like(disp_LOS) * -np.inf

        return disp_LOS

    def patch_slip(m, coords, look):
        """
        Compute LOS displacement for a single uniform-slip fault patch.
        """
        x_patch, y_patch, z, l, w, strike, dip, strike_slip, dip_slip = m
        x, y = coords
        slip = [strike_slip, dip_slip, 0]

        patch = pyffit.finite_fault.Patch()
        patch.add_self_geometry((x_patch, y_patch, z), strike, dip, l, w, slip=slip)

        disp = patch.disp(x, y, 0, alpha, slip=slip).reshape(x.size, 3, 1, 1)
        disp_LOS = pyffit.finite_fault.proj_greens_functions(disp, look)[:, :, 0].reshape(-1)

        if -np.inf in disp:
            print('Error')
            return np.ones_like(disp_LOS) * -np.inf

        return disp_LOS

    def patch_slip_xyz(m, x, y):
        """
        Compute 3-component displacement for a single uniform-slip patch.
        """
        x_patch, y_patch, z, l, w, strike, dip, strike_slip, dip_slip = m
        slip = [strike_slip, dip_slip, 0]

        patch = pyffit.finite_fault.Patch()
        patch.add_self_geometry((x_patch, y_patch, z), strike, dip, l, w, slip=slip)

        disp = patch.disp(x.flatten(), y.flatten(), 0, alpha, slip=slip)

        if -np.inf in disp:
            print('Error')
            return np.ones_like(disp) * -np.inf

        return disp

    # -------------------------------------------------------------------------
    # Load and prepare InSAR data
    # -------------------------------------------------------------------------

    if sampling_mode == 0:
        datasets = pyffit.insar.prepare_datasets(
            insar_files,
            look_dirs,
            weights,
            ref_point,
            EPSG=EPSG,
            rms_min=rms_min,
            nan_frac_max=nan_frac_max,
            width_min=width_min,
            width_max=width_max
        )

    elif sampling_mode == 1:
        datasets = pyffit.insar.prepare_datasets_smooth(
            insar_files,
            smooth_files,
            look_dirs,
            weights,
            ref_point,
            EPSG=EPSG,
            mean_low=mean_low,
            rms_min2=rms_min2,
            rms_min=rms_min,
            nan_frac_max=nan_frac_max,
            width_min=width_min,
            width_max=width_max
        )

    elif sampling_mode == 2:
        datasets = pyffit.insar.prepare_datasets_model(
            insar_files,
            model_dirs,
            look_dirs,
            weights,
            ref_point,
            EPSG=EPSG,
            mean_low=mean_low,
            rms_min2=rms_min2,
            rms_min=rms_min,
            nan_frac_max=nan_frac_max,
            width_min=width_min,
            width_max=width_max
        )

    # -------------------------------------------------------------------------
    # Plot quadtree subsampling for each dataset
    # -------------------------------------------------------------------------

    for dataset in datasets.keys():
        pyffit.figures.plot_quadtree(
            datasets[dataset]['data'],
            datasets[dataset]['extent'],
            (datasets[dataset]['x_samp'], datasets[dataset]['y_samp']),
            datasets[dataset]['data_samp'],
            vlim_disp=vlim_disp,
            file_name=f'{out_dir}/quadtree_init_{dataset}.png',
            cmap_disp='coolwarm',
            trace=[],
            original_data=[],
            cell_extents=[],
            show=False,
            dpi=500,
            figsize=(14, 8.2)
        )

    # -------------------------------------------------------------------------
    # Aggregate subsampled data from all tracks
    # -------------------------------------------------------------------------

    i_nans = np.concatenate([datasets[name]['i_nans'] for name in datasets.keys()])

    coords = (
        np.concatenate([datasets[name]['x_samp'] for name in datasets.keys()])[~i_nans],
        np.concatenate([datasets[name]['y_samp'] for name in datasets.keys()])[~i_nans],
    )

    look = np.concatenate([datasets[name]['look_samp'] for name in datasets.keys()])[~i_nans]
    data = np.concatenate([datasets[name]['data_samp'] for name in datasets.keys()])[~i_nans]
    data_std = np.concatenate([datasets[name]['data_samp_std'] for name in datasets.keys()])[~i_nans]

    weights = np.concatenate([
        np.ones_like(datasets[name]['data_samp']) * datasets[name]['weight']
        for name in datasets.keys()
    ])[~i_nans]

    B = np.diag(weights)
    S_inv = np.diag(data_std ** -2)

    # -------------------------------------------------------------------------
    # Quadtree-subsampled forward model using all tracks together
    # -------------------------------------------------------------------------

    if forward == 1:
        disp_fit = patch_slip(m_forward, coords, look)

        fig, axes = plt.subplots(1, 1, figsize=(14, 8.2))
        sc = axes.scatter(
            coords[0],
            coords[1],
            c=disp_fit,
            cmap='coolwarm',
            vmin=vlim_disp[0],
            vmax=vlim_disp[1]
        )

        cbar = fig.colorbar(sc, label='LOS Displacement (m)')
        axes.set_xlabel('X (km)', fontsize=20)
        axes.set_ylabel('Y (km)', fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=20)
        cbar.ax.tick_params(labelsize=20)

        axes.set_xlim(np.nanmin(coords[0]), np.nanmax(coords[0]))
        axes.set_ylim(np.nanmin(coords[1]), np.nanmax(coords[1]))

        plt.savefig(f'{out_dir}/forward_all_tracks.png', dpi=500)
        plt.close()

        print('Quadtree subsampled forward model computed')

    # -------------------------------------------------------------------------
    # Full-resolution outputs for each individual track
    # -------------------------------------------------------------------------

    for dataset in datasets.keys():

        x = datasets[dataset]['x']
        y = datasets[dataset]['y']
        look_full = datasets[dataset]['look']
        data_full = datasets[dataset]['data']

        coords_full = np.vstack((x.flatten(), y.flatten()))

        # ---------------------------------------------------------------------
        # Full-resolution vector field
        # ---------------------------------------------------------------------

        if vector == 1:
            disp_xyz = patch_slip_xyz(m_forward, x, y)

            n_rows, n_cols = x.shape
            skip = (slice(None, None, 15), slice(None, None, 15))

            u_x = disp_xyz[:, 0].reshape(n_rows, n_cols)
            u_y = disp_xyz[:, 1].reshape(n_rows, n_cols)
            u_z = disp_xyz[:, 2].reshape(n_rows, n_cols)

            fig, axes = plt.subplots(1, 1, figsize=(14, 8.2))
            axes.set_xlabel('X (km)', fontsize=20)
            axes.set_ylabel('Y (km)', fontsize=20)
            plt.quiver(x[skip], y[skip], u_x[skip], u_y[skip])

            axes.set_xlim(np.nanmin(coords_full[0]), np.nanmax(coords_full[0]))
            axes.set_ylim(np.nanmin(coords_full[1]), np.nanmax(coords_full[1]))

            plt.savefig(f'{out_dir}/forward_{dataset}_vec.png', dpi=500)
            plt.close()

            np.savetxt(f'{out_dir}/x_{dataset}.txt',  x,   delimiter=',')
            np.savetxt(f'{out_dir}/y_{dataset}.txt',  y,   delimiter=',')
            np.savetxt(f'{out_dir}/ux_{dataset}.txt', u_x, delimiter=',')
            np.savetxt(f'{out_dir}/uy_{dataset}.txt', u_y, delimiter=',')
            np.savetxt(f'{out_dir}/uz_{dataset}.txt', u_z, delimiter=',')

            print(f'Full-resolution vector field computed for {dataset}')

        # ---------------------------------------------------------------------
        # Full-resolution forward model
        # ---------------------------------------------------------------------

        if forward_full == 1 or mis == 1 or full_forward_data == 1:
            disp_full = patch_slip_ellip(m_forward, coords_full, look_full)

        if forward_full == 1:
            fig, axes = plt.subplots(1, 1, figsize=(14, 8.2))
            axes.set_xlabel('East (km)', fontsize=20)
            axes.set_ylabel('North (km)', fontsize=20)

            sc = axes.scatter(
                coords_full[0],
                coords_full[1],
                c=disp_full,
                cmap='coolwarm',
                vmin=vlim_disp[0],
                vmax=vlim_disp[1]
            )

            cbar = fig.colorbar(sc, label='LOS Displacement (m)')
            cbar.ax.tick_params(labelsize=20)
            plt.tick_params(axis='both', which='major', labelsize=20)

            axes.set_xlim(np.nanmin(coords_full[0]), np.nanmax(coords_full[0]))
            axes.set_ylim(np.nanmin(coords_full[1]), np.nanmax(coords_full[1]))

            plt.savefig(f'{out_dir}/fullforward_{dataset}.png', dpi=500)
            plt.close()

            print(f'Full-resolution forward model computed for {dataset}')

        # ---------------------------------------------------------------------
        # Full-resolution data plot
        # ---------------------------------------------------------------------

        if full_data == 1:
            fig, axes = plt.subplots(1, 1, figsize=(14, 8.2))
            axes.set_xlabel('East (km)', fontsize=20)
            axes.set_ylabel('North (km)', fontsize=20)

            sc = axes.scatter(
                coords_full[0],
                coords_full[1],
                c=data_full.flatten(),
                cmap='coolwarm',
                vmin=vlim_disp[0],
                vmax=vlim_disp[1]
            )

            cbar = fig.colorbar(sc, label='LOS Displacement (m)')
            cbar.ax.tick_params(labelsize=20)
            plt.tick_params(axis='both', which='major', labelsize=20)

            axes.set_xlim(np.nanmin(coords_full[0]), np.nanmax(coords_full[0]))
            axes.set_ylim(np.nanmin(coords_full[1]), np.nanmax(coords_full[1]))

            plt.savefig(f'{out_dir}/full_{dataset}.png', dpi=500)
            plt.close()

            print(f'Full-resolution data plotted for {dataset}')

        # ---------------------------------------------------------------------
        # Full-resolution misfit
        # ---------------------------------------------------------------------

        if mis == 1:
            misfit = data_full.flatten() - disp_full

            np.savetxt(f'{out_dir}/misfit_{dataset}.txt', misfit, delimiter=',')

            fig, axes = plt.subplots(1, 1, figsize=(14, 8.2))
            axes.set_xlabel('East (km)', fontsize=20)
            axes.set_ylabel('North (km)', fontsize=20)

            sc = axes.scatter(
                coords_full[0],
                coords_full[1],
                c=misfit,
                cmap='coolwarm',
                vmin=vlim_disp[0],
                vmax=vlim_disp[1]
            )

            cbar = fig.colorbar(sc, label='LOS Misfit (m)')
            cbar.ax.tick_params(labelsize=20)
            plt.tick_params(axis='both', which='major', labelsize=20)

            axes.set_xlim(np.nanmin(coords_full[0]), np.nanmax(coords_full[0]))
            axes.set_ylim(np.nanmin(coords_full[1]), np.nanmax(coords_full[1]))

            plt.savefig(f'{out_dir}/misfit_{dataset}.png', dpi=500)
            plt.close()

            print(f'Full-resolution misfit computed for {dataset}')

        # ---------------------------------------------------------------------
        # Save full forward model data
        # ---------------------------------------------------------------------

        if full_forward_data == 1:
            np.savetxt(f'{out_dir}/x_{dataset}.txt',    x,         delimiter=',')
            np.savetxt(f'{out_dir}/y_{dataset}.txt',    y,         delimiter=',')
            np.savetxt(f'{out_dir}/disp_{dataset}.txt', disp_full, delimiter=',')

        # ---------------------------------------------------------------------
        # Save subsampled data
        # ---------------------------------------------------------------------

        if subsampled_data == 1:
            np.savetxt(f'{out_dir}/coords_{dataset}.txt', coords, delimiter=',')
            np.savetxt(f'{out_dir}/look_{dataset}.txt',   look,   delimiter=',')
            np.savetxt(f'{out_dir}/data_{dataset}.txt',   data,   delimiter=',')


if __name__ == '__main__':
    main()
