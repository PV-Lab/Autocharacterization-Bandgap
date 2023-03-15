# import
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import numpy as np
from sklearn.linear_model import LinearRegression
from shapely.geometry import LineString, Point


def autoextract(data, savepath, ev_range=[1.4, 1.65], smoothing=20, intensity_scale=10000, verbose=False):
    '''
    ===================================
    == AUTOMATIC BAND GAP EXTRACTOR  ==
    ===================================

    Inputs:
        data:              An (n x m) pandas array with n reflectance data points and (m - 1) measured spectra, where m = 0 is the wavelength
        savepath:          Local file path to save Tauc plot images and band gap csv
        ev_range:          List of eVs to run band gap extractor between: [lower, upper]. The tighter the range, the more accurate the band gap calculation
        smoothing:         The degree to smooth the Tauc plots to remove noise that may impact band gap calculation
        intensity_scale:   The value to divide the raw reflectance intensity by to output a decimal E [0,1]
                           For the Resonon PikaL, intensity_scale = 10000
        verbose:           If True, plots the fitted band gap line and Tauc plot inline

    Outputs:
        EG:                An (e x (m-1)) array with e band gaps per spectra if more than one band gap is observed, where (m - 1) are the number of measured spectra
    '''
    # calcualte (F*E)^2 Tauc from reflectivity
    tauc = data.iloc[:, 1:]
    wl = data.iloc[:, 0]
    R = tauc / intensity_scale  # convert from 10,000 percentage points to decimal reflectivity
    k = (1. - R) ** 2  # k=(1-R)^2
    s = 2 * R  # s = 2*R
    F = k / s  # absorpotion coefficient
    ev = 1240. / wl  # calcualte eV from wavelength
    ev.name = 'eV'  # rename column
    tauc = F.mul(ev, axis=0) ** 2  # calculate tauc
    tauc = pd.concat([ev, tauc], axis=1)  # add ev column back to data
    tauc_ev = tauc[(ev >= ev_range[0]) & (ev <= ev_range[1])]  # bound ev range
    tauc_smooth_raw = tauc_ev.copy()
    smooth = scipy.signal.savgol_filter(tauc_ev.iloc[:, 1:], window_length=smoothing, polyorder=3,
                                        axis=0)  # savitsky-golay smoothing

    # upsample the number of datapoints from ~100 to 1000
    upsample = 1000  # number of points to upsample to
    f = scipy.interpolate.interp1d(tauc_smooth_raw.iloc[:, 0], smooth, axis=0)
    ev_upsample = np.linspace(np.max(tauc_smooth_raw.iloc[:, 0]), np.min(tauc_smooth_raw.iloc[:, 0]), upsample)
    tauc_smooth_1 = pd.DataFrame(np.hstack([ev_upsample.reshape(upsample, 1), f(ev_upsample)]),
                                 columns=tauc_smooth_raw.columns.values)
    tauc_smooth = tauc_smooth_1.iloc[::-1].reset_index(drop=True)  # sort ascending eV

    ############################################################################################
    ######################       Recursive Segmentation of Spectra       #######################
    ############################################################################################
    np.random.seed(0)
    bandgaps = []  # list of bandgaps for all spectra
    for i in range(tauc_smooth.shape[1] - 1):
        bandgaps_per_tauc = []
        current_tauc = tauc_smooth.iloc[:, i + 1].name  # name of current tauc spectra
        TAUC_X = np.array(tauc_smooth.iloc[:, 0]).reshape(-1, 1)
        TAUC_Y = np.array(tauc_smooth.iloc[:, i + 1]).reshape(-1, 1)
        X0 = [TAUC_X]  # initialize X values
        Y0 = [TAUC_Y]  # initialize Y values
        target_len = len(X0[0])  # target length to stop recursion
        R_tol = 0.990  # R^2 linear regression fit tolerance for line segments
        X_tol = []  # X segments above R_tol
        Y_tol = []  # Y segments above R_tol
        m = []  # list of slopes
        current_len = 0
        # run recursive segmentation
        while current_len < target_len:
            X = []
            Y = []
            for segX, segY in zip(X0, Y0):
                mid = len(segX) // 2
                # left segments
                X_L = segX[:mid + 1]  # left segment
                Y_L = segY[:mid + 1]  # left segment
                model_L = LinearRegression().fit(X_L, Y_L)
                if model_L.score(X_L, Y_L) >= R_tol:
                    X_tol.append(X_L)
                    Y_tol.append(Y_L)
                    m.append(model_L.coef_.item())
                else:
                    X.append(X_L)
                    Y.append(Y_L)
                # right segments
                X_R = segX[mid:]  # right segment
                Y_R = segY[mid:]  # right segment
                model_R = LinearRegression().fit(X_R, Y_R)
                if model_R.score(X_R, Y_R) >= R_tol:
                    X_tol.append(X_R)
                    Y_tol.append(Y_R)
                    m.append(model_R.coef_.item())
                else:
                    X.append(X_R)
                    Y.append(Y_R)
            X0 = X  # reinit
            Y0 = Y  # reinit
            # count num of element in X_tol. When X_tol == target_len, end recrusion.
            current_len = 0
            medians = []  # get list of all medians to sort list of lists later
            for l in X_tol:
                current_len += len(l)
                medians.append(np.median(l))
        # sort lists of lists based on X_tol order
        sort_mask = np.argsort(medians)
        X_tol_sort = np.array(X_tol, dtype=object)[sort_mask]
        Y_tol_sort = np.array(Y_tol, dtype=object)[sort_mask]
        thetas = np.rad2deg(
            np.arctan(np.array(m)[sort_mask]))  # calcualte inclination angles between segment slopes and x-axis

        # second smooth of tauc plot just to get clear signal for peaks
        peaks = scipy.signal.find_peaks(tauc_smooth.iloc[:, i + 1], height=200,
                                        width=(0, 999999))  # only find peaks with height TPs>20
        peak_height_idx = np.argsort(peaks[1]['peak_heights'])[::-1]  # sort from highest to lowest peak
        working_e = []  # list of all e-values that work to produce a line that intersects x=0 axis
        rmse_list = []  # list of all rmse
        for p_idx in peak_height_idx:
            for e in range(len(X_tol_sort)):
                try:
                    model = LinearRegression().fit(np.vstack([X_tol_sort[-1 - e], X_tol_sort[-2 - e]]), np.vstack(
                        [Y_tol_sort[-1 - e],
                         Y_tol_sort[-2 - e]]))  # take foruth last and fifth last segment to run regression
                    y_fit = model.predict(TAUC_X)
                    tngt = LineString([(np.min(TAUC_X), np.min(y_fit)), (np.max(TAUC_X), np.max(y_fit))])
                    xax = LineString([(np.min(TAUC_X), 0.), (np.max(TAUC_X), 0.)])
                    int_pt = tngt.intersection(xax)
                    Eg = int_pt.x  # bandgap, if no x-intercept, will throw an error

                    upper = int(
                        peaks[0][p_idx] - peaks[1]['widths'][p_idx] / 2)  # index location of peak minus half peak width
                    lower = int(np.abs(y_fit - 0).argmin())  # where y_fit intersects with x-axis

                    if len(peak_height_idx) > 1 and lower < peaks[0][p_idx - 1] and peaks[0][p_idx] > peaks[0][
                        p_idx - 1]:  # if there are multiple peaks
                        lower = int(peaks[0][p_idx - 1] - peaks[1]['widths'][p_idx - 1] / 2)

                    if model.coef_ > 0 and upper > lower:  # positive slope and lower < upper bound
                        rmse = np.sqrt(np.mean(
                            (TAUC_Y[lower:upper] - y_fit[lower:upper]) ** 2))  # find RMSE of fit and true tauc curve
                        rmse_list.append(rmse)
                        working_e.append(e)
                except:
                    pass  # keep looping if unsuccessful

            # only save band gap and tangent line where rmse is lowest
            try:
                best_e = working_e[np.array(rmse_list).argmin()]
                model = LinearRegression().fit(np.vstack([X_tol_sort[-1 - best_e], X_tol_sort[-2 - best_e]]), np.vstack(
                    [Y_tol_sort[-1 - best_e],
                     Y_tol_sort[-2 - best_e]]))  # take fourth last and fifth last segment to run regression
                y_fit = model.predict(TAUC_X)
                tngt = LineString([(np.min(TAUC_X), np.min(y_fit)), (np.max(TAUC_X), np.max(y_fit))])
                xax = LineString([(np.min(TAUC_X), 0.), (np.max(TAUC_X), 0.)])
                int_pt = tngt.intersection(xax)
                Eg = np.round(int_pt.x, 3)  # bandgap
                bandgaps_per_tauc.append(Eg)  # append bandgap

                fig, ax = plt.subplots(figsize=(3, 2))
                for n in range(len(X_tol_sort)):
                    ax.scatter(X_tol_sort[n][0], Y_tol_sort[n][0], marker='o', c='#040357', s=30, zorder=15)
                    ax.scatter(X_tol_sort[n][-1], Y_tol_sort[n][-1], marker='o', c='#040357', s=30, zorder=15)
                plt.plot(TAUC_X, TAUC_Y, c='#040357', lw=2.5, zorder=15)
                plt.plot(TAUC_X, y_fit, c='r', lw=3, alpha=0.75, zorder=20)
                plt.axvline(Eg, c='#e00000', ls='--', lw=1.5, zorder=5)
                plt.ylim([0, np.max(TAUC_Y)])
                plt.title(f'Material Deposit {current_tauc}')
                t = plt.text(Eg, np.max(TAUC_Y) * 0.95, f'E$_g$={Eg}eV', c='#e00000', ha='center', va='top',
                             weight='bold', zorder=10)
                t.set_bbox(dict(facecolor='w', alpha=0.9, lw=0))
                plt.xlabel(r'$hv$ (eV)')
                plt.ylabel(r'$(F(R)*hv)^2$ (a.u.)')
                ax.minorticks_on()
                ax.grid(which='minor', color='gray', linestyle='--', alpha=0.35)
                ax.grid(which='major', color='gray', linestyle='--', alpha=0.35)
                plt.tick_params(axis='both', labelleft=False)
                if len(bandgaps_per_tauc) > 1 and bandgaps_per_tauc[-1] != bandgaps_per_tauc[-2]:
                    plt.savefig(f'{savepath}\\{current_tauc}_bandgap{len(bandgaps_per_tauc)}.png', dpi=300,
                                bbox_inches='tight')
                else:
                    plt.savefig(f'{savepath}\\{current_tauc}.png', dpi=300, bbox_inches='tight')
                if verbose:  # plot
                    plt.show()

            except:  # if the tangent line does not intersect with x-axis => no band gap
                bandgaps_per_tauc.append(np.nan)

        bandgaps.append(np.array(bandgaps_per_tauc))

    for b in range(len(bandgaps)):
        bandgaps[b] = bandgaps[b][np.sort(np.unique(bandgaps[b], return_index=True)[1])]  # remove duplicates

    EG = pd.DataFrame(bandgaps).T.set_axis(tauc_smooth.columns.values[1:], axis=1)
    EG = EG.set_index(f'bandgap{n}' for n in range(EG.shape[0]))
    EG.to_csv(f'{savepath}\\Extracted_Band_Gaps.csv')
    return EG