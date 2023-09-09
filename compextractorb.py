# import
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as PathEffects
import cv2

def get_compositions(data, segmentation, motor_speeds, gcode, rasterSpeed, rasterScaleY, rasterScaleX, rasterOffsetX, 
                     rasterOffsetY, savepath, maxComp=100, minComp=0, bufferTimeStart=1, bufferTimeEnd=1):
    ####################
    # Motor Speed Data #
    ####################
    data_split = motor_speeds.split("_")
    data_cleaned = []
    for d in data_split:
        if int(d[1]) == 0:  # if T is 0
            pass
        else:
            data_cleaned.append(d)
    data_cleaned = data_cleaned[:-1]
    T = [0]  # time stamp
    A = [0]  # motor A speed
    B = [0]  # motor B speed
    for d in data_cleaned:  # format data to lists
        i = d.split("T")[1].split("A")
        T.append(int(i[0]))
        i = i[1].split("B")
        A.append(int(i[0]))
        B.append(int(i[1]))
    # fill in data
    motor_speeds = np.zeros((T[-1], 4))
    t0 = T[0]
    i = 0  # list iterator
    for t in range(motor_speeds.shape[0]):  # actual time iterator
        motor_speeds[t, 0] = t
        motor_speeds[t, 1] = A[i]
        motor_speeds[t, 2] = B[i]
        if motor_speeds[t, 0] == T[i]:
            i += 1  # add to list iterator
    # scale values to be % of total composition
    div = motor_speeds[:, 1].max()
    motor_speeds[:, 1] = motor_speeds[:, 1] / div * maxComp
    motor_speeds[:, 2] = motor_speeds[:, 2] / div * maxComp
    # update motor_speeds to add start up and ramp down buffer time
    motor_speeds[:, 2][0] = maxComp  # fix zero at start
    motorA = np.concatenate((np.array([0.0] * int(1000 * bufferTimeStart)), motor_speeds[:, 1]))
    motorA = np.concatenate((motorA, np.array([maxComp] * int(1000 * bufferTimeEnd))))
    motorB = np.concatenate((np.array([maxComp] * int(1000 * bufferTimeStart)), motor_speeds[:, 2]))
    motorB = np.concatenate((motorB, np.array([0.0] * int(1000 * bufferTimeEnd))))
    time = np.arange(0, len(motorA), 1)
    motor_speeds_full = np.array([time, motorA, motorB]).T
    # plot motor speed traces over time
    '''
    plt.figure(figsize=(3, 2))
    plt.plot(motor_speeds_full[:, 0], motor_speeds_full[:, 1], label='Motor B')
    plt.plot(motor_speeds_full[:, 0], motor_speeds_full[:, 2], label='Motor A')
    plt.xlabel('Time [ms]')
    plt.ylabel('Motor Speed Proportion [%]')
    plt.title('Motor Speed Traces')
    plt.grid()
    plt.legend()
    plt.show()
    ''' 

    #####################
    # Gcode Raster Data #
    #####################
    # get XY positions and timestamps from printer gcode
    gcode['sum'] = gcode.iloc[:, 0] + gcode.iloc[:, 1]
    gcode['time'] = np.zeros(gcode.shape[0])
    total_mm = 0
    total_time = 0
    x0 = gcode['X'][0]
    y0 = gcode['Y'][0]
    for n in range(gcode.shape[0] - 1):
        if gcode['X'][n] != gcode['X'][n + 1]:
            d = np.abs(gcode['X'][n + 1] - gcode['X'][n])
            total_mm += d
            total_time += d / rasterSpeed
            gcode['time'][n + 1] = total_time
        if gcode['Y'][n] != gcode['Y'][n + 1]:
            d = np.abs(gcode['Y'][n + 1] - gcode['Y'][n])
            total_mm += d
            total_time += d / rasterSpeed
            gcode['time'][n + 1] = total_time
    gcode['time'] = gcode['time'] * 1000  # convert from s to ms
    max_time_coords = int(gcode['time'][len(gcode['time']) - 1])  # interpolate position
    # upsample raster coordinates
    t = np.arange(0, max_time_coords, 1)  # new resolution
    fx = scipy.interpolate.interp1d(gcode['time'], gcode['X'])
    X_highres = fx(t)
    fy = scipy.interpolate.interp1d(gcode['time'], gcode['Y'])
    Y_highres = fy(t)
    # scale to sample image shape
    Xpath_scale = ((X_highres - X_highres.min()) / (X_highres.max() - X_highres.min()) * segmentation.shape[
        1]) * rasterScaleX + rasterOffsetX
    Ypath_scale = ((Y_highres - Y_highres.min()) / (Y_highres.max() - Y_highres.min()) * segmentation.shape[
        0]) * rasterScaleY + rasterOffsetY
    pixel_error = 15
    Xlower = Xpath_scale - pixel_error
    Xupper = Xpath_scale + pixel_error
    coords_full = np.array([t, Xpath_scale, Ypath_scale, Xlower, Xupper]).T
    # inner join the gcode coordinate location with the composition along the time axis
    df_coords = pd.DataFrame(coords_full, columns=['time', 'X', 'Y', 'Xlow', 'Xup'])
    df_motors = pd.DataFrame(motor_speeds_full, columns=['time', 'A', 'B'])
    motor_coords = pd.merge(df_motors, df_coords, on='time', how='inner')

    ###########################
    # Composition Calculation #
    ###########################
    drop_idx = np.unique(segmentation)[1:]
    drops_m = cv2.dilate(segmentation.copy().astype('uint8'), np.ones((7, 7), np.uint8))
    drops = np.ma.array(drops_m, mask=drops_m == 0)
    cmap = mpl.cm.viridis.copy()
    cmap.set_bad('white', 0.)
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.fill_between([-45, 800], -30, 450, hatch='//////', zorder=0, facecolor='w', ec='lightgray')
    ax.imshow(drops, cmap=cmap, zorder=15)
    ax.plot(Xpath_scale, Ypath_scale, c='k', lw=2, zorder=10, alpha=0.5)
    comp = []
    for d in drop_idx:
        # find x locations
        xpos = np.where(np.sum(segmentation == d, axis=0) != 0)[0]
        xlower = motor_coords['Xlow'] > np.inf  # init with all False
        xupper = xlower.copy()  # init with all False
        for x in xpos:
            xlower = xlower | (x >= motor_coords['Xlow'])  # OR
            xupper = xupper | (x <= motor_coords['Xup'])  # OR
        xloc = xlower & xupper  # AND
        # find y locations
        ypos = np.where(np.sum(segmentation == d, axis=1) != 0)[0]
        ylower = motor_coords['Y'] > np.inf  # init with all False
        yupper = ylower.copy()  # init with all False
        for y in ypos:
            ylower = ylower | (y >= motor_coords['Y'])  # OR
            yupper = yupper | (y <= motor_coords['Y'])  # OR
        yloc = ylower & yupper  # AND
        xstart, xend = np.where(np.sum(segmentation == d, axis=0) != 0)[0][[0, -1]]  # find x extents
        ystart, yend = np.where(np.sum(segmentation == d, axis=1) != 0)[0][[0, -1]]  # find y extents
        xyloc = xloc & yloc  # total boolean mask
        dAdt = np.average(motor_coords['A'][xyloc])  # motor A over time
        dBdt = np.average(motor_coords['B'][xyloc])  # motor B over time
        A = (dBdt / (dAdt + dBdt)).round(4)
        B = (dAdt / (dAdt + dBdt)).round(4)
        comp.append(B)
        composition = f'A{A}B{B}'
        t = ax.text((xstart + xend) / 2, (ystart + yend) / 2,
                    'A%s' % ((A * 100).round(1)) + '\nB%s' % ((B * 100).round(1)), ha='center', fontsize=8,
                    va='center', weight='bold', zorder=20)
        t.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground='w', alpha=0.7)])
    ax.text(motor_coords['X'].min(), motor_coords['Y'].min(), 'START', ha='center', va='top', fontsize=8, c='#e00000',
            alpha=0.5, weight='bold')
    ax.text(motor_coords['X'].max(), motor_coords['Y'].max() * 1.07, 'END', ha='center', va='bottom', fontsize=8,
            c='#e00000', alpha=0.5, weight='bold')
    ax.invert_yaxis()
    plt.axis('off')
    plt.savefig(savepath+'\\composition_map.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # return sorted idex
    sorted_idx = np.array([drop_idx, comp])
    sort = sorted_idx[0][np.argsort(sorted_idx[1])].astype(int)
    comp = np.sort(sorted_idx[1])
    
    # Sort droplets
    seg_idx = np.unique(segmentation)
    seg_idx = seg_idx[seg_idx != 0]
    seg_order = [np.where(seg_idx == i)[0].item() for i in sort]
    vals = data.iloc[:, 1:]
    vals = vals.iloc[:, seg_order]
    vals.columns=np.arange(0,vals.shape[1],1) # reset column names
    vals.insert(loc=0, column='wavelength', value=data.iloc[:,0])
    return sort, comp, vals
