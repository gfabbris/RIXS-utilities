from pylab import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from processCCD_image_new import *

def process_emap(scan_list, base_name, file_out, previous_Ein, e_step, elastic_guess,
                 curvature, clean = [10,6,4], fname_list_BG = [], fit_good = 5.0, guess_range = 10):
    
    """
    Process energy map at once. The most problematic section is fitting the elastic peak,
    a good elastic_guess is needed! (Consider using the "raw_1d_plot" function before)
    Once a reasonable elastic guess is given, the function will estimate the position
    for the next file based on the energy variation read from the header. The elastic
    peak fitting is more reliable in the sigma geometry, therefore it's recommended to
    arrange the scan_list to put the sigma scan before the pi scan for a given energy.
    """    
    
    plt.figure()
    if isinstance(scan_list, list) is False:
        tmp = scan_list
        scan_list = []
        scan_list.append(tmp)
        
    for scanno in scan_list:
        print '\n'
        filename = '{}{:04d}'.format(base_name, scanno)
        try:
            ims = CCD([filename], poly_order=2, binpix=1,
                      fname_list_BG = fname_list_BG, exclude=[1,9999])
            print "found {}".format(filename)
        except IOError:
            print "Didn't find {}".format(filename)
            continue
        
        Ein = ims.file_dictionary['beamline_monoenergy']        
        ims.photon_E = Ein
            
        ims.curvature = curvature
        ims.clean_std_new(clean)
        ims.sub_backgrounds()
                
        #Use if only 1 scan exists
        ims.get_specs_error()
        ims.spectrum = ims.specs[0]
        ims.background = ims.BGspecs[0]

        
        #ims.plot_spectrum()        
        
        ## These need to be used if there is more than 1 scan in filename.
        #ims.get_specs()
        #ims.sum_specs()
        
        print "Energy is {:.1f}".format(Ein)
        print "Energy changed by {:.1f}".format(Ein-previous_Ein)
        elastic_guess = elastic_guess - (Ein-previous_Ein)/e_step
        
        elastic = ims.fit_elastic(cen=None, sigma=2, x_window = [elastic_guess-guess_range, elastic_guess+guess_range])
        if np.abs(elastic - elastic_guess) < fit_good:
            print "Sucessful fit: guess = {:.1f}, value = {:.1f}".format(elastic_guess, elastic)
            elastic_guess = elastic
        else:
            print "Fitted value looks wrong. Using guess {:.1f}".format(elastic_guess)
            elastic = elastic_guess    
       
        ims.calibrate(elastic, e_step)
    
        ims.plot_spectrum(label=filename[-4:])
        
        previous_Ein = Ein
        
        ims.fileout = "{}{:04d}.dat".format(file_out,scanno)
        ims.write_file()
    
    plt.legend()
    plt.show()    
    
def avg_single_e(scan_list, base_name, previous_Ein, e_step, elastic_guess,
                 curvature, clean = [10,6,4], bkg = [], file_out = '', fit = True):
    """
    Processes all scan_list files as a single energy.The scans are aligned to each other
    using the correlate_specs, and the resulted spectrum energy loss scale is calibrated
    by fitting the elastic peak.
    """
    
    filename = []
    for scanno in scan_list:
        filename.append('{}{:04d}'.format(base_name, scanno))
        
    
    fname_list_BG = []
    if len(bkg) == 1:
        for i in range(len(filename)):
            fname_list_BG.append(bkg[0])
    else:
        fname_list_BG = bkg
    
    try:
        ims = CCD(filename, poly_order=2, binpix=1,
                  fname_list_BG = fname_list_BG, exclude=[1,9999])
    except IOError:
        print "Didn't find one or more files!"
        return
        
    ims.photon_E = ims.file_dictionary['beamline_monoenergy']
    Ein = ims.photon_E
            
    ims.curvature = curvature
    ims.clean_std_new(clean)
    ims.sub_backgrounds()
    ims.get_specs()

    if ims._verbose:
        plt.figure()
        ims.plot_specs()
        plt.title('Before align!')

    ims.correlate_specs()
    ims.sum_specs()
    
    if fit == True:
        elastic_guess = elastic_guess - (Ein-previous_Ein)/e_step
        
        elastic = ims.fit_elastic(cen=None, sigma=2, x_window = [elastic_guess-10, elastic_guess+10])
        if np.abs(elastic - elastic_guess) < 5.0:
            print "Sucessful fit: guess = {:.1f}, value = {:.1f}".format(elastic_guess, elastic)
        else:
            print "Fitted value looks wrong. Using guess = {}".format(elastic_guess)
            elastic = elastic_guess
    else:
        print "Using guess = {}".format(elastic_guess)
        elastic = elastic_guess
    
    if ims._verbose:
        plt.figure()
        ims.plot_specs()
        plt.title('After align!')
        
    ims.calibrate(elastic, e_step)        
    
    plt.figure()
    ims.plot_spectrum()
    plt.title('Summed spectrum')
    
    if file_out != '':
        ims.fileout = "{}".format(file_out)
        ims.write_file()
        if ims._verbose:
            print('Saved file {}').format(file_out)
        
    plt.show()          
    
    
def plot_emap(epoints, eloss, Mh, Mv, minh = -1e5, minv = -1e5, maxh = -1e5,
              maxv = -1e5, x_min = -1e5, x_max = -1e5, y_min = -1e5,
              y_max = -1e5, file_out = ''):
                  
    """
    Plot the energy maps for LH and LV polarizations. The epoints, eloos, Mh and Mv
    entries can be generated using the "build_M" process below. 
    """

    a = np.copy(eloss)
    eloss = -1*a
    
    x,y = np.meshgrid(epoints,eloss)


    if minh == -1e5:
        minh = np.min(Mh[np.isfinite(Mh)])
    if maxh == -1e5:
        maxh = np.max(Mh[np.isfinite(Mh)])

    if x_min == -1e5:
        x_min = np.min(epoints)
    if x_max == -1e5:
        x_max = np.max(epoints)
        
    if y_min == -1e5:
        y_min = np.min(eloss)
    if y_max == -1e5:
        y_max = np.max(eloss)
    
    plt.figure()
    plt.pcolor(x, y, Mh, vmin = minh, vmax = maxh)
    plt.colorbar()
    plt.ylim(y_min, y_max)
    plt.xlim(x_min, x_max)
    plt.title(r'LTNAO film $\pi$', fontsize = 16)
    plt.xlabel('Incident energy (eV)', fontsize = 14)
    plt.ylabel('Energy loss (eV)', fontsize = 14)
    if file_out != '':
        plt.savefig(file_out + '_LH_emap.png', format='png', dpi=1000)
    plt.show()
        
    if minv == -1e5:
        minv = np.min(Mv[np.isfinite(Mv)])
    if maxv == -1e5:
        maxv = np.max(Mv[np.isfinite(Mv)])
        
    plt.figure()
    plt.pcolor(x, y, Mv, vmin = minv, vmax = maxv)
    plt.colorbar()
    plt.ylim(y_min, y_max)
    plt.xlim(x_min, x_max)
    plt.title(r'LTNAO film $\sigma$', fontsize = 16)
    plt.xlabel('Incident energy (eV)', fontsize = 14)
    plt.ylabel('Energy loss (eV)', fontsize = 14)
    if file_out != '':
        plt.savefig(file_out + '_LV_emap.png', format='png', dpi=1000)
    plt.show()
    
def raw_1d_plot(scan, curvature, clean = [10,6,4]):
    """
    Plot the the raw spec of "scan". It's useful to obtain the first guess for
    the elastic peak when running the "process_emaps".
    """
    
    try:
        ims = CCD(scan, poly_order=2, binpix=1, exclude=[1,9999])
        print "File found!"
    except IOError:
        print "Didn't find {}.".format(scan)
        return
        
    ims.photon_E = ims.file_dictionary['beamline_monoenergy']
            
    ims.curvature = curvature
    ims.clean_std_new(clean)
    ims.get_specs()
    ims.plot_specs()
    
    
def plot_xas(epoints, Mh, Mv, ind, minh = -1e5, minv = -1e5, maxh = -1e5, 
             maxv = -1e5, x_min = -1e5, x_max = -1e5, file_out = ''):
                 
    """
    Calculates and plots the XAS obtained from RIXS
    """

    xasv = np.empty((len(epoints),2))
    xash = np.empty((len(epoints),2))
    
    for i in range(len(epoints)):
        xasv[i,1] = np.sum(Mv[ind,i])
        xash[i,1] = np.sum(Mh[ind,i])
        
    xash[:,0] = epoints
    xasv[:,0] = epoints
        
    mh = np.min(xash[np.isfinite(xash[:,1]),1])
    mv = np.min(xasv[np.isfinite(xasv[:,1]),1])

    plt.figure()
    plt.plot(epoints, xash[:,1] - mh,label = 'LH')
    plt.plot(epoints, xasv[:,1] - mv,label = 'LV')
    plt.legend()
    plt.xlabel('Energy (eV)')
    plt.ylabel('Intensty')

    if file_out != '':
        #savetxt(file_out + '_xas_lh.txt', xash)
        #savetxt(file_out + '_xas_lv.txt', xasv)
        plt.savefig(file_out + '_xas.png', format = 'png', dpi = 1000)
    
    plt.show()
    
def build_M(scan_numbers, base_name, eloss, file_out = ''):
    """
    Generates the epoints, eloss, Mh, and Mv files that contain the incident
    energy points, energy loss points, LH and LV incident x loss RIXS matrices,
    respectively. Saving the matrices to a file saves significant time when
    one intents to replot the RIXS maps multiple times.
    """

    lh = []
    lv = []
    epoints = []
    repeath = []
    repeatv = []
    
    Ein = 0.0
    for i in range(len(scan_numbers)):
        
        ee= []
        rixs = []
        
        scan = open(base_name + '%04d.dat' %(scan_numbers[i]))
            
        alltxt = scan.read()
        
        for line in alltxt.splitlines():
            try:
                variable, valuestr = line.split(' = ')
                
                if variable.lower() == '# monoenergy':
                    Ein_new = float(valuestr)
                
                if variable.lower() == '# idpolar':
                    polar = float(valuestr)
                
            except ValueError:
                try:
                    data = line.split()
                    if data[0] != '#':
                        ee.append(float(data[0]))
                        rixs.append(float(data[1]))
                except ValueError:
                    pass  
        
        f = interp1d(ee, rixs, kind='cubic', bounds_error=False, fill_value=np.nan) # last two added MPMD
        
        
        ind = []
        for i in range(len(epoints)):
            if abs(epoints[i] - Ein_new) < 0.05:
                ind.append(True)
            else:
                ind.append(False)
        
        if True in ind:
            if polar == 0:
                try:
                    i = ind.index(True)
                    lh[i] = lh[i] + f(eloss)
                    repeath[i] = repeath[i] + 1
                except Exception:
                    lh.append(f(eloss))
                    repeath.append(1)
            else:
                try:
                    i = ind.index(True)
                    lv[i] = lv[i] + f(eloss)
                    repeatv[i] = repeatv[i] + 1
                except Exception:
                    lv.append(f(eloss))
                    repeatv.append(1)
        else:
            Ein = Ein_new
            epoints.append(Ein)           
                
            if polar == 0:
                lh.append(f(eloss))
                repeath.append(1)
            else:
                lv.append(f(eloss))
                repeatv.append(1)
   
    for i in range(len(lh)):
        lh[i] = lh[i]/repeath[i]
    
    for i in range(len(lv)):
        lv[i] = lv[i]/repeatv[i]
            
            
    Mh = np.zeros((len(eloss), len(epoints)))
    Mv = np.zeros((len(eloss), len(epoints)))
    
    i=0
    for data in lh:
        Mh[:,i] = data
        i=i+1
        
    i=0
    for data in lv:
        Mv[:,i] = data
        i=i+1
        
        
    if file_out != '':
        savetxt(file_out + '_Mh.txt', Mh)
        savetxt(file_out + '_Mv.txt', Mv)
        savetxt(file_out + '_epoints.txt', epoints)
        savetxt(file_out + '_eloss.txt', eloss)
        
    return Mh, Mv, epoints