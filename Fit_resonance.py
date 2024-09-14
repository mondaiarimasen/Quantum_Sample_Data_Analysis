'''
Plots and fits the absorption dip, measures center wavelength (nm),
linewidth (GHz), Absorption %, and converts the x-axis to GHz.
we use remove_noise, sweep, & time_convert_wavelength 
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import remove_noise
import sweep
import time_convert_wavelength

c = 299792458.

timeseries = remove_noise.time
voltages = remove_noise.V
peakVoltage = max(voltages)
lowVoltage = min(voltages)
dipPercent = 100.*(peakVoltage-lowVoltage)/peakVoltage
voltageStep = remove_noise.resolution/50.

######### Fano Resonance fitting
def fano_fit(params, xvals):
    return -params[0]*((params[4]*params[2]/2 + xvals-params[1])**2 / ((params[2]/2)**2 + (xvals-params[1])**2)) + params[3]
def fano_residual(params, xvals, yvals, errors):
    return (fano_fit(params, xvals)-yvals)/errors
#########

######### Gaussian fitting
def gauss_fit(params, xvals):
    return -params[0]/(params[2]*np.sqrt(2.*np.pi))*np.exp(-(xvals-params[1])**2/(2.*params[2]**2))+params[3]
def gauss_residual(params, xvals, yvals, errors):
    return (gauss_fit(params, xvals)-yvals)/errors
#########



wavelengthList = []

def transformTimeToWavelength(t_series):
    maxSweep = sweep.maxMean
    minSweep = sweep.minMean
    waveMax = time_convert_wavelength.aveMax
    waveMin = time_convert_wavelength.aveMin
    slopeConversion = (waveMax - waveMin) / (maxSweep - minSweep)
    for t in timeseries:
        wavelengthList.append(slopeConversion * (t - maxSweep) + waveMax)

transformTimeToWavelength(timeseries)

######## Fano Resonance fitting parameter guessing
def guessFanoParams():
    baseLine = max(voltages)
    centralVal = np.mean(wavelengthList)
    dataLength = len(wavelengthList)
    #spread = wavelengthList[int(dataLength*0.6)] - wavelengthList[int(dataLength*0.4)]
    spread = 0.3
    areaEstimate = 2.*spread*(baseLine - min(voltages))
    fanoFactor = -0.5
    return [areaEstimate, centralVal, spread, baseLine, fanoFactor]


initial_guess_fano = guessFanoParams()
print("here is initial_guess_fano:")
print(initial_guess_fano)


######## Gaussian fitting parameter guessing
def guessGaussianParams():
    baseLine = max(voltages)
    centralVal = np.mean(wavelengthList)
    dataLength = len(wavelengthList)
    spread = wavelengthList[int(dataLength*0.6)] - wavelengthList[int(dataLength*0.4)]
    areaEstimate = 2.*spread*(baseLine - min(voltages))
    return [areaEstimate, centralVal, spread, baseLine]

initial_guess = guessGaussianParams()
print("here is initial_guess:")
print(initial_guess)
########


fittedParams, covariance, info, message, successCode = optimize.leastsq(fano_residual, initial_guess_fano,
                            args=(wavelengthList, voltages, voltageStep), full_output=1)



if covariance is None:
    print('Fit fail')
    print('Success code:', success)
    print(mesg)
else:
    print('Fit success')
    print(fittedParams)
    fittedParams_err = [np.sqrt(cov[i,i]) for i in range(len(fittedParams))]

    maxVoltage = max(V)
    minVoltage = min(V)
    dipPercentage = 100.*(maxVoltage-minVoltage)/maxVoltage

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x=np.linspace(min(scaledWavelengths), max(scaledWavelengths), 5000)
    xScaled = x-fittedParams[1]
    for i in range(len(xScaled)):
        xScaled[i]= - (c/fittedParams[1]) + (c/(fittedParams[1]-xScaled[i]))
    ax.plot(xScaled, fitfunc(fittedParams, x), 'r-', label='Fit')
    scaledWavelengths = scaledWavelengths - fittedParams[1]
    for i in range(len(scaledWavelengths)):
        scaledWavelengths[i]= - (c/fittedParams[1]) + (c/(fittedParams[1]-scaledWavelengths[i]))
    ax.errorbar(scaledWavelengths, V, dV, fmt='g.', label='Data', markersize=1.2)
    ax.plot(scaledWavelengths, V, 'g-', label='Data', markersize=1.2)

    HWHM = fittedParams[2]*np.sqrt(2.*np.log(2.))
    print('HWHM is: ' + str(fittedParams[2]))
    FWHM = (c/((fittedParams[1]-HWHM)*(10**(-9)))) - (c/((fittedParams[1]+HWHM)*(10**(-9))))
    print(('FWHM is: ' + str(FWHM)))
    textfit = 'λ = %.3f nm\n' 'Linewidth = %.3f GHz\n' 'Depth = %.1f%% \n' 'Q_factor = %.1f\n'\
        % (fittedParams[1],FWHM*(10**(-9)), dipPercentage, Q)
    textfit = 'λ = %.3f nm\n' 'Q = %.5g\n'\
        % (fittedParams[1], Q)
    ax.text(0.32, 0.95,textfit, transform = ax.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='left',fontweight='bold')

    legend_properties = {'weight':'bold'}
    plt.legend(loc='best',fontsize=10)
    plt.legend(prop=legend_properties)
    plt.xticks(fontsize=10,fontweight='bold')
    plt.yticks(fontsize=10,fontweight='bold')
    plt.ylim(0.2, 1.1)
    ax.set_xlabel('Detuning (GHz)',fontsize=10,fontweight='bold')
    ax.set_ylabel('Reflection (a.u.)',fontsize=10,fontweight='bold')
    ax.set_title(remove_noise.title)
    fig.savefig(remove_noise.title + '.eps',dpi=600)
    fig.savefig(remove_noise.title + '.png',dpi=600)
    plt.show()
