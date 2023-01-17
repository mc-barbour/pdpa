import numpy as np

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots


from source.loader import *
from source.experiment import *
from source.processing import *

pio.renderers.default = "browser"


colors = [px.colors.qualitative.Dark24[23],  px.colors.qualitative.Dark24[22],
          px.colors.qualitative.Dark2[7], px.colors.qualitative.Vivid[0]]


"""
Metrics that are computed:
    VFD - volume flux density of different particle classes at each ring location
        the total volume of particles (per class and measurement zone) divided by the aquisition time and the probe surface area
        It's currently called MF

"""


#%% Load validated_data

data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/20230103_Downstream/'
data_dir = "../Data/20230103_Downstream/"
experiment_names = sorted(get_experiment_names(data_dir))

experiments = []

for name in experiment_names:

    experiments.append(load_processed_data(name, data_dir))

n_exp = len(experiments)

exp = experiments[0]
#%% Experimental Constants

d_needle_liq = 0.016*2.54e-2
od_orrifice = 3./32. * 2.54e-2
d_needle_air = 0.028*2.54e-2

A_air = np.pi / 4. * (od_orrifice**2 - d_needle_air**2)
A_water = np.pi / 4. * d_needle_liq**2

Q_water = 6 #mL/min
V_water = Q_water/1e6/60/A_water

Q_air = exp.estimated_flowrate # LPM
V_air = Q_air / 1e3 / 60 / A_air

rho_air = 1.185
rho_water = 1000
M = rho_air * V_air**2 / (rho_water * V_water**2)
print("Momentum ratio = {:3f}".format(M))

slot = 150e-6 / np.sin(150. * np.pi/180) * (750/250) # is this eperiment dependent?

maxsize = 150
num_dbin = 51
median_size = [5, 25]

#%% median droplet size


#%% Calculate measurement annulii

measurement_area = slot*1e-3 #mm
radial_pos = []

for exp in experiments[0:10]:
    radial_pos.append(exp.radial_position * 1e-3)
    
num_rad_pos = len(radial_pos)
annulus_area = np.zeros(num_rad_pos)


delta_r = abs(radial_pos[1] - radial_pos[0])

for count, pos in enumerate(radial_pos):
    if count == 0:
        annulus_area[count] = np.pi * (radial_pos[count] + delta_r / 2.)**2
    else:
        annulus_area[count] = np.pi * ((radial_pos[count] + delta_r / 2.)**2 
            - (radial_pos[count] - delta_r / 2.)**2)

total_spray_area = np.pi*(radial_pos[-1] + delta_r/2)**2
     
print(np.sum(annulus_area), total_spray_area)   

#%% Determine beam dependent beam width
"""
Transcibing code from kee on. Need to upgdate variable names.

devide by zero at end - nans


"""

bins = np.linspace(0, maxsize, num_dbin)

num_particles_bin = np.zeros(len(bins))
trav_len_particles_bin = np.zeros(len(bins))
trav_time_bin = np.zeros(len(bins))
Lm = 0

ct = 0;
dm = 0; #total diameteres
nd = 0; #total number of particles
ttot = 0; #total acquisition time
nnn = 0; #number of particles at median size
ln = 0; #traverse length at median size
vxm = 0; #velocities at median size
txm = 0; #gate times at median size

for count, exp in enumerate(experiments[0:10]):
    
    valid_int = np.where(exp.valid == 1)[0]
    
    lx = abs(exp.gate_time[valid_int] * 1e-6 * exp.velocity_chan1[valid_int]) # gate time * velocity
    
    T = exp.time_chan1[valid_int][-1]    #T=rawdata.Data{FN}.T1(keepD(end)); % length of measurement time (s)
    d = exp.scaled_diameter[valid_int]
    vx = exp.velocity_chan1[valid_int]    
    tx = exp.gate_time[valid_int] * 1e-6

    # finding specifically for 95 and 115 um (m - median?)    
    ind = np.where((d > median_size[0]) & (d < median_size[1]))
    
    if len(ind)==0:
        print("Warning, not particles found")
    else:
        ln = ln + np.nansum(vx[ind] * tx[ind])
        vxm = vxm + np.nansum(vx)
        txm = txm + np.nansum(tx)
        nnn = nnn + len(ind)
    
    tx = exp.gate_time[valid_int] * 1e-6
    
    bin_count, bin_edge = np.histogram(d, bins)
    digitized_diameter = np.digitize(d, bins)
    
    for j in range(num_dbin - 1):
        bin_ind = np.where(digitized_diameter==j+1)[0]
        num_particles_bin[j] = num_particles_bin[j] + len(bin_ind) * annulus_area[count] / measurement_area
        trav_len_particles_bin[j] = trav_len_particles_bin[j] + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area
        trav_time_bin[j] = trav_time_bin[j] +  np.nansum(tx[bin_ind]) * annulus_area[count] / measurement_area
        
        dm = dm + np.nansum(d[bin_ind]) * annulus_area[count] / measurement_area
        Lm = Lm + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area
        nd = nd + len(bin_ind) * annulus_area[count] / measurement_area
        #ttot = ttot + T
        
dm = dm / nd
LN = ln / nnn
Lm = Lm/np.nansum(num_particles_bin) # single value
trav_len_particles_bin = trav_len_particles_bin / num_particles_bin
trav_time_bin = trav_time_bin / num_particles_bin



#%%    

    for i=1:numel(n)-1 % for each bin
        ind=find(nn==i); % find particles in bin
        Nnn(i)=Nnn(i)+numel(ind).*annuliA(ct)./measA; % summing numbers of particles per bin normalized by area                 
        L(i)=L(i)+nansum((lx(ind))).*annuliA(ct)./measA; % summing crossing distances per bin 
        Tx(i)=Tx(i)+nansum(tx(ind)).*annuliA(ct)./measA; % summing gate times per bin
        dm=dm+nansum(d(ind)).*annuliA(ct)./measA; % summing total diameters 
        Lm=Lm+nansum((lx(ind))).*annuliA(ct)./measA; % summing crossing distances 
        nd=nd+numel(ind).*annuliA(ct)./measA; % summing number of particles 
        ttot=ttot+T; %total sampling time
    end

end
dm=dm./nd; %dm is total diameter/number of particles   
LN=ln./nnn;%vxm.*txm./nnn.^2; % median traverse length / median number of particles 
Lm=(Lm./nansum(Nnn)); % summmed crossing distances / total number of particles found
L=L./Nnn; % summed crossing distances per bin / total number of particles found per bin
Tx=Tx./Nnn; %summed crossing time per bin / total number of particles found p


#%% Calculate Parameter "f"

bin_centers = bins + (bins[1] - bins[0]) / 2.
bin_centers = bin_centers[0:-1]

L1 = lambda bin_centers: (30 + 2.25 * bin_centers / (1 + 0.15 * np.exp(0.005*bin_centers)))

ind = np.where(abs(bin_centers - dm) == min(abs(bin_centers - dm)))[0]
f = (L1(bin_centers) - L1(bin_centers[ind])) / L1(bin_centers[ind])


#%% initialize Dxx arrays
D30 = np.zeros(num_rad_pos);
D10 = np.zeros(num_rad_pos);
D32 = np.zeros(num_rad_pos);
D43 = np.zeros(num_rad_pos);

L = np.zeros(num_dbin)
Tx = np.zeros(num_dbin)
NN = np.zeros(num_dbin)
Nf = np.zeros((len(annulus_area), len(bins)-1))
Nfreal = Nf


#%% Calculate Nf - what is Nf?? - don't actually use
"""
        % tally number of droplets found per position per diameter class
        % normalized by measurement time, probe area and diameter
        % correction function 'f'
        % Nfreal is the real count (you will see that sometimes 1 particle
        % in the data becomes 1e^4 in the Nf
        
"""

count = 0
slit = 150e-6  # LM should be a single value, not an array
A = Lm * slit * np.sin(70 * np.pi / 180) * 750. / 250.
V = A * Lm

for count, exp in enumerate(experiments[0:10]):

    
    valid_int = np.where(exp.valid == 1)[0]
    
    lx = abs(exp.gate_time[valid_int] * 1e-6 * exp.velocity_chan1[valid_int]) # gate time * velocity
    
    T = exp.time_chan1[valid_int][-1]   #T=rawdata.Data{FN}.T1(keepD(end)); % length of measurement time (s)
    d = exp.scaled_diameter[valid_int]
    vx = exp.velocity_chan1[valid_int]    
    tx = exp.gate_time[valid_int] * 1e-6
    

    bin_count, bin_edge = np.histogram(d, bins)
    digitized_diameter = np.digitize(d, bins)


    for j in range(num_dbin - 1):
        bin_ind = np.where(digitized_diameter==j+1)[0]
        Nfreal[count,j] = Nfreal[count,j] + len(bin_ind)
        Nf[count, j] = Nf[count, j] + len(bin_ind) / T/ A/ (1+f[j])
    
        if len(bin_ind) == 0:
            Nf[count, j] = 0
        

#%% Calculate NNft
# Need to double check the diameter / bin count

count = 0
slit = 150e-6 
A = Lm * slit * np.sin(150 * np.pi / 180) * 750. / 250.
V = A * Lm

VFD = np.zeros((len(annulus_area), len(bin_centers)))
VF = np.zeros((len(annulus_area), len(bin_centers)))
NNft = np.zeros(len(bins)-1)

for count, exp in enumerate(experiments[0:10]):
    
    valid_int = np.where(exp.valid == 1)[0]
    d = exp.scaled_diameter[valid_int]
    T = exp.time_chan1[valid_int][-1]
    tx = exp.gate_time[valid_int] * 1e-6
    

    bin_count, bin_edge = np.histogram(d, bins)
    digitized_diameter = np.digitize(d, bins)
    
    for d_bin in range(num_dbin - 1):
        bin_ind = np.where(digitized_diameter==d_bin+1)[0]
        
        print(bin_ind, d_bin)
        if len(bin_ind) == 0:
            tx_bin_mean = 0
        else:
            tx_bin_mean = trav_time_bin[d_bin]
        
        VFD[count, d_bin] = np.pi / 6 * (len(bin_ind) * (bins[d_bin]*1e-6)**3) / T / A / (1+f[d_bin]) 
        NNft[d_bin] = NNft[d_bin] + len(bin_ind) / T / A / (1+f[d_bin])*annulus_area[count] / total_spray_area # area might be wrong        
        
        VF[count, d_bin] = len(bin_ind) * (bins[d_bin]*1e-6)**3 * tx_bin_mean / T / V / (1.+f[d_bin])
        # Vf(ct,DBin)=(numel(indD).*(bins(DBin).*1e-6).^3).*Tx(DBin)'./T./V./(1+f(DBin)).^2
       
#%% pdf

pdfN = NNft / np.trapz(NNft, x=bin_centers)
fig = go.Figure()
fig.add_trace(go.Scatter(y=pdfN, x=bin_centers))
fig.show()



#%% integrate VFD across annulus
    
VF_annulus = np.zeros(len(bins)-1)
VFD_annulus = np.zeros(len(bins)-1)
for bin_count in range(num_dbin - 1):
    VF_annulus[bin_count] = np.sum(VF[:, bin_count] * annulus_area)/total_spray_area
    VFD_annulus[bin_count] = np.sum(VFD[:, bin_count] * annulus_area)/total_spray_area

#%% particle volumes

bin_volume = NNft * np.pi / 6 * (bins[0:-1]*1e-6)**3 # particle volume per unit area and time - integrated over spray


fig = go.Figure()
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(bin_volume)/np.sum(bin_volume)))
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VFD_annulus)/np.sum(VFD_annulus) ))
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VF_annulus)/np.sum(VF_annulus)))

fig.update_yaxes(title='CDF(VDF) - Volume density flux')
fig.update_xaxes(title='Particle Diameter (um)')
fig.update_layout(title="Tidal Volume = 750mL, Q_syr = 6ml/min")
fig.show()



#%% Create function

def full_spray_analysis(experiments):
    
    
    slot = 150e-6 / np.sin(150. * np.pi/180) * (750/250) # is this eperiment dependent?

    maxsize = 150
    num_dbin = 51
    median_size = [5, 25]

    measurement_area = slot*1e-3 #mm
    radial_pos = []


    ### Define spray radius and annulus area
    
    for exp in experiments:
        radial_pos.append(exp.radial_position * 1e-3)
        
    num_rad_pos = len(radial_pos)
    annulus_area = np.zeros(num_rad_pos)


    delta_r = abs(radial_pos[1] - radial_pos[0])

    for count, pos in enumerate(radial_pos):
        if count == 0:
            annulus_area[count] = np.pi * (radial_pos[count] + delta_r / 2.)**2
        else:
            annulus_area[count] = np.pi * ((radial_pos[count] + delta_r / 2.)**2 
                - (radial_pos[count] - delta_r / 2.)**2)

    total_spray_area = np.pi*(radial_pos[-1] + delta_r/2)**2
         
    print(np.sum(annulus_area), total_spray_area)   
    

        
    """
    Transcibing code from kee on. Need to upgdate variable names.
    
    devide by zero at end - nans
    
    
    """
    
    bins = np.linspace(0, maxsize, num_dbin)
    
    num_particles_bin = np.zeros(len(bins))
    trav_len_particles_bin = np.zeros(len(bins))
    trav_time_bin = np.zeros(len(bins))
    Lm = 0
    
    ct = 0;
    dm = 0; #total diameteres
    nd = 0; #total number of particles
    ttot = 0; #total acquisition time
    nnn = 0; #number of particles at median size
    ln = 0; #traverse length at median size
    vxm = 0; #velocities at median size
    txm = 0; #gate times at median size
    
    for count, exp in enumerate(experiments):
        
        valid_int = np.where(exp.valid == 1)[0]
        
        lx = abs(exp.gate_time[valid_int] * 1e-6 * exp.velocity_chan1[valid_int]) # gate time * velocity
        
        T = exp.time_chan1[valid_int][-1]    #T=rawdata.Data{FN}.T1(keepD(end)); % length of measurement time (s)
        d = exp.scaled_diameter[valid_int]
        vx = exp.velocity_chan1[valid_int]    
        tx = exp.gate_time[valid_int] * 1e-6
    
        # finding specifically for 95 and 115 um (m - median?)    
        ind = np.where((d > median_size[0]) & (d < median_size[1]))
        
        if len(ind)==0:
            print("Warning, not particles found")
        else:
            ln = ln + np.nansum(vx[ind] * tx[ind])
            vxm = vxm + np.nansum(vx)
            txm = txm + np.nansum(tx)
            nnn = nnn + len(ind)
        
        tx = exp.gate_time[valid_int] * 1e-6
        
        bin_count, bin_edge = np.histogram(d, bins)
        digitized_diameter = np.digitize(d, bins)
        
        for j in range(num_dbin - 1):
            bin_ind = np.where(digitized_diameter==j+1)[0]
            num_particles_bin[j] = num_particles_bin[j] + len(bin_ind) * annulus_area[count] / measurement_area
            trav_len_particles_bin[j] = trav_len_particles_bin[j] + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area
            trav_time_bin[j] = trav_time_bin[j] +  np.nansum(tx[bin_ind]) * annulus_area[count] / measurement_area
            
            dm = dm + np.nansum(d[bin_ind]) * annulus_area[count] / measurement_area
            Lm = Lm + np.nansum(lx[bin_ind]) * annulus_area[count] / measurement_area
            nd = nd + len(bin_ind) * annulus_area[count] / measurement_area
            #ttot = ttot + T
            
    dm = dm / nd
    LN = ln / nnn
    Lm = Lm/np.nansum(num_particles_bin) 
    trav_len_particles_bin = trav_len_particles_bin / num_particles_bin
    trav_time_bin = trav_time_bin / num_particles_bin
    
    ### Calculate F parameter - beam length adjustment
    bin_centers = bins + (bins[1] - bins[0]) / 2.
    bin_centers = bin_centers[0:-1]
    
    L1 = lambda bin_centers: (30 + 2.25 * bin_centers / (1 + 0.15 * np.exp(0.005*bin_centers)))
    
    ind = np.where(abs(bin_centers - dm) == min(abs(bin_centers - dm)))[0]
    f = (L1(bin_centers) - L1(bin_centers[ind])) / L1(bin_centers[ind])

    # compute spray parameters
    count = 0
    slit = 150e-6 
    A = Lm * slit * np.sin(150 * np.pi / 180) * 750. / 250.
    V = A * Lm
    
    VFD = np.zeros((len(annulus_area), len(bin_centers)))
    VF = np.zeros((len(annulus_area), len(bin_centers)))
    NNft = np.zeros(len(bins)-1)
    
    for count, exp in enumerate(experiments):
        
        valid_int = np.where(exp.valid == 1)[0]
        d = exp.scaled_diameter[valid_int]
        T = exp.time_chan1[valid_int][-1]
        tx = exp.gate_time[valid_int] * 1e-6
        
    
        bin_count, bin_edge = np.histogram(d, bins)
        digitized_diameter = np.digitize(d, bins)
        
        for d_bin in range(num_dbin - 1):
            bin_ind = np.where(digitized_diameter==d_bin+1)[0]
            
            print(bin_ind, d_bin)
            if len(bin_ind) == 0:
                tx_bin_mean = 0
            else:
                tx_bin_mean = trav_time_bin[d_bin]
            
            VFD[count, d_bin] = np.pi / 6 * (len(bin_ind) * (bins[d_bin]*1e-6)**3) / T / A / (1+f[d_bin]) 
            NNft[d_bin] = NNft[d_bin] + len(bin_ind) / T / A / (1+f[d_bin])*annulus_area[count] / total_spray_area # area might be wrong        
            
            VF[count, d_bin] = len(bin_ind) * (bins[d_bin]*1e-6)**3 * tx_bin_mean / T / V / (1.+f[d_bin])
    


    ### compute metrics to return
    
    pdfN = NNft / np.trapz(NNft, x=bin_centers)
    bin_volume = NNft * np.pi / 6 * (bins[0:-1]*1e-6)**3

    ### integrate VFD across annulus
        
    VF_annulus = np.zeros(len(bins)-1)
    VFD_annulus = np.zeros(len(bins)-1)
    for bin_count in range(num_dbin - 1):
        VF_annulus[bin_count] = np.sum(VF[:, bin_count] * annulus_area)/total_spray_area
        VFD_annulus[bin_count] = np.sum(VFD[:, bin_count] * annulus_area)/total_spray_area

    return pdfN, VF_annulus, VFD_annulus, bin_volume, total_spray_area



#%% Test function
rad_pos = []
#experiments_60 = 
for count, exp in enumerate(experiments[10::]):
    rad_pos.append(exp.radial_position)
    print(exp._name, count)
    

experiments_60 = np.array(experiments[10::])[np.argsort(rad_pos)]

for count, exp in enumerate(experiments_60):

    print(exp._name, count)
    

#%% Test Function

pdf_30, VF_30, VFD_30, bin_volume_30, total_spray_area_30 = full_spray_analysis(experiments[0:10])
pdf_60, VF_60, VFD_60, bin_volume_60, total_spray_area_60 = full_spray_analysis(experiments_60)

#%%
fig = go.Figure()
# fig.add_trace(go.Scatter(x=bins, y=np.cumsum(bin_volume_30)/np.sum(bin_volume_30)))
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VFD_30)/np.sum(VFD_30) , name="y=30mm",
                         mode='markers+lines', marker_color=colors[0]))
# fig.add_trace(go.Scatter(x=bins, y=np.cumsum(bin_volume_60)/np.sum(bin_volume_60)))
fig.add_trace(go.Scatter(x=bins, y=np.cumsum(VFD_60)/np.sum(VFD_60) , name="y=60mm",
                         mode='markers+lines', marker_color=colors[2]))

fig.add_hline(y=0.9, line_color='black', line_dash='dashdot')
fig.update_yaxes(title='CDF(VDF) - Volume density flux')
fig.update_xaxes(title='Particle Diameter (um)')
fig.update_layout(title="Tidal Volume = 750mL, Q_syr = 6ml/min")
fig.show()


#%% plot VFD

fig = go.Figure()
# fig.add_trace(go.Scatter(x=bins, y=np.cumsum(bin_volume_30)/np.sum(bin_volume_30)))
fig.add_trace(go.Scatter(x=bins, y=(VFD_30) , name="y=30mm",
                         mode='markers+lines', marker_color=colors[0]))
# fig.add_trace(go.Scatter(x=bins, y=np.cumsum(bin_volume_60)/np.sum(bin_volume_60)))
fig.add_trace(go.Scatter(x=bins, y=VFD_60 , name="y=60mm",
                         mode='markers+lines', marker_color=colors[2]))


fig.update_yaxes(title='VFD - Volume flux density')
fig.update_xaxes(title='Particle Diameter (um)')
fig.update_layout(title="Tidal Volume = 750mL, Q_syr = 6ml/min")
fig.show()


#%% Integrate VDF- droplet flowrate

VFD_full_spray_30 = np.sum(VFD_30)
VFD_full_spray_60 = np.sum(VFD_60)

print(VFD_full_spray_30*total_spray_area_30*60/1e-6*10, VFD_full_spray_60*total_spray_area_60*60/1e-6*10)

injection_flow = np.zeros(200)
injection_flow[0:25] = 6
print(np.trapz(injection_flow, dx=1e-2))



#%% transistion function to seperate module
data_dir = '/Users/mbarbour/OneDrive - UW/Downstream/Data/20230103_Downstream/'
data_dir = "../Data/20230103_Downstream/"
experiment_names = sorted(get_experiment_names(data_dir))

experiments = []

for name in experiment_names:

    experiments.append(load_processed_data(name, data_dir))

n_exp = len(experiments)

exp = experiments[0]

#%%
from source.processing import *

params = {"max_diameter": 150,
          "slot_area": 150e-6 / np.sin(150. * np.pi/180) * (750/250),
          "median_particle_diameter": [10,20],
          "number of diameter bins": 51}


pdf_30, VF_30, VFD_30, bin_volume_30, total_spray_area_30, bins, diam_dict_30 = process_full_spray(experiments[0:10], params)
pdf_60, VF_60, VFD_60, bin_volume_60, total_spray_area_60, bins, diam_dict_60 = process_full_spray(experiments_60, params)

#%%
VFD_full_spray_30 = np.sum(VFD_30)
VFD_full_spray_60 = np.sum(VFD_60)

print(VFD_full_spray_30*total_spray_area_30*60/1e-6*5, VFD_full_spray_60*total_spray_area_60*60/1e-6*5)

#%% plot the pdf

fig = go.Figure()
# fig.add_trace(go.Scatter(x=bins, y=np.cumsum(bin_volume_30)/np.sum(bin_volume_30)))
fig.add_trace(go.Scatter(x=bins, y=pdf_30, name="y=30mm",
                         mode='markers+lines', marker_color=colors[0]))
# fig.add_trace(go.Scatter(x=bins, y=np.cumsum(bin_volume_60)/np.sum(bin_volume_60)))
fig.add_trace(go.Scatter(x=bins, y=pdf_60 , name="y=60mm",
                         mode='markers+lines', marker_color=colors[2]))


fig.update_yaxes(title='PDF(d)')
fig.update_xaxes(title='Particle Diameter (um)')
fig.update_layout(title="Tidal Volume = 750mL, Q_syr = 6ml/min")
fig.show()
